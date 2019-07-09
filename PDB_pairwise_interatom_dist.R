.libPaths(c(
  .libPaths(),
  "/home/bsc51/bsc51927/R/x86_64-pc-linux-gnu-library"
))
##### Load packages
library(stringr)
library(bio3d)
library(tidyr)
library(parallel)
library(veriNA3d)

###################################################################
# Note: to install veriNA3d, since it is in a private repository, #
# follow the instructions explained in this link.                 #
# (R version >= 3.5 is needed):                                   #
#                                                                 #
# ---> Optional: if you have a GitHub lab account, then you can   #
# use the following command in R:                                 #
# 
# devtools::install_git("http://mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d.git",
#                                credentials = git2r::cred_user_pass("USER", getPass::getPass()))
# In USER: type your GitLab user. Password will be                #
# asked afterwards                                                #
# Note: to avoid this problems a docker container will            #
# be availabe soon                                                #
###################################################################

##### Load necessary functions
source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/PDB_pairwise_interaction.R")
source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/calc_interatom_dist.R")


#' Title: CALCULATION INTER-ATOMIC DISTANCES
#'
#' @param biounit
#' @param download
#' @param PDB_ID
#' @param type_of_interaction
#' @param atom1
#' @param atom2
#' @param dist_threshold
#'
#' @return
#' @export
#'
#' @examples
#'
PDB_iter_atom_distances <- function(biounit_filename,
                                    pdb_file,
                                    type_of_interaction = "protein",
                                    atom_temp = "all",
                                    atom_int = "all",
                                    dist_threshold = 5,
                                    output_dir) {
  # Get all chains to make pairs type_of_interactions
  ChainsList <- PDB_pairwise_interaction(pdb_file, type_of_interaction)
  TemplateChains <-   ChainsList[[1]]
  InteractionChains <- ChainsList[[2]]
  names(InteractionChains)[names(InteractionChains)=="real.pos"] = "real.pos.1"
  
  #Select type of atom to calculate distances
  if (atom_temp != "all") {
    TemplateChains    <- subset(TemplateChains, elety == atom_temp)
  }
  if (atom_int != "all") {
    InteractionChains <-
      subset(InteractionChains, elety == atom_int)
  }
  
  # put dataframe into separate list of data frames according to the different chains
  TemplateChains_list <-
    split.data.frame(TemplateChains, f = TemplateChains$chain)
  InteractionChains_list <-
    split.data.frame(InteractionChains, f = InteractionChains$chain)
  
  # Determine pairwise type_of_interaction between chains 
  # depending on the type of interaction
  # If protein interaction we skip the comparison of one chain with itself
  # If ligand or acid nucleic interaction then all chains compared against all.
  if (type_of_interaction == "protein") {
    PairwiseComb <- combn(1:length(TemplateChains_list), 2)
    
  } else if (type_of_interaction == "nucleic" |
             type_of_interaction == "ligand") {
    PairwiseComb <-
      t(expand.grid(
        1:length(TemplateChains_list),
        1:length(InteractionChains_list)
      ))
    
  } else {
    stop("type_of_interaction not found")
  }
  
  # Select the chains to be comapred
  CompareTemplateChains <- TemplateChains_list[PairwiseComb[1, ]]
  CompareInteractionChains <-
    InteractionChains_list[PairwiseComb[2, ]]
  
  CompareTemplateChains.Vector <- lapply(CompareTemplateChains,
                                         function(x) {
                                           x <- x[, c("x", "y", "z")]
                                           x <- as.vector(t(x))
                                           return(x)
                                         })
  CompareInteractionChains.Vector <-
    lapply(CompareInteractionChains,
           function(x) {
             x <- x[, c("x", "y", "z")]
             x <- as.vector(t(x))
             return(x)
           })
  
  distList <- list()
  i = 1
  while (i <= length(CompareInteractionChains.Vector)) {
    try(distList[[paste("iter", i, sep = "")]] <-
          calc_PDB_dist(
            pdb_file,
            CompareTemplateChains[[i]],
            CompareInteractionChains[[i]],
            CompareTemplateChains.Vector[i],
            CompareInteractionChains.Vector[i],
            dist_threshold = 5,
            type_of_interaction = type_of_interaction
          ),
        silent = TRUE)
    # Filter by selected distance threshold
    i = i + 1
  }
  
  
  pdb_atom_inter_df <-
    ldply(distList, function(x)
      data.frame(x))
  
  #add identifier of the PDB structure
  PDB_ID <- sub(".gz+", "", biounit_filename)
  pdb_atom_inter_df$pdb.id <-  PDB_ID
  # Eliminate dummy columns
  if(type_of_interaction =="protein"){
  pdb_atom_inter_df <- pdb_atom_inter_df[, c(
    "type",
    "eleno",
    "elety",
    "chain",
    "resid",
    "resno",
    "type.1",
    "eleno.1",
    "elety.1",
    "chain.1",
    "resid.1",
    "resno.1",
    "distance",
    "interaction",
    "pdb.id",
    "real.pos",
    "real.pos.1"
  )]
  } else {
    pdb_atom_inter_df <- pdb_atom_inter_df[, c(
      "type",
      "eleno",
      "elety",
      "chain",
      "resid",
      "resno",
      "type.1",
      "eleno.1",
      "elety.1",
      "chain.1",
      "resid.1",
      "resno.1",
      "distance",
      "interaction",
      "pdb.id",
      "real.pos"
    )]
  }
  
  # Append results to txt file
  output_filePath <-
    file.path(
      output_dir,
      paste(
        PDB_ID,
        "-",
        type_of_interaction,
        "_",
        "predicted_interfaces.txt",
        sep = ""
      )
    )
  
  if (!file.exists(output_filePath)) {
    write.table(
      pdb_atom_inter_df,
      file = output_filePath,
      append = TRUE,
      quote = FALSE,
      sep = "\t",
      row.names = F,
      col.names = !file.exists(output_filePath)
    )
  } else {
    stop("Calculation already done.")
  }
  
}
