#! /usr/bin/Rscript
##### Load necessary functions

options(echo = FALSE, verbose = F,warn = -1) 

#' Title: CALCULATION INTER-ATOMIC DISTANCES
#'
#' @param type_of_interaction
#' @param dist_threshold
#' @param pdb_filename
#' @param pdb_file
#' @param atom_select
#' @param output_dir
#'
#' @return
#' @export
#'
#' @examples
#'
PDB_iter_atom_distances <- function(pdb_filename,
                                    pdb_file,
                                    type_of_interaction,
                                    atom_select, 
                                    # backbone heavy atoms are only considered as default
                                    dist_threshold,
                                    output_dir, 
                                    ROOT, 
                                    biolip) {
  biounit_filename = basename(pdb_filename)
  if (biolip == "True" & type_of_interaction == "ligand"){
    bp = fread(file.path(ROOT,'biolip_list.txt'), header = F)
    pdb_file$atom = subset(pdb_file$atom, !resid %in% bp$V1 )
    if (nrow(pdb_file$atom)<1) {

      stop("This PDB file does not interact with any ligand.")
    } 
  }
  # Get all chains to make pairs type_of_interactions
  ChainsList <-
    PDB_pairwise_interaction(pdb_file, type_of_interaction, atom_select)
  TemplateChains <-   ChainsList[[1]]
  InteractionChains <- ChainsList[[2]]
  names(InteractionChains)[names(InteractionChains) == "real.pos"] = "real.pos.1"

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
  CompareTemplateChains <- TemplateChains_list[PairwiseComb[1,]]
  CompareInteractionChains <-
    InteractionChains_list[PairwiseComb[2,]]
  
  rm(TemplateChains_list, InteractionChains_list)
  gc()
  
  #CompareTemplateChains.Vector <- lapply(CompareTemplateChains,
  #                                       function(x) {
  #                                         x <- x[, c("x", "y", "z")]
  #                                         x <- as.vector(t(x))
  #                                         return(x)
  #                                       })
  #CompareInteractionChains.Vector <-
  #  lapply(CompareInteractionChains,
  #         function(x) {
  #           x <- x[, c("x", "y", "z")]
  #           x <- as.vector(t(x))
  #           return(x)
  #         })
  distList <- list()
  struct = list()
  i = 1
  while (i <= length(PairwiseComb)/2) {

    #distList[[paste("iter", i, sep = "")]] <-
      pdb_atom_inter_df <-
        calc_PDB_dist(
          CompareTemplateChains[[i]],
          CompareInteractionChains[[i]],
          as.vector(t(CompareTemplateChains[[i]][, c("x", "y", "z")])),
          as.vector(t(CompareInteractionChains[[i]][, c("x", "y", "z")])),
          dist_threshold ,
          type_of_interaction = type_of_interaction
        )
    # Filter by selected distance threshold
    

    
   # pdb_atom_inter_df <-
  #    ldply(distList, function(x)
  #      data.frame(x))
    if (nrow(pdb_atom_inter_df) > 0) {
    #  print(head(pdb_atom_inter_df))
      #add identifier of the PDB structure
      PDB_ID <- sub("\\.gz+", "", biounit_filename)
      pdb_atom_inter_df$pdb.id <-  PDB_ID

      # Eliminate dummy columns
      if (type_of_interaction == "protein") {
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
          "b",
          "b.1",
          "interaction",
          "pdb.id",
          "real.pos",
          "real.pos.1"
        )]
        pdb_atom_inter_df$resid.1 = aa321(pdb_atom_inter_df$resid.1)
        struct[[i]] = unique(pdb_atom_inter_df[,c ("chain", "resno")])
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
          "b",
          "b.1",
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
            "_chain",
            unique(pdb_atom_inter_df$chain),
            "_",
            type_of_interaction,
            "_",
            "predicted_interfaces.txt",
            sep = ""
          )
        )
      #print(head(pdb_atom_inter_df))
      
      pdb_atom_inter_df$resid = aa321(pdb_atom_inter_df$resid)
      write.table(
        pdb_atom_inter_df,
        file = output_filePath,
        append = TRUE,
        quote = FALSE,
        sep = "\t",
        row.names = F,
        col.names = !file.exists(output_filePath)
      )
      rm(pdb_atom_inter_df)
      distList <- list()
      gc()
      
    }
    i = i + 1
  }

  if (type_of_interaction == "protein"){
    return(rbindlist(struct))  
  }
  
}
