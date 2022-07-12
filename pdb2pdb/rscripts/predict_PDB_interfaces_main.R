#! /usr/bin/Rscript
################################################################
warn.conflicts = FALSE
library(veriNA3d)
requiredPackages = c('tidyr', 'stringr', 'bio3d',  'plyr','dplyr', 'data.table') #'parallel',
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  suppressMessages(library(p, character.only = TRUE))
}


#########################
# Input arguments       #
#########################
ROOT = as.character(commandArgs(TRUE)[1])
PDB_FILENAME = as.character(commandArgs(TRUE)[2])
OUTPUT_DIR = as.character(commandArgs(TRUE)[3])
DIST_THRESHOLD = as.numeric(commandArgs(TRUE)[4])
ATOMS_INTERACTION = as.character(commandArgs(TRUE)[5])
INT = as.character(commandArgs(TRUE)[6])
BIOLIP = commandArgs(TRUE)[7]

# Load necessary functions
source(file.path(ROOT, "PDB_pairwise_interatom_dist.R"))
source(file.path(ROOT, "readPDB.R"))
source(file.path(ROOT,"PDB_pairwise_interaction.R"))
source(file.path(ROOT,"calc_interatom_dist.R"))

######################################################################
# Compute distance calculation for each type of possible interaction #
######################################################################

if(INT == 'all'){
  INT = c("protein", "nucleic", "ligand")
} else {
  INT = str_split(INT, ' ')[[1]]
}

tryCatch({
  pdb_file <- readPDB(PDB_FILENAME,
                      download = "NO")
}, error = function(e) {
  cat("ERROR :", conditionMessage(e), "\n")
})

########################
## TO MAP TO STRUCTURE
#########################
# Select information of all the atoms of the pdb file
all_chains <- pdb_file$atom
# Select position of protein chains in PDB file
pdb_protChains <-
  combine.select(
    atom.select(pdb_file, ATOMS_INTERACTION, verbose = F),
    atom.select(pdb_file, "protein", verbose = F)
  )
# Select protein chains
protein_chains <- all_chains[pdb_protChains$atom ,]
if (ATOMS_INTERACTION == "calpha") {
  protein_chains = subset(protein_chains , elety == "CA")
}
if (ATOMS_INTERACTION == "cbeta") {
  df = rbind(
    subset(protein_chains, resid == "GLY" & elety == "CA"),
    subset(protein_chains , elety == "CB")
  )
  protein_chains = df[order(df$eleno),]
}
protein_chains <-
  split(protein_chains, f = protein_chains$chain)
protein_chains <-
  lapply(protein_chains, function(x) {
    respos <-
      unique(x[, c("resno", "insert")])  # take into account insertions of residues
    df <-
      data.frame(respos, real.pos = 1:nrow(respos))
    left_join(x, df, by = c("resno", "insert"))
  })
protein_chains <- ldply(protein_chains, data.frame)

PDB_ID <- basename(sub("\\.gz+", "", PDB_FILENAME))


for (j in 1:length(INT)) {
  
  skip_to_next <- FALSE
  tryCatch({
    interfaces = PDB_iter_atom_distances(
      PDB_FILENAME,
      pdb_file,
      atom_select = ATOMS_INTERACTION,
      type_of_interaction = INT[j],
      dist_threshold = DIST_THRESHOLD,
      output_dir = OUTPUT_DIR,
      ROOT,
      biolip = BIOLIP)
    if (INT[j] == "protein"){
    struct <<- anti_join(unique(protein_chains[, c(
      "chain",
      "resid",
      "resno",
      "real.pos",
      "b")]), interfaces, by=c("chain", "resno"))
    } else { skip_to_next <<- TRUE}
    # test if works
   },error = function(e) { 
      
    if (INT[j] == "protein"){
      struct <<- unique(protein_chains[, c(
       "chain",
       "resid",
       "resno",
       "real.pos",
       "b")])
    } else{ skip_to_next <<- TRUE }
  })
  
  if(skip_to_next) { next }  
  
   struct = unique(
     data.frame(type = "ATOM", 
                eleno = "NA", 
                elety = "NA",
                struct[, c(
                  "chain",
                  "resid",
                  "resno")], 
                          type.1 = 'NA',
                          eleno.1 = 'NA',
                          elety.1 = 'NA',
                          chain.1 = 'NA',
                          resid.1 = "NA",
                          resno.1 = "NA",
                          distance = "NA",
                          b = struct[,"b"],
                          b.1 = 'NA',
                          interaction = "NA",
                          pdb.id = PDB_ID,
                          real.pos = struct[,"real.pos"],
                          real.pos.1 ="NA"))
    output_filePath <-
      file.path(
        OUTPUT_DIR,
        paste(
          PDB_ID,
          "_",
          INT[j],
          "_",
          "predicted_interfaces.txt",
          sep = ""
        )
      ) 
    write.table(
      struct,
      file = output_filePath,
      append = TRUE,
      quote = FALSE,
      sep = "\t",
      row.names = F,
      col.names = !file.exists(output_filePath)
    )  
    

}
  



#file.remove(sub("\\.gz", "", PDB_FILENAME))
