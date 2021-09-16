#! /usr/bin/Rscript
################################################################
warn.conflicts = FALSE
library(veriNA3d)
requiredPackages = c('tidyr', 'stringr', 'bio3d',  'plyr','dplyr') #'parallel',
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

for (j in 1:length(INT)) {
  skip_to_next <- FALSE
    tryCatch({
    PDB_iter_atom_distances(
      PDB_FILENAME,
      pdb_file,
      atom_select = ATOMS_INTERACTION,
      type_of_interaction = INT[j],
      dist_threshold = DIST_THRESHOLD,
      output_dir = OUTPUT_DIR,
      ROOT,
      biolip = BIOLIP
    )},error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next } 
}
  



#file.remove(sub("\\.gz", "", PDB_FILENAME))
