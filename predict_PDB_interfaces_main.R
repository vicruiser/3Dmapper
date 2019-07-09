################################################################
##  Execute function in a BASH script to parallelize the task. 

#indicate library path
.libPaths(c(.libPaths(),"/gpfs/home/bsc51/bsc51927/R/x86_64-pc-linux-gnu-library")) 
library(plyr)
library(dplyr)
# Load necessary functions
source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/PDB_pairwise_interatom_dist.R")  
source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/readPDB.R") 

#########################
# Pass input arguments  #
#########################

PDB_FILENAME = as.character(commandArgs(TRUE)[1])  
INPUT_DIR = as.character(commandArgs(TRUE)[2])
OUTPUT_DIR = as.character(commandArgs(TRUE)[3])
DIST_THRESHOLD = as.numeric(commandArgs(TRUE)[4]) 

######################################################################
# Compute distance calculation for each type of possible interaction #
######################################################################
int <- c("protein", "ligand", "nucleic")
for (j in 1:length(int)) {
  tryCatch({
    pdb_file <- readPDB(PDB_FILENAME,
                       download = "NO",
                       input_dir = INPUT_DIR
                       )
    
    PDB_iter_atom_distances(
      PDB_FILENAME,
      pdb_file,
      type_of_interaction = int[j],
      dist_threshold = DIST_THRESHOLD,
      output_dir = OUTPUT_DIR
    )
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

file.remove(file.path(INPUT_DIR,sub("\\.gz", "", PDB_FILENAME)))
