library(bio3d)
library(veriNA3d)
#pass input arguments (compressed pdb file)
INPUT_PDB = as.character(commandArgs(TRUE)[1])  
INPUT_DIR = as.character(commandArgs(TRUE)[2]) 

list_pdb_files = list.files("/home/vruizser/PhD/2018-2019/pdb_db/complete_PDB_9Jan2019")

# read PDB file function
source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/readPDB.R") 
source("/home/vruizser/Dropbox/Interfaces_mapping/task3_Predict_prot_interfaces/readPDB.R")

INPUT_DIR="/home/vruizser/PhD/2018-2019/pdb_db/complete_PDB_9Jan2019"
for (i in 57863:length(list_pdb_files)){
  print(i)
  INPUT_PDB = list_pdb_files[i]
  pdb_file = readPDB(INPUT_PDB, download = "NO", input_dir = INPUT_DIR)
  inserts <- any(!is.na(pdb_file$atom$insert))
  
  if(inserts == TRUE){
    write.table( data.frame(pdb.id = paste(INPUT_PDB)),"pdb_files_with_inserts.csv", append = T, col.names = !file.exists("pdb_files_with_inserts.csv"))
  } else {
    write.table( data.frame(pdb.id =paste(INPUT_PDB)), "pdb_files_without_inserts.csv", append = T, col.names = !file.exists("pdb_files_without_inserts.csv")) 
  }
}
