#### Prepare the pdb ids greasy

# greasy_selected_pdb.ids <- paste("Rscript /gpfs/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/execute_PDB_calc_dist_main.R",
#                                  IDs_formatted_name, 
#                                  "/gpfs/scratch/bsc51/bsc51927/pdb_db",
#                                  "/gpfs/home/bsc51/bsc51927/PDB_map_interfaces/output/out_calc_dist/ 5")
# 
# write.table(greasy_selected_pdb.ids, "/home/vruizser/Escritorio/pdb_db/greasy_calc_pdb_dist.txt", quote=F, row.names = F, col.names = F)
program                  = as.character(commandArgs(TRUE)[1])  
path_to_excutable_script = as.character(commandArgs(TRUE)[2])
tasks_list               = read.csv((commandArgs(TRUE)[3]))
input_arguments          = as.character(commandArgs(TRUE)[4])
out_path                 = as.character(commandArgs(TRUE)[5])
ouput_filename           = as.character(commandArgs(TRUE)[6])


greasy_tasks <- paste(program,
                      path_to_excutable_script,
                      tasks_list,
                      input_arguments,
                      sep = " ")

write.table(greasy_taks,
            paste(out_path,
                  ouput_filename,
                  sep =""),
            quote=F,
            row.names = F,
            col.names = F)

  
 
  
