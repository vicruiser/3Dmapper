#! /usr/bin/Rscript
###############################################################################
# Load necessary packages
###############################################################################
warn.conflicts = FALSE
options(echo = FALSE, verbose = F,warn = -1) 
#indicate library path
requiredPackages = c('plyr', 'dplyr','data.table','stringr','reshape2','seqinr','bio3d','flock')
suppressMessages(
  for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}
)


###############################################################################
# Input arguments
###############################################################################
args <- commandArgs(trailingOnly = T)
argsLen <- length(args)

ROOT = as.character(commandArgs(TRUE)[1])
pdb_id                = as.character(args[2])
blast_outfiles_dir    = as.character(args[3])
interfaces_dir        = as.character(args[4])
output_dir            = as.character(args[5])


#######################################
# Get files with blast and interfaces
#######################################
blast_outfiles_paths <- c(Sys.glob(file.path(
    blast_outfiles_dir,
   paste(pdb_id, "*filtered*", sep =
           "")
  )))

#chain_ids = str_split_fixed(blast_outfiles_paths, '_',3)[,2]
for( i in blast_outfiles_paths){
  blast_file = i
  
  #predicted_int_file = Sys.glob(paste(
  #  pdb_id, '_chain',chain_id "4udf-assembly1.cif_chain39_protein_predicted_interfaces"))

###############################################################################
# BLAST "PDB - Ensembl" results
###############################################################################
source(file.path(ROOT,"read_Blast.R"))
blast_filtered <- read.Blast(blast_file)

###############################################################################
# Calculated pi results
###############################################################################
source(file.path(ROOT,"new_mapped_interfaces.R"))

chain_id  <- str_match(basename(i), "chain\\s*(.*?)\\s*\\.fasta")[,2]

#new_mpi <-
new.Mapped.Interfaces(pdb_id, chain_id, interfaces_dir, blast_filtered)

# store the results if the are any results associated to the pdb_id
# if (nrow(new_mpi) > 0) {
#   new_mpi = new_mpi[order(new_mpi$spos),]
#   
#   individual_new_mpi  =  unique(new_mpi[, c(
#     "prot.id",
#     "slen",
#     "spos",
#     "resid_sseq",
#     "pdb.id.x",
#     "chain",
#     "qlen",
#     "resno",
#     "qpos",
#     "resid_qseq",
#     "evalue",
#     "pident",
#     "qcov",
#     "length",
#     "interaction",
#     "chain.1",
#     "resno.1",
#     "resid.1",
#     "distance",
#     "b",
#     "b.1",
#     "qstart",
#     "qend",
#     "sstart",
#     "send"
#   )])
#   
#   # individual_new_mpi  <- unique(suppressMessages(
#   #   new_mpi3 %>%
#   #     group_by_at(
#   #       c(
#   #         "prot.id",
#   #         "slen",
#   #         "pdb.id.x",
#   #         "chain",
#   #         "qlen",
#   #         "evalue",
#   #         "pident",
#   #         "qcov",
#   #         "length",
#   #         "interaction",
#   #         "chain.1",
#   #         "qstart",
#   #         "qend",
#   #         "sstart",
#   #         "send"
#   #       )
#   #     ) %>%
#   #     summarise(
#   #       spos = paste(spos, collapse = "/"),
#   #       resid_sseq = paste(resid_sseq, collapse = "/"),
#   #       resno = paste(resno, collapse = "/"),
#   #       resid_qseq = paste(resid_qseq, collapse = "/"),
#   #       qpos = paste(qpos, collapse = "/"),
#   #       resno.1 = paste(resno.1, collapse = "/"),
#   #       resid.1 = paste(ifelse(
#   #         interaction == "protein", aa321(resid.1), resid.1
#   #       ),
#   #       collapse = "/"),
#   #       distance = paste(round(distance, 2), collapse = "/"),
#   #       b = paste(b, collapse = "/"),
#   #       b.1 = paste(b.1,  collapse = "/")
#   #     )
#   # ))#,
#   
#   individual_new_mpi = individual_new_mpi[, c(
#     "prot.id",
#     "slen",
#     "spos",
#     "resid_sseq",
#     "pdb.id.x",
#     "chain",
#     "qlen",
#     "resno",
#     "qpos",
#     "resid_qseq",
#     "evalue",
#     "pident",
#     "qcov",
#     "length",
#     "interaction",
#     "chain.1",
#     "resno.1",
#     "resid.1",
#     "distance",
#     "b",
#     "b.1",
#     "sstart",
#     "send",
#     "qstart",
#     "qend"
#   )]
#   
#   colnames(individual_new_mpi) =
#     c(
#       "Protein_accession",
#       "Protein_length",
#       "Protein_position",
#       "Protein_aa",
#       "PDB_code",
#       "PDB_chain",
#       "PDB_chain_length",
#       "PDB_3D_position",
#       "PDB_seq_position",
#       "PDB_aa",
#       "Evalue",
#       "Pident",
#       "Protein_coverage",
#       "Length_alignment",
#       "Interaction_type",
#       "PDB_interacting_chain",
#       "PDB_interacting_3D_position",
#       "PDB_interacting_aa",
#       "Interface_min_distance",
#       "PDB_B_factor",
#       "PDB_interacting_B_factor",
#       "Protein_alignment_start",
#       "Protein_alignment_end",
#       "PDB_alignment_start",
#       "PDB_alignment_end"
#     )
#   
#   # individual_new_mpi$Interface_interacting_positions = interface_new_mpi$PDB_interacting_position
#   # individual_new_mpi$Chimera_interacting_position = interface_new_mpi$PDB_interacting_3D_position
#   # individual_new_mpi$Chimera_3D_position = interface_new_mpi$PDB_3D_position
#   
#   #individual_new_mpi  <-individual_new_mpi %>%
#   #     group_by_at(
#   #       c("Protein_accession",
#   #         "Protein_length",
#   #       "Protein_position",
#   #       "Protein_aa",
#   # "PDB_code",
#   # "PDB_chain",
#   # "PDB_chain_length",
#   # "PDB_3D_position",
#   # "PDB_seq_position",
#   # "PDB_aa",
#   # "Evalue",
#   # "Pident",
#   # "Protein_coverage",
#   # "Length_alignment",
#   # "Interaction_type",
#   # "PDB_interacting_chain",
#   # "PDB_interacting_3D_position",
#   # "PDB_interacting_aa",
#   # "Interface_min_distance",
#   # "PDB_B_factor",
#   # "PDB_interacting_B_factor",
#   # "Protein_alignment_start",
#   # "Protein_alignment_end",
#   # "PDB_alignment_start",
#   # "PDB_alignment_end")%>%
#   #     summarise(
#   #       Interface_interacting_positions = paste(Interface_interacting_positions, collapse = "/"),
#   #       Chimera_interacting_position = paste(Chimera_interacting_position, collapse = "/"),
#   #       Chimera_3D_position = paste(Chimera_3D_position, collapse = "/")
#   #     )
#   
#   
#   # add protfeatureid
#   individual_new_mpi$Structure_feature_id = paste(
#     individual_new_mpi$PDB_code,
#     individual_new_mpi$Protein_accession,
#     individual_new_mpi$PDB_chain,
#     individual_new_mpi$PDB_interacting_chain,
#     individual_new_mpi$Interaction_type,
#     sep = "_"
#   )
# 
#   
#   individual_new_mpi_list = split(individual_new_mpi, individual_new_mpi$Protein_accession)
#   for (prot in names(individual_new_mpi_list)){
#     fn = file.path(output_dir, paste(prot, '.txt.gz', sep = ""))
#     file.lock = lock(fn)
#     fwrite(
#       individual_new_mpi_list[[prot]],
#       fn,
#       quote = F,
#       row.names = F,
#       col.names = file.info(fn)$size==0,
#       append = T,
#       sep = "\t",
#       compress = "gzip"
#       )
#     unlock(file.lock)
#   }
#   
# } else{
#   #####################################################################################
#   # Store pdb ids associated to interfaces that did not map with any Ensembl prot seq #
#   #####################################################################################
#   
#   write.table(
#     paste(
#       "The predicted interfaces of",
#       pdb_id ,
#       "did not map with any protein sequence.",
#       sep = " "
#     ),
#     file.path(script.dir, "unmapped_interfaces_pdbids.txt"),
#     quote = F,
#     row.names = F,
#     col.names = F,
#     sep = "\t",
#     append = T
#   )
#   rm(list = ls())
#   suppressMessages(gc(reset = TRUE))
#   stop(
#     paste(
#       "The predicted interfaces of",
#       pdb_id ,
#       "did not map with any protein sequence.",
#       sep = " "
#     )
#   )
# }

}
