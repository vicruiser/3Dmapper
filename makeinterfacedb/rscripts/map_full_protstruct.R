###############################################################################
# Load necessary packages
###############################################################################
#indicate library path
.libPaths(c(.libPaths(),"/gpfs/home/bsc08/bsc08927/R/x86_64-pc-linux-gnu-library")) 
suppressMessages(library(plyr, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(stringr, quietly = T))
suppressMessages(library(reshape2, quietly = T))
suppressMessages(library(seqinr))
suppressMessages(library(bio3d))
suppressMessages(library(here,lib.loc ="/gpfs/home/bsc08/bsc08927/R/x86_64-pc-linux-gnu-library" ))

script.dir <- here()
###############################################################################
# Input arguments
###############################################################################
args <- commandArgs(trailingOnly=T)
argsLen <- length(args)

pdb_id                = as.character(args[1])
blast_outfiles_dir    = as.character(args[2])
interfaces_dir        = as.character(args[3])
#dist_threshold        = as.numeric(args[4])
output_dir            = as.character(args[4])
protdb_version        = as.character(args[5])

if (argsLen == 6){
  biomart_filepath      = as.character(args[6])
} else {
  biomart_filepath      = NULL
}

## NOTE:  (m)pi makes reference to (mapped) protein interface
# pdb_id = "1a2x.pdb1"
# output_dir = interfaces_dir = blast_outfiles_dir = "/home/vruizser/PhD/2018-2019/git/Interfacer/task4_Map_prot_interfaces/test"
# output_dir = interfaces_dir = blast_outfiles_dir = "/home/victoria/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/TEST"

###############################################################################
# BLAST "PDB - Ensembl" results
###############################################################################

source(file.path(script.dir,"read_Blast.R"))
#source("/home/vruizser/PhD/2018-2019/git/Interfacer/task4_Map_prot_interfaces/read_Blast.R")
blast_output <- read.Blast(pdb_id, blast_outfiles_dir)
suppressMessages(gc())
###############################################################################
# Calculated pi results
###############################################################################
source(file.path(script.dir,"new_mapped_full_protstruct.R"))
#source("/home/vruizser/PhD/2018-2019/git/Interfacer/new_mapped_interfaces.R")
new_mpi <-
  new.Mapped.Interfaces(pdb_id , interfaces_dir, blast_output)

#new_mpi = subset(new_mpi, distance <= dist_threshold)
suppressMessages(gc())



# store the results if the are any results associated to the pdb_id
if (nrow(new_mpi) > 0) {
  #filter out repeated rows
  # new_mpi_compare <-
  #   unique(new_mpi[, c(
  #     "chain",
  #     "pdb.id.x",
  #     "prot.id",
  #     "sstart",
  #     "send",
  #     "chain.1",
  #     "interaction",
  #     "mapped.real.pos"
  #   )])
  # new_mpi_compare$pdb.id.x <-
  #   as.character(levels(new_mpi_compare$pdb.id.x))[new_mpi_compare$pdb.id.x]
  # new_mpi_compare$interaction <-
  #   as.character(levels(new_mpi_compare$interaction))[new_mpi_compare$interaction]
  # 
  # if (nrow(new_mpi_compare) > 0) {
  #   new_mpi_compare <- data.frame(new_mpi_compare, id = "new")
  # } else{
  #   new_mpi_compare <- data.frame(new_mpi_compare, id = character())
  # }
  # 
  ############################################################################
  # Reformat Eduard'S results file to a R friendly format.
  ############################################################################
#  source(file.path(script.dir,"ref_mapped_interfaces.R"))
 # source("/home/vruizser/PhD/2018-2019/git/Interfacer/task4_Map_prot_interfaces/ref_mapped_interfaces.R")
#  ref_mpi <- ref.Mapped.Interfaces(pdb_id)
#  suppressMessages(gc())
  
  ############################################################################
  # Comparison reference with updated version
  ############################################################################
 # ref_new_mpi <- rbindlist(list(unique(ref_mpi),
#                                new_mpi_compare[!colnames(new_mpi_compare) == "interaction"]),
#                           use.names = TRUE)
#  mpi_compared <-
#    setDT(ref_new_mpi)[, .N, by = c(names(ref_new_mpi)[-8])]
#  names(mpi_compared)[8] <- "freq"
  
  ############################################################################
  # COMPARISON between reference results and updated
  ############################################################################
#  source(file.path(script.dir,"ref_new_mpi_comparison.R"))
  #source("/home/vruizser/PhD/2018-2019/git/Interfacer//task4_Map_prot_interfaces/ref_new_mpi_comparison.R")
  # comparison_ref_new_table <-
  #   Ref.New.Mpi.Comparison(ref_mpi, new_mpi, mpi_compared, pdb_id)
  # suppressMessages(gc())
  # # Save file
  # write.table(
  #   comparison_ref_new_table,
  #   file.path(output_dir,paste( "mpi_v89_5_vs_",protdb_version,".txt", sep = "")),
  #   quote = FALSE,
  #   append = T,
  #   col.names = !file.exists(file.path(output_dir,paste("mpi_v89_5_vs_",protdb_version,".txt", sep = ""
  #   ))), sep ="\t",
  #   row.names = F
  # )

  ############################################################################
  # Associate PFAM domains with each Ensembl ID if possible
  ############################################################################
  source(file.path(script.dir,"mapped_interfaces_pfam.R"))

  #source("/home/vruizser/PhD/2018-2019/git/Interfacer/task4_Map_prot_interfaces/mapped_interfaces_pfam.R")
  if (!is.null(biomart_filepath)){
  new_mpi_biomart <-  Mapped.Interfaces.Biomart(new_mpi, biomart_filepath)
  }else{
    new_mpi_biomart = new_mpi
  }
  

 # new_mpi_biomart = new_mpi_biomart[order(new_mpi_biomart$mapped.real.pos),]
  
  new_mpi_biomart3=  unique(new_mpi_biomart[,c("prot.id",
                                               "slen",
                                               "spos",
                                               "resid_sseq",
                                               "pdb.id.x",
                                               "chain",
                                               "qlen",
                                               "resno", 
                                               "qpos",
                                               "resid_qseq",
                                               "evalue",
                                               "pident",
                                               "length",
                                              # "interaction",
                                             #  "chain.1",
                                             #  "resno.1",
                                             #  "resid.1",
                                             #  "distance",
                                             
                                               "sstart",
                                               "send",
                                             "qstart",
                                             "qend"
  )])
  
  individual_new_mpi_biomart  <- unique(suppressMessages(new_mpi_biomart3%>%
                                                      group_by_at( c("prot.id",
                                                                     "slen",
                                                                     "pdb.id.x",
                                                                     "chain",
                                                                     "qlen",
                                                                     "evalue",
                                                                     "pident",
                                                                     "length",
                                                                  #   "interaction",
                                                                 #    "chain.1", 
                                                                 "qstart",
                                                                 "qend",
                                                                     "sstart", 
                                                                     "send") ) %>%
                                                           #group_by(pdb.id, prot.id,temp.chain,int.chain, interaction)%>%
                                                           summarise(
                                                             spos = paste(spos, collapse = "-"),
                                                             resid_sseq = paste(resid_sseq, collapse = "-"),
                                                             resno = paste(resno, collapse = "-"),
                                                             resid_qseq = paste(resid_qseq, collapse = "-"),
                                                             qpos = paste(qpos, collapse = "-")#,
                                                         #    resno.1 = paste(resno.1, collapse = "-"),
                                                         #    resid.1 = paste(ifelse(interaction == "protein", aa321(resid.1), resid.1),
                                                         #                    collapse = "-"),
                                                          #   distance = paste(round(distance, 2), collapse = "-")))
                                        )))#,

  individual_new_mpi_biomart = individual_new_mpi_biomart[,c("prot.id",
                                                             "slen",
                                                             "spos",
                                                             "resid_sseq",
                                                             "pdb.id.x",
                                                             "chain",
                                                             "qlen",
                                                             "resno", 
                                                             "qpos",
                                                             "resid_qseq",
                                                             "evalue",
                                                             "pident",
                                                             "length",
                                                             # "interaction",
                                                             #  "chain.1",
                                                             #  "resno.1",
                                                             #  "resid.1",
                                                             #  "distance",
                                                             
                                                             "sstart",
                                                             "send",
                                                             "qstart",
                                                             "qend"
  )]
 
  colnames(individual_new_mpi_biomart) = 
    c("Protein_accession",
      "Protein_length",
      "Protein_position",
      "Protein_aa",
      "PDB_code",
      "PDB_chain",
      "PDB_chain_length",
      "PDB_3D_position",
      "PDB_seq_position",
      "PDB_aa", 
      "Evalue",
      "Pident", 
      "Length_alignment",
      "Protein_alignment_start",
      "Protein_alignment_end",
      "PDB_alignment_start",
      "PDB_alignment_end"
      )#,
      #"Interaction_type",
      #"PDB_interacting_chain",
      #"PDB_interacting_position", 
      #"PDB_interacting_aa", 
      #"Interface_min_distance")
  
  
  # add protfeatureid
  individual_new_mpi_biomart$Structure_feature_id = paste(individual_new_mpi_biomart$PDB_code, 
                                                          individual_new_mpi_biomart$Protein_accession, 
                                                          individual_new_mpi_biomart$PDB_chain,
                                                          #individual_new_mpi_biomart$PDB_interacting_chain,
                                                          #individual_new_mpi_biomart$Interaction_type,
                                                          sep="_")
  
  write.table(
    individual_new_mpi_biomart,
    file.path(output_dir,paste( "mpi_",protdb_version,"_", pdb_id,"_.txt", sep = "")),
    quote = F,
    row.names = F,
    sep = "\t"
  )
  
  
  
  
  
  # setnames(
  #   new_mpi_biomart,
  #   c(
  #     "chain",
  #     "pdb.id.x",
  #     "slen",
  #     "sstart",
  #     "send",
  #     "length",
  #     "resno",
  #     "resno.1",
  #     "chain.1",
  #     "resid.1",
  #     "Gene stable ID",
  #     "Transcript stable ID",
  #     "Gene start (bp)",
  #     "Gene end (bp)",
  #     "Chromosome/scaffold name"
  #   ),
  #   c(
  #     "temp.chain",
  #     "pdb.id",
  #     "temp.length",
  #     "temp.start",
  #     "temp.end",
  #     "length.ali",
  #     "temp.resno",
  #     "int.resno",
  #     "int.chain",
  #     "int.resid",
  #     "gene.id",
  #     "transcript.id",
  #     "gene.start.bp",
  #     "gene.end.bp",
  #     "Chr"
  #   )
  # )
  # 
  # new_mpi_biomart3  <-suppressMessages( new_mpi_biomart%>%
  #                                         group_by_at( c( "pdb.id",
  #                                                         "prot.id",
  #                                                         "gene.id",
  #                                                         "transcript.id",
  #                                                         "gene.start.bp",
  #                                                         "gene.end.bp",
  #                                                         "Chr",
  #                                                         "pfam.dom",
  #                                                         "dom.start",
  #                                                         "dom.end",
  #                                                         "pident" )) %>%
  #                                         #group_by(pdb.id, prot.id,temp.chain,int.chain, interaction)%>%
  #                                         summarise(
  #                                           mapped.real.pos = paste(unique(mapped.real.pos), collapse = "-"))) #,
  # #max.dist = paste(round(max.dist, 2), collapse = "-"),
  # #mean.dist = paste(round(mean.dist,2), collapse = "-"))
  # 
  # new_mpi_biomart_all <- unique(new_mpi_biomart3[c(
  #   "pdb.id",
  #   "prot.id",
  #   "gene.id",
  #   "transcript.id",
  #   "gene.start.bp",
  #   "gene.end.bp",
  #   "Chr",
  #   "pfam.dom",
  #   "dom.start",
  #   "dom.end",
  #   "pident",
  #   "mapped.real.pos"
  # )])
  # 
  # 
  # 
  # write.table(
  #   new_mpi_biomart_all,
  #   file.path(output_dir,paste( "mpi_",protdb_version,"_genepos_all.txt", sep = "")),
  #   quote = F,
  #   row.names = F,
  #   append = T,
  #   sep = "\t",
  #   col.names = !file.exists(file.path(
  #     output_dir,paste( "mpi_",protdb_version,"_genepos_all.txt", sep = ""
  #     )))
  # )
  
  
  # 
  #####################################################################################
  # Store pdb ids associated to interfaces that did not map with any Ensembl prot seq #
  #####################################################################################
 } else{
  #pdb_id  = as.character(commandArgs(TRUE)[1])
  write.table(
    paste(
      "The predicted interfaces of",
      pdb_id ,
      "did not map with any Ensembl protein sequence.",
      sep = " "
    ),
    file.path(script.dir,"log/unmapped_interfaces_pdbids.txt"),
    quote = F,
    row.names = F,
    col.names = F,
    sep = "\t",
    append = T
  )
  rm(list = ls())
  suppressMessages(gc(reset = TRUE))
  stop(
    paste(
      "The predicted interfaces of",
      pdb_id ,
      "did not map with any Ensembl protein sequence.",
      sep = " "
    )
  )
}
