###############################################################################
# Load necessary packages
###############################################################################
#indicate library path
.libPaths(c(.libPaths(),"/gpfs/home/bsc51/bsc51927/R/x86_64-pc-linux-gnu-library")) 
library(plyr, quietly = T)
library(dplyr, quietly = T)
library(data.table, quietly = T)
library(stringr, quietly = T)
library(reshape2, quietly = T)
library(seqinr)
library(bio3d)

###############################################################################
# Input arguments
###############################################################################
pdb_id                = as.character(commandArgs(TRUE)[1])
blast_outfiles_dir    = as.character(commandArgs(TRUE)[2])
interfaces_dir        = as.character(commandArgs(TRUE)[3])
output_dir            = as.character(commandArgs(TRUE)[4])

## NOTE:  (m)pi makes reference to (mapped) protein interface
# pdb_id = "1pzu.pdb2"
# output_dir = interfaces_dir = blast_outfiles_dir = "/home/vruizser/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/TEST"
# output_dir = interfaces_dir = blast_outfiles_dir = "/home/victoria/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/TEST"

###############################################################################
# BLAST "PDB - Ensembl" results
###############################################################################

source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/read_Blast.R")
#source("~/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/read_Blast.R")
blast_output <- read.Blast(pdb_id, blast_outfiles_dir)
gc()
###############################################################################
# Calculated pi results
###############################################################################
source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/new_mapped_interfaces.R")
#source("~/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/new_mapped_interfaces.R")
new_mpi <-
  new.Mapped.Interfaces(pdb_id , interfaces_dir, blast_output)
gc()

# store the results if the are any results associated to the pdb_id
if (nrow(new_mpi) > 0) {
  #filter out repeated rows
  new_mpi_compare <-
    unique(new_mpi[, c(
      "chain",
      "pdb.id.x",
      "ensembl.prot.id",
      "sstart",
      "send",
      "chain.1",
      "interaction",
      "mapped.real.pos"
    )])
  new_mpi_compare$pdb.id.x <-
    as.character(levels(new_mpi_compare$pdb.id.x))[new_mpi_compare$pdb.id.x]
  new_mpi_compare$interaction <-
    as.character(levels(new_mpi_compare$interaction))[new_mpi_compare$interaction]
  
  if (nrow(new_mpi_compare) > 0) {
    new_mpi_compare <- data.frame(new_mpi_compare, id = "new")
  } else{
    new_mpi_compare <- data.frame(new_mpi_compare, id = character())
  }
  
  ############################################################################
  # Reformat Eduard'S results file to a R friendly format.
  ############################################################################
  source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/ref_mapped_interfaces.R")
  #source("/home/vruizser/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/ref_mapped_interfaces.R")
  ref_mpi <- ref.Mapped.Interfaces(pdb_id)
  gc()
  
  ############################################################################
  # Comparison reference with updated version
  ############################################################################
  ref_new_mpi <- rbindlist(list(unique(ref_mpi),
                                new_mpi_compare[!colnames(new_mpi_compare) == "interaction"]),
                           use.names = TRUE)
  mpi_compared <-
    setDT(ref_new_mpi)[, .N, by = c(names(ref_new_mpi)[-8])]
  names(mpi_compared)[8] <- "freq"
  
  ############################################################################
  # COMPARISON between reference results and updated
  ############################################################################
  source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/ref_new_mpi_comparison.R")
  #source("/home/vruizser/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/ref_new_mpi_comparison.R")
  comparison_ref_new_table <-
    Ref.New.Mpi.Comparison(ref_mpi, new_mpi, mpi_compared, pdb_id)
  gc()
  # Save file
  write.table(
    comparison_ref_new_table,
    paste(output_dir, "mpi_v89_5_vs_v94.csv", sep = ""),
    quote = FALSE,
    append = T,
    col.names = !file.exists(paste(
      output_dir, "mpi_v89_5_vs_v94.csv", sep = ""
    )),
    row.names = F
  )
  ############################################################################
  # Associate PFAM domains with each Ensembl ID if possible
  ############################################################################
  source("/home/bsc51/bsc51927/PDB_map_interfaces/r_scripts/mapped_interfaces_pfam.R")
  #source("~/Dropbox/Interfaces_mapping/task4_Map_prot_interfaces/mapped_interfaces_pfam.R")
  new_mpi_biomart <-  Mapped.Interfaces.Biomart(new_mpi)
  gc()
  # Save file
  setnames(
    new_mpi_biomart,
    c(
      "chain",
      "pdb.id.x",
      "slen",
      "sstart",
      "send",
      "length",
      "resno",
      "resno.1",
      "chain.1",
      "resid.1",
      "Gene stable ID",
      "Transcript stable ID",
      "Gene start (bp)",
      "Gene end (bp)",
      "Chromosome/scaffold name"
    ),
    c(
      "temp.chain",
      "pdb.id",
      "temp.length",
      "temp.start",
      "temp.end",
      "length.ali",
      "temp.resno",
      "int.resno",
      "int.chain",
      "int.resid",
      "ensembl.gene.id",
      "ensembl.transcript.id",
      "gene.start.bp",
      "gene.end.bp",
      "Chr"
    )
  )
  
  drops <-
    c(
      "type.1",
      "pdb.id.y",
      "nident",
      "evalue",
      "type",
      "qlen",
      "qstart",
      "qend",
      "temp.resno",
      "int.resno",
      "mapped.pos",
      "is.gap",
      "n_gap",
      "sseq",
      "qseq",
      "sgaps",
      "qgaps",
      "resid"
    )
  
  new_mpi_biomart <-
    unique(new_mpi_biomart[,!(names(new_mpi_biomart) %in% drops)])
  
  new_mpi_biomart3  <- ddply(new_mpi_biomart,.(pdb.id,
                          ensembl.prot.id,
                          ensembl.gene.id,
                          ensembl.transcript.id,
                          gene.start.bp,
                          gene.end.bp,
                          Chr,
                          temp.chain,
                          int.chain,
                          temp.length,
                          temp.start,
                          temp.end,
                          length.ali,
                          pident,
                          interaction,
                          pfam.dom,
                          dom.start,
                          dom.end),
        summarise,
        mapped.real.pos = paste(mapped.real.pos, collapse = "-"),
        resid_qseq = paste(resid_qseq, collapse = "-"),
        resid_sseq = paste(resid_sseq, collapse = "-"),
        qpos = paste(qpos, collapse = "-"),
        spos = paste(spos, collapse = "-"),
        q_ali_pos = paste(q_ali_pos, collapse = "-"),
        pdb.pos = paste(real.pos, collapse = "-"))
  
  # write individual files per PDB id to keep all the distance calculations
  individual_new_mpi_biomart <- unique(new_mpi_biomart3[,c("pdb.id", "ensembl.prot.id",
                          "temp.chain",
                          "int.chain",
                          "temp.length",
                          "temp.start",
                          "temp.end",
                          "length.ali",
                          "pident",
                          "interaction",
                          "resid_qseq",
                          "resid_sseq",
                          "qpos",
                          "spos",
                          "q_ali_pos",
                          "mapped.real.pos",
                          "pdb.pos"
                          
                          )])
  write.table(
    individual_new_mpi_biomart,
    paste(output_dir, "mpi_v94_", pdb_id, ".csv", sep = ""),
    quote = F,
    row.names = F
  )
  
  # generate one file containing only PFAM, gene position, ensembl id and pdb ids
  new_mpi_biomart_all <- unique(new_mpi_biomart3[c(
    "pdb.id",
    "ensembl.prot.id",
    "ensembl.gene.id",
    "ensembl.transcript.id",
    "gene.start.bp",
    "gene.end.bp",
    "Chr",
    "pfam.dom",
    "dom.start",
    "dom.end",
    "pident",
    "mapped.real.pos",
    "resid"
  )])
  
  write.table(
    new_mpi_biomart_all,
    paste(output_dir, "mpi_v94_genepos_all.csv", sep = ""),
    quote = F,
    row.names = F,
    append = T,
    col.names = !file.exists(paste(
      output_dir, "mpi_v94_genepos_all.csv", sep = ""
    ))
  )
  
  
  rm(list = ls())
  gc(reset = TRUE)
  
  #####################################################################################
  # Store pdb ids associated to interfaces that did not map with any Ensembl prot seq #
  #####################################################################################
} else{
  pdb_id = as.character(commandArgs(TRUE)[1])
  write.table(
    paste(
      "The predicted interfaces of",
      pdb_id ,
      "did not map with any Ensembl protein sequence.",
      sep = " "
    ),
    "/gpfs/home/bsc51/bsc51927/PDB_map_interfaces/debug/unmapped_interfaces_pdbids.txt",
    quote = F,
    row.names = F,
    col.names = F,
    append = T
  )
  rm(list = ls())
  gc(reset = TRUE)
  stop(
    paste(
      "The predicted interfaces of",
      pdb_id ,
      "did not map with any Ensembl protein sequence.",
      sep = " "
    )
  )
}
