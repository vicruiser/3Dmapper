#! /usr/bin/Rscript
###############################################################################
# BLAST "PDB - Ensembl" results
###############################################################################
#' Title
#'
#' @param pdb_id
#' @param blast_outfiles_dir
#'
#' @return
#' @export
#'
#' @examples
#' 
options(echo = FALSE, verbose = F,warn = -1) 
read.Blast <- function( blast_file) {
  # read all blast output files of all chains corresponding to the selected pdb_id
  # blast_outfiles_paths <- c(Sys.glob(file.path(
  #   blast_outfiles_dir,
  #   paste(pdb_id, "*filtered*", sep =
  #           "")
  # )))
  # blast_output_list <-
  #   lapply(blast_outfiles_paths, function(x)
  #     fread(x, header = F, sep =" "))
  # blast_output <- as.data.frame(rbindlist(blast_output_list))
  
  blast_output = fread(blast_file, header = F, sep =" ")
  #print(blast_output)
  # change column names conveniently to merge data frames afterwards by colname
  colnames(blast_output) <-
    c(
      "pdb.id",
      "chain",
      "qlen",
      "prot.id" ,
      "slen",
      "qstart",
      "qend",
      "sstart",
      "send" ,
      "evalue" ,
      "length" ,
      "pident" ,
      "nident",
      "qseq",
      "sseq",
      "gaps",
      "qcov"
    )
  
  # remove word "chain" from chain column for comparative purposes
  blast_output$chain <- sub("chain", "", blast_output$chain)
  blast_output$prot.id <-
    sub("\\..*", "", blast_output$prot.id)
  
  blast_output$qseq = as.character(blast_output$qseq)
  blast_output = subset(as.data.frame(blast_output), pdb.id!="pdbid")
  blast_output$qstart = as.numeric(blast_output$qstart)
  blast_output$sstart = as.numeric(blast_output$sstart)
  return(blast_output)
}
