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
read.Blast <- function(pdb_id, blast_outfiles_dir) {
  # read all blast output files of all chains corresponding to the selected pdb_id
  blast_outfiles_paths <- c(Sys.glob(file.path(
    blast_outfiles_dir,
    paste(pdb_id, "*filtered*", sep =
            "")
  )))
  blast_output_list <-
    lapply(blast_outfiles_paths, function(x)
      read.table(x, header = T))
  blast_output <- as.data.frame(rbindlist(blast_output_list))
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
  
  return(blast_output)
}
