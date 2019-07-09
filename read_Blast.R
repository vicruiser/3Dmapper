###############################################################################
# BLAST "PDB - Ensembl" results
###############################################################################
read.Blast <- function(pdb_id, blast_outfiles_dir) {
  # read all blast output files of all chains corresponding to the selected pdb_id
  blast_outfiles_paths <- c(Sys.glob(file.path(
    blast_outfiles_dir,
    paste(pdb_id, "*.blast50pident", sep =
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
      "ensembl.prot.id" ,
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
      "sseq"
    )
  
  # remove word "chain" from chain column for comparative purposes
  blast_output$chain <- sub("chain", "", blast_output$chain)
  blast_output$ensembl.prot.id <-
    sub("\\..*", "", blast_output$ensembl.prot.id)
  
  # select emsembl prot ids with pident >= 50%
  blast_output <- subset(blast_output, pident >= 50)
  
  return(blast_output)
}

