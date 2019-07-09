# path to packages
.libPaths(c( .libPaths(), "/home/bsc51/bsc51927/R/x86_64-pc-linux-gnu-library"))
# load packages
library(tidyr)
library(stringr)

# input arguments
INPUT_FILE = as.character(commandArgs(TRUE)[1])
OUTPUT_DIR = as.character(commandArgs(TRUE)[2])

# read the output file if not empty. 
if (file.size(INPUT_FILE) > 0) {
  # read the output files as data frame
  blast_out <- read.table(INPUT_FILE)
  # add column names
  colnames(blast_out) <-
    c(
      "qseqid",
      "qlen",
      "sseqid",
      "slen",
      "qstart",
      "qend",
      "sstart",
      "send",
      "evalue",
      "length",
      "pident",
      "nident",
      "qseq",
      "sseq"
    )
  # filter by percentage of sequence identity
  if (any(blast_out$pident >= 50)) {
    
    blast_out_50pident <- subset(blast_out, pident >= 50)
    # split pdb.id and chain variables into two columsn
    blast_out_50pident <-
      separate(blast_out_50pident, qseqid, c("pdbid", "chain"), sep = "_")
    # write new file 
    file_outputName <- file.path(OUTPUT_DIR,
                                 paste(
                                   unique(blast_out_50pident$pdbid),
                                   "_",
                                   unique(blast_out_50pident$chain),
                                   ".blast50pident",
                                   sep = ""
                                 ))
    write.table(
      blast_out_50pident,
      file_outputName,
      col.names = !file.exists(file_outputName),
      quote = FALSE,
      row.names = F
    )
    
    # keep a record of all the results together
    write.table(
      blast_out_50pident[, c(
        "pdbid",
        "chain",
        "qlen",
        "sseqid",
        "slen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "length",
        "pident",
        "nident"
      )],
      file.path(OUTPUT_DIR, "hits_pdb_ensembl_v94_blast_jan_2019_pident50.txt"),
      col.names = !file.exists(file.path(OUTPUT_DIR, "hits_pdb_ensembl_v94_blast_jan_2019_pident50.txt")),
      quote = FALSE,
      row.names = F,
      append = T
    )
  }
  
} else {
  # keep a record of which pdb id the corresponding chains that dit not score a single blast hit
  no_hits_file <- tools::file_path_sans_ext(basename(Sys.glob(INPUT_FILE)))
  no_hits_table<- str_split_fixed(no_hits_file, "_", 2)
  colnames(no_hits_table) <- c("pdb.id", "chain")
  file_outputName <- paste(OUTPUT_DIR, "nohits_BLAST.out", sep ="")
  
  write.table(
    no_hits_table,
    file_outputName,
    col.names = !file.exists(file_outputName),
    quote = FALSE,
    row.names = F,
    append=T
  )
  stop("No BLAST hits were found for this chain.")
}
