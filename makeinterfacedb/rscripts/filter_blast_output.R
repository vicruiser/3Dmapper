#! /usr/bin/Rscript
# path to packages
requiredPackages = c('tidyr', 'stringr')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  suppressMessages(library(p, character.only = TRUE))
}

# input arguments
INPUT_FILE = as.character(commandArgs(TRUE)[1])
SEQIDENT_FILTER = as.numeric(commandArgs(TRUE)[2])
COVERAGE_FILTER = as.numeric(commandArgs(TRUE)[3]) # percent of query length  that is included in the aligned segments
EVALUE_FILTER = as.numeric(commandArgs(TRUE)[4])
OUTPUT_DIR = as.character(commandArgs(TRUE)[5])
OUTPUT_FNAME = as.character(commandArgs(TRUE)[6])


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
      "sseq",
      "gaps"
    )
  #add query coverage colum
  blast_out$qcov = ((blast_out$length - blast_out$gaps) /  blast_out$slen)  * 100
  
  # filter by percentage of sequence identity
  if (any(blast_out$pident >= SEQIDENT_FILTER) |
      any(blast_out$qcov >= COVERAGE_FILTER) |
      any(blast_out$evalue <= EVALUE_FILTER)) {
    blast_out_filtered <- subset(
      blast_out,
      pident >= SEQIDENT_FILTER &
        qcov >= COVERAGE_FILTER &
        evalue <= EVALUE_FILTER
    )
    # split pdb.id and chain variables into two columsn
    blast_out_filtered <-
      separate(blast_out_filtered, qseqid, c("pdbid", "chain"), sep = "chain")
    # write new file
    write.table(
      blast_out_filtered,
      file.path(OUTPUT_DIR, OUTPUT_FNAME),
      col.names = !file.exists(file.path(OUTPUT_DIR, OUTPUT_FNAME)),
      quote = FALSE,
      row.names = F
    )
    
    
  }
  
} else {
  # keep a record of which pdb id the corresponding chains that dit not score a single blast hit
  no_hits_file <-
    tools::file_path_sans_ext(basename(Sys.glob(INPUT_FILE)))
  no_hits_table <- str_split_fixed(no_hits_file, "_", 2)
  colnames(no_hits_table) <- c("pdb.id", "chain")
  file_outputName <- file.path(OUTPUT_DIR, "nohits_BLAST.out")
  
  write.table(
    no_hits_table,
    file_outputName,
    col.names = !file.exists(file_outputName),
    quote = FALSE,
    row.names = F,
    append = T
  )
  stop("No BLAST hits were found for this chain.")
}
