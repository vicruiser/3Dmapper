#!/usr/bin/env Rscript
#####################################################################
# BLAST OF EVERY CHAIN OF EVERY PDB STRUCTURE
# Extract the chain sequences of every PDB structures and
# save them separatedly in FASTA format with header.
#####################################################################
# Load necessary packages
# # Open connection to black hole
# con=file(open="/dev/null")
# # Don't print anything to screen
# sink(file=con, type="output")
# # Don't print messages (e.g. errors/warnings)
# sink(file=con, type="message")
options(echo = FALSE, verbose = F,warn = -1) 

requiredPackages = c("bio3d", "stringr", "data.table") #'parallel',
suppressMessages(
  for (p in requiredPackages) {
    if (!require(p, character.only = TRUE))
      install.packages(p)
    library(p, character.only = TRUE)
  }
)

###################################################################
# Note: to install veriNA3d, since it is in a private repository, #
# follow the instructions explained in this link.                 #
# (R version >= 3.5 is needed):                                   #
#                                                                 #
# ---> Optional: if you have a GitHub lab account, then you can   #
# use the following command in R:                                 #
#
# devtools::install_git("http://mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d.git",
#                                credentials = git2r::cred_user_pass("USER", getPass::getPass()))
# In USER: type your GitLab user. Password will be                #
# asked afterwards                                                #
# Note: to avoid this problems a docker container will            #
# be availabe soon                                                #
###################################################################

suppressMessages(
  if(!require("veriNA3d")){
  system('wget mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d-dev/repository/archive.zip?ref=master -O ~/veriNA3d_0.99.0.zip')
  system('unzip ~/veriNA3d_0.99.0.zip')
  system('mv ~/veriNA3d-*master* ~/veriNA3d_0.99.0')
  system('R CMD build ~/veriNA3d_0.99.0 --no-build-vignettes')
  system('R CMD INSTALL ~/veriNA3d*.tar.gz')
  
})

suppressMessages(library("veriNA3d", character.only = TRUE))


args = commandArgs(trailingOnly = TRUE)

PDB_file = args[1]
output_dir = args[2]

# decompress the downloaded pdb file
#R.utils::gunzip(PDB_file, remove = FALSE,  overwrite = TRUE)
# pdb filename (without the ".gz")
pdb_filename = PDB_file #substr(PDB_file, 1, nchar(PDB_file) - 3)


# read pdb file (watch out, this wont work for cif files, ckeck next for loop)
tryCatch({
  if (str_detect(pdb_filename, ".cif")) {
    cif_file <- cifParser(pdb_filename)
    pdb_file <- cifAsPDB(cif_file)
    
  } else if (str_detect(pdb_filename, ".pdb")) {
    pdb_file <- read.pdb(pdb_filename, verbose = F, rm.insert = F, rm.alt = F)
    
    
  } else{
    stop("Error!, wrong PDB file format. Accepted formats are .pdb or .cif.")
    
  }
  
  
  
},
error = function(err) {
  chains_seqs_table <- data.frame(
    pdbid = basename(pdb_filename),
    n_prot_chains = 0,
    n_read_prot_chains = 0,
    name_prot_chains = NA,
    name_read_prot_chains = NA,
    len_prot_chain = NA,
    n_unique_ligands = 0,
    ligands = NA,
    n_ligands = 0,
    n_nuc_chains = 0,
    nuc_chains =  NA,
    len_nuc_chain = NA
  )
  
  write.table(
    chains_seqs_table,
    file.path(output_dir, "PDB_prot_nuc_lig_db.txt"),
    col.names = !file.exists(file.path(output_dir, "PDB_prot_nuc_lig_db.txt")),
    quote = FALSE,
    row.names = F,
    append = T
  )
  gc()
})


# select only protein chains
prot_info <- atom.select(pdb_file, string = "protein", verbose = F)
lig_info <- atom.select(pdb_file, string = "ligand", verbose = F)
nuc_info <- atom.select(pdb_file, string = "nucleic", verbose = F)
chains_info <- pdb_file$atom[prot_info$atom,]
# count number of chains
nChains <- length(unique(chains_info$chain))
#count number of ligands
lig_chains <-
  unique(pdb_file$atom[lig_info$atom, c("resid", "insert", "resno")])$resid
n_lig_chains <- length(lig_chains)
n_unique_lig_chains <- length(unique(lig_chains))

#count number of nucleic chains
nuc_chains <- pdb_file$atom[nuc_info$atom ,]

if (length(unique(nuc_chains$resno)) > 4) {
  len_nuc_chains <-
    as.vector(sapply(split(nuc_chains, f = nuc_chains$chain), function(x)
      length(unique(x$resno))))
  nuc_chains <- unique(nuc_chains$chain)
  n_nuc_chains <- length(nuc_chains)
} else {
  #Update ligand chains and include nucleotide ligand
  lig_chains <-
    c(lig_chains, unique(nuc_chains[, c("resid", "insert", "resno")])$resid)
  n_lig_chains <- length(lig_chains)
  n_unique_lig_chains <- length(unique(lig_chains))
  #update nucleic chains information
  n_nuc_chains <- 0
  nuc_chains <- NA
  len_nuc_chains <- NA
}

if (n_lig_chains == 0) {
  lig_chains <- NA
}

# get every chain names
chainNames <- unique(chains_info$chain)
n_read_chains <- 0
name_read_chains <- c()
len_chain <- c()
# iterate over every chain to retrieve its sequence
for (j in 1:nChains) {
  # retrieve the chain sequence
  chain_info <- subset(chains_info, chain == chainNames[j])
  ChainSeq <-
    aa321(unique(chain_info[c("resno", "insert", "resid")])[, "resid"])
  len_chain <- c(len_chain, length(ChainSeq))
  
  if (length(ChainSeq) != 0) {
    # write sequence into fasta format
    seqinr::write.fasta(
      ChainSeq,
      names = paste(basename(pdb_filename), "_chain", chainNames[j], sep =
                      "") ,
      file.out = paste(
        output_dir,
        "/",
        basename(pdb_filename),
        "_chain",
        chainNames[j],
        ".fasta",
        sep = ""
      )
    )
    n_read_chains <- n_read_chains + 1
    name_read_chains <- c(name_read_chains, chainNames[j])
  } else{
    n_read_chains <- n_read_chains + 0
    next
  }
}
# remove pdb file (keep only the compressed one)
# tryCatch(
#   file.remove(pdb_filename),
#   error = function()
#     next
# )

if (nChains == 0) {
  chainNames <- NA
  name_read_chains <- NA
}

chains_seqs_table <- data.frame(
  pdbid = basename(pdb_filename),
  n_prot_chains = nChains,
  n_read_prot_chains = n_read_chains,
  name_prot_chains = paste(as.vector(chainNames), collapse =
                             "-"),
  name_read_prot_chains = paste(name_read_chains, collapse =
                                  "-"),
  len_prot_chain = paste(len_chain, collapse =
                           "-"),
  n_unique_ligands = n_unique_lig_chains,
  ligands = paste(lig_chains, collapse =
                    "-"),
  n_ligands = n_lig_chains,
  n_nuc_chains = n_nuc_chains,
  nuc_chains =  paste(nuc_chains, collapse =
                        "-"),
  len_nuc_chain = paste(len_nuc_chains, collapse =
                          "-")
  
)
write.table(
  chains_seqs_table,
  file.path(output_dir, "PDB_prot_nuc_lig_db.txt"),
  col.names = !file.exists(file.path(output_dir, "PDB_prot_nuc_lig_db.txt")),
  quote = FALSE,
  row.names = F,
  append = T
)
if (nChains == 0) {
  stop('Error: No protein chains found in this PDB file. ')
}


# Turn off output sink
# sink()
# # Turn off message sink
# sink(type="message")
# # Close connection to black hole
# close(con)