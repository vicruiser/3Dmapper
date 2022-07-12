#! /usr/bin/Rscript
#' Title  READ PDB FILE
#'
#' @param pdb_file
#' @param download
#'
#' @return
#' @export
#'
#' @examples
options(echo = FALSE, verbose = F,warn = -1) 
readPDB <- function(pdb_file) {
  biounit_filename = basename(pdb_file)
  pdb_filename <- sub("\\.gz", "", biounit_filename)
  
  #Check whether the PDB exists or not.
  if (!file.exists(pdb_file)) {
      tryCatch({
        pdb <- read.pdb(pdb_file, verbose =F, rm.alt=F, rm.insert=F)
      }, error = function(e) {
        cat("ERROR : The input PDB does not exists. \n")
      })
  
    }else {
    # If the file exists or after it is downloade we can read it.
    #if (!file.exists(file.path(dirname(pdb_file), pdb_filename))) {
    #  R.utils::gunzip(pdb_file, remove = FALSE, overwrite = TRUE)
    #}
    
    if (any(str_detect(pdb_filename, "\\.cif"))) {
      cif_filepath <- file.path(dirname(pdb_file), pdb_filename)
      cif <- cifParser(pdb_file)#cif_filepath)
      pdb <- cifAsPDB(cif)
      
    } else if (any(str_detect(pdb_filename, "\\.pdb"))) {
      #pdb_filepath <- file.path(dirname(pdb_file), pdb_filename)
      pdb <- read.pdb(pdb_file, verbose =F, rm.alt=F, rm.insert=F)#path)
      
      } else {
      stop("Sorry, wrong PDB structure format. The file must have '.pdb' or '.cif' extensions.")
      
    }
  }
  return(pdb)
}

