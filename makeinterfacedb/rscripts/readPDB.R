#' Title  READ PDB FILE
#'
#' @param pdb_file
#' @param download
#'
#' @return
#' @export
#'
#' @examples

readPDB <- function(pdb_file,
                    download = "NO") {
  biounit_filename = basename(pdb_file)
  pdb_filename <- sub("\\.gz", "", biounit_filename)
  
  #Check whether the PDB exists or not.
  if (!file.exists(pdb_file)) {
    if (download == "YES") {
      # Read files accordingly to their format (mmCIF or PDB)
      if (any(str_detect(biounit_filename, "\\.cif"))) {
        # Choose between downloading the files or use an already downloaded local data base.
        pdb_url <-
          file.path("ftp://ftp.pdbj.org/pub/pdb/data/biounit/mmCIF/all",
                    biounit_filename)
      } else if (any(str_detect(biounit_filename, "\\.pdb"))) {
        pdb_url <-
          file.path("ftp://ftp.pdbj.org/pub/pdb/data/biounit/PDB/all",
                    biounit_filename)
      }
      #Download file
      tmpdir <- tempdir()
      file <- basename(pdb_url)
      download.file(pdb_url, file)
      R.utils::gunzip(file)
      
    } else if (download == "NO") {
      stop ("Sorry, it is not possible to download the PDB file. Please, provive a PDB file.")
      
    }
  } else {
    # If the file exists or after it is downloade we can read it.
    #if (!file.exists(file.path(dirname(pdb_file), pdb_filename))) {
    #  R.utils::gunzip(pdb_file, remove = FALSE, overwrite = TRUE)
    #}
    
    if (any(str_detect(pdb_filename, "\\.cif"))) {
      cif_filepath <- file.path(dirname(pdb_file), pdb_filename)
      cif_file <- cifParser(pdb_file)#cif_filepath)
      pdb_file <- cifAsPDB(cif_file)
      
    } else if (any(str_detect(pdb_filename, "\\.pdb"))) {
      #pdb_filepath <- file.path(dirname(pdb_file), pdb_filename)
      pdb_file <- read.pdb(pdb_file)#path)
      
    } else {
      stop("Sorry, wrong PDB structure format.")
      
    }
  }
  return(pdb_file)
}
