#' Downloads the list of ID of ALL current PDB entries
#' 
#' Function to get the list of ALL PDB IDs in the Protein Data Bank at the 
#' moment.
#'
#' @param justIDs A logical to return only pdbIDs without other information.
#' 
#' @return A vector with all the PDB ID entries (updated weekly).
#' 
#' @examples 
#' # pdblist <- queryEntryList()
#' 
#' @author Diego Gallego
#' 

queryEntryList <-
function(justIDs=TRUE){
    ## Send query
    URL <- "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt"
    out <- .launchquery(URL, FUN=read.table)

    ## Extract PDB IDs and sort them
    if (justIDs) {
        out <- toupper(as.character(out[, 1]))
    } else {
        names(out) <- c("pdbID", "type", "technique")
        out$pdbID <- toupper(as.character(out$pdbID))
        out$type <- as.character(out$type)
        out$technique <- as.character(out$technique)
    }
    return(out)
}

## Usually behind the most updated list
queryEntryList2 <-
function(){
    URL <- "http://mmb.pcb.ub.es/api/pdb/?fields=ids&noheaders"
    return(.launchquery(URL, FUN=..launchquery, JSON=FALSE))
}
