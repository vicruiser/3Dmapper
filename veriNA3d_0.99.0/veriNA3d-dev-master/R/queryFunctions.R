#' General functions to query PDB (Protein Data Bank) data 
#' 
#' Strightforward way to access structural data by making queries through the
#' EBI or MMB mirrors of the PDB.
#'
#' @param pdbID A 4-character string that matches a structure in the Protein 
#'     Data Bank. 
#' @param ... For advanced usage, arguments to be passed to subfunction 
#'     [queryAPI()].
#' @param chain A string with the chain identifier (in case you are only
#'     interested in a particular chain). If NULL, the info about all the
#'     chains is returned.
#' @param force A logical to force the query to the API (TRUE) or allow 
#'     checking presaved data.
#' @param subset Optional argument indicating "type", "length" or 
#'     "description". If NULL, all the columns in the data.frame are returned.
#' @param NAtoNa A logical. If TRUE, sodium ion (NA) is modified as "Na".
#' @param onlyligands A logical. If TRUE, the function only returns the list
#'     of unique ligands.
#' @param onlymodres A logical. If TRUE, only the modified residues are 
#'     returned.
#' 
#' @return A character vector or data.frame with the desired information:
#'     * {queryAuthors} List of authors.
#'     * {queryChains} Data frame with list of chains and properties.
#'     * {queryDescription} Author description of the entry.
#'     * {queryCompType} Type of entry as defined in PDB (e.g. Prot-nuc).
#'     * {queryDepdate} Deposition date.
#'     * {queryEntities} Data frame with list of entities and properties.
#'     * {queryFormats} Files available for the entry (e.g. to check if pdb 
#'                     format is available for the structure).
#'     * {queryHeader} Classification of the structure as it appears in the 
#'                     header (PDB format) or in the
#'                     "_struct_keywords.pdbx_keywords" field 
#'                     (mmCIF format).
#'     * {queryHetAtms} List of HETATM (modified residues and ligands).
#'     * {queryModres} List of modified residues.
#'     * {queryNDBId} NDB ID for Nucleic Acids.
#'     * {queryLigands} Retrieves ligands.
#'     * {queryOrgLigands} Retrieves just the organic ligands (not ions).
#'     * {queryReldate} Release date.
#'     * {queryResol} Resolution.
#'     * {queryRevdate} Revision date.
#'     * {queryStatus} PDB status.
#'     * {queryTechnique} Experimental technique.
#'
#' @examples
#'     queryTechnique("4y4o")
#'     queryAuthors("1s72")
#'     queryNDBId("1bau")
#'
#' @author Diego Gallego
#' @references Official PDBe REST API site:
#'     http://www.ebi.ac.uk/pdbe/pdbe-rest-api\cr
#'     Official MMB API site:
#'     http://mmb.irbbarcelona.org/api/
#' @name queryFunctions
#'
NULL
##############################################################################

#' @export
#' @rdname queryFunctions
queryAuthors <-
function(pdbID, ...) {
    info <- "autsList"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryChains <-
function(pdbID, chain=NULL, subset=NULL, ...) {
    pdbID <- tolower(pdbID)
    data <- queryAPI(pdbID, info="chains/header", ...=...)

    if (!is.null(chain)) {
        data <- data[data$chain == chain,]
    }
    if (!is.null(subset)) {
        data <- data[, subset]
    }
    return(data)
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryDescription <-
function(pdbID, ...) {
    info <- "compound"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryCompType <-
function(pdbID, ...) {
    info <- "compType"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryDepdate <-
function(pdbID, ...) {
    info <- "ascDate"
    output <- queryAPI(ID=pdbID, info=info, ...=...)
    output <- gsub(pattern="\\/", replacement="-", output, fixed=TRUE)
    return(output)
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryEntities <-
function(pdbID, ...) {
    info <- "entities"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryFormats <-
function(pdbID, ...) {
    info <- "formats"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryHeader <-
function(pdbID, ...) {
    info <- "header"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryHetAtms <-
function(pdbID, NAtoNa=TRUE, ...) {
    info <- "hetAtms"
    out <- queryAPI(ID=pdbID, info=info, ...=...)
    if (is.null(out)) return(NULL)
    if (NAtoNa && any(is.na(out))) out[which(is.na(out))] <- "Na"
    return(out)
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryModres <-
function(pdbID, onlymodres=FALSE, ...) {
    out <- queryAPI(pdbID, info="modres", ...=...) 
    if (is.null(out)) return(NULL)
    if (onlymodres && !is.na(out)) out <- out$chem_comp_id
    return(out)
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryNDBId <-
function(pdbID, ...) {
    info <- "NDBId"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryLigands <-
function(pdbID, onlyligands=FALSE, NAtoNa=TRUE, ...) {
    string1 <- "pdb/entry/ligand_monomers/"
    string2 <- ""
    out <- queryAPI(ID=pdbID, API="ebi", 
                    string1=string1, string2=string2)[[1]]

    if (NAtoNa) {
        if (any(is.na(out$chem_comp_id) | out$chem_comp_id == "NA")) {
            ind <- which(is.na(out$chem_comp_id) | out$chem_comp_id == "NA")
            out$chem_comp_id[ind] <- "Na"
        }
    }

    if (onlyligands) {
        return(unique(out$chem_comp_id))
    } else {
        return(out)
    }
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryOrgLigands <-
function(pdbID, ...) {
    ## Query info about all ligands (includes ions) --------------------------
    ligands <- queryLigands(pdbID, onlyligands=TRUE, NAtoNa=TRUE, ...=...)

    ## Check if there are ligands --------------------------------------------
    if (is.null(ligands)) return(NULL)

    ## If all ligands are in the list of ions, there are no organic ligands --
    if (sum(ligands %in% .ions) == length(ligands)) {
        return(NULL)
    }

    ## If reached this point, return the organic ligands ---------------------
    indices <- which(ligands %in% .ions)

    if (length(indices) == 0) {
        orgligands <- ligands
    } else {
        orgligands <- ligands[-indices]
    }

    return(orgligands)
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryReldate <-
function(pdbID, ...) {
    info <- "relDate"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions 
queryResol <- 
function(pdbID, force=FALSE, ...) {
    ## Check if presaved
    if (!force) {
        fast <- .fast_check(pdbID, "Resol")
        if (fast[[1]])
            return(fast[[2]])
    }

    info <- "resol"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryRevdate <-
function(pdbID, ...) {
    info <- "revDate"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryStatus <- 
function(pdbID, ...) {
    info <- "status"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################

#' @export
#' @rdname queryFunctions
queryTechnique <- 
function(pdbID, force=FALSE, ...) {
    ## Check if presaved
    if (!force) {
        fast <- .fast_check(pdbID, "Technique")
        if (fast[[1]])
            return(fast[[2]])
    }

    info <- "expType"
    return(queryAPI(ID=pdbID, info=info, ...=...))
}
##############################################################################
## Internal objects
## ===========================================================================

## Define a 'dictionary' of ions
.ions <- c(
    "2HP", "3CO", "ACT", "AG", "ALF", "AU3", "BA", "BEF", "BO4",
    "BR", "CA", "CAC", "CD", "CL", "CO", "CS", "CU", "F", "FE2", "FLC",
    "HG", "IOD", "IR3", "IRI", "IUM", "K", "LU", "MG", "MLI", "MMC", "MN",
    "Na", "NH4", "NI", "NO3", "OS", "PB", "PO4", "PT", "PT4", "RB", "RHD",
    "RU", "SE4", "SF4", "SO4", "SR", "TB", "TL", "UNX", "VO4", "ZN")
