#' Check if a given PDB contains the ligand/modbase of interest
#'
#' Given a 4-character string (PDB ID) and a ligand/modbase ID, the function
#' checks the presence of the ligand/modres in the given PDB ID. To check for
#' the presence of sodium ions use hetAtms="Na" instead of NA. If you are
#' interested on the whole list of heterogeneous atoms see 
#' [queryHetAtms()].
#'
#' @param pdbID A 4-character string.
#' @param hetAtms A string with the ligand/modbase ID.
#'
#' @return A logical. TRUE if the given hetAtms is present in the structure.
#'
#' @examples
#' hasHetAtm("1s72", "MG") # Check if structure has Magnesium ion
#'
#' @author Diego Gallego
#'

hasHetAtm <-
function(pdbID, hetAtms) {
    ## Make sure input is correct
    if (length(hetAtms) == 0) {
        stop("Introduce hetAtms")
    }
    ## If there's a NA (Sodium ion) in input, change it by a "Na" string
    if (any(hetAtms == "NA" || is.na(hetAtms))) { 
        ind <- which(hetAtms == "NA" || is.na(hetAtms))
        hetAtms[ind] <- "Na"
    }

    ## Query to find heteroatoms
    lig <- queryHetAtms(pdbID, NAtoNa=TRUE)

    if (length(lig) == 0) {
        ## If the structure has no heteroatoms, return as many FALSE as needed
        out <- rep(FALSE, length(hetAtms))
    } else {
        ## Iterate over list of hetAtms and save logical TRUE/FALSE
        out <- c()
        for (i in seq_along(hetAtms)) {
            out[i] <- any(grepl(lig, 
                                pattern=paste("^", hetAtms[i], "$", sep="")))
        }
    }

    ## Add names to the output vector
    names(out) <- hetAtms
    return(out)
}
