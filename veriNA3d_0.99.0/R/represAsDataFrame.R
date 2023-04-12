#' Coerce Representative list to a data.frame
#'
#' Takes the output of getRNAList or getAltRepres, which
#' represent molecules with the format "XXXX|M|C+XXXX|M|C" (XXXX: PDB ID; 
#' M: Model; C: Chain) and returns a data.frame with a more friendly 
#' structure:
#' \itemize{ 
#'     \item Col 1: Equivalence Class.
#'     \item Col 2: PDB ID.
#'     \item Col 3: Model.
#'     \item Col 4: Chain.
#' }
#' Columns 2 to 4 can be the direct input of 
#' [pipeNucData()].
#' 
#' @param nrlist The output of 
#'     [getRNAList()] or 
#'     [getAltRepres()].
#'
#' @return A data frame with the data of the representative structures.
#' 
#' @examples 
#'  data <- getRNAList(release=3.2, threshold="1.5A")
#'  reps <- represAsDataFrame(nrlist=data)
#'
#' @author Diego Gallego
#'
represAsDataFrame <-
function(nrlist) {
    ## Get repreesntative list -----------------------------------------------
    rep        <- nrlist$Representative
    names(rep) <- nrlist$Equivalence_class
    rep        <- sort(rep[!is.na(rep)])

    ## Manage "XXXX|M|C+XXXX|M|C" cases --------------------------------------
    rep <- unlist(strsplit(rep, split="+", fixed=TRUE))
    eq_classes <- names(rep)

    ## Generate new data.frame -----------------------------------------------
    rep <- as.data.frame(matrix(
                                unlist(strsplit(rep, split="|", fixed=TRUE)),
                                ncol=3,
                                byrow=TRUE),
                            stringsAsFactors=FALSE)
    out <- cbind(eq_classes, rep)
    names(out) <- c("Equivalence_class", "pdb", "model", "chain")
    if (any(nchar(out$pdb) > 4)) {
        out$pdb <- gsub(" ", "", out$pdb)
    }
    return(out)
}
