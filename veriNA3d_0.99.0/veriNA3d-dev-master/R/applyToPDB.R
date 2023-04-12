#' Applies a function over a list of PDB ID entries.
#'
#' Given a function of interest, it is applied to all the PDB entries. 
#' See supported functions in ?queryFunctions.
#'
#' @param listpdb A list/vector containing the PDB IDs of interest. If NULL,
#'     the complete list of PDB entries is downloaded and used.
#' @param FUNCTION A function of interest.
#' @param as.df A logical that stands for "as.data.frame". If TRUE, the output
#'     will be returned in the form of a data.frame, otherwise a list.
#' @param cores Number of CPU cores to be used.
#' @param progressbar A logical to print in screen a progress bar.
#' @param ... optional arguments to FUNCTION.
#'
#' @return A data.frame with the PDB IDs (first colunm) and the output of the
#'     function of interest (second column) or a list with the results.
#'
#' @examples
#'     listpdb <- c("1s72", "1bau", "1rna")
#'     applyToPDB(listpdb, queryTechnique)
#'
#' @author Diego Gallego
#'

applyToPDB <-
function(listpdb=NULL, FUNCTION, as.df=TRUE, cores=1, progressbar=TRUE, ...) {

    ## Match function --------------------------------------------------------
    FUNCTION <- match.fun(FUNCTION)

    ## Download full PDB list if necessary -----------------------------------
    if (is.null(listpdb))
        listpdb <- queryEntryList()

    ## Print progress bar ----------------------------------------------------
    total <- length(listpdb)
    if (progressbar) {
        pbar <- txtProgressBar(min=0, max=total, style=3)
    } else {
        pbar <- NULL
    }

    ## Apply function over the list ------------------------------------------
    output_list <- .xlapply(seq_along(listpdb),
        FUN=function(i, FUNCTION, listpdb, pbar, progressbar, ...) {

            if (progressbar) {
                ## Print progress bar
                setTxtProgressBar(pbar, i)
            }

            tryCatch({
                return(FUNCTION(listpdb[i], ...))
            }, error=function(e) {
                return(NA)
            })
        }, FUNCTION=FUNCTION, listpdb=listpdb, pbar=pbar, 
            progressbar=progressbar, ...=..., 
            mc.cores=cores, mc.preschedule=TRUE)

    if (progressbar) {
        cat("\n")
    }

    ## Any query might receive NA with a certain frequency, even when it 
    ## shoudn't, thus the NA in the list are double-checked ------------------
    torepeat <- which(is.na(output_list))
    if (length(torepeat) != 0) {

        for (i in torepeat) {
            tryCatch({
                output_list[[i]] <- FUNCTION(listpdb[i], reuse=FALSE, ...)
            }, error=function(e) {
                output_list[[i]] <- NA
            })
        }
    }

    ## Give format to the output ---------------------------------------------
    if (as.df) {
        ## If a list entry has more than one element, concatenate them
        if (any(unlist(lapply(output_list, length)) != 1)) {
            inds <- which(unlist(lapply(output_list, length)) != 1)
            for (j in inds) {
                output_list[[j]] <- paste(output_list[[j]], collapse="|")
            }
        }

        ## If a list entry is NULL, replace by empty string
        if (any(unlist(lapply(output_list, is.null)))) {
            inds <- which(unlist(lapply(output_list, is.null)))
            for (j in inds) {
                output_list[[j]] <- ""
            }
        }

        ## Create data.frame
        output <- as.data.frame(
                    matrix(c(listpdb, unlist(output_list)),
                            byrow=FALSE,
                            ncol=2),
                    stringsAsFactors=FALSE)
    } else {
        names(output_list) <- listpdb
        output <- output_list
    }
    return(output)
}
