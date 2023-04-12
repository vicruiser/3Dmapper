#' Subset nucleotide data according with puckering
#'
#' Function to clean raw data after the pipeline 
#' [pipeNucData()]. 
#' It takes a data.frame that should have the columns "pu_phase", "delta" and 
#' "Dp", and returns the nucleotides matching the desired puckering state.
#'
#' @param ntinfo A data.frame. The output of 
#'     [pipeNucData()].
#' @param surenorth A logical to return nucleotides in north with restrictions
#'     in delta and Dp.
#' @param suresouth A logical to return nucleotides in south with restrictions
#'     in delta and Dp.
#' @param pucker A string with the puckering state of interest. Only necessary
#'     if surenorth and suresouth are FALSE and range is NULL. When using this
#'     option, only the phase is used to subset the data.
#' @param range A numeric vector with the desired phaserange. Only used if no
#'     other argument above could be applyed.
#' @param verbose A logical to print details of the process.
#'
#' @return An integer vector with the nucleotides (ntID) matching the desired 
#'     puckering state.
#'
#' @examples
#'     ntinfo <- pipeNucData("1bau")
#'     north_ntID <- cleanByPucker(ntinfo, surenorth=TRUE)
#'     north <- ntinfo[ntinfo$ntID %in% north_ntID,]
#' 
#' @author Diego Gallego
#'
cleanByPucker <-
function(ntinfo, surenorth=FALSE, suresouth=FALSE, 
            pucker="C3'endo", range=NULL, verbose=TRUE) {

    ## Check input -----------------------------------------------------------
    if (surenorth && suresouth) {
        stop("Please, chose only one puckering state")
    }

    ## Return north without outliers -----------------------------------------
    if (surenorth) {
        if (verbose) {
            cat("Returning nucleotides in north with phase between 342-54,",
                    " delta between 54-114 and Dp distance > 2.9 A\n", 
                    sep="")
        }
        return(ntinfo[which(complete.cases(ntinfo$pu_phase) &
                        (ntinfo$pu_phase > 342 | ntinfo$pu_phase < 54) &
                        (ntinfo$delta > 54 | ntinfo$delta < 114) &
                        ntinfo$Dp > 2.9), "ntID"])

    ## Return south without outliers -----------------------------------------
    } else if (suresouth) {
        if (verbose) {
            cat("Returning nucleotides in south with phase between 126-198,",
                    " delta between 117-177 and Dp distance < 2.9 A\n", 
                    sep="")
        }
        return(ntinfo[which(complete.cases(ntinfo$pu_phase) & 
                        ntinfo$pu_phase > 126 & ntinfo$pu_phase < 198 & 
                        (ntinfo$delta>117 | ntinfo$delta < 177) &
                        ntinfo$Dp < 2.9), "ntID"])
    }

    ## Check input -----------------------------------------------------------
    if (is.null(pucker) && is.null(range)) {
        stop("Introduce a valid pucker or range")
    } 

    ## If a pucker conformation is provided, the range is set automatically --
    if (!is.null(pucker)) {
        if (!is.character(pucker)) {
            stop("introduce valid pucker state")
        }

        if (!pucker %in% .puckerdata$pucker) {
            stop(paste("Introduce valid pucker state: ", 
                    paste(.puckerdata$pucker, collapse="; "), sep=""))
        }

        range <- as.integer(.puckerdata[.puckerdata$pucker == pucker, 
                            c("from", "to")])
    }
    
    ## Check the range is correct and use it to subset -----------------------
    if (range[2] > range[1]) {
        output <- ntinfo[which(ntinfo$pu_phase > range[1] &
                            ntinfo$pu_phase <= range[2]), "ntID"]
        if (verbose) {
            cat("Returning nucleotides in phase between: ", 
                range[1], "-",range[2], "\n", sep="")
        }
    } else {
        output <- ntinfo[which(ntinfo$pu_phase > range[1] |
                            ntinfo$pu_phase <= range[2]), "ntID"]
        if (verbose) {
            cat("Returning nucleotides in phase between: ", 
                range[2], "-",range[1], "\n", sep="")
        }
    }
    return(output)
}
##############################################################################
## Internal objects
## ===========================================================================
.puckerdata <- data.frame(pucker=c("north", "south", "east", "west", 
                                    "C3'endo", "C4'exo", "O4'endo",
                                    "C1'exo", "C2'endo", "C3'exo", 
                                    "C4'endo", "O4'exo", "C1'endo", 
                                    "C2'exo"),
                            from=c(315, 135, 45, 225, 0, 36, 72, 108,
                                    144, 180, 216, 252, 288, 324),
                            to=c(45, 225, 135, 315, 36, 72, 108, 144,
                                    180, 216, 252, 288, 324, 360),
                            stringsAsFactors=FALSE)
