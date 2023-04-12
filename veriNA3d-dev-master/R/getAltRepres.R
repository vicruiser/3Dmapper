#' Get Alternative Representants
#'
#' This function is closely related with getRNAList(). From its 
#' output, the family members of each equivalence class are checked for
#' desired features. The first member that matches all the desired features
#' is returned.
#'
#' @param rnalist The output of getRNAList.
#' @param technique One or more techniques of interest (For correct use, see 
#'  example below). For the list of techniques, see "veriNA3d:::.allowedtechs".
#' @param resol A positive real number to specify a desired resolution.
#' @param type A string indicating the type of desired RNA, according with 
#'  the {classifyRNA} function.
#' @param length To be passed to {classifyRNA}.
#' @param progressbar A logical to print in screen a progress bar.
#' @param verbose A logical to print details of the process.
#' @param as.df A logical to return the output as a data.frame.
#' @param na.rm A logical to clean missing values out of the output.
#'
#' @return A data.frame with info about all the "Equivalence Classes" and
#'  the selected Representants according with the specified conditions. 
#'
#' @examples 
#'     #rnalist <- getRNAList(release=3.2, threshold="1.5A")
#'     #alternative <- getAltRepres(rnalist=rnalist, 
#'     #                            type="nakedRNA")
#'
#' @author Diego Gallego
#'

getAltRepres <-
function(rnalist, technique=NULL, resol=NULL, type=NULL, length=3, 
            progressbar=TRUE, verbose=FALSE, as.df=FALSE, na.rm=TRUE) {
    ## Make sure the inputs make sense ---------------------------------------
    if (is.null(c(technique, resol, type))) {
        stop("Which features should the alternative representants have?")
    }

    if (!is.null(technique)) { 
        if (!all(technique %in% .allowedtechs)) {
            stop(paste("Introduce one or more techniques: ", 
                paste(.allowedtechs, collapse="; ") , sep=""))
        }
    } 

    if (!is.null(type) && !type %in% 
            c("protRNA", "DNARNA", "ligandRNA", "nakedRNA")) {

        stop(paste("Introduce a valid type (according with", 
            " the RNAclassifier): 'protRNA', 'nakedRNA', ",
            "'ligandRNA' or 'DNARNA'", sep=""))
    }

    if (!is.null(resol) && resol <= 0) {
        stop("'resol' must be a positive value")
    }

    ## Get necessary presaved data to speed up the process -------------------
    fastquery <- NULL
    data("fastquery", envir=environment())
    obsolete <- queryObsoleteList()

    ## Print progress bar ----------------------------------------------------
    total <- nrow(rnalist)
    if (progressbar) {
        if (verbose)
            print("Set progressbar to FALSE to get the verbose option")

        verbose <- FALSE
        pbar <- txtProgressBar(min=0, max=total, style=3)
    } else {
        pbar <- NULL
    }

    ## Do the real work ------------------------------------------------------
    rnalist$Representative <- invisible(mapply(
                                    FUN=.get_alternative_representant,
                                        seq_len(nrow(rnalist)),
                                        MoreArgs=list(data=rnalist,
                                                    technique=technique,
                                                    resol=resol, 
                                                    type=type,
                                                    length=length,
                                                    progressbar=progressbar,
                                                    pbar=pbar,
                                                    verbose=verbose,
                                                    fastquery=fastquery,
                                                    obsolete=obsolete)))
    if (progressbar)
        cat("\n")

    if (na.rm) {
        if (any(is.na(rnalist$Representative))) {
            rnalist <- rnalist[which(!is.na(rnalist$Representative)), ]
        }
    }

    if (as.df) {
        return(represAsDataFrame(rnalist[, 2:1]))
    } else {
        return(rnalist[, 2:1])
    }
}

###############################################################################
## Subfunctions
## ============================================================================
.get_alternative_representant <-
function(index, data, #eqclass, members,
            technique, resol, type, length,
            verbose, pbar, progressbar, fastquery, obsolete) {

    eqclass <- data[index, 1]
    members <- data[index, 3]

    if (verbose) 
        cat("\n", eqclass, "\n")

    members <- strsplit(members, split=";")[[1]]
    Members <- substr(members, 1, 4)

    if (is.null(technique)) { 
        Tech <- ""; technique <- "" 
    }
    if (is.null(resol)) { 
        Resol <- ""; resol <- " " 
    }
    if (is.null(type)) { 
        Type <- ""; type <- "" 
    }
    ## If the script does not find a proper representative, will return NA ---
    out <- NA

    for (i in seq_along(Members)) {
        ## Save pdbID instead of Leontis format ------------------------------
        pdbID <- substr(Members[i], 1, 4)

        ## Check status
        if (pdbID %in% obsolete) {
            status <- queryStatus(pdbID)
            if (verbose) {
                cat(paste(pdbID, " superceded by ", 
                    status$superceded_by, "... ", sep=""))
            }
            members[i] <- gsub(pdbID, x=members[i],
                                replacement=toupper(status$superceded_by))
            pdbID <- toupper(status$superceded_by)
            if (pdbID == "NA" | is.na(pdbID)) {
                next()
            }
        }

        if (verbose) 
            cat(paste(pdbID, " ... ", sep=""))

        ## Optimization with presaved data -----------------------------------
        if (pdbID %in% fastquery$pdbID) {
            fast <- TRUE
            fast_ind <- which(fastquery$pdbID == pdbID)
        } else {
            fast <- FALSE
        }

        ## If interested in any technique, query technique and cache result --
        if (technique[1] != "") {
            if (fast) {
                Tech <- fastquery[fast_ind, "Technique"]
            } else {
                Tech <- queryTechnique(pdbID, reuse=TRUE, 
                                        #verbose=verbose,
                                        envir=parent.frame(n=1))
            }
        }

        ## If necessary check resol and cache result -------------------------
        if (resol != " ") {
            if (fast) {
                Resol <- suppressWarnings(as.numeric(
                                fastquery[fast_ind, "Resol"]))
            } else {
                Resol <- as.numeric(queryResol(pdbID, reuse=TRUE,
                                                #verbose=verbose,
                                                envir=parent.frame(n=1)))
            }

            ## If Resol is NA, check whether the query actually made sense
            if (is.na(Resol) || length(Resol) == 0) {
                ## Make sure to know the experimental technique
                if (Tech == "") {
                    if (fast) {
                        Tech <- fastquery[fast_ind, "Technique"]
                    } else {
                        Tech <- queryTechnique(pdbID, reuse=TRUE,
                                                #verbose=verbose,
                                                envir=parent.frame(n=1))
                    }
                }
                ## If the technique is NMR, set resol to 0, so that any NMR 
                ## structure is accepted
                if (Tech %in% .nmrtechs) {
                    Resol <- 0
                } else {
                ## If the technique is not NMR and we don't have data about 
                ## the resolution, avoid using that structure
                    Resol <- 250 #Arbitrary value
                }
            }
        }

        ## If a particular type of RNA is specified, query and cache ---------
        if (type != "") {
            Type <- classifyRNA(pdbID, reuse=TRUE, length=length, 
                    verbose=verbose, 
                    envir=parent.frame(n=1))
        }

        if (tolower(Tech) %in% tolower(technique) && 
                Resol <= resol && Type %in% type) {

            out <- members[i]
            break()
        }
    }

    ## Print progress bar
    if (progressbar) {
        setTxtProgressBar(pbar, index)
    }

    if (verbose) 
        cat("\n")
    return(out)
}
############################################################################## 
## Subfunctions
## ===========================================================================
## None

## ===========================================================================
## Internal objects

.allowedtechs <- c(
    "X-RAY DIFFRACTION",        #X-RAY DIFFRACTION
    "SOLUTION NMR",             #SOLUTION NMR
    "SOLID-STATE NMR",          #SOLID-STATE NMR
    "ELECTRON MICROSCOPY",      #ELECTRON MICROSCOPY
    "FIBER DIFFRACTION",        #FIBER DIFFRACTION
    "FLUORESCENCE TRANSFER",    #FLUORESCENCE TRANSFER
    "POWDER DIFFRACTION",       #POWDER DIFFRACTION
    "ELECTRON CRYSTALLOGRAPHY", #ELECTRON CRYSTALLOGRAPHY
    "SOLUTION SCATTERING",      #SOLUTION SCATTERING
    "NEUTRON DIFFRACTION",      #NEUTRON DIFFRACTION
    "INFRARED SPECTROSCOPY")    #INFRARED SPECTROSCOPY

.nmrtechs <- c(
                "SOLUTION NMR", "SOLID-STATE NMR", 
                "Solution NMR", "Solid-state NMR")
