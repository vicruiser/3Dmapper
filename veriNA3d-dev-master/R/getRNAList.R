#' Download Representative List of RNA structures
#'
#' According with Leontis & Zirbel (references below), the PDB contains 
#' structures that are redundant. Their server provides weekly releases 
#' of Representative Sets of RNA structures (also called non-redundant 
#' lists). This function access their website and returns the desired list.
#'
#' @param release A number indicating the list of interest.
#' @param threshold A string that matches one of the lists in the BGSU RNA 
#'     site ("1.5A", "2.0A", "2.5A", "3.0A", "3.5A", "4.0A", "20.0A", "all").
#'     Note that "all" returns structures solved by any technique.
#' @param as.df A logical to return the output as a data.frame
#'
#' @return A data frame with the list of Equivalence Classes and the 
#'     Representant and Members of each Equivalence Class. Note that the 
#'     output is formated according with Leontis&Zirbel nomenclature 
#'     (AAAA|M|C), where "AAAA" is the PDB accession code, "M" is the model
#'     and "C" is the Chain to be used.
#'
#' @examples 
#'     data <- getRNAList(release=3.2, threshold="1.5A")
#'
#' @author Diego Gallego
#' @references 
#'     Official site: http://rna.bgsu.edu/rna3dhub/nrlist/
#'     Publication: Leontis, N.B., and C.L. Zirbel. 2012. “Nonredundant 3D 
#'     Structure Datasets for RNA Knowledge Extraction and Benchmarking.” 
#'     In RNA 3D Structure Analysis and Prediction, edited by N. Leontis and
#'     E. Westhof, 27:281–98. Springer Berlin Heidelberg
#'
getRNAList <- 
function(release="current", threshold="all", as.df=FALSE) {

    ## Check argument threshold is correct -----------------------------------
    if (!threshold %in% thresholds) { 
        stop(paste('threshold can only be into one of the following ',
                    'categories: "1.5A", "2.0A", "2.5A", "3.0A", "3.5A", ',
                    ' "4.0A", "20.0A" or "all"', sep=""))
    }

    ## Check Leontis server --------------------------------------------------
    URL <- "http://rna.bgsu.edu/rna3dhub/nrlist"
    if (!.check_internet(url = URL)) {
        stop('Server not responding, try again later.')
    }

    ## Find last release -----------------------------------------------------
    if (release == "current") {
        release <- .getLastReleaseOLD()[[1]]
    }

    ## Dowload data ----------------------------------------------------------
    URL <- paste(URL, "/release/", release, "/", threshold, sep="")
    suppressWarnings(text <- .launchquery(URL=URL, FUN=readLines, N.TRIES=2))

    ## Prepare string to grep lines ------------------------------------------
    if (strsplit(threshold, split="")[[1]][
        length(strsplit(threshold, split="")[[1]])] == "A") {

        threshold <- substr(threshold, 1, nchar(threshold) - 1)
    }
    indices <- grep(pattern=paste("NR_", threshold, "_", sep=""), text)

    ## Give format to the data -----------------------------------------------
    data <- unlist(lapply(seq_along(indices), 
                    FUN=.see_equivalence_class,
                        text=text, indices=indices,
                        release=release, threshold=threshold))

    ## Updated: 2017-Jul-21 CORNER CASE, 4R3I does not have model "0" any more
    ind <- grep(pattern="4R3I|0", data, fixed=TRUE) 
    if (length(ind) > 0) {
        data[ind] <- gsub(pattern="0", replacement="1", x=data[ind])
    }

    output <- as.data.frame(matrix(data, ncol=3, byrow=TRUE),
                            stringsAsFactors=FALSE)
    names(output) <- c("Equivalence_class", "Representative", "Class_members")

    if (as.df) {
        output <- represAsDataFrame(output)
    }
    return(output)
}

###############################################################################
## Subfunctions
## ============================================================================

## Manage format of data for a given entry (equivalence class)
.see_equivalence_class <- function(x, text, indices, release, threshold) {
    splited_terms <- strsplit(text[indices[x]], split="/")[[1]]

    ## Obtain the name of the equivalence class ------------------------------
    ind1 <- grep(pattern=paste("NR_", threshold, "_", sep=""), splited_terms)
    if (length(ind1)!=1) {
    stop("!")
    }
    eq_class <- strsplit(splited_terms[ind1], split="\"")[[1]][1]

    ## Obtain the representant of the equivalence class ----------------------
    ind2 <- grep(pattern="strong class=\"pdb\"", splited_terms)
    if (length(ind2)!=1) {
        stop("!")
    }
    repstr <- substring(strsplit(splited_terms[ind2], 
                                    split="\ ")[[1]][1], first=8)

    ## Obtain the list of members of the equivalence class -------------------
    ind3 <- grep(pattern="class='pdb'", splited_terms)
    all_members <- unlist(lapply(ind3, FUN=.kk, splited_terms))

    return(c(eq_class, repstr, paste(all_members, collapse=";")))
}

.kk <- function(.ind3, splited_terms) {
    temp <- strsplit(splited_terms[.ind3], split=">")[[1]]
    return(substr(temp[length(temp)], 1, nchar(temp[length(temp)])-1))
}

thresholds <- 
    c("1.5A", "2.0A", "2.5A", "3.0A", "3.5A", "4.0A", "20.0A", "all")

## ============================================================================
## Get Leontis List last release number and date
.getLastRelease <-
function() {
    ## Find last release -----------------------------------------------------
    URL <- "http://mmb.irbbarcelona.org/api/RNANRList/info"
    info <- .launchquery(URL=URL, FUN=..launchquery, JSON=T, N.TRIES=2L)

    ## Save release info -----------------------------------------------------
    return(list(info$Data$release, info$Data$date))
}
.getLastReleaseOLD <-
function() {
    ## Find last release -----------------------------------------------------
    URL <- "http://rna.bgsu.edu/rna3dhub/nrlist/release/current/1.5A"
    text <- .launchquery(URL=URL, FUN=readLines, N.TRIES=2L, n=200)
    info <- text[grep("<small>Release", text)]
    
    info <- strsplit(info, "</small>")[[1]]
    info <- strsplit(info, "<small>")[[1]][2]
    
    info <- strsplit(info, ", ")[[1]]
    
    ## Save release info -----------------------------------------------------
    release <- info[1]
    releasenum <- gsub("Release ", "", release)
    date <- info[2]

    return(list(releasenum, date))
}
## ============================================================================
