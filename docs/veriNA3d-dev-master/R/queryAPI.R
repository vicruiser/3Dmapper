#' Launch queries to the MMB or EBI APIs
#'
#' Given a 4-character string (PDB ID) and the desired "info", it sends a
#' query to the desired API and returns the output. This is an intermediate
#' wrapper called by most of the queryFunctions (for documentation see
#' ?queryFunctions).
#'
#' @param ID A 4 character string that matches a structure in the Protein 
#'     Data Bank, or a 3 character string matching a compound.
#' @param info A string with the desired query name.
#' @param API A string that matches "ebi" or "mmb".
#' @param string1 A string to configure the query. See example below.
#' @param string2 A string to configure the query. See example below.
#' @param reuse A logical. Set to TRUE if the same query is going to be send
#'     repeatedly, so that the result is saved in RAM (ir provides faster 
#'     user access and avoids unnecessary work in the servers).
#' @param envir Environment to save&retrieve the data if reuse is TRUE.
#' @param verbose A logical. TRUE to print details of the process.
#' @param noerrors A logical. TRUE to suppress errors.
#' @param https A logical. TRUE to use https instead of http.
#'
#' @return A vector or data.frame with the desired data.
#'
#' @examples
#'     ## Imagine you want to programmatically access the EBI API contents
#'     ## through "http://www.ebi.ac.uk/pdbe/api/topology/entry/1s72/chain/H".
#'     ## 'queryAPI' understands it with four intructions:
#'     ## 'API="ebi"' stands for the root of the website name ("http.../api/").
#'     ## 'string1' is the string from the root to the pdb ID.
#'     ## 'ID' is just the PDB code.
#'     ## 'string2' is the string after the pdb ID.
#'     ## Thus, the call would be:
#'     #data <- queryAPI(ID="1s72", API="ebi", 
#'     #                    string1="topology/entry/", string2="chain/H/")
#'
#' @author Diego Gallego
#'

## Higher level common function to make API calls
queryAPI <-
function(ID, info=NULL, API="default", string1=NULL, string2=NULL,
            reuse=TRUE, envir=parent.frame(n=2), verbose=FALSE,
            noerrors=FALSE, https=FALSE) {

    ## Check that the input ID is not over 4 character string ----------------
    if (nchar(ID) > 4) {
        stop("Please provide a correct PDB ID")
    }

    ## If a query name (info) is provided, get strings to make query ---------
    if (!is.null(info)) {
        strings <- .get_strings(info=info, API=API)
        API <- strings$API
        string1 <- strings$string1
        string2 <- strings$string2
        process <- TRUE
    } else if (is.null(string1) || is.null(string2) || API == "default") {
        stop("Provide a valid 'info' string or check documentation ",
                "for advanced usage of queryAPI", sep="")
    } else {
        API <- tolower(API)
        process <- FALSE
    }

    ## Find website name root
    if (API == "mmb") {
        webroot <- "http://mmb.pcb.ub.es/api/pdb/"
    } else if (API == "mmb_internal") {
        webroot <- "http://web.mmb.pcb.ub.es/MMBApi/web/pdb/"
    } else if (API == "ebi") {
        webroot <- "http://www.ebi.ac.uk/pdbe/api/"
    }

    ## Let the user decide if to use secure conection
    if (https) {
        webroot <- gsub("http", "https", webroot)
    }

    ## Generate string with the website name
    URL <- paste(webroot, string1, ID, sep="")
    if (string2 != "") {
        URL <- paste(URL, "/", string2, sep="")
    }
    if (verbose)
        print(paste("Querying: ", URL, sep=""))

    ## If data is in RAM, just retrieve it
    if (reuse) {
        ## Generate string to save/reuse data in RAM -------------------------
        if (is.null(info)) {
            info <- URL
        }
        infoname <- paste(".", toupper(ID), info, API, sep="")

        if (exists(infoname, envir=envir)) {
            if (verbose) 
                print(paste("Getting ", info, " from RAM", sep=""))
            return(get(infoname, envir=envir))
        }
    }

    ## Otherwise, send query to the API
    if (verbose)
        print(paste("Getting ", info, " from API", sep=""))

    ## Query to MMB API ------------------------------------------------------
    if (API %in% c("mmb", "mmb_internal")) {
        text <- tryCatch({
                    suppressWarnings(
                        .launchquery(URL, FUN=..launchquery, JSON=FALSE))
                    }, error = function(e) {
                        if (noerrors) {
                            return(NULL)
                        } else {
                            stop(e)
                        }
                    })
        #.launchquery(URL, FUN=..launchquery, JSON=FALSE)
        if (process) {
            output <- tryCatch({
                            .process_mmb_call(text, info, ID)
                        }, error = function(e) {
                            if (noerrors) {
                                return(NULL)
                            } else {
                                stop(e)
                            }
                        })
        } else {
            output <- text
        }
    }

    ## Query to EBI API ------------------------------------------------------
    if (API == "ebi") {
        text <- tryCatch({
                    suppressWarnings(
                        .launchquery(URL, FUN=..launchquery, JSON=TRUE))
                    }, error = function(e) {
                        if (noerrors) {
                            return(NULL)
                        } else {
                            stop(e)
                        }
                    })
#        text <- .launchquery(URL, FUN=..launchquery, JSON=TRUE)
        if (process) {
            output <- tryCatch({
                            .process_ebi_call(text, info)
                        }, error = function(e) {
                            if (noerrors) {
                                return(NULL)
                            } else {
                                stop(e)
                            }
                        })
        } else {
            output <- text
        }
    }

    ## Save in RAM if desired, so that a later (same) call of the function 
    ## will be faster --------------------------------------------------------
    if (reuse) {
        if (verbose)
            print(paste("Saving ", info , " in RAM", sep=""))
        assign(infoname, output, envir=envir)
    }

    return(output)
}
onlyebiqueries <- c("relDate", "revDate", "entities", "modres", "status")
onlymmbqueries <- c("header", "compType", "NDBId", "hetAtms",
                        "chains/header")
