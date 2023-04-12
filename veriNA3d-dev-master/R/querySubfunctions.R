## Query Subfunctions
## ===========================================================================
## Ping function inspired in:
## "http://stackoverflow.com/questions/5076593/how-to-determine-if-you-
## have-an-internet-connection-in-r"
.check_internet <- 
function(url="http://mmb.irbbarcelona.org/www/") {
    # test the http capabilities of the current R build
    if (!capabilities(what="http/ftp")) 
        return(FALSE)
    # test connection by trying to read first line of url
    test <- try(suppressWarnings(readLines(url, n=1)), silent=TRUE)
    # return FALSE if test inherits 'try-error' class
    !inherits(test, "try-error")
}

## ===========================================================================
## Send query function and handle errors, adapted from bioconductor template
.launchquery <-
function(URL, FUN, ..., N.TRIES=3L, SLEEP=0.05) {
    ## Match function
    FUN <- match.fun(FUN)
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

    while (N.TRIES > 0L) {
        result <- tryCatch(FUN(URL, ...=...), error=identity)
        if (!inherits(result, "error"))
            break
        N.TRIES <- N.TRIES - 1L
        Sys.sleep(SLEEP)
        SLEEP <- SLEEP * 1.5
    }

    if (N.TRIES == 0L) {
        stop(paste("query failed:",
                "  URL: ", URL,
                "  error: ", conditionMessage(result), sep=""))
    }

    return(result)
}

## ===========================================================================
## Basic query
#' @importFrom jsonlite fromJSON
..launchquery <-
function(URL, JSON=FALSE) {
    ## Open connection and scan site
    con <- url(URL)
    if (JSON) {
        text <- fromJSON(con)
    } else {
        ## Make sure the connection is closed even if an error occurs
        on.exit(close.connection(con))
        text <- scan(con, "character", quiet=TRUE)
    }
    return(text)
}

## ===========================================================================
## From a type of query (given by info) and API, returns the necessary data to
## make the call
.get_strings <-
function(info, API) {
    ## If no API is selected, set API ----------------------------------------
    if (API == "default") {
        if (info %in% onlymmbqueries) {
            API <- "mmb"
        } else {
            API <- "ebi"
        }
    } else {
        ## Check that the selected API is in the list below ------------------
        API <- tolower(API)
        .check_api(API, supported=c("mmb", "ebi", "mmb_internal"))
    }

    ## Warn the user that some data is only taken from mmb/ebi api -----------
    if (info %in% onlymmbqueries && API == "ebi") {
        .APIwarning("MMB")
        API <- "mmb"

    } else if (info %in% onlyebiqueries && 
                API %in% c("mmb", "mmb_internal")) {
        .APIwarning("EBI")
        API <- "ebi"
    }

    ## Return the strings to make the correct call to the API ----------------
    if (API %in% c("mmb", "mmb_internal")) {
        string1 <- ""
        string2 <- paste("entry/", info, "/", sep="")
    } else if (API == "ebi") {
        string1 <- .translate_string(info)
        string2 <- ""
    }

    return(list(API=API, string1=string1, string2=string2))
}

## ===========================================================================
## Translate the info string into the correct query names for the EBI API
.translate_string <-
function(info) {
    if (info %in% c("expType", "compound", "autsList", 
            "ascDate", "relDate", "revDate")) {
        return("pdb/entry/summary/")

    } else if (info == "formats") {
        return("pdb/entry/files/")

    } else if (info == "resol") {
        return("pdb/entry/experiment/")

    } else if (info == "entities") {
        return("pdb/entry/molecules/")

    } else if (info == "modres") {
        return("pdb/entry/modified_AA_or_NA/")

    } else if (info == "status") {
        return("pdb/entry/status/")

    } else {
        stop("Query not supported")
    }
}

## ===========================================================================
## Function to manage the data retrieved from the mmb API to return to the 
## user the desired output
.process_mmb_call <-
function(text, info, pdbID) {
    if (info %in% c("hetAtms", "formats")) {
        start <- grep("[", text, fixed=TRUE)
        end <- grep("]", text, fixed=TRUE)

        if (length(start) == 0 | length(end) == 0)
            return(NULL)

        if ((end-start) == 1) 
            return(NULL)

        text <- text[(start+1):(end-1)]
        if (!all(is.na(text)) && any(text == ",")) {
            text <- text[-which(text == ",")]
        }
        return(text)
    } else if (info == "chains/header") {
        ind <- grep(pdbID, text)
        ind <- ind[-1]
        text <- as.data.frame(
            matrix(unlist(strsplit(text[ind], split="  ")), 
                    ncol=3, byrow=TRUE), stringsAsFactors=FALSE)
        text <- as.data.frame(
                    cbind(do.call(rbind, strsplit(text$V1, " ")),
                            do.call(rbind, strsplit(text$V2, " ")),
                            text$V3),
                    stringsAsFactors=FALSE)
        names(text) <- c("pdbID", "chain", "type", "length", "description")
    return(text)
    } else {
        return(text[grep(pattern=info,text)+2])
    }
}
## ===========================================================================
## Function to manage the data retrieved from the ebi API to return to the
## user the desired output
.process_ebi_call <-
function(text, info) {
    if (is.null(text)) 
        return(NULL)
    if (info == "expType") {
        return(text[[1]]$experimental_method[[1]])
    } else if (info == "formats") {
        return(text[[1]]$PDB$downloads)
    } else if (info == "compound") {
        return(text[[1]]$title)
    } else if (info ==  "autsList") {
        return(text[[1]]$entry_authors[[1]])
    } else if (info == "ascDate") {
        return(text[[1]]$deposition_date)
    } else if (info == "relDate") {
        return(text[[1]]$release_date)
    } else if (info == "revDate") {
        return(text[[1]]$revision_date)
    } else if (info == "resol") {
        return(text[[1]]$resolution)
    } else if (info == "entities") {
        text <- text[[1]][order(text[[1]]$entity_id),]
        return(text)
    } else if (info == "status") {
        return(text[[1]])
    } else if (info == "modres") {
        return(text[[1]])
    }    
}

## ===========================================================================
## Function to check if the input string for the API is valid
.check_api <-
function(API, supported = c("mmb", "ebi")) {
    if (!API %in% supported) {
        stop("Introduce valid API (mmb or ebi)")
    }
}

## ===========================================================================
.APIwarning <-
function(apiname) {
    warning(paste("The desired info can only be retrieved",
                "from the ", apiname ," API", sep=""))
}
## ===========================================================================
