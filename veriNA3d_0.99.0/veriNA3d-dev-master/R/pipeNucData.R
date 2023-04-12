#' Obtain nucleotide details from a data set of RNA structures
#' 
#' Pipeline to generate a data.frame with the desired info for a list of PDB. 
#' Nucleotides are labeled with a unique identifier (column ntID).
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param model A vector with same length of pdbID containing the
#'     desired model for each pdbID. If all models are desired, use "all".
#'     If no models are specified, the first one will be used for each pdbID.
#' @param chain A vector with same length of pdbID containing the
#'     desired chain for each pdbID. If no chain is specified, all chains will
#'     be analysed by default. Non-nucleic acid chains will be ignored.
#' @param range A numeric vector with two values to indicate the desired
#'     length range for the Nucleic Acid chains. If a chain falls outside the
#'     range, it is not analysed.
#' @param path Directory in which the PDB/CIF files can be found (if NULL, the
#'     function will download them). If you provide a "path", make sure the
#'     file names are the PDB IDs followed by ".cif" or "pdb". The function
#'     will find them using the strings in pdbID, so make sure you use the 
#'     same case.
#' @param extension A string matching the files extension (e.g. ".pdb", 
#'     ".cif", "pdb.gz", "cif.gz").
#'     Only necessary if the PDB files are to be read from disk and a path is
#'     provided.
#' @param cores Number of CPU cores to be used.
#' @param progressbar A logical to print in screen a progress bar.
#' @param ... Arguments to be passed to [measureNuc()].
#'
#' @return A data.frame with data about every nucleotide in the input set.
#'
#' @examples 
#'     ## This is a toy example, see vignettes for real-world usages.
#'     pdblist <- list("1bau", "2rn1")
#'     model <- list("1", "2")
#'     chain <- list("all", "all")
#'     ntinfo <- pipeNucData(pdbID=pdblist, model=model, chain=chain)
#'
#' @author Diego Gallego
#'
pipeNucData <-
function(pdbID, model=NULL, chain=NULL, range=c(3, 100000),
            path=".veriNA3d_mmCIF_files", extension=".cif.gz", cores=1, 
            progressbar=TRUE, ...) {

    ## Make sure the input pdbID is a list -----------------------------------
    if (.isCIF(pdbID))
        pdbID <- list(pdbID)
    if (is.pdb(pdbID))
        pdbID <- list(pdbID)
    if (!is.list(pdbID))
        pdbID <- as.list(pdbID)

    ## Checking input vectors are equal in length ----------------------------
    if (is.null(model)) {
        model <- rep(1, length(pdbID))
    } else if (length(pdbID) != length(model)) {
        stop("pdbID and model vectors should have the same length")
    }
    if (is.null(chain)) {
        chain <- rep("all", length(pdbID))
    } else if (length(pdbID) != length(chain)) {
        stop("pdbID and chain vectors should have the same length")
    }

    ## Determine whether to read CIF/pdb objects from file, internet or input
    read <- .whereToRead(pdbID=pdbID, path=path, extension=extension,
                            verbose=progressbar)

    ## Download files if necessary
    if (any(read == "download.RAM")) {
        down <- unique(unlist(pdbID[which(read == "download.RAM")]))
        applyToPDB(listpdb=down, FUNCTION=cifDownload, 
                    cores=cores, progressbar=progressbar)
        print(paste("Download completed, saved in: ", tempdir()))
    }

    ## Print progress bar ----------------------------------------------------
    total <- length(pdbID)
    if (progressbar) {
        pbar <- txtProgressBar(min=0, max=total, style=3)
    } else {
        pbar <- NULL
    }

    ## Iterate over the list of entries to obtain the desired information ---- 
    ntinfo <- .xmapply(FUN=.manage_PDB,
                        index=seq_len(total),
                        pdbID=pdbID,
                        model=model,
                        chain=chain,
                        read=read,
                        mc.cores=cores,
                        MoreArgs=list(...=...,
                                        FUN=.make_chain_ntinfo,
                                        range=range,
                                        path=path,
                                        extension=extension,
                                        progressbar=progressbar,
                                        pbar=pbar), 
                        SIMPLIFY=FALSE)

    ## Print new line after progress bar -------------------------------------
    if (progressbar) {
        cat("\n")
    }

    ## Prepare output format -------------------------------------------------
    ntinfo <- ntinfo[which(lapply(ntinfo, length)>0)]
    if (length(ntinfo) == 0) {
        stop("Are you sure your input data is correct?")

    } else {
        ## Coerce list to data.frame
        ntinfo <- do.call(rbind, ntinfo)
        ntinfo <- cbind(ntID=seq_len(nrow(ntinfo)), ntinfo)
#        ntinfo$ntID <- seq_len(nrow(ntinfo))
        return(ntinfo)
    }
}
##############################################################################
## Subfunctions
## ===========================================================================
## Where should the input be read from?
.whereToRead <-
function(pdbID, path=NULL, extension=NULL, verbose=TRUE) {
    ## Preallocate and fill 'read' object ------------------------------------
    read <- vector("character", length(pdbID))

    ## If the user provides a path and extension -----------------------------
    if (!is.null(path) & !is.null(extension)) {
        ## Find files in path that match the extension
        files <- dir(path, pattern=extension)

        ## Check if input is the complete file name or without extension
        file_bool <- all(grepl(extension, pdbID, perl=TRUE))

        ## Find which of the input pdbID has its correspondent file
        if (file_bool) {
            inds <- which(pdbID %in% files)
        } else {
            inds <- which(paste(pdbID, extension, sep="") %in% files)
        }

        ## Fill 'read' vector
        read[inds] <- paste("read", extension, sep="")
        pdbID[inds] <- ""
    }

    ## If the pdbID input contains pdb objects -------------------------------
    if (any(unlist(lapply(pdbID, function(x) { 
                        return(class(x)[1] == "pdb") 
                    })))) {

        inds <- which(unlist(lapply(pdbID, function(x) {
                            return(class(x)[1] == "pdb")
                        })))
        ## Fill 'read' vector with string 'read.list'
        read[inds] <- "read.list"
        pdbID[inds] <- ""
    } 
    ## If the pdbID input contains CIF objects -------------------------------
    if (any(unlist(lapply(pdbID, function(x) { 
                        return(class(x)[1] == "CIF") 
                    })))) {

        inds <- which(unlist(lapply(pdbID, function(x) {
                                return(class(x)[1] == "CIF")
                            })))
        ## Fill 'read' vector with string 'read.list'
        read[inds] <- "read.list.cif"
        pdbID[inds] <- ""
    } 
    ## If the pdbID input contains 4 char PDB ID -----------------------------
    if (any(unlist(lapply(pdbID, function(x) {
                                            return(nchar(x) == 4)
                                        })))) {
        inds <- which(unlist(lapply(pdbID, function(x) {
                                            return(nchar(x) == 4)
                                        })))
        if (verbose) {
            print(paste(
                "The PDB IDs: ", 
                substr(paste(unique(pdbID[inds]), collapse="; "), 1, 40), 
                "... are going to be downloaded", 
                sep=""))
        }

        read[inds] <- "download.RAM"
        pdbID[inds] <- ""
    }
    
    ## If anything in input pdbID is not recognized, stop --------------------
    if (any(read == "")) { 
        inds <- which(read == "")
        stop("Check input pdbID:\n", paste(pdbID[inds], collapse="; "), 
                "\nThey should be 4 character PDB IDs or ",
                "files in path or ",
                "pdb/CIF objects", sep="")
    }
    return(read)
}

## ===========================================================================
## Intermediate wrapper that finds the pdb/CIF object and generates all the 
## possible model&chain combinations.
.manage_PDB <-
function(pdbID, model, chain, read, ..., 
            path=NULL, extension=NULL, index, pbar, FUN, progressbar=TRUE) {

    ## Find if the given pdb is multi model ----------------------------------
    if (length(model) == 1 && model == 1) {
        multi <- FALSE
    } else {
        multi <- TRUE
    }

    ## Save pdb ID if possible -----------------------------------------------
    if (read == "read.list") {
        name <- as.character(pdbID$call)[[1]]
    } else if (read == "read.list.cif") {
        name <- as.character(cifEntry(pdbID))
    } else {
        name <- pdbID
    }

    ## To cope with ALT records later set variables accordingly --------------
    ALT <- NULL
    rm.alt <- FALSE

    ## Save the atom coordinates in the form of pdb object -------------------
    if (read == "read.list") {
        temp_PDB <- pdbID
    } else if (read == "read.list.cif") {
        temp_PDB <- cifAsPDB(pdbID, alt=ALT)
    } else if (read %in% c("read.pdb", "read.pdb.gz", 
                            "read.ent", "read.ent.gz")) {
        temp_PDB <- suppressWarnings(read.pdb(
                                        paste(path, name, extension, sep=""), 
                                        multi=multi, 
                                        rm.alt=rm.alt, 
                                        verbose=FALSE))
    } else if (read %in% c("read.cif", "read.cif.gz")) {
        temp_PDB <- cifAsPDB(paste(path, name, extension, sep=""), alt=ALT)
    } else if (read == "download.RAM") {
        temp_PDB <- cifAsPDB(name, alt=ALT)
    }

    ## Make sure the pdb object has the necessary format ---------------------
    temp_PDB <- .perfect_input_format(temp_PDB)

    ## Save the different model numbers --------------------------------------
    if (model == "all" | model == 0) {
        model <- seq_len(nrow(temp_PDB$xyz))
    }
    ## Save the different chain ids ------------------------------------------
    if (chain == "all") {
        chain <- as.character(unique(temp_PDB$atom$chain))
    }
    ## Fins all the necessary combinations of models and chains --------------
    .combinations <- expand.grid(model, chain, stringsAsFactors=FALSE )
    names(.combinations) <- c("model", "chain")

    ## Iterate over every combination of chain and model to get data ---------
    FUN <- match.fun(FUN) # .make_chain_ntinfo
    ntinfo <- mapply(FUN=FUN,
                        model=.combinations[, "model"],
                        chain=.combinations[, "chain"],
                        MoreArgs=list(pdb=temp_PDB,
                                        name=name,
                                        ...=...),
                        SIMPLIFY=FALSE)

    ## Print progress bar
    if (progressbar) {
        setTxtProgressBar(pbar, index)
    }

    ## Remove empty entries of the list --------------------------------------
    ntinfo <- ntinfo[which(lapply(ntinfo, length) > 0)]

    ## Free up memory --------------------------------------------------------
    gc()

    ## Return output for every chain and model as given by input -------------
    if (length(ntinfo) == 0) {
        return()
    } else {
        ntinfo <- do.call(rbind, ntinfo)
        return(ntinfo)
    }
}

## ===========================================================================
## Takes a chain and model and calls the functions to get the desired data
.make_chain_ntinfo <-
function(pdb, model, chain, range, ..., name) {

    ## Selection of Model of interest ----------------------------------------
    pdb <- selectModel(pdb=pdb, model=model, verbose=FALSE)

    ## Select chain of interest ----------------------------------------------
    selection <- atom.select(pdb, chain=chain)

    ## pdb contains the PDB object ONLY with the selected model and chain ----
    pdb_ch <- trim(pdb, selection)

    ## Check that input contains a Nucleic Acid ------------------------------
    resid <- unique(pdb_ch$atom$resid)
    if (!any(resid %in% .nucleotides)) {
        string <- paste("Nothing to analyse in ", 
                    name, "|", model, "|", chain, " ", sep="")
        cat("\r   |", string)
        return()
    }

    ## Obtain number of (R/D)NA residues -------------------------------------
    refatm <- "C4'"
    reslist <- pdb_ch$atom$resno[which(pdb_ch$atom$elety == refatm)]
    total <- length(reslist)

    ## When reached this point, if the model does only contain the backbone
    ## phosphates, it's necessary to change the reference atom
    if (total == 0) {
        refatm <- "P"
        reslist <- pdb_ch$atom$resno[which(pdb_ch$atom$elety == refatm)]
        total <- length(reslist)
    }

    ## Check that the given chain is the desired length range ----------------
    if (total == 0 | total < range[1] | total > range[2]) {
        string <- paste("Nothing to analyse in ", 
                    name, "|", model, "|", chain, " ", sep="")
        cat("\r   |", string)
        return()
    }

    ## Final check for ALT records, in case of doubt use A -------------------
    if (any(is.na(pdb_ch$atom$alt))) {
        ind <- which(is.na(pdb_ch$atom$alt))
        pdb_ch$atom$alt[ind] <- "."
    }
    if (any(unique(pdb_ch$atom$alt) != ".")) {
        alt <- unique(pdb_ch$atom$alt)
        ind <- which(alt != ".")
        ## Find ALT strings and save to a vector
        validalt <- alt[ind]
        if (length(validalt) > 1) {
            ## If possible use A, otherwise just the first one
            if (any(validalt == "A")) {
                ALT <- c(".", "A")
            } else {
                ALT <- c(".", validalt[1])
            }
            ## Select the atoms using element numbers
            eleno <- pdb_ch$atom$eleno[pdb_ch$atom$alt %in% ALT]
            selection <- atom.select(pdb_ch, eleno=eleno)
            pdb_ch <- trim(pdb_ch, selection)
        }
    }

    pdb_ch <- .adddummy(pdb_ch, refatm=refatm)
    ## Check and measure the chain and make common data.frame ----------------
    ntinfo1 <- .checkNuc(pdb_ch, model=model, chain=chain, id=name, 
                            refatm=refatm, select=FALSE)
    ntinfo2 <- .measureNuc(pdb_ch, model=model, chain=chain, 
                            refatm=refatm, select=FALSE, ...)

    ntinfo <- cbind(ntinfo1, ntinfo2[, 
        which(!names(ntinfo2) %in% names(ntinfo1))])

    return(ntinfo)
}
