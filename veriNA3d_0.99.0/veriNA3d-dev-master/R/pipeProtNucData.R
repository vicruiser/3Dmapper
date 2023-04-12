#' Obtain nucleotide-protein interactions from a data set of structures
#'
#' Pipeline to generate a data.frame with the data about the closests 
#' nucleotides to the protein for a list of PDB. 
#' The data can be related to unique nucleotide
#' indentifiers (ntID) by providing the output of the independent
#' pipeline [pipeNucData()].
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param model A vector with same length of pdbID containing the
#'     desired model for each pdbID. If all models are desired, use "all".
#'     If no models are specified, the first one will be used for each pdbID.
#' @param chain A vector with same length of pdbID containing the
#'     desired chain for each pdbID. If no chain is specified, all chains will
#'     be analysed by default. Non-nucleic acid chains will be ignored.
#' @param ntinfo Optional. A data.frame obtained from 
#'     [pipeNucData()] 
#'     for the same dataset (or bigger), but not smaller.
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
#' @param cutoff A numeric with the maximum distance to return. To be passed 
#'     to [findBindingSite()].
#' @param ... Additional arguments to be passed to 
#'     [findBindingSite()].
#'
#' @return A data.frame with data about the atomic distances in the
#'     interacting sites of every structure in the input set.
#'
#' @examples 
#'     ## This is a toy example, see vignettes for more usages.
#'     pdblist <- list("1nyb", "2ms1")
#'     aantinfo <- pipeProtNucData(pdbID=pdblist)
#'
#' @author Diego Gallego
#'
pipeProtNucData <-
function(pdbID, model=NULL, chain=NULL, ntinfo=NULL,
            path=NULL, extension=NULL, cores=1, 
            progressbar=TRUE, cutoff=15, ...) {

    ## Make sure the input pdbID is a list -----------------------------------
    if (class(pdbID) == "CIF")
        pdbID <- list(pdbID)
    if (class(pdbID) == "pdb")
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
    read <- .whereToRead(pdbID=pdbID, path=path,
                            extension=extension, verbose=FALSE)

    ## Make sure that all pdbID are actually prot-nuc complexes --------------
    to_remove <- .find_wrong_input(pdbID, read, cores)
    if (length(to_remove) > 0) {
        pdbID <- pdbID[-to_remove]
        model <- model[-to_remove]
        chain <- chain[-to_remove]
        read  <- read[-to_remove]
    }
    if (length(pdbID) == 0) {
        stop("You pdbID input is not a valid set of prot-nuc structures")
    }

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
    interactionsdata <- .xmapply(FUN=.manage_PDB,
                                    index=seq_len(total),
                                    pdbID=pdbID,
                                    model=model,
                                    chain=chain,
                                    read=read,
                                    mc.cores=cores,
                                    MoreArgs=list(...=...,
                                                    FUN=.WrapperBindingSite,
                                                    ntinfo=ntinfo,
                                                    path=path,
                                                    extension=extension,
                                                    pbar=pbar,
                                                    progressbar=progressbar,
                                                    cutoff=cutoff),
                                    SIMPLIFY=FALSE)

    ## Print new line after progress bar -------------------------------------
    if (progressbar) {
        cat("\n")
    }

    ## Return output for every chain and model as given by input -------------
    interactionsdata <- interactionsdata[
                                which(lapply(interactionsdata, length) > 0)]

    ## Coerce list to data.frame
    output <- do.call(rbind, interactionsdata)

    return(output)
}
##############################################################################
## Subfunctions
## ===========================================================================
.WrapperBindingSite <-
function(pdb, model, chain, ..., name, ntinfo=NULL) {

    residues <- unique(pdb$atom[pdb$atom$chain == chain, "resid"])
    if (any(residues %in% .nucleotides)) {
        ## Check and measure the chain and make common data.frame ------------
        aantinfo <- findBindingSite(pdb=pdb,
                                            model=model,
                                            nchain=chain,
                                            ...)

        ## Add the nucleotide identifier if ntinfo was provided --------------
        if (!is.null(ntinfo)) {
            residues <- paste(name, model,
                                aantinfo$resno_A, aantinfo$insert_A, 
                                aantinfo$chain_A, sep=".")
            ntinfo_res <- paste(ntinfo$pdbID, ntinfo$model,
                                ntinfo$resno, ntinfo$insert, 
                                ntinfo$chain, sep=".")
    
            nt_id <- lapply(residues, FUN=function(x, ntinfo_res) {
                                            which(ntinfo_res == x)
                                        }, ntinfo_res=ntinfo_res)
    
            ## Not all contacts are going to be with a nucleotide
            if (any(lapply(nt_id, length) == 0)) {
                ind <- which(lapply(nt_id, length) == 0)
                nt_id[ind] <- ""
            }

            ntID <- ntinfo[unlist(nt_id), "ntID"]
            
            aantinfo <- cbind(ntID=ntID, 
                                pdbID=rep(name, nrow(aantinfo)), 
                                model=rep(model, nrow(aantinfo)),
                                aantinfo)

        ## In any case, add the pdbID and model to the output data.frame -----
        } else {
            aantinfo <- cbind(
                            pdbID=rep(name, nrow(aantinfo)),
                            model=rep(model, nrow(aantinfo)),
                            aantinfo)
        }
        return(aantinfo)

    ## If there is no nucleotide in the given chain, return nothing ----------
    } else {
        return()
    }
}

## Check that all pdbID are actually prot-nuc complexes
.find_wrong_input <-
function(pdbID, read, cores) {
    inds <- .xmapply(FUN=.is_wrong_input,
                        pdbID=pdbID,
                        read=read,
                        mc.cores=cores)
    return(which(inds))
}

.is_wrong_input <-
function(pdbID, read) {

    if (read == "read.list") {
        if (any(pdbID$calpha)) {
            return(FALSE)
        } else {
            return(TRUE)
        }

    } else if (read == "read.list.cif") {
        resid <- unique(cifAtom_site(pdbID)$label_comp_id)
        
        if (any(resid %in% .aa)) {
            return(FALSE)
        } else {
            return(TRUE)
        }

    } else {
        #name <- pdbID

        #if (sum(countEntities(name)[c("Prot", "Dprot")]) > 0) {
            return(FALSE)
        #} else {
        #    return(TRUE)
        #}
    }
}
