#' Trim a pdb/cif object to obtain a nucleotide/s of interest and the
#' surrounding area.
#'
#' From a pdb/CIF object, the nucleotide of interest and a radius, the 
#' function finds all the atoms in the given area and returns a pdb object
#' that only includes the nearest atoms.
#'
#' @param cif A cif/pdb object obtained from cifParser/read.pdb respectively
#'     or a pdb ID so that the function can download the data.
#' @param model The model of interest to use in the calculations. The first 
#'     model is always the default.
#' @param ntindex A numeric index/indices for the position of the desired
#'     nucleotides in the given chain. Not necessary if you provide sel (see 
#'     below).
#' @param chain A string with the chain of interest. Not necessary if you 
#'     provide sel (see below).
#' @param sel A "select" object as obtained from atom.select (bio3d). Note 
#'     that if you are using this option, cif must be the same input object
#'     you used for the atom.select function.
#' @param cutoff A numeric indicating the radius in angstroms to select around 
#'     the desired nucleotides. If 0 only the nucleotides are returned.
#' @param cutres A logical. TRUE to return only what it is found in the cutoff
#'     (residues in the boundaries of the cutoff are usually truncated) or 
#'     FALSE to return whole residues even if further than the cutoff.
#' @param file A string to save the output in a pdb formated file. If NULL the
#'     fucntions just returns the pdb object.
#' @param verbose A logical to print details of the process.
#' @param alt Value of alternative records to keep in output object.
#' @param rename A logical to rename the chains of the oucput pdb.
#' @param ... Arguments to be passed to internal functions.
#'
#' @return A smaller pdb object or a pdb file. 
#'
#' @examples
#'     ## Toy example:
#'     cif <- cifParser("1s72")
#'
#'     ## Generate a smaller pdb with the residues 55 to 58 of the RNA chain
#'     ## "9" with a sorrounding sphere of 5 Angstroms:
#'     smallerpdb <- trimSphere(cif, ntindex=seq(55, 58, 1), chain="9", 
#'                                 cutoff=5, verbose=FALSE)
#'
#'     ## Same process saving the output in a file:
#'     smallerpdb <- trimSphere(cif, ntindex=seq(55, 58, 1), chain="9", 
#'                                 cutoff=5, verbose=FALSE, file="output.pdb")
#'
#'     ## Second example:
#'     ## Obtain a PDB with just the interacting region between RNA and prot
#'     #pdb <- cifAsPDB("1nyb")
#'     #data <- findBindingSite(pdb, select="RNA", byres=TRUE)
#'     #sel <- bio3d::atom.select(pdb,
#'     #                            eleno=append(data$eleno_A, data$eleno_B))
#'     #trimSphere(pdb, sel=sel, file="interacting_site.pdb", verbose=FALSE)
#'
#' @author Diego Gallego
#'
trimSphere <-
function(cif, model=NULL, ntindex=NULL, chain=NULL, sel=NULL, cutoff=8, 
            cutres=FALSE, file=NULL, verbose=TRUE, alt="uniq", rename=TRUE, 
            ...) {

    ## Make sure the object is a S3 pdb object with the desired model --------
    cif <- .input_to_pdb(cif=cif, model=model, verbose=verbose, alt=NULL, ...)

    ## Make sure the pdb object has the necessary format ---------------------
    cif <- .perfect_input_format(cif)

    ## Find eleno numbers ----------------------------------------------------
    data <- paste(cif$atom$resno, cif$atom$insert, cif$atom$chain, sep="|")
    if (is.null(sel)) {
        if (is.null(chain)) {
            stop("Please, provide a sel or chain argument")
        }
        inds <- which(cif$atom$elety == "C4'" & cif$atom$chain == chain)
        inds2 <- which(duplicated(data[inds]))
        if (length(inds2) != 0) {
            inds <- inds[-inds2]
        }
        if (length(inds) == 0) {
            inds <- which(cif$atom$elety == "CA" & cif$atom$chain == chain)
        }
        if (length(inds) == 0) {
            inds <- which(cif$atom$chain == chain)
        }
        if (is.null(ntindex)) {
            ntindex <- seq_along(inds)
        }
        resno <- cif$atom$resno[inds][ntindex]
        insert <- cif$atom$insert[inds][ntindex]
        query <- paste(resno, insert, chain, sep="|")
        refeleno <- cif$atom$eleno[data %in% query]
    } else {
        refeleno <- cif$atom$eleno[sel$atom]
        chain <- unique(cif$atom$chain[sel$atom])
        query <- unique(data[refeleno])
    }
    eleno <- cif$atom$eleno

    ## Find closest neighbors in radius=cutoff -------------------------------
    if (cutoff > 0) {
        dis_map <- measureElenoDist(pdb=cif, 
                    refeleno=refeleno, 
                    eleno=eleno, 
                    n=NULL, 
                    cutoff=cutoff,
                    detailedoutput=TRUE,
                    data_of_interest=c("resno",
                                        "insert",
                                        "chain"),
                    verbose=verbose)

        ## Select whole residues if necessary --------------------------------
        if (cutres) {
            outeleno <- unique(dis_map$eleno_B)
        } else {
            query2 <- unique(paste(dis_map$resno_B, 
                                    dis_map$insert_B, 
                                    dis_map$chain_B, sep="|"))
            outeleno <- cif$atom$eleno[data %in% query2]
        }

    ## If cutoff is 0 only the selection is returned -------------------------
    } else {
        outeleno <- refeleno
    }

    ## Trim the input pdb to prepare a smaller one ---------------------------
    pdb <- trim.pdb(cif, eleno=outeleno)
    if (rename) {
        if (any(nchar(pdb$atom$chain) > 1) | any(pdb$atom$chain == "?")) {
            query3 <- pdb$atom[as.character(outeleno), "chain"]
            Unique <- unique(query3)
            for (i in seq_along(Unique)) {
                pdb$atom$chain[query3 == Unique[i]] <- toupper(letters)[i]
                if (chain == Unique[i]) 
                    chain <- toupper(letters)[i]
            }
        }
    }

    ## Get just desired alt records ------------------------------------------
    if (alt[1] == "uniq") {
        alts <- alts <- sort(unique(pdb$atom$alt))
        alt <- alts[!alts == c(".")][1]
        eleno <- pdb$atom$eleno[pdb$atom$alt %in% c(".", alt)]
        sel <- atom.select(pdb, eleno=eleno)
        pdb <- trim.pdb(pdb, inds=sel)
    }

    ## Ensure that the output pdb has the correct format ---------------------
    tmp <- .perfect_output_format(pdb, outeleno, query, rename)
    resno2 <- tmp$resno2
    pdb <- tmp$pdb

    ## Save the output to a file if a name is specified ----------------------
    if (is.null(file)) {
        return(pdb)
    } else {
        if (length(grep("chain", file))>0) {
            file <- sub("chain.+_", paste("chain", chain, "_", sep=""), file)
        }
        if (length(grep("resno", file))>0 & exists("resno2")) {
            for (i in seq_along(resno2)) {
                file <- sub(paste("resno", resno[i], "\\.", sep=""), 
                paste("resno", resno2[i], ".", sep=""), file)
            }
        }
    tryCatch(
            {
                write.pdb(pdb, file=file, segid=pdb$atom$segid)
            }, error=function(e) {
                write.pdb(pdb, file=file, segid=pdb$atom$entid)
            })
    }
}
##############################################################################
## Subfunctions
## ===========================================================================
## Input to pdb
.input_to_pdb <-
function(cif, model=NULL, verbose, alt=NULL) {
    ## Check if input cif argument is a PDB ID or a "cif" object -------------
    if (length(class(cif) == 1) && class(cif) == "character") {
        ## If the input is a PDB ID, the data is downloaded
        if (nchar(cif) == 4) {
            if (verbose)
                print(cif)
            cif <- cifParser(cif)
        } else {
            stop("Your input string is not a pdb ID")
        }

    } else if (!.isCIF(cif) & !is.pdb(cif)) {
        stop(paste("Your input data is not a cif or pdb object, ",
                    "please refer to the cifParser or read.pdb functions",
                    sep=""))
    }

    ## Select model of interest ----------------------------------------------
    if (!is.null(model)) {
        if(.isCIF(cif)) {
            cif <- selectModel(cif=cif, model=model, verbose=verbose)
        } else if(is.pdb(cif)) {
            cif <- selectModel(pdb=cif, model=model, verbose=verbose)
        }
    }

    ## Make sure the object is a S3 pdb object -------------------------------
    if (.isCIF(cif)) {
        if(is.null(alt))
            alt <- unique(cifAtom_site(cif)$label_alt_id)

        cif <- cifAsPDB(cif, alt=alt)
    }
    return(cif)
}

.perfect_input_format <-
function(pdb) {
    if (any(is.na(pdb$atom$insert)) | 
            any(pdb$atom$insert == "") |
            any(pdb$atom$insert == " ")) {

        pdb$atom$insert[is.na(pdb$atom$insert)] <- "?"
        pdb$atom$insert[pdb$atom$insert == ""]  <- "?"
        pdb$atom$insert[pdb$atom$insert == " "] <- "?"
    }

    if (any(is.na(pdb$atom$chain))) {
        pdb$atom$chain[which(is.na(pdb$atom$chain))] <- "?"
    }

    if (any(is.na(pdb$atom$alt))) {
        pdb$atom$alt[is.na(pdb$atom$alt)] <- "."
    }

    if (any(is.na(pdb$atom$b))) {
        pdb$atom$b[which(is.na(pdb$atom$b))] <- 0
    }

    row.names(pdb$atom) <- pdb$atom$eleno
    return(pdb)
}

.perfect_output_format <-
function(pdb, outeleno, query, rename=TRUE) {
    if (any(is.na(pdb$atom$alt))) {
        pdb$atom$alt <- ""
    }

    if (any(pdb$atom$alt == ".")) {
        pdb$atom$alt <- ""
    }

    pdb$atom$charge <- ""
    pdb$atom$entid <- ""
    if (rename) {
        if (any(outeleno > 99999))
            pdb$atom$eleno <- seq_len(nrow(pdb$atom))
        if (any(pdb$atom$resno > 9999)) {
            query3 <- paste(pdb$atom[as.character(outeleno), "resno"],
                            pdb$atom[as.character(outeleno), "insert"],
                            pdb$atom[as.character(outeleno), "chain"],
                            sep="|")
            Unique <- unique(query3)
            resno2 <- c()
            for (i in seq_along(Unique)) {
                pdb$atom$resno[query3 == Unique[i]] <- i
                if (any(query == Unique[i])) resno2[query == Unique[i]] <- i
            }
        } else {
            resno2 <- NULL
        }
    } else {
        resno2 <- NULL
    }

    if (any(is.na(pdb$atom$insert))) { 
        pdb$atom$insert[which(is.na(pdb$atom$insert))] <- ""
    }

    if (any(pdb$atom$insert == "?")) {
        pdb$atom$insert[which(pdb$atom$insert == "?")] <- ""
    }

    return(list(pdb=pdb, resno2=resno2))
}
