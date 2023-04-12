## Methods for CIF objects

##############################################################################
## cif-accessors

#' @rdname cif_accessors
setMethod("cifEntry",
    signature="CIF",
    function(x) x@entry)

#' @rdname cif_accessors
setMethod("cifAudit_conform",
    signature="CIF",
    function(x) x@audit_conform)

#' @rdname cif_accessors
setMethod("cifDatabase_2",
    signature="CIF",
    function(x) x@database_2)

#' @rdname cif_accessors
setMethod("cifPdbx_database_status",
    signature="CIF",
    function(x) x@pdbx_database_status)

#' @rdname cif_accessors
setMethod("cifAudit_author",
    signature="CIF",
    function(x) x@audit_author)

#' @rdname cif_accessors
setMethod("cifEntity",
    signature="CIF",
    function(x) x@entity)

#' @rdname cif_accessors
setMethod("cifChem_comp",
    signature="CIF",
    function(x) x@chem_comp)

#' @rdname cif_accessors
setMethod("cifExptl",
    signature="CIF",
    function(x) x@exptl)

#' @rdname cif_accessors
setMethod("cifStruct",
    signature="CIF",
    function(x) x@struct)

#' @rdname cif_accessors
setMethod("cifStruct_keywords",
    signature="CIF",
    function(x) x@struct_keywords)

#' @rdname cif_accessors
setMethod("cifStruct_asym",
    signature="CIF",
    function(x) x@struct_asym)

#' @rdname cif_accessors
setMethod("cifAtom_sites",
    signature="CIF",
    function(x) x@atom_sites)

#' @rdname cif_accessors
setMethod("cifAtom_type",
    signature="CIF",
    function(x) x@atom_type)

#' @rdname cif_accessors
setMethod("cifAtom_site",
    signature="CIF",
    function(x) x@atom_site)

## End of section cif-accessors
##############################################################################

##############################################################################
## CIF S4 constructor

## cifParser 
#' @rdname cifParser
setMethod("cifParser",
    definition=function(pdbID, verbose=FALSE, cache=TRUE) {

        ## Print ID if verbose is TRUE
        if (verbose)
            print(pdbID)

        ## Retrieve structure from RAM if it's there
        if (cache) {
            ID <- toupper(substr(pdbID, 1, 4))
            if (exists(".cacheCIF", envir=.GlobalEnv) && 
                    cifEntry(get(".cacheCIF", envir=.GlobalEnv)) == ID) {

                return(get(".cacheCIF", envir=.GlobalEnv))            
            }
        }

        ## Read CIF block ----------------------------------------------------
        ## Save extension, in case its a file
        ext <- substr(pdbID, nchar(pdbID) - 3, nchar(pdbID))

        if (file.exists(pdbID) && (ext == ".cif" || ext == "f.gz")) { # Read
            pdb <- readLines(pdbID)

        } else if (nchar(pdbID) == 4) { # Otherwise download by pdb ID
            destfile <- cifDownload(pdbID=pdbID, verbose=verbose)
            pdb <- readLines(destfile)

        } else { # Otherwise it is just an error
            stop("Please, provide a valid pdbID or file")
        }
    
        ## Find #, they indicate the beggining/end of each section
        hash_inds <- grep("^#", pdb, perl=T)
        loop_inds <- which(pdb == "loop_")

        ## Define a list of indices for sections of interest
        sections <- lapply(cifAttr,
                            function(x) {
                                x  <- paste("^_", x, "\\.", sep="")
                                st <- grep(x, pdb, perl=TRUE)[1] - 1
                                if (st %in% loop_inds) {
                                    st = st - 1
                                }
                                return(which(hash_inds == st))
                            })

        ## Parse the CIF sections of interest
        cif <- lapply(sections,
                        FUN=.cifParser,
                        pdb=pdb, hash_inds=hash_inds)
        names(cif) <- cifAttr

        ## Create CIF S4 object and return output ----------------------------
        out <- CIF( entry                = cif$entry,
                    audit_conform        = cif$audit_conform,
                    database_2           = as.data.frame(cif$database_2),
                    pdbx_database_status = cif$pdbx_database_status,
                    audit_author         = as.data.frame(cif$audit_author),
                    entity               = as.data.frame(cif$entity),
                    chem_comp            = as.data.frame(cif$chem_comp),
                    exptl                = as.data.frame(cif$exptl),
                    struct               = cif$struct,
                    struct_keywords      = cif$struct_keywords,
                    struct_asym          = as.data.frame(cif$struct_asym),
                    atom_sites           = as.character(cif$atom_sites),
                    atom_type            = as.data.frame(cif$atom_type),
                    atom_site            = cif$atom_site)

        ## Save structure in RAM if cache is TRUE
        if (cache) {
            assign(".cacheCIF", out, envir=.GlobalEnv)
        }

        return(out)
    })

## End of section CIF S4 constructor
##############################################################################

##############################################################################

#' Download Protein Data Bank structures
#'
#' Given a 4-character string (PDB ID), download structure.
#'
#' @param pdbID A 4 character string that matches a structure in the Protein 
#'     Data Bank.
#' @param path The path to save the downloaded file. In basic use of the
#'     function, the "default" path is a temporal directory, which is removed
#'     after the R session is terminated. 
#'     Some users might prefer to save permanently the files in a local path.
#'     This advanced use of the function can be adquired creating (manually) a 
#'     local hidden directory called ".veriNA3d_mmCIF_files" with 
#'     `dir.create("~/.veriNA3d_mmCIF_files")`. This function will 
#'     automatically recognize the new "default" path and use it from there on.
#'     Note that the function will never overwrite existing files.
#' @param destfile File name to save the downloaded structure. If NULL, the 
#'     file name is constructed based on the PDB ID and extension.
#' @param extension A string with the desired file extension.
#' @param URL A string with the URL of use. If NULL, the RCSB is used as
#'     default.
#' @param verbose A logical to print process details.
#'
#' @return The file name.
#'
#' @examples
#'     ## cifDownload("1bau")
#'
#' @author Diego Gallego
#'
#' @rdname cifDownload
cifDownload <- 
function(pdbID, path="default", destfile=NULL, extension=".cif.gz", URL=NULL, 
            verbose=FALSE) {

    ## Make sure the input is lower case
    pdbID <- tolower(pdbID)

    ## If path is not provided, use temp directory
    if (path == "default") {
        if (dir.exists("~/.veriNA3d_mmCIF_files")) {
            path <- "~/.veriNA3d_mmCIF_files/"
        } else {
            path <- tempdir()
        }
    }

    ## If file name is not provided, and filename
    if (is.null(destfile)) {
        destfile <- paste(path, "/", pdbID, extension, sep="")
    }

    ## If file is already there but has size 0, remove it and download again
    if (file.exists(destfile) && file.info(destfile)$size == 0) {
        file.remove(destfile)
    }

    ## If file is not there, downlaod it
    if (!file.exists(destfile)) {
        if (verbose)
            cat("Downloading file from Internet\n")

        ## if URL is not provided, use the RCSB
        if (is.null(URL)) {
            URL <- paste("https://files.rcsb.org/download/",
            #URL <- paste("http://mmb.pcb.ub.es/api/pdb/", 
            ## For development tests I rather use the internal call
            #URL <- paste("http://web.mmb.pcb.ub.es/MMBApi/web/pdb/", 
                            pdbID, extension, sep ="")
        } else {
            URL <- paste(URL, pdbID, extension, sep ="")
        }

        ## Use lanchquery internal function to have error-handling
        .launchquery(URL, 
                        FUN=download.file,
                        destfile=destfile,
                        method="auto", 
                        quiet=!verbose)
    }
    return(destfile)
}
##############################################################################

##############################################################################
## Function to check if an object is CIF and related

# Is an Object of Class CIF?
#
# Checks whether an object is of Class CIF.
#
# @rdname .isCIF
#
# @param x An R object.
#
# @return A logical.
#
# @examples
# cif <- cifParser("1bau")
# .isCIF(cif)
#
# @author Diego Gallego
#
.isCIF <-
function(x) {
    inherits(x, "CIF")
}


# Is it a CIF object? Make it be!
#
# Internal function to check if an input is actually a CIF object
# If not, the cif file is read from the MMB API.
#
# @rdname cifMakeSure
#
# @param cif A cif object obtained from cifParser or a pdb ID so that the
#    function can download the data.
# @param verbose A logical indicating whether to print details of the process.
# @param check A string with the name of the function to use. It has been
#    thought to be used with '.isCIF' function.
#
# @return A cif object, which might be the same input or the downloaded data.
#
# @examples 
# cif <- veriNA3d:::.cifMakeSure("1bau")
#
# @author Diego Gallego
#
## .cifMakeSure
.cifMakeSure <-
function(cif, verbose=FALSE, check=".isCIF") {
    ## Check if input cif argument is a PDB ID -------------------------------
    if (length(class(cif) == 1) && class(cif) == "character") {

        ## If the input is a PDB ID or a file, cifParser is called -----------
        if (nchar(cif) == 4 || file.exists(cif)){
            if(verbose)
                print(cif)

            cif <- cifParser(cif)

        } else {
            stop("Your input string is not a pdb ID")
        }

    ## Check if it is a CIF --------------------------------------------------
    } else if( !do.call(check, list(cif)) ) {

        stop(paste(" Your input is not a 'CIF' object (i.e. from ",
                    "'cifParser') nor a pdb ID",
                    sep=""))
    }
    return(cif)
}
## End of section .isCIF 
##############################################################################

##############################################################################
## Method to coerce CIF S4 object to pdb S3 object as found in bio3d package

## cifAsPDB
#' @rdname cifAsPDB
setMethod("cifAsPDB",
    signature(cif="CIF"),
    definition=function(cif, model=NULL, chain=NULL, 
                        alt=c("A"), verbose=FALSE) {

        return(.cifAsPDB(cif, model=model, chain=chain, 
                            alt=alt, verbose=verbose))
    })

#' @rdname cifAsPDB
setMethod("cifAsPDB",
    signature(cif="character"),
    definition=function(cif, model=NULL, chain=NULL, 
                        alt=c("A"), verbose=FALSE) {
        ## Make sure it is a CIF ---------------------------------------------
        cif <- .cifMakeSure(cif)
        return(.cifAsPDB(cif, model=model, chain=chain,
                            alt=alt, verbose=verbose))
    })

## End of section cifAsPDB method
##############################################################################

##############################################################################
## Model selection methods

#' @rdname selectModel
setMethod("selectModel",
    signature(cif="CIF"),
    definition=function(cif, model, verbose=FALSE) {
        atom   <- cifAtom_site(cif)
        models <- unique(atom$pdbx_PDB_model_num)
        if (!model %in% models)
            stop("The model selected does not exist")
        cif@atom_site <- atom[atom$pdbx_PDB_model_num == model, ]
        return(cif)
    })

#' @rdname selectModel
setMethod("selectModel",
    definition=function(pdb, model, verbose=FALSE) {

        if (length(grep("trim", pdb$call)) > 0) {
            stop(paste(
                "Please, select the model you desire before applying other ",
                "functions", sep=""))
        }

        model <- as.numeric(model)
        if("model" %in% attributes(pdb)$names &&
                length(pdb$model) == 1 && 
                pdb$model == model) {

            if(verbose) 
                print(paste("The input is already the desired model, thus ",
                        "output = input", sep=""))
            return(pdb)
        }

        ## "flag" is an attirbute given by cifAsPDB. If TRUE, the pdb has
        ## models with different number of atoms, thus they are treated in a 
        ## special way
        if ("flag" %in% attributes(pdb)$names && pdb$flag) {
            pdb$atom <- pdb$model[[model]]
            pdb$flag <- FALSE
            pdb$xyz  <- as.xyz(matrix(c(t(pdb$atom[, c("x", "y", "z")])),
                                        nrow=1))
        } else {
            xyz <- pdb$xyz
            if (model > nrow(xyz)) {
                stop("The model selected does not exist")
            }
            pdb$xyz <- trim(xyz, row.inds=model)
            coords <- matrix(pdb$xyz, ncol=3, byrow=TRUE)
            pdb$atom[, "x"] <- coords[, 1]
            pdb$atom[, "y"] <- coords[, 2]
            pdb$atom[, "z"] <- coords[, 3]
        }

        pdb$model <- model
        return(pdb)
    })

## End of model selection methods
##############################################################################

##############################################################################
## rVector. Implemented to reproduce barnaba (Bottaro et al. NAR, 2014)

#' @rdname rVector
setMethod("rVector",
    signature(cif="CIF"),
    definition=function(cif, outformat="rvector", simple_out=TRUE) {
        pdb <- cifAsPDB(cif)
        out <- tryCatch({
                    .rVector(pdb1=pdb, outformat=outformat, 
                                simple_out=simple_out, elety="C4'")
                }, error=function(e) {
                    .rVector(pdb1=pdb, outformat=outformat, 
                                simple_out=simple_out, elety="C4")
                })
        return(out)
    })

#' @rdname rVector
setMethod("rVector",
    signature(cif="character"),
    definition=function(cif, outformat="rvector", simple_out=TRUE) {
        pdb <- cifAsPDB(cif)
        out <- tryCatch({
                    .rVector(pdb1=pdb, outformat=outformat, 
                                simple_out=simple_out, elety="C4'")
                }, error=function(e) {
                    .rVector(pdb1=pdb, outformat=outformat, 
                                simple_out=simple_out, elety="C4")
                })
        return(out)
    })

#' @rdname rVector
setMethod("rVector",
    definition=function(pdb, outformat="rvector", simple_out=TRUE) {
        out <- tryCatch({
                    .rVector(pdb1=pdb, outformat=outformat, 
                                simple_out=simple_out, elety="C4'")
                }, error=function(e) {
                    .rVector(pdb1=pdb, outformat=outformat, 
                                simple_out=simple_out, elety="C4")
                })
        return(out)
    })

## End of section for rVector methods
##############################################################################

##############################################################################
## eRMSD. Implemented to reproduce barnaba (Bottaro et al. NAR, 2014)

#' @rdname eRMSD
setMethod("eRMSD",
    signature(cif1="CIF", cif2="CIF"),
    definition=function(cif1=NULL, cif2=NULL) {

        rvectors1 <- rVector(cif1, outformat="rvector", simple_out=TRUE)
        rvectors2 <- rVector(cif2, outformat="rvector", simple_out=TRUE)

        len <- sqrt(nrow(rvectors1))
        deltaG <- t(apply(cbind(rvectors1[, c(1, 2, 3)],
                                rvectors2[, c(1, 2, 3)]),
                            MARGIN=1, FUN=.deltaGmodule))
        return(sqrt(sum(deltaG^2) / len))
    })

#' @rdname eRMSD
setMethod("eRMSD",
    definition=function(pdb1=NULL, pdb2=NULL,
                        rvectors1=NULL, rvectors2=NULL) {

        if (!is.null(rvectors1) && !is.null(rvectors2)) {
            if (!nrow(rvectors1) == nrow(rvectors2)) {
                stop("Different number of rvectors.", 
                        " The original PDB had a different length!", sep="")
            }
        } else if (!is.null(pdb1) && !is.null(pdb2)) {

            if (!sum(pdb1$atom$elety == "C4'") == 
                    sum(pdb2$atom$elety == "C4'")) {

                stop("Different lengths in input PDB objects")
            }
            rvectors1 <- rVector(pdb=pdb1, outformat="rvector",
                                            simple_out=TRUE)
            rvectors2 <- rVector(pdb=pdb2, outformat="rvector",
                                            simple_out=TRUE)
        } else {
            stop("Introduce two PDB objects or two set of rvectors")
        }

        len <- sqrt(nrow(rvectors1))
        deltaG <- t(apply(cbind(rvectors1[, c(1, 2, 3)], 
                                rvectors2[, c(1, 2, 3)]),
                            MARGIN=1, FUN=.deltaGmodule))
        return(sqrt(sum(deltaG^2) / len))
    })

## End of section for epsilon RMSD methods
##############################################################################

##############################################################################

#' @rdname RMSD
setMethod("RMSD",
    signature(cif1="CIF", cif2="CIF"),
    definition=function(cif1=NULL, cif2=NULL, sel1=NULL, sel2=NULL, ...) {

        pdb1 <- cifAsPDB(cif1)
        pdb2 <- cifAsPDB(cif2)

        if (is.null(sel1)) {
            sel1 <- atom.select(pdb1, ...)
        }
        if (is.null(sel2)) {
            sel2 <- atom.select(pdb2, ...)
        }

        fit <- fit.xyz(pdb1$xyz, pdb2$xyz,
                        fixed.inds=sel1$xyz, mobile.inds=sel2$xyz)

        return(rmsd(pdb1$xyz[,sel1$xyz],fit[,sel2$xyz]))

    })

#' @rdname RMSD
setMethod("RMSD",
    definition=function(cif1=NULL, cif2=NULL, sel1=NULL, sel2=NULL, ...) {

        pdb1 <- cif1
        pdb2 <- cif2

        if (is.null(sel1)) {
            sel1 <- atom.select(pdb1, ...)
        }
        if (is.null(sel2)) {
            sel2 <- atom.select(pdb2, ...)
        }

        fit <- fit.xyz(pdb1$xyz, pdb2$xyz,
                        fixed.inds=sel1$xyz, mobile.inds=sel2$xyz)

        return(rmsd(pdb1$xyz[,sel1$xyz],fit[,sel2$xyz]))

    })

## End of section for RMSD methods
##############################################################################
