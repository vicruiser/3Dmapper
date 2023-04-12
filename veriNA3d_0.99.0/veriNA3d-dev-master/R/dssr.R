#' Dissecting the Spatial Structure of RNA with DSSR
#'
#' Wrapper function to execute DSSR (see reference below) on a DNA or RNA 
#' structure and parse the result.
#'
#' @param pdb It can be: 
#'     \itemize{
#'         \item A 4 character string corresponding to a PDB ID.
#'         \item A pdb/mmcif file.
#'         \item A pdb object as provided by [cifAsPDB()] or 
#'             [bio3d::read.pdb()].
#'     }
#'
#' @param exefile A string with the program name.
#' @param dssrargs A vector of strings with the desired arguments to feed DSSR.
#' @param verbose A logical indicating whether to print details of the process.
#'
#' @return A list with the json output of DSSR.
#'
#' @examples
#'     # Not run
#'     # dssr_1bau <- dssr("1bau")
#'
#' @author Diego Gallego
#'
#' @references 
#'     Lu, X.J., H.J. Bussemaker, and W.K. Olson. 2015. “DSSR: An Integrated 
#'     Software Tool for Dissection the Spatial Structure of RNA.” Nucleic 
#'     Acids Research 43 (21): e142
#'
dssr <- 
function(pdb, exefile="x3dna-dssr",
            dssrargs=c("--nmr", "--torsion360", "--more"), verbose=FALSE) {

    ## Check if the program is executable ------------------------------------
    os1 <- .Platform$OS.type
    status <- system(paste(exefile, "--version"),
                        ignore.stderr=TRUE, ignore.stdout=TRUE)

    if (status != 0) {
        stop(paste("x3dna-dssr execution failed\n",
                    " Make sure '", exefile, "' is in your search path.",
                    "Alternatively, fill the 'exefile' argument with the",
                    " whole path to the executable", sep=""))
    }

    ## Check if user provides file or pdb ID ---------------------------------
    if (length(class(pdb)) == 1 && class(pdb) == "character") {

        ## When file exists, infile is defined as input pdb
        if (file.exists(pdb)) {
            infile <- pdb
            #flag <- FALSE

        ## Otherwise, download structure if possible
        } else if (nchar(pdb) == 4) {
            #flag <- TRUE
            url <- "http://www.ebi.ac.uk/pdbe/entry-files/download/"
            infile <- cifDownload(pdbID=pdb, extension=".cif", 
                                    URL=url, verbose=verbose)

        ## If the string is not a file or ID, stop
        } else {
            stop("Provide a valid pdb")
        }

    } else {
        ## Write file from pdb object
        if ("pdb" %in% class(pdb)) {
            infile <- tempfile()
            #flag <- TRUE
            tryCatch(
                    { 
                        write.pdb(pdb, file=infile, segid=pdb$atom$segid)
                    }, error=function(e) {
                        write.pdb(pdb, file=infile, segid=pdb$atom$entid)
                    })

        ## If the input is not a S3 bio3d pdb object
        } else {
            stop("Provide a valid pdb")
        }
    }

    ## Create output file name
    outfile <- tempfile(fileext=".json")
    ## Generate command instruction
    ## --auxfile=no avoids the generation of additional files
    ## --json ensures an easy parsing of the json formated data
    dssrargs <- paste(" --auxfile=no --json ", dssrargs, collapse="", sep="")
    cmd <- paste(exefile, dssrargs,
                    " -i=", infile, " -o=", outfile, sep="")

    ## Execute command -------------------------------------------------------
    if (os1 == "windows") {
        FUN <- match.fun(shell)
    } else {
        FUN <- match.fun(system)
    }
    success <- FUN(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)

    ## Check whether the command was executed properly or not
    if(success!=0) {
        stop(paste("An error occurred while running command\n '",
                    cmd, "'", sep=""))
    }
    
    ## Remove input file if necessary and leave things as they were before
    #if (flag) {
    #    file.remove(infile)
    #}

    ## Parse the result, remove outfile and return output
    out <- fromJSON(outfile)
    file.remove(outfile)

    return(out)
}
