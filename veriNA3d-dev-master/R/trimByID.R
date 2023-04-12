#' Calls trimSphere to generate smaller pdb files
#'
#' Given a data frame with nucleotide info (as obtained from pipeNucData)
#' and the desired nucleotide index (ntID), the function returns a pdb 
#' object or file allowing the user to select a number of 5' and 3' neighbors
#' in sequence and non-conected residues in a cutoff radius.
#'
#' @param cif A cif/pdb object obtained from cifParser/read.pdb respectively
#'     or a pdb ID so that the function can download the data. If NULL, the 
#'     function will extract the pdb ID from the ntinfo data frame 
#'     (pdbID col).
#' @param ntID An integer/string with the desired nucleotide ID for
#'     analysis.
#' @param ntinfo a data.frame with the data. It should contain at least the
#'     columns "pdbID", "chain", "model", "resno", "insert" and "ntID" (as the
#'     output of pipeNucData function).
#' @param prev Number of desired 5' neigbours to be returned.
#' @param post Number of desired 3' neigbours to be returned.
#' @param verbose A logical to print details of the process.
#' @param file A string to save the output in a pdb formated file. If NULL the
#'     functions just returns the pdb object.
#' @param justbackbone A logical to keep only the bakcbone of the output pdb.
#' @param ... Arguments to be passed to trimSphere (type ?trimSphere for 
#'     details).
#' 
#' @return A smaller pdb object or a pdb file. 
#'
#' @examples 
#'     cif <- cifParser("1bau")
#'     ntinfo <- pipeNucData(cif, torsionals=NULL, distances=NULL, angles=NULL)
#'
#'     ## Obtain a smaller pdb of the 4th nucleotide +-2 neigbours and a 
#'     ## sorrounding sphere of 5 Angstroms
#'     pdb <- trimByID(cif=cif, ntinfo=ntinfo, ntID=4, prev=2, post=2, 
#'                             cutoff=5)
#'
#'     ## Same process saving the output in a file:
#'     trimByID(cif=cif, ntinfo=ntinfo, ntID=4, prev=2, post=2, 
#'                             cutoff=5, file="output.pdb")
#'
#' @author Diego Gallego
#'
trimByID <-
function(cif=NULL, ntID, ntinfo, prev=2, post=2,
            verbose=TRUE, file=NULL, justbackbone=FALSE, ...) {

    id <- ntinfo$ntID == ntID
    #row.names(ntinfo) <- ntinfo$ntID
    if (is.null(cif)) 
        cif <- ntinfo[id, "pdbID"]
        #cif <- ntinfo[as.character(ntID), "pdbID"]

    desired <- .select_ntID_neighbours(ntID=ntID, ntinfo=ntinfo,
                                        prev=prev, post=post, 
                                        verbose=verbose)

    ntindex <- ntinfo[ntinfo$ntID %in% desired, "ntindex"]
    chain   <- ntinfo[id, "chain"]
    model   <- ntinfo[id, "model"]

    pdb <- trimSphere(cif=cif, ntindex=ntindex, model=model,
                        chain=chain, verbose=verbose, ...=...)

    if (justbackbone) {
        ## Select backbone atoms
        ind <- which(pdb$atom$elety == "C4'")
        rlist <- pdb$atom$resno[ind]
        ilist <- pdb$atom$insert[ind]

        ## For first residue of the chain, keep atoms from C4' onwards
        inds1 <- which(pdb$atom$resno == rlist[1] & 
                        pdb$atom$insert == ilist[1] &
                        pdb$atom$elety %in% c("C4'", "C3'", "O3'"))
        inds2 <- c()
        for (j in 2:(length(rlist) - 1)) {
            inds <- which(pdb$atom$resno == rlist[j] &
                            pdb$atom$insert == ilist[j] &
                            pdb$atom$elety %in% c("P", "O5'", "C5'", 
                                                    "C4'", "C3'", "O3'"))
            inds2 <- append(inds2, inds)
        }
        ## For the last residue of the chain, keep atoms from P to C4'
        inds3 <- which(pdb$atom$resno == rlist[length(rlist)] &
                        pdb$atom$insert == ilist[length(rlist)] &
                        pdb$atom$elety %in% c("P", "O5'", "C5'", "C4'"))

        ele <- pdb$atom$eleno[c(inds1, inds2, inds3)]
        sel = atom.select(pdb, eleno=ele)

        ## Trim pdb
        pdb <- trim.pdb(pdb, inds=sel)
    }

    if (!is.null(file)) {
        tryCatch(
            {
                write.pdb(pdb, file=file, segid=pdb$atom$segid)
            }, error=function(e) {
                write.pdb(pdb, file=file, segid=pdb$atom$entid)
            })
    } else {
        return(pdb)
    }
}
##############################################################################
## Subfunctions
## ===========================================================================
# Select neighboring nucleotides
#
# Given a ntID and a number of neighbors returns the ntIDs of the whole 
# polinulceotide. The interesting point of this function is that in the case
# of asking for too much neighbors, a vector containing as many NA will be
# returned, so the output vector will always have the desired length.
#
# @param ntID an obejct of class vector with the desired nucleotide of 
#    analysis.
# @param ntinfo a data.frame with the data. It should contain at least the
#    columns "pdbID", "chain", "model", "resno", "insert" and "ntID" (as the
#    output of [pipeNucData()] function.
# @param prev Number of desired 5' neigbours to be returned.
# @param post Number of desired 3' neigbours to be returned.
# @param info Column name of the desired data to be returned.
# @param verbose A logical to print details of the process.
#
# @return A vector with the desired data, extracted from the input data.frame
#
# @author Diego Gallego
#
.select_ntID_neighbours <-
function(ntID, ntinfo, prev=2, post=2,
            info="ntID", verbose=TRUE) {

    pdbID <- ntinfo[ntinfo$ntID == ntID, "pdbID"]
    chain <- ntinfo[ntinfo$ntID == ntID, "chain"]
    model <- ntinfo[ntinfo$ntID == ntID, "model"]

    resno_chain <- ntinfo[ntinfo$pdbID == pdbID &
                            ntinfo$chain == chain &
                            ntinfo$model == model, "resno"]
    length_chain <- length(resno_chain)
    insert_chain <- ntinfo[ntinfo$pdbID == pdbID &
                            ntinfo$chain == chain &
                            ntinfo$model == model, "insert"]
    chain_pos <- which(resno_chain == ntinfo[ntinfo$ntID == ntID, "resno"] &
                        insert_chain == ntinfo[ntinfo$ntID == ntID, "insert"])

    out.1 <- c()
    if (chain_pos-prev <= 0) {
        if (verbose) {
            print(paste("The nucleotide ", ntID,
                        " (ntID) doesn't have as many neighbours at",
                        " 5' as specified", sep=""))
        }
        while (chain_pos-prev <= 0) {
            out.1 <- append(out.1, NA)
            prev <- prev-1
        }
    }

    out.2 <- c()
    if (chain_pos+post > length_chain) {
        if (verbose)
            print(paste("The nucleotide ", ntID,
                        " (ntID) doesn't have as many neighbours at",
                        " 3' as specified", sep=""))
        while (chain_pos + post > length_chain) {
            out.2 <- append(NA, out.2)
            post <- post-1
        }
    }

    inds <- seq((ntID - prev), (ntID + post), 1)
    out <- c(out.1,
                ntinfo[ntinfo$ntID %in% inds, info],
                out.2)

    return(out)
}
