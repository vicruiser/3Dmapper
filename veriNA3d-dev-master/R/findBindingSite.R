#' Function to get data about the atoms in interacting site.
#'
#' For pdb structures with protein-nucleic acid complexes, the function finds
#' the atoms in the interacting site. It allows the user to set as reference
#' the nucleic acid, the protein, or particular desired chains.
#'
#' @param pdb A cif/pdb object obtained from cifParser/read.pdb respectively
#'     or a pdb ID so that the function can download the data.
#' @param cutoff A numeric to set the maximum distance for atoms to be 
#'     returned.
#' @param select A string that should match "Nuc", "Prot", "DNA" or "RNA", to
#'     be used as reference.
#' @param nchain A string with the nucleic acid chain to get data about.
#'     If NULL, all of them are selected (according with select argument).
#' @param pchain A string with the protein chain to get data about.
#'     If NULL, all of them are selected.
#' @param hydrogens A logical to use the hydrogens in the structure or remove 
#'     them.
#' @param byres A logical to indicate if the output should be referred to the 
#'     residues rather than atoms.
#' @param verbose A logical to print details of the process.
#' @param ... Arguments to selectModel and/or alt records.
#'
#' @return A data.frame with the atomic distances in the interacting site.
#'
#' @examples
#'     pdb <- cifParser("1b3t") # A protein-DNA complex
#'     data <- findBindingSite(pdb, select="DNA", byres=TRUE)
#'
#' @author Diego Gallego
#'
findBindingSite <-
function(pdb, cutoff=5, select="Nuc", nchain=NULL,
            pchain=NULL, hydrogens=FALSE, byres=FALSE, verbose=FALSE, ...) {

    ## Check input
    if (!select %in% c("Nuc", "RNA", "DNA", "Prot"))
        stop("Introduce a valid 'select' argument")

    ## Make sure the object is a S3 pdb object with the desired model --------
    pdb <- .input_to_pdb(cif=pdb, verbose=verbose, ...=...)

    ## Make sure the pdb object has the necessary format ---------------------
    pdb <- .perfect_input_format(pdb)

    ## Remove hydrogens from coordinates if ncessary -------------------------
    tmp_pdb <- .treatH(pdb, hydrogens, verbose)

    ## Which nucleotides should be taken into account? -----------------------
    nucleotides <- .select_nucleotides(select)

    ## If no nucleic chain is provided, all are used -------------------------
    if (is.null(nchain)) {
        nchain <- as.character(unique(tmp_pdb$atom[
                            tmp_pdb$atom$resid %in% nucleotides, "chain"]))
    }
    ## Double check
    resid <- unique(tmp_pdb$atom[tmp_pdb$atom$chain %in% nchain, "resid"])
    if (!any(resid %in% nucleotides)) {
        stop("Insufficent ", select, " atoms in structure. Chain:",
                    paste(nchain, collapse=", "), sep="")
    }

    ## Same for protein chains -----------------------------------------------
    if (is.null(pchain)) {
        pchain <- as.character(unique(
                        tmp_pdb$atom[tmp_pdb$atom$resid %in% .aa, "chain"]))
    }

    ## Select element numbers (eleno) ----------------------------------------
    ## Selection by entity molecule is more desirable, but not always possible
    if ("entid" %in% names(tmp_pdb$atom)) {
        entity <- TRUE
        ## Save nucleic and protein entity IDs
        nentid <- unique(tmp_pdb$atom[tmp_pdb$atom$chain %in% nchain & 
                            (tmp_pdb$atom$resid %in% nucleotides), "entid"])
        pentid <- unique(tmp_pdb$atom[tmp_pdb$atom$chain %in% pchain & 
                            (tmp_pdb$atom$resid %in% .aa), "entid"])

        ## Save nucleic and protein atoms
        neleno <- tmp_pdb$atom[tmp_pdb$atom$entid %in% nentid, "eleno"]
        peleno <- tmp_pdb$atom[tmp_pdb$atom$entid %in% pentid, "eleno"]

        ## Duble-check
        if (length(neleno) == 0 | length(peleno) == 0) {
            entity <- FALSE
        }
    } else {
        entity <- FALSE
    }

    if (!entity) {
        ## Save nucleic and protein indices
        ninds  <- atom.select(tmp_pdb, chain=nchain)
        pinds <- atom.select(tmp_pdb, string='protein', chain=pchain)

        ## Save nucleic and protein atoms
        neleno <- tmp_pdb$atom[ninds$atom, "eleno"]
        peleno <- tmp_pdb$atom[pinds$atom, "eleno"]
    }

    ## Check atoms were selected ---------------------------------------------
    if (length(neleno) == 0)
        stop("Insufficent 'nucleic' atoms in structure")
    if (length(peleno) == 0)
        stop("Insufficent 'protein' atoms in structure")

    ## According with desired analysis, set reference ------------------------
    if (select %in% c("RNA", "DNA", "Nuc")) {
        refeleno <- neleno
        eleno <- peleno
    } else if (select == "Prot") {
        refeleno <- peleno
        eleno <- neleno
    }

    ## Obtain data.frame of distances ----------------------------------------
    out <- measureElenoDist(pdb, 
                            refeleno=refeleno,
                            eleno=eleno,
                            cutoff=cutoff,
                            verbose=verbose,
                            data_of_interest=c("elety", "resid", "resno",
                                                "chain", "insert", "alt", 
                                                "b"))

    ## A different selection is made if the desired output is by residue -----
    if (byres) {
        residues <- paste(out$resno_A, out$insert_A, out$chain_A, sep=".")
        unique_res <- unique(residues)

        ## Find smallest distance for every residue
        inds <- lapply(unique_res,
                        FUN=function(x, residues, data) {
                            ind <- which(residues == x)
                            return(ind[which.min(data[ind, "distance"])])
                        }, residues=residues, data=out)
        out <- out[unlist(inds),]
    }

    return(out)
}
##############################################################################
## Internal objects 
## ===========================================================================
.nucleotides <- c("A", "G", "C", "U", "T", "I", 
                    "DA", "DG", "DC", "DT", "DU", "DI")
.Rnucleotides <- c("A", "G", "C", "U", "T", "I")
.Dnucleotides <- c("DA", "DG", "DC", "DT", "DU", "DI")
.aa <- c("GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "CYS",
            "PRO", "THR", "SER", "TYR", "GLN", "ASN",
            "ASP", "GLU", "HIS", "LYS", "ARG")

## ===========================================================================
## Subfunctions

.select_nucleotides <-
function(select) {
    if (select == "RNA") {
        return(.Rnucleotides)
    } else if (select == "DNA") {
        return(.Dnucleotides)
    } else if (select == "Nuc") {
        return(.nucleotides)
    }
}

.treatH <-
function(pdb, hydrogens, verbose) {
    if (!hydrogens) {
        pdb.inds <- atom.select(pdb, string="noh", verbose=verbose)
        if (length(pdb.inds$atom) == 0) {
            stop("Wrong hydrogen removal, try with hydrogens=TRUE")
        }
        return(trim.pdb(pdb, pdb.inds, verbose=verbose))
    } else {
        return(pdb)
    }
}

