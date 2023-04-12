#' Classify RNA or DNA structures
#'
#' `classifyRNA` classifes a structure in the following groups:\cr
#' * {NoRNA}: the structure does not contain RNA or it is shorter
#' than a threshold set by "length".\cr
#' * {nakedRNA}: the only molecule(s) in the structure is RNA. It can
#' include ionic ligands, but not organic ligands.\cr
#' * {protRNA}: the PDB contains a protein-RNA complex.\cr
#' * {DprotRNA}: the PDB contains a protein-RNA complex and the 
#' protein has D aminoacids.\cr
#' * {DNARNA}: the PDB contains DNA-RNA.\cr
#' * {PNARNA}: the PDB contains PNA-RNA.\cr
#' * {ANARNA}: the PDB contains ANA-RNA.\cr
#' * {LNARNA}: the PDB contains LNA-RNA.\cr
#' * {ligandRNA}: the RNA is interacting with an organic 
#' ligand. Ions are not considered as ligands in this class.\cr
#'
#' `classifyDNA` classifes a structure in the following groups:\cr
#' * {NoDNA}: the structure does not contain DNA.\cr
#' * {nakedDNA}: the only molecule(s) in the structure is DNA. It can
#' include ionic ligands, but not organic ligands.\cr
#' * {protDNA}: the PDB contains a protein-DNA complex.\cr
#' * {DprotDNA}: the PDB contains a protein-DNA complex and the 
#' protein has D aminoacids.\cr
#' * {DNARNA}: the PDB contains DNA-RNA.\cr
#' * {PNADNA}: the PDB contains PNA-DNA.\cr
#' * {ligandDNA}: the DNA is interacting with an organic 
#' ligand. Ions are not considered as ligands in this class.\cr
#'
#' In `classifyRNA`, nucleic acid hybrids are considered RNA, while in
#' in `classifyDNA` they are considered DNA (e.g. pdb ID 2HVR).
#'
#' @param pdbID A 4-character string that matches a structure ID in the
#' Protein Data Bank.
#' @param length A positive integer. Minimum numer of nucleotides to consider
#' RNA as a polymer. An RNA shorter than this threshold is classified in
#' the NoRNA group.
#' @param force A logical to force the query instead of getting presaved data.
#' @param ... Arguments to be passed to query function (see ?queryFunctions).
#'
#' @return A string with the type of RNA.
#'
#' @examples
#' classifyRNA("1S72")
#'
#' @author Diego Gallego
#'
#' @name classifyNA
NULL
##############################################################################

#' @export
#' @rdname classifyNA
classifyRNA <-
function(pdbID, length=3, force=FALSE, ...) {

    pdbID <- toupper(pdbID)

    ## Check if the desired data is already presaved in the package ----------
    if (!force & (length <= 1 | length == 3)) {
        if (length <= 1) {
            fast <- .fast_check(pdbID, "RNAclassOver0")
        } else {
            fast <- .fast_check(pdbID, "RNAclassOver2")
        }
        if (fast[[1]]) 
            return(fast[[2]])
    }

    ## Check for some corner cases manually annotated ------------------------
    check <- .corner_cases(pdbID)
    if (check[[1]]) 
        return(check[[2]])

    ## Get entity data -------------------------------------------------------
    MM <- countEntities(pdbID, force=force, ...=...)
    MM <- as.data.frame(rbind(MM), stringsAsFactors=FALSE)

    ## If the PDB entry does not contain RNA it is classified as "NoRNA" -----
    if (MM$RNA + MM$Hybrid == 0)
        return("NoRNA")

    ## I a length >0 is provided, check that the RNA is longer than that -----
    if (length > 1) {
        MM2 <- queryEntities(pdbID, ...=...)
        ## Index for RNA in the data.frame -----------------------------------
        RNA_ind <- which(MM2$molecule_type %in% c("polyribonucleotide",
                        "polydeoxyribonucleotide/polyribonucleotide hybrid"))
    
        ## RNA that does not surpass a threshold is also classified as "NoRNA"
        if (all(MM2[RNA_ind, "length"] < length)) 
            return("NoRNA")
    }

    ## Logical, are there organic ligands? -----------------------------------
    ## Ions do not categorize a structure as ligandRNA since they are always
    ## in buffers
    if (MM$Ligands > 0) {
        ligands <- length(queryOrgLigands(pdbID, ...=...)) > 0
    } else {
        ligands <- FALSE
    }

    ## If there are proteins, the PDB entry is classified as "protRNA" -------
    if (MM$Prot > 0)
        return("protRNA")
                
    ## If there are D proteins, the PDB entry is classified as "DprotRNA" ----
    if (MM$Dprot > 0)
        return("DprotRNA")
                
    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (MM$DNA > 0)
        return("DNARNA")

    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (MM$PNA > 0)
        return("PNARNA")

    ## If there are ligands, the PDB entry is classified as "ligandRNA" ------
    if (ligands)
        return("ligandRNA")

    ## If the only molecule is RNA, then the PDB entry is classified as 
    ## "nakedRNA" ------------------------------------------------------------
    return("nakedRNA")
}
##############################################################################
#' @export
#' @rdname classifyNA
classifyDNA <-
function(pdbID, force=FALSE, ...) {

    pdbID <- toupper(pdbID)

    ## Check if the desired data is already presaved in the package ----------
    if (!force) {
        fast <- .fast_check(pdbID, "DNAclass")
        if (fast[[1]])
            return(fast[[2]])
    }

    ## Check for some corner cases manually annotated ------------------------
    check <- .corner_cases(pdbID)
    if (check[[1]]) 
        return("NoDNA")

    ## Get entity data -------------------------------------------------------
    MM <- countEntities(pdbID, force=force, ...=...)
    MM <- as.data.frame(rbind(MM), stringsAsFactors=FALSE)

    ## If the PDB entry does not contain RNA it is classified as "NoRNA" -----
    if (MM$DNA + MM$Hybrid == 0)
        return("NoDNA")

    ## Logical, are there organic ligands? -----------------------------------
    ## Ions do not categorize a structure as ligandDNA since they are always
    ## in buffers
    if (MM$Ligands > 0) {
        ligands <- length(queryOrgLigands(pdbID, ...=...)) > 0
    } else {
        ligands <- FALSE
    }

    ## If there are proteins, the PDB entry is classified as "protRNA" -------
    if (MM$Prot > 0)
        return("protDNA")

    ## If there are D proteins, the PDB entry is classified as "DprotRNA" ----
    if (MM$Dprot > 0)
        return("DprotDNA")

    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (MM$RNA > 0)
        return("DNARNA")

    ## If there are DNA molecules, the PDB entry is classified as "PNADNA" ---
    if (MM$PNA > 0)
        return("PNADNA")

    ## If there are ligands, the PDB entry is classified as "ligandRNA" ------
    if (ligands)
        return("ligandDNA")

    ## If the only molecule is DNA, then the PDB entry is classified as 
    ## "nakedRNA" ------------------------------------------------------------
    return("nakedDNA")
}

##############################################################################
## Subfunctions
## ===========================================================================

## Wrong or incomplete data in the API might generate a wrong classification,
## here I fix the detected ones
.corner_cases <-
function(pdbID) {
    pdbID <- toupper(pdbID)
    if (pdbID %in% c("2P7E",
                "3CR1")) {
    return(list(TRUE, "nakedRNA"))
    } else if (pdbID %in% c("3OK2",
                "3OK4")) {
    return(list(TRUE, "ANARNA"))
    } else if (pdbID %in% c("1HHW",
                "1HHX")) {
    return(list(TRUE, "LNARNA"))
    }
    return(list(FALSE, ""))
}

.fast_check <-
function(pdbID, info, verbose=FALSE) {
    fastquery <- NULL
    if (!exists("fastquery")) {
        data("fastquery", envir=environment())
    }
    pdbID <- toupper(pdbID)
    if (!info %in% names(fastquery)) {
        return(list(FALSE, ""))
    }
    if (pdbID %in% fastquery$pdbID) {
        if (verbose)
            print(paste("Presaved data for ", pdbID, " ", info, sep=""))
        ind <- which(fastquery$pdbID == pdbID)
        return(list(TRUE, fastquery[ind, info]))
    } else {
        return(list(FALSE, ""))
    }
}
