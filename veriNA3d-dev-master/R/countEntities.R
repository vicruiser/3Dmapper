#' Count entities
#'
#' For a given pdbID, the function gets the Entity data and counts the number
#' of instances of the different entities (e.g. the number of different RNA,
#' the number of different proteins...).
#'
#' @param pdbID A 4-character string that matches a structure ID in the
#' Protein Data Bank.
#' @param force A logical to force the query instead of getting presaved data.
#'
#' @param ... Arguments to be passed to query function (see ?queryFunctions).
#'
#' @return A list with the number of instances of each entity.
#'
#' @examples
#' countEntities("1S72")
#'
#' @author Diego Gallego
#'
countEntities <- 
function(pdbID, force=FALSE, ...) {

    ## Check if presaved data is available or force the query
    pdbID <- toupper(pdbID)
    if (!force) {
        entities <- NULL
        data("entities", envir=environment())
        if (pdbID %in% entities$pdbID) {
            ind <- which(entities$pdbID == pdbID)
            MM <- as.integer(entities[ind, 2:ncol(entities)])
            names(MM) <- c("RNA", "DNA", "Hybrid", "PNA", "Prot", "Dprot", 
                            "Ligands", "Water", "Other")
            return(MM)
        }
    }

    ## Download info about entities, chains and length -----------------------
    Ent <- queryEntities(pdbID, ...=...)

    ## Solve corner cases (e.g. 2ICY) ----------------------------------------
    if (any(is.na(Ent$molecule_type))) {
        ind <- which(is.na(Ent$molecule_type))
        MM <- Ent[-ind, ]
    } else {
        MM <- Ent
    }

    ## Count number of each possible entity ----------------------------------
    RNA <- sum(MM$molecule_type == "polyribonucleotide")
    DNA <- sum(MM$molecule_type == "polydeoxyribonucleotide")
    Hybrid <- sum(MM$molecule_type == 
                    "polydeoxyribonucleotide/polyribonucleotide hybrid")
    PNA <- sum(MM$molecule_type == "peptide nucleic acid")
    Prot <- sum(MM$molecule_type %in% c("polypeptide(L)", 
                                        "cyclic-pseudo-peptide"))
    Dprot <- sum(MM$molecule_type == "polypeptide(D)")

    Ligands <- sum(MM$molecule_type == "Bound")
    Water <- sum(MM$molecule_type == "Water")
    Other <- sum(MM$molecule_type == "other")

    ## Prepare output --------------------------------------------------------
    out <- c(
                RNA=RNA,
                DNA=DNA,
                Hybrid=Hybrid,
                PNA=PNA,
                Prot=Prot,
                Dprot=Dprot,
                Ligands=Ligands,
                Water=Water,
                Other=Other)

    ## Double check if all entities were counted -----------------------------
    if(sum(unlist(out)) != nrow(Ent)) {
        warning(paste(pdbID, ": At least one entity not counted!", sep=""))
    }

    return(out)
}
