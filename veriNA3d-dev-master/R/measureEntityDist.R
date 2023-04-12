#' Computes distances between all the atoms of selected entities in a mmCIF 
#' structure
#' 
#' Given a cif object (or a pdb ID), the function computes the distances 
#' between atoms of the selected entity IDs. For each atom/residue of the 
#' reference entity the function returns the closest atoms of the other 
#' entities. This function is a wrapper of
#' [measureElenoDist()] 
#' and simplifies its use. If you are unfamiliar with the concept of 
#' entities in a mmCIF structure see example below.
#'
#' @param cif A cif object obtained from cifParser or a pdb ID.
#' @param model The model of interest to use in the calculations. The first 
#'     model is always the default.
#' @param refent A string with the entity ID of reference. The distance output
#'     will be referred to the atoms/residues of this entity.
#' @param entities A character vector with the entities of interest. The 
#'     default "all" will select all of them except the refent.
#' @param ... Additional arguments to be passed to
#'     [measureElenoDist()].
#'
#' @return A data.frame with the nearest atoms neighbour information.
#'
#' @examples 
#'     ## To see the entities of a given structure use:
#'     cif <- cifParser("1enn")
#'     cifEntity(cif)
#'
#'     ## Supose you are interested on the interactions of water and DNA
#'     water_entity <- 5
#'     dna_entity <- 1
#'     
#'     ## Find which DNA atoms are in 5 Angstroms distance from the water
#'     data <- measureEntityDist(cif, refent=water_entity, 
#'                 entities=dna_entity, n=10, cutoff=5)
#'     ## An equivalent run without downloading the cif file previously
#'     data <- measureEntityDist("1enn", refent=water_entity, 
#'                 entities=dna_entity, n=10, cutoff=5)
#'
#'     ## This option is better than using the example in ?measureElenoDist,
#'     ## since this way it would also take into account modified residues, if
#'     ## any.
#'
#' @author Diego Gallego
#'

measureEntityDist <- 
function(cif, model=NULL, refent, entities=c("all"), ...) { 

    ## Check if input cif argument is a PDB ID or a "cif" object -------------
    cif <- .cifMakeSure(cif)

    ## Select model of interest ----------------------------------------------
    if (!is.null(model)) {
        cif <- selectModel(cif, model)
    }

    ## If the reference entity is to be compared with all the rest, they must
    ## be specified ----------------------------------------------------------
    if (length(entities) == 1 && entities == "all") {
        ent_data <- cifEntity(cif)
        entities <- ent_data$id[-which(ent_data$id == as.character(refent))]
    }

    ## Now coerce the CIF S4 object to a pdb S3 object -----------------------
    alt <- unique(cifAtom_site(cif)$label_alt_id)
    pdb <- cifAsPDB(cif, alt=alt)

    ## Find element numbers (eleno) ------------------------------------------
    refent_ind <- which(pdb$atom$entid == refent)
    entities_ind <- which(pdb$atom$entid %in% entities)

    eleno <- pdb$atom[, "eleno"]
    A_eleno <- eleno[refent_ind]
    B_eleno <- pdb$atom[entities_ind, "eleno"]

    ## Use the eleno numbers to call the measureElenoDist function ------------
    return(measureElenoDist(
        pdb=pdb,
        refeleno=A_eleno,
        eleno=B_eleno,
        ...=...))
}
