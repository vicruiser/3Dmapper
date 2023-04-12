#' Computes distances between the atoms of interest in a mmCIF structure
#' 
#' Given a pdb object (or a pdb ID), the function computes the distances 
#' between the desired atoms and returns the closest ones. Note that eleno 
#' numbers might be different in the PDB vs mmCIF formats and this may lead to
#' errors. 
#'
#' @param pdb A pdb object obtained as from 
#'     [cifAsPDB()] 
#'     or read.pdb/read.cif (bio3d functions).
#' @param model The model of interest to use in the calculations. The first
#'     model is always the default.
#' @param refeleno A vector of eleno (element number) to take as reference.
#' @param eleno A vector of eleno to measure the distances.
#' @param n An integer indicating how many closests atoms to return. The
#'     default n=1 returns only the closest atom; n=2 would return the two
#'     closest atoms and so on. If NULL, any number of atoms within the 
#'     cutoff will be returned.
#' @param cutoff A numeric vector indicating the distance range to consider in
#'     angstroms. Atoms further than the cutoff won't be returned.
#' @param verbose A logical indicating whether to print details of the process.
#' @param detailedoutput A logical indicating whether to include additional 
#'     information for each atom (see data_of_interest below). If FALSE, only
#'     the eleno (element number) and distances are returned.
#' @param data_of_interest A vector of strings. Only used if detailedoutput is 
#'     TRUE. The vector should only contain the strings between the following:
#'     "type", "elety", "alt", "resid", "chain", "resno", "insert", "x", "y",
#'     "z", "o", "b", "entid", "elesy", "charge", "asym_id", "seq_id", 
#'     "comp_id", "atom_id", "model".
#'     The selected fields will be returned for both atoms.
#'
#' @return A data.frame with the nearest atom neighbours information. Fields
#'     suffixed with '_A' refer to the atoms used as reference. Fields 
#'     suffixed with '_B' refer to the 'contacting'/closest atoms.
#'
#' @examples 
#'     ## Dowload cif file and save coordinates data
#'     cif <- cifParser("1enn")
#'     coordinates <- cifAtom_site(cif)
#'
#'     ## Find atom numbers for desired entities (e.g. water and DNA)
#'     water_eleno <- coordinates[coordinates$label_atom_id == "O", "id"]
#'     dna_eleno <- coordinates[coordinates$label_comp_id %in% 
#'                                         c("DA", "DT", "DG", "DU"), "id"]
#'
#'     ## Find which DNA atoms are in 5 Angstroms distance from the water
#'     data <- measureElenoDist(cif, refeleno=water_eleno, eleno=dna_eleno, 
#'                                 n=NULL, cutoff=5)
#'
#'     ## To see the data
#'     head(data)
#'
#' @author Diego Gallego
#'
measureElenoDist <- 
function(pdb, model=NULL, refeleno, eleno, n=1, cutoff=c(0, 5), verbose=FALSE,
            detailedoutput=TRUE, data_of_interest=NULL) {

    ## Make sure the object is a S3 pdb object -------------------------------
    if (any(class(pdb) == "CIF")) {
        alt <- unique(cifAtom_site(pdb)$label_alt_id)
        pdb <- cifAsPDB(pdb, alt=alt)
    }

    ## Select model of interest ----------------------------------------------
    if (!is.null(model))
        pdb <- selectModel(pdb, model, verbose=verbose)

    ## Save the coordinates of the desired atoms -----------------------------
    A_eleno <- refeleno
    B_eleno <- eleno

    A_eleno_ind <- which(pdb$atom$eleno %in% A_eleno)
    B_eleno_ind <- which(pdb$atom$eleno %in% B_eleno)

    A <- pdb$atom[A_eleno_ind, c("x", "y", "z")]
    B <- pdb$atom[B_eleno_ind, c("x", "y", "z")]
    if (!all(row.names(B) == B_eleno))
        row.names(B) <- B_eleno

    ## Make sure cutoff format is correct ------------------------------------
    if (is.null(cutoff) || cutoff < 0) {
        cutoff <- c(0, 5)
        warning(paste("If you don't select a positive cutoff, ",
                        "5A is set as default", sep=""))
    }

    if (length(cutoff) == 1) 
        cutoff <- c(0, cutoff)

    ## Estimate number of atoms to expect within the cutoff ------------------
    if (is.null(n)) {
        n <- cutoff[2] * cutoff[2] * 3
        n <- min(nrow(B), n)
    }

    ## Compute the distances -------------------------------------------------
    if (verbose) 
        print("Computing distances ...")

    dis_map <- nn2(query=A, 
                    data=B, 
                    searchtype="standard", 
                    radius=cutoff[2], 
                    k=n)

    df_map <- lapply(seq_len(nrow(A)),
                        FUN=function(i, A_eleno, dis_map, B) {
                            elenoA <- A_eleno[i]
                            indsB <- dis_map$nn.idx[i,]
                            elenoB <- as.integer(row.names(B[indsB,]))
                            distances <- dis_map$nn.dists[i,]
                            out <- cbind(elenoA=rep(elenoA, length(elenoB)),
                                            elenoB=elenoB, 
                                            distances=distances)
                            out <- out[which(out[, 3] >= cutoff[1] & 
                                                    out[, 3] <= cutoff[2]), ]
                            return(c(t(out)))
                        }, A_eleno=A_eleno, dis_map=dis_map, B=B)

    ## Coerce from list to data.frame ----------------------------------------
    out <- as.data.frame(matrix(unlist(df_map), ncol=3, byrow=TRUE),
                            stringsAsFactors=FALSE)
    names(out) <- c("eleno_A", "eleno_B", "distance")

    ## Round distances -------------------------------------------------------
    out$distance <- round(out$distance, 3)

    if (verbose) 
        print(" ... done")

    ## Should the output include detailed fields about the atoms involved? ---
    if (detailedoutput) {
        if (is.null(data_of_interest)) {
            data_of_interest <- c("elety", "resid", "resno",
                                    "chain", "insert", "alt", "b")
        }
        if (nrow(out) == 0) {
            out2 <- rep(NA, length(data_of_interest))
            names(out2) <- data_of_interest
            return(c(out, out2))
        }

        if (verbose) 
            print("Finding the atom details ...")

        row.names(pdb$atom) <- pdb$atom$eleno
        df_A <- pdb$atom[as.character(as.integer(out$eleno_A)),data_of_interest]
        df_B <- pdb$atom[as.character(as.integer(out$eleno_B)),data_of_interest]
        names(df_A) <- paste(data_of_interest, "_A", sep="")
        names(df_B) <- paste(data_of_interest, "_B", sep="")

        if (verbose) 
            print(" ... done, the output is coming")

        out <- cbind(out, df_A, df_B)
    }

    ## For pair-wise distances, keep only one combination
    ## If distance of atom 2 to 3 is included, 3 to 2 is not necessary
    if (length(eleno) == length(refeleno) && all(eleno == refeleno)) {
        ## Filter step
        m <- t(combn(eleno, 2))
        kk_m <- paste(m[,1], m[,2], sep="_")
        pastedins <- paste(out$eleno_A, out$eleno_B, sep="_")
        out <- out[which(pastedins %in% kk_m), ]
    }
    return(out)
}
