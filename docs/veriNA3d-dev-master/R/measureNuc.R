#' Obtain desired nucleotide measurements 
#' 
#' From a nucleic acid structure (pdb object), it computes the desired 
#' atomic distances, angles, dihedral angles, puckering conformation and Dp
#' distance (See definition of Dp in MolProbity paper by Chen et al. 
#' 2010).
#' 
#' @param pdb A pdb object as obtained from cifAsPDB or read.cif/read.pdb 
#'     (bio3d package).
#' @param model A string with the desired model number.
#' @param chain A string with the desired chain id.
#' @param v_shifted A logical. If TRUE, puckering angles (nu0 to nu4) are
#'     returned in the range 0 to 360 degrees. Otherwise, -180 to +180.
#' @param b_shifted A logical. If TRUE, backbone angles, chi and kappa are
#'     returned in the range 0 to 360 degrees. Otherwise, -180 to +180.
#' @param distances A data.frame indicating all the intra and inter-nucleotide
#'     atomic distances of interest. See details section. A default option is
#'     preconfigured to simplify the use of the function and can be seen 
#'     typing 'veriNA3d::.distances'.
#' @param angles A data.frame indicating all the intra and inter-nucleotide
#'     angles of interest. See details section. A default option is 
#'     preconfigured to simplify the use of the function and can be seen
#'     typing 'veriNA3d::.angles'.
#' @param torsionals A data.frame indicating all the intra and inter-
#'     nucleotide torsional angles of interest. See details section. A default
#'     option is preconfigured to simplify the use of the function and can be
#'     seen typing 'veriNA3d::.torsionals'.
#' @param pucker A logical indicating whether to compute the puckering.
#' @param Dp A logical indicating whether to compute the Dp distance.
#' @param refatm A string with the atom to use to identify the nucelotides.
#'     Important to analyse models with just (in example) phosphate atoms,
#'     in which refatm should be set to "P" (it was thought when analysing
#'     the structure with PDB code: 1Y1Y).
#' @param force A logical to force the analysis. Useful when the function 
#'     does not recognise a nucleic acid in the structure (e.g. because all
#'     bases are non-canonical: 1PBL, 1XV6, 1DV4, etc).
#'
#' @details The format of 'distances', 'angles' and 'torsionals' is:
#'     First column should indicate the first atom, second column second
#'     atom (and so on in the case of angles and torsional angles). An extra 
#'     last column is optional and should contain the names to identify each
#'     measurement in the output. Plane atom names are interpreted as intra-
#'     nucleotide measurments. For inter-nucleotide measurments use the prefix
#'     "pre_" or "post_" before the atom name. In example, to compute all 
#'     inter-phosphate distances, use as argument: \cr
#'     distances=data.frame(atomA=c("P"), atomB=c("post_P"), 
#'                             labels=c("interphosphate"), 
#'                             stringsAsFactors=FALSE).\cr
#'
#' @return A data.frame with the measurements for every nucleotide.
#'
#' @examples
#'     distances <- data.frame(atomA=c("P"), atomB=c("post_P"), 
#'                             labels=c("interphosphate"), 
#'                             stringsAsFactors=FALSE)
#'     measureNuc(cifAsPDB("1bna"), distances=distances, angles=NULL, 
#'                     torsionals=NULL, Dp=NULL)
#'
#' @author Diego Gallego 
#' 
measureNuc <-
function(pdb, model=1, chain="all", v_shifted=TRUE, b_shifted=TRUE, 
            distances="default", angles="default", torsionals="default", 
            pucker=TRUE, Dp=TRUE, refatm="C4'", force=FALSE) {

    ## Save desired model if necessary ---------------------------------------
    if (model == "all") {
        model <- seq_len(nrow(pdb$xyz))
    }
    ## Save desired chain if necessary ---------------------------------------
    if (chain == "all") {
        chain <- as.character(unique(pdb$atom$chain))
    }

    if (!force) {
        ## Check that the input has a nucleic acid ---------------------------
        resid <- unique(pdb$atom$resid[pdb$atom$chain %in% chain])
        if (!any(resid %in% .nucleotides)) {
            stop("Does the input pdb object contain a nucleic acid?")
        }
    }

    ## Make sure input is correct --------------------------------------------
    distances  <- .check_distances(distances)
    angles     <- .check_angles(angles)
    torsionals <- .check_torsionals(torsionals)

    ## Check if Dp dostance should be computed -------------------------------
    if (!is.null(Dp) && !is.na(Dp) && (Dp == TRUE | Dp == "default")) {
        Dp <- TRUE
    } else {
        Dp <- NULL
    }

    ## Check the user actually wants to compute something --------------------
    if (is.null(c(distances, angles, torsionals, Dp, pucker))) {
        stop("What do you want to measure? Arguments are NULL")
    }

    ## Make sure the pdb object has the necessary format ---------------------
    pdb <- .perfect_input_format(pdb)
    ## Find all combinations of models and chains to be computed -------------
    combinations <- expand.grid(model, chain, stringsAsFactors=FALSE)
    names(combinations) <- c("model", "chain")

    ## Time to measure -------------------------------------------------------
    ntinfo <- mapply(FUN=.measureNuc,
                        model=combinations[, "model"],
                        chain=combinations[, "chain"],
                        MoreArgs=list(pdb=pdb,
                                        v_shifted=v_shifted,
                                        b_shifted=b_shifted,
                                        distances=distances,
                                        angles=angles,
                                        torsionals=torsionals,
                                        pucker=pucker,
                                        Dp=Dp,
                                        refatm=refatm,
                                        force=force),
                        SIMPLIFY=FALSE)

    ## Give format to the output ---------------------------------------------
    ntinfo <- ntinfo[which(lapply(ntinfo, length)>0)]
    ntinfo <- do.call(rbind, ntinfo)
    ntinfo <- cbind(ntID=seq_len(nrow(ntinfo)), ntinfo)
    return(ntinfo)
}

###############################################################################
## Subfuncitons 
## ============================================================================
## Necessary internal functions to check input objects

.check_distances <-
function(distances) {
    ## Check user input for atomic distances ---------------------------------
    if (!is.null(distances) && !is.na(distances) && 
            class(distances) == "character" && distances == "default") {

        distances <- .distances

    ## Finish user data.frame if necessary -----------------------------------
    } else if ((is.matrix(distances) | is.data.frame(distances))) {

        if (ncol(distances) == 2 & 
                    length(grep("atom", colnames(distances))) == 2) {

            labels <- gsub("'", "p", paste("dist", distances[, 1],
                            distances[, 2], sep="."))
            distances <- as.data.frame(cbind(distances, labels), 
                            stringsAsFactors=FALSE)

        } else if (ncol(distances) > 3) {
            stop("Wrong format of input 'distances': too many columns")

        } else if (!all(colnames(distances) %in% 
                    c("atomA", "atomB", "labels"))) {
            stop("Wrong format of input 'distances': 
                    colnames should be 'atomA', 'atomB' and 'labels'")
        }

    ## Do not compute distances ----------------------------------------------
    } else {
        distances <- NULL
    }
    return(distances)
}

.check_angles <-
function(angles) {
    ## Check user input for angles -------------------------------------------
    if (!is.null(angles) && !is.na(angles) && 
            class(angles) == "character" && angles == "default") {

        angles <- .angles

    ## Finish user data.frame if necessary -----------------------------------
    } else if ((is.matrix(angles) | is.data.frame(angles))) {

        if (ncol(angles) == 3 &
            length(grep("atom", colnames(angles))) == 3) {

            labels <- gsub("'", "p", paste("angle", angles[, 1], angles[, 2],
                                angles[, 3], sep="."))
            angles <- as.data.frame(cbind(angles, labels), 
                                    stringsAsFactors=FALSE)
        } else if (ncol(angles) > 4) {
            stop("Wrong format of input 'angles': too many columns")

        } else if (!all(colnames(angles) %in% 
                                c("atomA", "atomB", "atomC", "labels"))) {
            stop("Wrong format of input 'angles': 
                colnames should be 'atomA', 'atomB', 'atomC' and 'labels'")
        }

    ## Do not compute angles -------------------------------------------------
    } else {
        angles <- NULL
    }
    return(angles)
}

.check_torsionals <-
function(torsionals) {
    ## Check user input for torsional angles ---------------------------------
    if (!is.null(torsionals) && !is.na(torsionals) && 
            class(torsionals) == "character" && torsionals == "default") {

        torsionals <- .torsionals

    ## Finish user data.frame if necessary -----------------------------------
    } else if ((is.matrix(torsionals) | is.data.frame(torsionals))) {

        if (ncol(torsionals) == 4) {
            labels <- gsub("'", "p", paste("tors", torsionals[, 1],
                    torsionals[, 2], torsionals[, 3], torsionals[, 4], 
            sep="."))
            torsionals <- as.data.frame(cbind(torsionals, labels),
                                        stringsAsFactors=FALSE) 
        } else if (ncol(torsionals) > 5) {
            stop("Wrong format of input 'torsionals': too many columns")

        } else if (!colnames(torsionals) %in%
                            c("atomA", "atomB", "atomC", "atomD", "labels")) {
            stop("Wrong format of input 'torsionals': colnames should be 
                    'atomA', 'atomB', 'atomC', 'atomD', and 'labels'")
        }

    ## Do not compute torsionals ---------------------------------------------
    } else {
        torsionals <- NULL
    }
    return(torsionals)
}

## ============================================================================
## Necessary internal data.frames for measure usage
.distances <- as.data.frame(cbind(
                    c("P",   "O5'", "C5'", "C4'", "C3'", "O3'",
                        "C1'", "C2'", "C3'", "O4'", "C1'", "C1'"),
                    c("O5'", "C5'", "C4'", "C3'", "O3'", "post_P",
                        "O4'", "C1'", "C2'", "C4'", "N9",  "N1"),
                    c("dist.P.O5p", "dist.O5p.C5p", "dist.C5p.C4p",
                        "dist.C4p.C3p", "dist.C3p.O3p", "dist.O3p.post_P",
                        "dist.C1p.O4p", "dist.C2p.C1p", "dist.C3p.C2p",
                        "dist.O4p.C4p", "dist.C1p.N9", "dist.C1p.N1")),
                stringsAsFactors=FALSE)
colnames(.distances) <- c("atomA", "atomB", "labels")

.angles <- as.data.frame(matrix(c(
                        "P",    "O5'",  "C5'",          "angle.P.O5p.C5p",
                        "O5'",  "C5'",  "C4'",          "angle.O5p.C5p.C4p",
                        "C5'",  "C4'",  "C3'",          "angle.C5p.C4p.C3p",
                        "C5'",  "C4'",  "O4'",          "angle.C5p.C4p.O4p",
                        "C4'",  "C3'",  "O3'",          "angle.C4p.C3p.O3p",
                        "C4'",  "C3'",  "C2'",          "angle.C4p.C3p.C2p",
                        "C4'",  "O4'",  "C1'",          "angle.C4p.O4p.C1p",
                        "C3'",  "O3'",  "post_P",       "angle.C3p.O3p.post_P",
                        "C3'",  "C2'",  "C1'",          "angle.C3p.C2p.C1p",
                        "C2'",  "C1'",  "O4'",          "angle.C2p.C1p.O4p",
                        "O3'",  "C3'",  "C2'",          "angle.O3p.C3p.C2p",
                        "C3'",  "C2'",  "O2'",          "angle.C3p.C2p.O2p",
                        "C1'",  "C2'",  "O2'",          "angle.C1p.C2p.O2p"), 
                        ncol=4, byrow=TRUE), 
                stringsAsFactors=FALSE)
colnames(.angles) <- c("atomA", "atomB", "atomC", "labels")

.torsionals <- as.data.frame(matrix(c(
                    "pre_O3'",  "P",     "O5'",      "C5'",         "alpha",
                    "P",        "O5'",   "C5'",      "C4'",         "beta",
                    "O5'",      "C5'",   "C4'",      "C3'",         "gamma",
                    "C5'",      "C4'",   "C3'",      "O3'",         "delta",
                    "C4'",      "C3'",   "O3'",      "post_P",      "epsilon",
                    "C3'",      "O3'",   "post_P",   "post_O5'",    "zeta",
                    "C4'",      "O4'",   "C1'",      "C2'",         "nu0",
                    "O4'",      "C1'",   "C2'",      "C3'",         "nu1",
                    "C1'",      "C2'",   "C3'",      "C4'",         "nu2",
                    "C2'",      "C3'",   "C4'",      "O4'",         "nu3",
                    "C3'",      "C4'",   "O4'",      "C1'",         "nu4",
                    "H2'",      "C2'",   "O2'",      "HO2'",        "kappa",
                    "pre_C4'",  "P",     "C4'",      "post_P",      "eta",
                    "P",        "C4'",   "post_P",   "post_C4'",    "theta",
                    "pre_C1'",  "P",     "C1'",      "post_P",      "eta_prime",
                    "P",        "C1'",   "post_P",   "post_C1'",    "theta_prime",
                    "pre_borg", "P",     "borg",     "post_P",      "eta_prime2",
                    "P",        "borg",  "post_P",   "post_borg",   "theta_prime2",
                    "O4'",      "C1'",   "N_base",   "C_base",      "chi"
                    ), ncol=5, byrow=TRUE),
                stringsAsFactors=FALSE)
colnames(.torsionals) <- c("atomA", "atomB", "atomC", "atomD", "labels")

## ============================================================================
## Intermediate wrapper that generates all the possible combinations of
## models and chains and calls the function to really make the measurments

.measureNuc <-
function(pdb, model, chain, v_shifted=TRUE, b_shifted=TRUE,
            distances=.distances[c(6, 11, 12), ], angles=.angles, 
            torsionals=.torsionals, pucker=TRUE, 
            Dp=TRUE, refatm="C4'", force=FALSE, select=TRUE) {

    if (select) {
        ## Selection of Model of interest ------------------------------------
        pdb <- selectModel(pdb=pdb, model=model, verbose=FALSE)
        if (any(grepl("eta_prime2", torsionals$labels))) {
            pdb <- .adddummy(pdb, refatm=refatm)
        }

        ## Selection of Chain of interest ------------------------------------
        selection <- atom.select(pdb, chain=chain)

        ## pdb contains the PDB object ONLY with the selected model and chain
        pdb <- trim(pdb, selection)
    }

    ## Check if pucker should be computed ------------------------------------
    if (!is.null(pucker) && !is.na(pucker) && 
            (pucker == TRUE | pucker == "default")) {

        add_torsionals <- .torsionals[grep("nu", 
                                .torsionals$labels, perl=TRUE),]

        ## If torsionals were not going to be computed, now they are
        if (is.null(torsionals)) {
            torsionals <- add_torsionals

        ## Check that puckering torsionals are included in the data.frame
        } else {
            pucker_tor <- apply(add_torsionals[, seq_len(5)], 1, 
                function(x) paste(x, collapse= "."))
            all_tor <- apply(torsionals[, seq_len(5)], 1, 
                                function(x) paste(x, collapse= "."))
            if (!all(pucker_tor %in% all_tor)) {
                torsionals <- rbind(torsionals, add_torsionals)
            }
        }

    ## Else do not compute puckering
    } else {
        pucker <- NULL
    }

    ## Make sure the chain selected is a nucleic acid ------------------------
    if (!any(.is.nucleic(pdb)) && !force) {
        return()
    }

    ## ridlist contains the sequence
    ## reslist contains the number of each nucleotide
    ## inslist contains the insertion code, necessary to differentiate some
    ## nucleotides that appear with the same number
    ridlist <- pdb$atom$resid[which(pdb$atom$elety == c(refatm))]
    reslist <- pdb$atom$resno[which(pdb$atom$elety == c(refatm))]
    inslist <- pdb$atom$insert[which(pdb$atom$elety == c(refatm))]

    total <- length(reslist)
    if (total == 0) {
        return()
    }
    indices <- seq_len(total)

    ## Call to do the maths for the given chain ------------------------------
    ntinfo <- lapply(indices, 
                        FUN=.new_measure,
                        reslist=reslist, 
                        inslist=inslist, 
                        ridlist=ridlist, 

                        pdb=pdb, 
            
                        distances=distances, 
                        angles=angles, 
                        torsionals=torsionals,

                        v_shifted=v_shifted,
                        b_shifted=b_shifted,
                        Dp=Dp,
                        pucker=pucker)

    ## Prepare the output ----------------------------------------------------
    ntinfo <- do.call(rbind, ntinfo)
    ntinfo <- cbind(model=rep(model, total),
                    chain=as.character(rep(chain, total)),
                    resid=as.character(ridlist), 
                    resno=as.character(reslist),
                    insert=as.character(inslist),
                    ntinfo, 
                    stringsAsFactors=FALSE)
    return(ntinfo)
}

## ============================================================================
## Function that actually does the maths
.new_measure <-
function(index, reslist, inslist, ridlist, pdb, 
            distances, angles, torsionals, v_shifted, b_shifted, Dp, pucker) {

    ## Extract info about residue number, insert records and 
    ## identifier (A, G, C, U, DA, DG, DC, DT...) ----------------------------
    resno <- reslist[index]
    insert <- inslist[index]
    resid <- ridlist[index]

    ## Extract info about residue number and insert records for the previous 
    ## and following nucleotides for inter-nucleotide measures (e.g.eta-theta)
    if (index == 1) {
        preresno <- ""
        postresno <- reslist[index + 1]
        preinsert <- ""
        postinsert <- inslist[index + 1]
    }else if (index == length(reslist)) {
        preresno <- reslist[index-1]
        postresno <- "" 
        preinsert <- inslist[index-1]
        postinsert <- ""
    } else {
        preresno <- reslist[index-1]
        postresno <- reslist[index + 1]
        preinsert <- inslist[index-1]
        postinsert <- inslist[index + 1]
    }

    ## Find base type --------------------------------------------------------
    if (resid %in% c("A", "G", "DA", "DG")) {
        base_type <- "pu"
    } else if (resid %in% c("C", "U", "DC", "DT")) {
        base_type <- "py"
    } else {
        ## For modified bases check the existence of the atom N9
        if (any(pdb$atom$resno == resno & 
                pdb$atom$insert == insert & 
                pdb$atom$elety == "N9")) {
            base_type <- "pu"
        } else {
            base_type <- "py"
        }
    }

    ## Since the input may contain some atom names called N_base and C_base,
    ## they will be replaced by N9/N1 and C4/C2 in function of base type -----
    if (base_type == "pu") {
        N_base <- "N9"
        C_base <- "C4"
    } else if (base_type == "py") {
        N_base <- "N1"
        C_base <- "C2"
    }

    if (!is.null(distances)) {
        distances[distances == "N_base"] <- N_base
        distances[distances == "C_base"] <- C_base
    }
    if (!is.null(angles)) {
        angles[angles == "N_base"] <- N_base
        angles[angles == "C_base"] <- C_base
    }
    if (!is.null(torsionals)) {
        torsionals[torsionals == "N_base"] <- N_base
        torsionals[torsionals == "C_base"] <- C_base
    }

    ##########################################################################
    ## This is an old but useful block comment. Now the selection is automated
    ##########################################################################
    ## Selection of all atoms required for torsions
    ## alpha   O3'(j-1) P(j) O5'(j) C5'(j)
    ## beta         P(j) O5'(j) C5'(j) C4'(j)
    ## gamma         O5'(j) C5'(j) C4'(j) C3'(j)
    ## delta            C5'(j) C4'(j) C3'(j) O3'(j)
    ## epsilon                 C4'(j) C3'(j) O3'(j) P(j + 1)
    ## zeta                       C3'(j) O3'(j) P(j + 1) O5'(j + 1)
    ## chi(Pu) C4(j) N9(j) C1'(j) O4'(j)
    ## chi(Py) C2(j) N1(j) C1'(j) O4'(j)
    ## v0      C4'(j) O4'(j) C1'(j) C2'(j)                  
    ## v1         O4'(j) C1'(j) C2'(j) C3'(j)
    ## v2            C1'(j) C2'(j) C3'(j) C4'(j)
    ## v3               C2'(j) C3'(j) C4'(j) O4'(j)
    ## v4                  C3'(j) C4'(j) O4'(j) C1'(j)
    ## kappa       H2'(j) C2'(j) O2'(j) HO2'(j) 
    ## Selection for pseudo torsions eta-theta
    ## eta     C4'(j-1) P(j) C4'(j) P(j + 1)
    ## theta        P(j) C4'(j) P(j + 1) C4'(j + 1)
    ##########################################################################
    ## End of block ##
    ##########################################################################

    ## Proces to know the atoms necessary for the measures -------------------
    distatoms <- unique(unlist(
                            distances[, grep("atom", names(distances))],
                            use.names=FALSE))
    angatoms <- unique(unlist(
                            angles[, grep("atom", names(angles))],
                            use.names=FALSE))
    toratoms <- unique(unlist(
                            torsionals[, grep("atom", names(torsionals))],
                            use.names=FALSE))

    ## still to work!!!
    if (!is.null(Dp) && Dp == TRUE) { 
        moreatoms <- c(N_base, "C1'", "post_P")
    } else {
        moreatoms <- NA
    }
    ##

    ## Generate vector with all necessary atom names that will be used for the
    ## calculations ----------------------------------------------------------
    atomlist <- sort(unique(c(distatoms, 
                                angatoms, 
                                toratoms, 
                                moreatoms)
                            )) #add moreatoms object in the c() function!!!!

    ## Generate vector with the name of the objects that will contain the atom
    ## selection -------------------------------------------------------------
    atomlist_sel <- paste( gsub("\'", "p", atomlist), "_sel", sep="")

    ## Generate vector with real elety names (remove prefix "pre" and "post")
    atomelety <- lapply(strsplit(atomlist, split="_"), 
                        function(x) { 
                            return(x[length(x)]) 
                        })

    ## Generate vector with apropiate strings to call the correct object, 
    ## in case there are atoms of the previous or following nucleotides to be
    ## used ------------------------------------------------------------------
    prefix <- vector("character", length(atomlist))
    prefix[grep("pre", atomlist)] <- "pre"
    prefix[grep("post", atomlist)] <- "post"

    ## Generate vectors indicating which resno and insert objects should be 
    ## used ("", "pre" or "post" nucleotide records) -------------------------
    resnocall <- paste( prefix, "resno", sep="")
    insertcall <- paste( prefix, "insert", sep="")

    ## Use the previous vectors to generate as many objects as atoms, 
    ## containing the selection, the atom and xyz indices, necessary to
    ## compute distances, angles and torsionals ------------------------------
    invisible(mapply(FUN=.select_many,
                        out_object=atomlist_sel,
                        elety=atomelety,
                        resno=resnocall,
                        insert=insertcall,
                        MoreArgs=list(pdb=pdb)))

    ## Measure interatomic distances -----------------------------------------
    if (!is.null(distances)) {
        distances$atomA <- paste( 
                                    gsub("\'", "p", distances$atomA), "_sel", 
                                    sep="")
        distances$atomB <- paste( 
                                    gsub("\'", "p", distances$atomB), "_sel", 
                                    sep="")

        distout <- apply(distances, 
                        MARGIN=1, 
                        FUN=function(i, pdb) {
                            atomAsel <- get(i[1], envir=parent.frame(2))
                            atomBsel <- get(i[2], envir=parent.frame(2))
                            label <- i[3]
                            substraction <- pdb$xyz[atomAsel$xyz] - 
                                            pdb$xyz[atomBsel$xyz]
                            if (length(substraction) == 0) {
                                output <- NA
                            } else {
                                output <- round(
                                sqrt(sum((substraction^2))), 3)
                            }
                            return(output)
                        }, pdb=pdb)
        names(distout) <- distances$labels
    } else {
        distout <- NULL
    }

    ## Measure angles --------------------------------------------------------
    if (!is.null(angles)) {
        angles$atomA <- paste( 
                                gsub("\'", "p", angles$atomA), "_sel", 
                                sep="")
        angles$atomB <- paste( 
                                gsub("\'", "p", angles$atomB), "_sel", 
                                sep="")
        angles$atomC <- paste( 
                                gsub("\'", "p", angles$atomC), "_sel", 
                                sep="")

        angles_list <- substring(angles$labels, 7)

        angles_sel <- c(paste(angles_list, "_sel", sep=""))

        invisible(apply(cbind(  angles_sel, 
                                angles$atomA,
                                angles$atomB,
                                angles$atomC),
                        MARGIN=1, 
                        FUN=.append_selections))

        invisible(mapply(   FUN=.angles_many, 
                            out_object=angles_list,
                            selection=angles_sel, 
                            MoreArgs=list(pdb=pdb)))

        invisible(lapply(   angles_list, 
                            function(ang) {
                                if (is.null(get(ang))) {
                                    assign(ang, value=NA, 
                                    envir=parent.frame(n=2))
                                }
                            }))

        names(angles_list) <- angles$labels
        angles_list <- lapply(angles_list, function(ang) get(ang))
    } else {
        angles_list <- NULL
    }

    ## Measure torsionals ----------------------------------------------------
    if (!is.null(torsionals)) {
        torsionals$atomA <- paste(
                                    gsub("\'", "p", torsionals$atomA), 
                                            "_sel", sep="")
        torsionals$atomB <- paste(
                                    gsub("\'", "p", torsionals$atomB),
                                            "_sel", sep="")
        torsionals$atomC <- paste(
                                    gsub("\'", "p", torsionals$atomC),
                                            "_sel", sep="")
        torsionals$atomD <- paste(
                                    gsub("\'", "p", torsionals$atomD), 
                                            "_sel", sep="")

        torsions_list <- paste(torsionals$labels, sep ="")
        torsions_sel <- c(paste(torsions_list, "_sel", sep=""))

        invisible(apply(
                        cbind(  torsions_sel,
                                torsionals$atomA,
                                torsionals$atomB,
                                torsionals$atomC,
                                torsionals$atomD), 
                        MARGIN=1, FUN=.append_selections))

        invisible(mapply(FUN=.torsions_many, 
                            out_object=torsions_list,
                            selection=torsions_sel,
                            MoreArgs=list(pdb=pdb)))

        if (!is.null(pucker) && pucker == TRUE) {
            puc <- .measure_pucker( get("nu0"), 
                                    get("nu1"),
                                    get("nu2"),
                                    get("nu3"),
                                    get("nu4"))
            pu_phase <- puc$pu_phase
            pu_amp <- puc$pu_amp
        } else {
            pu_phase <- NA
            pu_amp <- NA
        }

        ## Find which torsionals should be shifted from -180>x>180 to 0>x>360
        toshift <- NULL
        if (b_shifted) {
            toshift <- torsions_list[-grep(pattern="nu", torsions_list)]
        }
        if (v_shifted) {
            if(!is.null(toshift)) {
                toshift <- torsions_list
            } else {
                toshift <- torsions_list[grep(pattern="nu", torsions_list)]
            }
        }
        not_toshift <- torsions_list[!(torsions_list %in% toshift)]

        if (length(toshift) > 0) {
            invisible(lapply(toshift, 
                            FUN=function(tor) {
                                assign(tor, value=.shift360(get(tor)),
                                        envir=parent.frame(n=2))
                            }))
        }
        ## Make sure no angle is NULL
        if (length(not_toshift) > 0) {
            invisible(lapply(not_toshift, 
                                function(tor) {
                                    if (is.null(get(tor, 
                                                envir=parent.frame(n=2)))) { 

                                        assign(tor, value=NA, 
                                                envir=parent.frame(n=2))
                                    }
                                }))
        }

        names(torsions_list) <- torsionals$labels
        torsions_list <- lapply(torsions_list, function(tor) get(tor))
    } else {
        torsions_list <- NULL
    }

    ## Measure Richardson distance -------------------------------------------
    if (!is.null(Dp) && Dp == TRUE) {
        N_chi <- get(paste( N_base, "_sel", sep=""))
        post_P_sel <- get("post_P_sel")
        C1p_sel <- get("C1p_sel")

        if (length(pdb$xyz[post_P_sel$xyz]) !=0 &&
                length(pdb$xyz[C1p_sel$xyz]) !=0 &&
                exists("N_chi") && length(pdb$xyz[N_chi$xyz]) != 0) {

            glyc_vector <- pdb$xyz[N_chi$xyz] - pdb$xyz[C1p_sel$xyz]
            glyc_vector[glyc_vector == 0] <- 0.00000000001
            PLANE <- c(glyc_vector,
                        sum(glyc_vector * -pdb$xyz[post_P_sel$xyz]))
    
            line_eq1 <- c(glyc_vector[2], -glyc_vector[1], 0,
                                -pdb$xyz[C1p_sel$xyz][1] * glyc_vector[2] +
                                pdb$xyz[C1p_sel$xyz][2] * glyc_vector[1])
            line_eq2 <- c(0, glyc_vector[3], -glyc_vector[2],
                                -pdb$xyz[C1p_sel$xyz][2] * glyc_vector[3] +
                                pdb$xyz[C1p_sel$xyz][3] * glyc_vector[2])
        
            equation_system <- matrix(c(PLANE[c(1, 2, 3)], 
                                            line_eq1[c(1, 2, 3)],
                                            line_eq2[c(1, 2, 3)]),
                                        nrow=3, byrow=TRUE)
            solution_system <- c(-PLANE[4], -line_eq1[4], -line_eq2[4])
            point <- solve(equation_system, solution_system)
            Rich_distance <- round(sqrt(sum(
                                        (pdb$xyz[post_P_sel$xyz] - point)^2
                                    )), 3)
        } else {
            Rich_distance <- NA
        }
    } else {
        Rich_distance <- NA
    }

    ## Give format to the output ---------------------------------------------
    output <- list()
    if (!(length(distout) == 1 && is.null(distout))) {
        output <- append(output, distout)
    }
    if (!(length(angles_list) == 1 && is.null(angles_list))) {
        output <- append(output, angles_list)
    }
    if (!(length(torsions_list) == 1 && is.null(torsions_list))) {
        output <- append(output, torsions_list)
#        if (all(c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", 
#                    "chi", "eta", "theta") %in% names(torsions_list))) {
#            ind <- which(names(torsions_list) %in% c("alpha", "beta", "gamma",
#                                                        "delta", "epsilon", 
#                                                        "zeta", "chi", "eta", 
#                                                        "theta"))
#            torsions <- unlist(torsions_list[ind])
#            is.helical <- .is.helical(torsions)
#            output <- append(output, list(is.helical=is.helical))
#        }
    }

    if (!is.null(pucker) && pucker == TRUE) {
        output <- append(output, 
                            unlist(list(pu_amp=pu_amp, pu_phase=pu_phase)))
    }

    if (!is.null(Dp) && Dp == TRUE) {
        output <- append(output, unlist(list(Dp=Rich_distance)))
    }

    return(data.frame(output, row.names=NULL))
}

## ============================================================================
## Additional subfunctions to .new_measure
.select_many <- function(out_object, elety, pdb, resno, insert) {
    resno <- get(resno, envir=parent.frame(n=2))
    insert <- get(insert, envir=parent.frame(n=2))
    assign(out_object,
            value=atom.select(pdb, 
                                eleno=pdb$atom[which(
                                                pdb$atom$resno == resno &
                                                pdb$atom$insert == insert &
                                                pdb$atom$elety == elety), 
                                                "eleno"], 
                                verbose=FALSE),
            envir=parent.frame(n=2))
}

.torsions_many <- function(out_object, selection, pdb) {
    sel <- get(selection, envir=parent.frame(n=2))
    if (length(sel) == 12) {
        assign(out_object, round(torsion.xyz(pdb$xyz[sel]), 3),
        envir=parent.frame(n=2))
    } else {
        assign(out_object, NULL, envir=parent.frame(n=2))
    }
}

.angles_many <- function(out_object, selection, pdb) {
    sel <- get(selection, envir=parent.frame(n=2))
    if (length(sel) == 9) {
        assign(out_object, round(angle.xyz(pdb$xyz[sel]), 3),
        envir=parent.frame(n=2))
    } else {
        assign(out_object, NULL, envir=parent.frame(n=2))
    }
}

.append_selections <- function(selections) {
    output <- invisible(lapply(seq(2, length(selections), 1), FUN=function(i) {
        return(get(selections[i], envir=parent.frame(n=4))$xyz)
    }))
    assign(selections[1], value=unlist(output), envir=parent.frame(n=2))
}

.measure_pucker <-
function(nu0, nu1, nu2, nu3, nu4) {
    ## Make vector -----------------------------------------------------------
    pu_vec <- c(nu2, nu3, nu4, nu0, nu1)
    ## Compute pucker --------------------------------------------------------
        sumA <- 0
        sumB <- 0
        for (rt in seq_len(5)) {
                sumA <- sumA + (pu_vec[rt] * cos((4 / 5) * pi * (rt - 1)))
                sumB <- sumB + (pu_vec[rt] * sin((4 / 5) * pi * (rt - 1)))
        }
        A <- (2 / 5) * sumA
        B <- -(2 / 5) * sumB
        pu_amp <- round(sqrt((A^2) + (B^2)) ,3)
        pu_phase <- round(atan2(B, A) * (180 / pi), 3)

    ## Shift 360 -------------------------------------------------------------
    pu_amp <- .shift360(pu_amp)
    pu_phase <- .shift360(pu_phase)
    return(list(pu_phase=pu_phase, pu_amp=pu_amp))
}

## Function to shift 360 degrees torsion angles
.shift360 <-
function(tor) {
    if (!is.null(tor) && !is.na(tor) && length(tor) > 0) {
        if (tor < 0) {
            tor_shifted <- tor + 360
        } else {
            tor_shifted <- tor
        }
        return(tor_shifted)

    } else {
        return(NA)
    }
}

## ============================================================================
## Code addapted from bio3d
.is.nucleic <- function(pdb) {
    nuc.aa <- c("A",   "U",  "G",  "C",   "T",  "I",
                "DA", "DU", "DG", "DC",  "DT", "DI")
    return(pdb$atom$resid %in% nuc.aa)
}

## ============================================================================
## Function to find helical nucleotides
.is.helical <- function(torsions) {
    angles <- names(torsions)
    ind <- which(angles %in% c("alpha", "beta", "gamma", "delta", 
                                        #"epsilon", "zeta", "chi"))
                                        "epsilon", "zeta"))

    angl <- angles[ind]
    tors <- torsions[ind]

    if (all(is.na(tors))) {
        return(FALSE)
    }

    if (any(is.na(tors))) {
        ind <- which(!is.na(tors))
        angl <- angl[ind]
        tors <- tors[ind]
    }

    ref <- .helical_df[angl, "mean"]


    min_diffs <- min(abs(abs(tors - ref) - 360), abs(tors - ref))
    rmsd_tor <- sqrt(sum(min_diffs^2)/length(tors))

    return(rmsd_tor)

#    for (i in seq_along(torsions)) {
#        if (is.na(torsions[i])) {
#            next()
#        }
#        if (!(torsions[i] < .helical_df[.helical_df$angle == angles[i], "max"]
#            &torsions[i] > .helical_df[.helical_df$angle == angles[i], "min"])
#            ) {
#
#            return(FALSE)
#        }
#    }
#    return(TRUE)
}

.helical_df <- 
data.frame(
    angle=c("alpha", "beta", "gamma", "delta", "epsilon", 
            "zeta", "chi", "eta", "theta"),
    mean=c(292.033, 172.424, 56.673, 81.644, 208.973, 
            288.315, 199.776, 168.044, 216.321),
    sd=c(11.782, 9.873, 10.195, 5.394, 10.522, 
            9.311, 9.346, 7.067, 9.726),
    max=c(327.377, 202.044, 87.258, 97.825, 240.538, 316.248, 
            227.813, 189.244, 245.499),
    min=c(256.688, 142.805, 26.089, 65.463, 177.408, 260.382, 
            171.739, 146.844, 187.143), stringsAsFactors=FALSE, 
            row.names="angle")

.adddummy <- function(pdb, refatm="C4'", cores=1) {
    #pdb <- selectModel(pdb=pdb, model=1)
    pdb$atom$id <- getID(ntinfo=pdb$atom)$id_dssr
    inds <- which(pdb$atom$elety == refatm)
    atom <- pdb$atom

    data("refframes")
    #Abase <- read.pdb("../A_ref.pdb")
    #Cbase <- read.pdb("../C_ref.pdb")
    #Gbase <- read.pdb("../G_ref.pdb")
    #Tbase <- read.pdb("../T_ref.pdb")
    #Ubase <- read.pdb("../U_ref.pdb")
    #sos <- read.pdb("../sos.pdb")
    #Asel2 <- atom.select(Abase, 
    #elety=c("C1'", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"))
    #Csel2 <- atom.select(Cbase, 
    #elety=c("C1'", "N1", "C2", "N3", "C4", "C5", "C6"))
    #Gsel2 <- atom.select(Gbase, 
    #elety=c("C1'", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"))
    #Tsel2 <- atom.select(Tbase, 
    #elety=c("C1'", "N1", "C2", "N3", "C4", "C5", "C6", "C5M"))
    #Usel2 <- atom.select(Ubase, 
    #elety=c("C1'", "N1", "C2", "N3", "C4", "C5", "C6"))
    #ssel2 <- atom.select(sos, elety=c("C1'", "N"))
    #refframes <- list(Abase=Abase, Cbase=Cbase, Gbase=Gbase, Tbase=Tbase, 
    #                  Ubase=Ubase, sos=sos, Asel2=Asel2, Csel2=Csel2, 
    #                  Gsel2=Gsel2, Tsel2=Tsel2, Usel2=Usel2, ssel2=ssel2)
    #save(refframes, file="data/refframes.rda")
    vec <- c("C1'", "N1", "C2", "N3", "C4", "C5", 
              "C6", "N7", "C8", "N9")
    #atom2 <- list()
    #counter <- 1
    #for (i in inds) {
    atom2 <- .xlapply(inds, mc.cores=cores, mc.preschedule=TRUE,
                    function(i, atom, vec) {
        id <- pdb$atom$id[i]
        #print(id)
        atominds <-  which(pdb$atom$id == id & 
                            pdb$atom$elety %in% vec)
        resinds <- which(atom$id == id)
        if (length(atominds) > 1) {
            ele <- pdb$atom$eleno[atominds]
            sel1 <- atom.select(pdb, eleno=ele)
            sel1$xyz <- matrix(sel1$xyz, ncol=3, byrow=T)
            if (pdb$atom[atominds, "resid"][1] %in% c("A", "DA")) {
                base <- refframes$Abase
                sel2 <- refframes$Asel2
            } else if (pdb$atom[atominds, "resid"][1] %in% c("C", "DC")) {
                base <- refframes$Cbase
                sel2 <- refframes$Csel2
            } else if (pdb$atom[atominds, "resid"][1] %in% c("G", "DG")) {
                base <- refframes$Gbase
                sel2 <- refframes$Gsel2
            } else if (pdb$atom[atominds, "resid"][1] %in% c("T", "DT")) {
                base <- refframes$Tbase
                sel2 <- refframes$Tsel2
            } else if (pdb$atom[atominds, "resid"][1] %in% c("U", "DU")) {
                base <- refframes$Ubase
                sel2 <- refframes$Usel2
            } else if (all(c("N7", "C8", "N9") %in% pdb$atom$elety[sel1$atom])) {
                if ("N6" %in% atom$elety[which(atom$id == id)]) {
                    base <- refframes$Abase
                    sel2 <- refframes$Asel2
                } else {
                    base <- refframes$Gbase
                    sel2 <- refframes$Gsel2
                }
            } else {
                if (all(c("O4", "C5M") %in% atom$elety[which(atom$id == id)])) {
                    base <- refframes$Tbase
                    sel2 <- refframes$Tsel2
                } else if ("O4" %in% atom$elety[which(atom$id == id)]) {
                    base <- refframes$Ubase
                    sel2 <- refframes$Usel2
                } else {
                    base <- refframes$Cbase
                    sel2 <- refframes$Csel2
                }
            }
            bas <- base$atom$elety[base$atom$elety %in% vec]
            if (!all(bas %in% pdb$atom$elety[atominds])) {
                keep <- which(bas %in% pdb$atom$elety[atominds])
                bas <- bas[keep]
                sel2 <- atom.select(base, elety=pdb$atom$elety[atominds])
            }

            mor <- order(match(pdb$atom$elety[sel1$atom], bas))
            sel1$atom <- sel1$atom[mor]
            sel1$xyz <- c(t(sel1$xyz[mor, ]))

            if (length(sel1$atom) != length(sel2$atom)) {
                #atom2[[counter]] <- atom[resinds,]
                #counter <- counter + 1
                #next()
                return(atom[resinds,])
            }

            xyz <- fit.xyz(pdb$xyz, base$xyz, fixed.inds=sel1$xyz, mobile.inds=sel2$xyz)

            a <- matrix(pdb$xyz[sel1$xyz], ncol=3, byrow=T)
            #c <- matrix(base$xyz[sel2$xyz], ncol=3, byrow=T)
            b <- matrix(xyz[sel2$xyz], ncol=3, byrow=T)
            #sqrt(rowSums((a - b)^2))
            if (sum((a - b)^2) > 0.1) {
                #atom2[[counter]] <- atom[resinds,]
                #counter <- counter + 1
                #next()
                return(atom[resinds,])
            }

            row <- pdb$atom[i, ]
            row$elety <- "borg"
            row$elesy <- "D"
            if ("atom_id" %in% names(pdb)) {
                row$atom_id <- "borg"
            }
            #row$x <- round(mean(pdb$atom$x[atominds]), 3)
            #row$y <- round(mean(pdb$atom$y[atominds]), 3)
            #row$z <- round(mean(pdb$atom$z[atominds]), 3)
            row[, c("x", "y", "z")] <- xyz[(length(xyz)-2):length(xyz)]

            #latest <- resinds[length(resinds)]
            #atom <- rbind(atom[1:latest,], row, atom[(latest + 1):nrow(atom),])
            return(rbind(atom[resinds,], row))
            #atom2[[counter]] <- atom[resinds,]
            #atom2[[counter + 1]] <- row
            #counter <- counter + 2
        } else {
            #atom2[[counter]] <- atom[resinds,]
            #counter <- counter + 1
            #next()
            return(atom[resinds,])
        }
    }, atom=atom, vec=vec)
    if (length(atom2) == 0) {
        #cat("No borg points added!\n")
        return(pdb)
    }
    atom <- do.call(rbind, atom2)
    atom$eleno <- 1:nrow(atom)

    pdb$atom <- atom
    pdb$xyz <- as.xyz(matrix(as.numeric(c(t(atom[,c("x","y","z")]))), nrow=1))
    return(pdb)
}
