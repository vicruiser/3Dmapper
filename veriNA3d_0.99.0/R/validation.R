#' Get validation data from official EBI-API
#' 
#' Get validation data for a given PDB on a nucleotide basis, based on existing
#' information on EBI-PDBe REST-API. Validation performed on the following 
#' parameters: clashes, suite outliers, pucker outliers, bond lengths, bond 
#' angles, chirals and rsrz
#'
#' @param pdbID A string containing the desired PDB ID
#' @param ntinfo Optional. A data.frame obtained from 
#'     [pipeNucData()].
#' @param model A string/integer with the model of interest.
#' @param force A logical to avoid returning errors.
#'
#' @return A data.frame with TRUE/FALSE data for the different validation 
#' parameters for each nucleotide, reported  with unique identifier.
#'
#' @examples 
#'     pdbID <- "1bau"
#'     vdata <- validation(pdbID)
#'
#' @author Eric Matamoros & Diego Gallego

validation <- function(pdbID, ntinfo=NULL, model="all", force=FALSE) {

    #Check for the entities that belong to RNA
    ent_query <- queryEntities(pdbID)
    ala <- c("polyribonucleotide", "polydeoxyribonucleotide", 
                "polydeoxyribonucleotide/polyribonucleotide hybrid")
    check_rib <- which(ent_query$molecule_type %in% ala)

    if (length(check_rib) == 0) {
        if (force) {
            return(NULL)
        } else {
            stop(paste("Are you sure ", pdbID, " has a nucleic acid?", sep=""))
        }
    }

    ## Make sure to have a reference data.frame to work with
    if (is.null(ntinfo)) {
        ntinfo <- .extract_ntinfo(pdbID=pdbID, model=model)
    }
    table3 <- getID(pdbID, ntinfo=ntinfo)
    seq <- table3$seq
    id_dssr <- table3$id_dssr
    total <- paste(ntinfo$resno, ntinfo$chain, sep=".")
    ## Modified id_dssr
    resid <- ntinfo$resid
    ntinfo$resid <- "nan"
    id_dssr_mod <- getID(ntinfo=ntinfo)$id_dssr
  
    #-------------------------------------------------------------------------
    ## Check for validation data of suite_outliers and pucker_outliers
    IDsummary <- queryAPI(ID=pdbID, API="ebi", string1=
                            "validation/RNA_pucker_suite_outliers/entry/", 
                            string2="")
    ## suite_outliers
    suite_out <- IDsummary[[1]]$suite_outliers
    if (nrow(suite_out) == 0 || is.null(suite_out) || length(suite_out) == 0) {
        table3$suite_outlier <- NA
    } else {
        suite_out <- .manage_it(suite_out)
        id_suite <- getID(ntinfo=suite_out)$id_dssr
    
        table3$suite_outlier <- id_dssr_mod %in% id_suite
    }
    ## pucker_outliers
    pucker_out <- IDsummary[[1]]$pucker_outliers
    if (nrow(pucker_out) == 0 || is.null(pucker_out) || length(pucker_out) == 0) {
        table3$pucker_outlier <- NA
    } else {
        pucker_out <- .manage_it(pucker_out)
        id_puc <- getID(ntinfo=pucker_out)$id_dssr
    
        table3$pucker_outlier <- id_dssr_mod %in% id_puc
    }

    #-------------------------------------------------------------------------
    #Check for validation data of geometrical outliers
    String1 <- "validation/protein-RNA-DNA-geometry-outlier-residues/entry/"
    IDsummary_geom <- queryAPI(ID=pdbID, API="ebi",
                               string1=String1, string2 = "")
    
    ## Initiate columns
    table3$clashes <- NA
    table3$bond_lengths <- NA
    table3$bond_angles <- NA
    table3$chirals <- NA
    table3$planes <- NA

    outlier_types <- names(IDsummary_geom[[1]]$molecules$chains[[1]]$models[[1]]$outlier_types)

    if (!(is.null(outlier_types) | length(outlier_types) == 0)) {
        ## Fill with FALSE the columns that will actually have data
        for (m in outlier_types) {
            table3[, m] <- FALSE
        }

        #Exctract the different types of outliers
        ## Loop over entities
        for (i in 1:length(check_rib)){
            # pick entity data
            entid <- check_rib[i]
            ent <- IDsummary_geom[[1]]$molecules[IDsummary_geom[[1]]$molecules$entity_id == entid, "chains"]
            # loop over chains in the given entity
            if(length(ent) != 0){
                for (row_ch in 1:nrow(ent[[1]])) {
                    # pick data of given chain
                    ent_ch <- ent[[1]][row_ch, ]
                    # save chain
                    chain <- ent_ch$chain_id
                    for (m in outlier_types){
                        ## Append data of different models (for NMR)
                        te <- ent_ch$models[[1]]$outlier_types[, m]
                        te <- lapply(1:length(te), function(x) {return(cbind(model=x, te[[x]]))})
                        te <- do.call(rbind, te)
                        te$chain <- chain
                        te <- .manage_it(te)
                        id_local <- getID(ntinfo=te)$id_dssr

                        # save data about clashes ... or whatever
                        table3[, m][which(id_dssr_mod %in% id_local)] <- TRUE
                    }
                }
            }
        }
    }
    #--------------------------------------------------------------------------
    
    #Check for validation data of rsrz
    IDsummary_rsrz <- queryAPI(ID = pdbID, API  = "ebi",
                               string1 = "/validation/outliers/all/", string2 = "")
    
    table3$rsrz <- FALSE
    
    rsrz <- IDsummary_rsrz[[1]]$rsrz$outliers
    rsrz_out <- IDsummary_rsrz[[1]]$rsrz$outliers$units
    #strsplit(as.character(rsrz_out), "\\|")
    
    find_rsrz <- IDsummary_rsrz[[1]]$types_of_outliers$outliers$types
    look <- which(grepl("rsrz", find_rsrz, perl=T))
    
    if (is.na(look[1]) == TRUE){
      #print("no rsrz") #nothing, means that there is no rsrz value in there
    }else{
        #Segmentate strings and get the chain_id, residue and resno
        maps <- list()
        for (i in rsrz_out){
            maps[i] <- list(strsplit(as.character(i), "\\|")[[1]][c(3,4, 5, 8)])
        }
        
        map2 <- as.data.frame(maps)
        
        
        #Creation of a data frame with the single nucleotides, remove any AA fragment
        c <- c("A", "C", "G", "U")
        
        m <- list()
        #l[2,90] %in% c
        
        for(i in 1:length(map2)){
            m[i] <- as.character(map2[2, i]) 
        }
        
        # Creation of map_ntd (contains units) ans value_ntd(contains values) only for the nucleotides. Remove all proteic residues
        map_ntd <- map2[which(m %in% c)]
        value_ntd <- rsrz$value[which(m %in% c)]
        
        #Extract the column that contains the insert values for each nucleotide
        map_insert <- list()
        
        if(length(map_ntd) != 0){
            for(i in 1:length(map_ntd)){
                map_insert[i] <- as.character(droplevels(map_ntd[4, i ]))
            }
            
            #Extrct insert indexes from the vector
            index7 <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", map_insert, perl=T))
            
            #Determine length of the whole nucleotide list and crete a list
            l <- length(as.list(droplevels(map_ntd[3,])))
            l <- list(1:l)
            l <- as.list(l[[1]])
            index7 <- as.list(index7)
            
            #Determine with TRUE the index which are not in l and with FALSE the ones which are in both lists
            `%notin%` <- Negate(`%in%`)
            clar <- l %notin% index7
            
            map_ntd_insert <- list()
            
            #Run through the index and copy the insert when needed
            for (i in 1:length(l)){
                if(clar[i] == TRUE){
                    map_ntd_insert[i] <- as.character(droplevels(map_ntd[3, i]))
                }else{
                    map_ntd_insert[i] <- paste(as.character(droplevels(map_ntd[3, i])), as.character(droplevels(map_ntd[4, i])), sep="^")
                }
            }
            
            # Create an empty list and operate over map_ntd to pase the index of the chain_id with the index of the resno
            chain_resno <- list()
            
            #Create the resno.chain identifier for the rsrz
            for (i in 1:length(l)){
                chain_resno[i] <- paste(as.character(map_ntd_insert[i]), as.character(droplevels(map_ntd[1, i])), sep = ".")
            }
            
            table3$rsrz <- FALSE
            
            #Creation of the column and set values to TRUE when there is an rsrz
            chain_resno <- as.character(chain_resno)
            table3$rsrz[which(total %in% chain_resno)] <- TRUE
            
            #----------------------------------------------------
            
            #Now we are going to create another row with the value of rsrz
            li <- list()
            table3$value_rsrz <- FALSE
            
            #Create a list that appends the reno and the value_ntd (which is the value of the rsrz)
            for (i in 1:length(chain_resno)){
                li[[i]] <- append(chain_resno[i], value_ntd[i])
            }
            
            #Change those chain_resno which are set to FALSE to their correct value. Leave the others as FALSE
            for(i in 1:length(li)){
                #check[[i]] <- which(total %in% li[[i]][1])
                table3$value_rsrz[which(total %in% li[[i]][1])] <- li[[i]][2]
        
            }
        }
        return(table3)
        #---------------------------------------------------------------------
    }
    return(table3)
}
##############################################################################
## Subfunctions
## ===========================================================================

#' Extract the unique identifier for each of the nucleotides
#' 
#' Generate the ID for every nucleotide according with dssr format: 
#' model:chain.residresno^insert.
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param ntinfo A data.frame containing the information retreived from the
#'     PipeNucData function from the veriNa3d package containing unique
#'     parameters for each nucleotide
#'     
#'
#' @return A string with the unique identifier for each nucleotide.
#'
#' @examples 
#'     pdblist <- list("1bau", "2rn1")
#'     ntinfo <- .extract_ntinfo(pdbID=pdblist)
#'     id_dssr <- id_dssr_extract(ntinfo)
#'     ntinfo$id_dssr <- id_dssr #Both id_dssr and ntinfo have the same length
#'
#' @author Eric Matamoros & Diego Gallego
#'
getID <- function(pdbID=NULL, ntinfo=NULL){

    if(is.null(ntinfo)){
        ntinfo <- .extract_ntinfo(pdbID=pdbID)
    }
    #Obtain sequence
    seq <- ntinfo$resid
    
    # Creation of resid_resno
    resid_resno <- paste(ntinfo$resid, ntinfo$resno, sep="")
    
    # Creation of chain.residresno
    chain_resid_resno <- paste(ntinfo$chain, resid_resno, sep=".")
    
    #Creation of id_dssr <- model:chain.residresno
    id_dssr <- paste(ntinfo$model, chain_resid_resno, sep=":")
    
    inds <- which(ntinfo$insert != "?")
    if (length(inds) != 0) {
      id_dssr[inds] <- paste(id_dssr[inds], ntinfo$insert[inds], sep='^')
    }

    data_table <- data.frame(seq=seq, id_dssr=id_dssr, stringsAsFactors=FALSE)
    return(data_table)
}

#' Obtain nucleotide details from a data set of RNA structures
#' 
#' Pipeline to generate a data.frame with the desired info for a list of PDB. 
#' Nucleotides are labeled with a unique identifier (column ntID).
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param ... Arguments to be passed to the checkNuc function
#'
#' @return A data.frame with data about every nucleotide in the input set.
#'
#' @examples 
#'     pdblist <- list("1bau", "2rn1")
#'     ntinfo <- extract_ntinfo(pdbID=pdblist)
#'
#' @author Eric Matamoros & Diego Gallego
#'
.extract_ntinfo <- function(pdbID, ...){
    #Download data from Protein Data Bank or retreive it from file
    dir <- "/orozco/projects/PDB.datamining/veriNA/"
    if (dir.exists(dir)) {
        txt_files2 <- dir(dir, pattern = "ntinfo.txt", full.names = TRUE)
        txt_files <- dir(dir, pattern = "ntinfo.txt")

        #Subset all PDB ID from the file names "PDB"_ntinfo.txt
        pdbs <-  strsplit(txt_files, "_")
        only_pdb <- list()

        for (i in 1:length(pdbs)){
            only_pdb[i] <- pdbs[[i]][1]
        }

        only_pdb <- toupper(as.character(only_pdb))

        if (any(toupper(pdbID) %in% only_pdb)){
            a <- which(only_pdb %in% toupper(pdbID))
            ntinfo <- lapply(txt_files2[a], function(x) {
                return(read.table(x, header=TRUE, stringsAsFactors=FALSE))})
            ntinfo <- do.call(rbind, ntinfo)
            ntinfo$ntID <- seq(1, nrow(ntinfo), 1)
        }
        return(ntinfo)
    }

    ntinfo <- lapply(pdbID, function(x) {return(checkNuc(cifAsPDB(x), ...=...))})
    ntinfo <- do.call(rbind, ntinfo)
    ntinfo$ntID <- seq(1, nrow(ntinfo), 1)

    return(ntinfo)
}

.manage_it <- function(df) {
    df$resid <- "nan"
    names(df) <- gsub("author_residue_number", "resno", names(df))
    names(df) <- gsub("chain_id", "chain", names(df))
    names(df) <- gsub("model_id", "model", names(df))
    names(df) <- gsub("author_insertion_code", "insert", names(df))
    df$insert[df$insert == ""] <- "?"
    return(df)
}
