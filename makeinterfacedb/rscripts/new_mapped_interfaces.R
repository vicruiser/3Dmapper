###############################################################################
# Reformat predicted interfaces results
###############################################################################

#' Title
#'
#' @param pdb_id
#' @param interfaces_dir
#' @param blast_output
#'
#' @return
#' @export
#'
#' @examples
new.Mapped.Interfaces <-
  function(pdb_id , interfaces_dir, blast_output) {
    # read output files of all calculated distances (prot, lig and ac. nuc)
    # corresponding to the selected pdb_id
    # pi stands for "predicted interfaces"
    pi_paths   <- Sys.glob(file.path(
      interfaces_dir,
      paste("*", pdb_id, "*predicted_interfaces*", sep =
              "")
    ))
    pi_df_list <-
      lapply(pi_paths, function(x)
        fread(x, header = T))
    pi         <-
      as.data.frame(rbindlist(pi_df_list, use.names = TRUE, fill = TRUE))

    
    # subset ligand interactions and change chain id and use residue id instead
    if (nrow(pi) > 0) {
      pi_ligand <- subset(pi, interaction == "ligand")
      pi_ligand$chain.1 <- pi_ligand$resid.1
      
      # Check for fake nucleic interactions that actually are nucleotide ligand
      # and change the type of interaction from "nucleic" to "ligand".
      pi_nuc <- subset(pi, interaction == "nucleic")
      if (any(pi_nuc$chain.1 %in% pi$chain)) {
        nucleotide_ligand <- which(pi_nuc$chain.1 %in% pi$chain)
        pi_nuc[nucleotide_ligand, "interaction"] <- "ligand"
        pi_nuc$chain.1 <- pi_nuc$resid.1
      }
      
      # rearrange dataframe so column "chain" contains info of all possible
      # combination of all existing chains
      pi_protein <- subset(pi, interaction == "protein" )
      pi_struct <- subset(pi, is.na(interaction) )
      pi_struct = as.data.frame(unique(pi_struct%>%group_by(resno)%>%mutate(b = min(b), b.1 = min(b.1))))

      if (nrow(pi_protein) > 0) {
        pi_reversed <-
          pi_protein[, c(
            "type.1" ,
            "eleno.1" ,
            "elety.1" ,
            "chain.1" ,
            "resid.1" ,
            "resno.1" ,
            "type" ,
            "eleno" ,
            "elety" ,
            "chain" ,
            "resid" ,
            "resno" ,
            "distance" ,
            "b",
            "b.1",
            "interaction" ,
            "pdb.id" ,
            "real.pos.1",
            "real.pos"
          )]
        names(pi_reversed) <- c(
          "type" ,
          "eleno" ,
          "elety" ,
          "chain" ,
          "resid" ,
          "resno" ,
          "type.1" ,
          "eleno.1" ,
          "elety.1" ,
          "chain.1" ,
          "resid.1" ,
          "resno.1" ,
          "distance" ,
          "b.1",
          "b",
          "interaction" ,
          "pdb.id" ,
          "real.pos",
          "real.pos.1"
        )
        
        pi_all <- rbindlist(
          list(
            pi_protein,
            pi_ligand,
            as.data.frame(pi_reversed),
            pi_nuc
          ),
          use.names = T,
          fill = T
        )
      } else{
        pi_reversed <- data.frame(
          type.1 = character() ,
          eleno.1 = integer() ,
          elety.1 = character() ,
          chain.1 = character() ,
          resid.1 = character() ,
          resno.1 = integer(),
          type = character() ,
          eleno = integer() ,
          elety = character() ,
          chain = character(),
          resid = character() ,
          resno = integer() ,
          distance = numeric() ,
          b.1= numeric(),
          b= numeric(),
          interaction = character() ,
          pdb.id = character() ,
          real.pos.1 = integer(),
          real.pos = integer
        )
        pi_all <-    data.table(rbind(pi_protein,
                                      pi_ligand,
                                      pi_reversed,
                                      pi_nuc))
        
        
      }
      
      pi_all$chain.1 <- as.character(pi_all$chain.1)
      
      # Calculate the min distance between each interacting residue pair
      dist <- pi_all[, c("resno", "resno.1", "distance", "b", "b.1")]
      min_dist <-
        lapply(split(dist, by = c("resno", "resno.1")), function(x)
          min(x$distance))
      min_b<-
        lapply(split(dist, by = c("resno", "resno.1")), function(x)
          min(x$b))
      min_b.1<-
        lapply(split(dist, by = c("resno", "resno.1")), function(x)
          min(x$b.1))
      
      
      min_max_mean_dist <-
        data.frame(
          resno = as.numeric(sub("\\..*", "", names(min_dist))),
          resno.1 = as.numeric(sub(".*\\.", "", names(min_dist))),
          distance = unlist(min_dist),
          b =unlist(min_b),
          b.1= unlist(min_b.1)
        )
      rownames(min_max_mean_dist) <- NULL
      
      # unique by number of residues
      pi_all_unique <-
        unique(pi_all[, c(
          "type" ,
          "chain" ,
          "resid" ,
          "resno" ,
          "type.1" ,
          "chain.1" ,
          "resid.1" ,
          "resno.1" ,
          "interaction" ,
          "pdb.id",
          "real.pos"
        )])
      
      #add the distance calculation
      pi_all_unique <-
        merge(pi_all_unique,
              unique(min_max_mean_dist),
              by = c("resno", "resno.1"))
      # add struct info
      if (nrow(pi_all_unique) >0){
      pi_all_unique = rbind(pi_all_unique, unique(pi_struct[, c(
        "type" ,
        "chain" ,
        "resid" ,
        "resno" ,
        "type.1" ,
        "chain.1" ,
        "resid.1" ,
        "resno.1" ,
        "interaction" ,
        "pdb.id",
        "real.pos",
        "distance",
        "b",
        "b.1"
      )]))
      } else {
        pi_all_unique = pi_struct
      }
      # merge BLAST info and new mapped pi
      new_mpi <- merge(blast_output, pi_all_unique, by = "chain")

      ###############################################################################################
      # Eliminate those residues that do not map to the protein of reference (Ensembl)              #
      ###############################################################################################
      # REMINDER: in PDB_pairwise_interactions.R (step 3), the "real position", i.e.,
      # the one that starts from 1, termed "real.pos" of the chains' residues was included.
      # Thus, only those residues containing coordinates will be considered for this index
      # position.
      
      # 1) Check if real.pos of a residue on a predicted interface is included in
      #    the alignment, i.e., maps to an Ensembl protein sequence. Here we have to
      #    take into account that the indexes "qstart", "qend", "sstar" and "ssend" are
      #    the residues position corresponding to the begining and end of the alignment
      #    between the pdb chain sequence (query) and the ensembl protein sequence (subject)
      #    without taking into account the presence of gaps.
      
      new_mpi <-
        subset(new_mpi, real.pos >= qstart & real.pos <= qend)
      
      # 2) Get mapped position
      # split new_mpi dataframe by chain and prot.id since the same blast hit is repeated
      new_mpi_template_split <-
        split(new_mpi, new_mpi[, c('chain', 'prot.id', 'qstart', 'sstart')], drop =
                TRUE)
      
      # for every data frame previously splitted
      new_mpi_template_split <- lapply(new_mpi_template_split,
                                       function(x) {
                                         # Get the query (PDB) and the subject (Ensembl) aligned sequences and put
                                         # them into a data.frame (every row is a residue position, including gaps!)
                                         df <-
                                           data.table(
                                             resid_qseq = s2c(as.character(x[1, "qseq"])),
                                             resid_sseq = s2c(as.character(x[1, "sseq"]))
                                           )
                                         df <- data.table(df,
                                                          sgaps =  as.numeric(!is.gap(df$resid_sseq)),
                                                          qgaps = as.numeric(!is.gap(df$resid_qseq)))
                                         # Create two columns, one for resid_sseq and one for resid_qseq with identifiers for the
                                         # residues and the gaps. If gap == 0 , if residue == 1.
                                         # Then calculate cumulative sum of the identifiers to figure out the real position of every residue
                                         # in each of the sequences.
                                         df <-
                                           df %>% group_by(qgaps) %>% mutate(qpos = cumsum(qgaps) + unique(x$qstart) -
                                                                               1)
                                         df <-
                                           df %>% group_by(sgaps) %>% mutate(spos = cumsum(sgaps) + unique(x$sstart) -
                                                                               1)
                                         # map the predicted interfaces ("real.pos") to the ensembl sequence
                                         # check the real.pos that match the position in the query sequence (qpos). Then
                                         # extract the real position in the subject sequence = s_ali_pos
                                         #y <- df[match(x$real.pos, df$q_ali_pos),]
                                         y <- df[match(x$real.pos, df$qpos), ]
                                         # columns descriptions:
                                         #  resid_qseq: residues of aligned query sequence (i.e.: only aligned and includes gaps.)
                                         #  resid_sseq: residues of aligned subject sequence (i.e.: only aligned and includes gaps.)
                                         # qgaps : gaps identifiers of query sequence (0 = gap, 1= resiude)
                                         # sgaps : gaps identifiers of subject sequence (0 = gap, 1= resiude)
                                         # qpos : index position of each residue of the aligned query sequence starting from 1
                                         # spos : index position of each residue of the aligned subject sequence starting from 1
                                         # q_ali_pos: real index position of each residue of the aligned query sequence  (i.e, + qstart)
                                         # s_ali_pos: real index position of each residue of the aligned subject sequence (i.e, + sstart, i.e.= mapped real position!)
                                         x <- cbind(x, as.data.frame(y))
                                         return(x)
                                       })
      new_mpi <- rbindlist(new_mpi_template_split)
      
      new_mpi = subset(new_mpi, resid_qseq != '-' &
                         resid_sseq != '-')
      return(new_mpi)
    } else{
      stop("no pi")
    }
  }
