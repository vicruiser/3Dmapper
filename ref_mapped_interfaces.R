###################################################################
# Reformat Eduard'S results file to a R friendly format.          #
###################################################################

###################################################################
# Preliminary steps                                               #
###################################################################
  # Ed_mapped <- data.table::fread("pi_mapped_to_v89_5.txt")
  # head(Ed_mapped)
  # newcols <- as.data.frame(str_split_fixed(Ed_mapped$V2, "-", 5))
  # names(newcols)[1:4] <- c("sstart", "send", "V2", "chain.1")
  # Ed_mapped_reformatted<- unique(cbind(Ed_mapped[,1], newcols[,1:4], Ed_mapped[,3]))
  # newcols <- as.data.frame(str_split_fixed(Ed_mapped_reformatted$V2, "\\|", 2))
  # names(newcols) <- c("pdb.id.x", "chain")
  # Ed_mapped_reformatted2<- cbind(Ed_mapped_reformatted[,1:3],
  #                                newcols, Ed_mapped_reformatted[,5:6])
  # names(Ed_mapped_reformatted2)[1]<- "ensembl.prot.id"
  # names(Ed_mapped_reformatted2)[7]<- "resno"
  # write.table(Ed_mapped_reformatted2,
  #             "~/Dropbox/PhD/2018-2019/Eduard/pi_mapped_to_v89_5_reformatted.txt", quote=F)

###################################################################  
# Check wether all the pdbids are in Eduard's results or not      #
# and viceversa.                                                  #
################################################################### 

  # list_calc <- fread("~/Escritorio/PDB_CALC_DIST_finalRes/calc_dist_output_list.txt", header = F)
  # s<- unique(str_split_fixed(string = list_calc$V1,"-", 2)[,1])
  #
  # length(which(!(s %in% unique(Ed_mapped_reformatted2$pdb.id.x))))
  # length(which(!( unique(Ed_mapped_reformatted2$pdb.id.x) %in% s)))
  # Ed_mapped_reformatted2$pdb.id.x[which(!( unique(Ed_mapped_reformatted2$pdb.id.x) %in% s))]
  #
  # length(setdiff(unique(Ed_mapped_reformatted2$pdb.id.x), s))
  # length(setdiff( s, unique(Ed_mapped_reformatted2$pdb.id.x)))
  
  # read only subset of reference dataset (only pdb id of interest)
  # using awk and put it in the rigth format

###########################################################################
# Once we have the desired format we build the ref.Mapped.Interfaces().   #
# This function takes as input the current pdb structure we are analyzing #
# and scans the reference mapped interfaces database of reference in      #
# search of existing results related to it.                               #
###########################################################################

ref.Mapped.Interfaces <- function(pdb_id) {
  # define the command line to execute. 
  command_line = paste(
    "awk 'NR==1; $5==\"",
    pdb_id,
    "\"' /gpfs/home/bsc51/bsc51927/PDB_map_interfaces/ensembl_db/interfaces_mapped_to_v89_5_reformatted_without_interdomains.txt",
    #"\"' /home/vruizser/PhD/2018-2019/Eduard/interfaces_mapped_to_v89_5_reformatted.txt",
    sep = ""
  )
  
  # retrieved data from ref database
  ref_mpi <- tryCatch(
    fread(
      text = system(command_line,
                    intern = TRUE),
      sep = " ",
      drop=1
    ),
    error = function(e)
      data.frame()
  )
  
  # Return a data frame (ref_mpi) independently of whether is 
  # empty or not, i.e., empty if no results were retrieved or
  # the corresponding results.
  if (nrow(ref_mpi) > 0) {
    ref_mpi <- rbindlist(lapply(1:nrow(ref_mpi), function(i)
      data.frame(ref_mpi[i, 1:6],
                 unlist(strsplit(
                   as.character(ref_mpi$resno[i]), "-"
                 )), row.names = NULL)))
    
    ref_mpi <- as.data.frame(ref_mpi, stringsAsFactors = FALSE)
    
    colnames(ref_mpi)[7] <- "mapped.real.pos"
    
    ref_mpi$mapped.real.pos <-
      as.integer(levels(ref_mpi$mapped.real.pos))[ref_mpi$mapped.real.pos]
    
    ref_mpi$id <- "ref"
    
  } else {
    
    ref_mpi <- data.frame(
      chain = character(),
      pdb.id.x =  character(),
      ensembl.prot.id = character(),
      sstart = integer(),
      send = integer(),
      chain.1 = character(),
      mapped.real.pos = character(),
      id = character()
    )
  }
  
  return(ref_mpi)
}