###############################################################################
# REFORMAT PFAM DOMAINS DATABASE                                              #
###############################################################################

# nf <-
#   max(count.fields(
#     "~/Dropbox/PhD/2018-2019/Eduard/iur_and_pfam_ensembl_v85.txt"
#   ))
# pfamscan <-
#   read.table(
#     "~/Dropbox/PhD/2018-2019/Eduard/iur_and_pfam_ensembl_v85.txt",
#     fill = TRUE,
#     col.names = 1:nf
#   )
# pfamscan_melted <-
#   melt(pfamscan,
#        id = "X1",
#        value.name = "pfam.dom",
#        na.rm = TRUE)
# pfamscan_melted2 <- pfamscan_melted[-2]
# pfamscan_melted3 <-
#   pfamscan_melted2[!(is.na(pfamscan_melted2$pfam.dom) |
#                        pfamscan_melted2$pfam.dom == ""),]
# sub_pfamscan_melted3 <-
#   subset(pfamscan_melted3, X1 == "ENSP00000251020")
# newcols <-
#   as.data.frame(str_split_fixed(pfamscan_melted3$pfam.dom, "-", 3))
# names(newcols) <- c("pfam.dom", "dom.start", "dom.end")
# pfamscan_reformatted <- cbind(pfamscan_melted3[-2], newcols)
# names(pfamscan_reformatted)[1] <- "ensembl.prot.id"
# write.table(
#   pfamscan_reformatted,
#   "~/Dropbox/PhD/2018-2019/Eduard/iur_and_pfam_ensembl_v85_reformatted.txt",
#   quote = F
# )

###############################################################################
# Associate PFAM domains with each Ensembl ID if possible                     #
###############################################################################

#pfamscan_reformatted <- read.table("~/PhD/2018-2019/Eduard/iur_and_pfam_ensembl_v85_reformatted.txt")
#biomart <- fread("~/PhD/2018-2019/PDB_CALC_DIST_finalRes/ensembl94_gene_prot_pfam_id.txt")

Mapped.Interfaces.Biomart <- function(new_mpi) {
  # read pfam and biomart databases
  #pfamscan_reformatted <-
  #  read.table(
  #    "/gpfs/home/bsc51/bsc51927/PDB_map_interfaces/ensembl_db/iur_and_pfam_ensembl_v85_reformatted.txt"
      #"~/PhD/2018-2019/Eduard/iur_and_pfam_ensembl_v85_reformatted.txt"
  #  )
  biomart <-
    fread(
      "/gpfs/home/bsc51/bsc51927/PDB_map_interfaces/ensembl_db/ensembl94_gene_prot_pfam_id.txt"
      #"~/PhD/2018-2019/pdb_db/ensembl94_gene_prot_pfam_id.txt"
    )
  # change some column names to match them when merging with other dfs
  names(biomart)[3] <- "ensembl.prot.id"
  names(biomart)[c(7,8,9)] <- c("pfam.dom", "dom.start", "dom.end")
  # merge new resulst with pfam domain information
  new_mpi_biomart <-
    left_join(new_mpi, biomart, by = "ensembl.prot.id")
  # new_mpi_pfam_ensembl1 <- left_join(new_mpi_pfam ,
  #                                    biomart[, 1:6],
  #                                    by = c("ensembl.prot.id"))
  
  # merge (new resulst + pfam) with ensembl info
  # new_mpi_pfam_ensembl <-
  #   left_join(new_mpi_pfam_ensembl1 ,
  #             biomart[, c(3, 7:9)],
  #             by = c("ensembl.prot.id", "pfam.dom", "dom.start", "dom.end"))
  
  # return data frame with the updated information.
  return(new_mpi_biomart)
}
