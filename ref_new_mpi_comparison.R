###############################################################################
# Comparison between reference results and updated by pdb id                  #
###############################################################################
library(dplyr)

Ref.New.Mpi.Comparison <-
  # input: ref_mpi -> reference results of a single pdb biounit
  #        new_mpi -> new results of a single pdb biounit
  #        mpi_compared -> df with ref and new data merged of a single pdb biounit
  function(ref_mpi, new_mpi, mpi_compared, pdb_id) {
    
    # dataframe summarizing results of reference
    ref_results <- data.frame(
      # total number of template chains with mapped interfaces
      chain =  length(setdiff(
        unique(ref_mpi$chain) ,
        unique(new_mpi$chain)
      )),
      # total number of interacting chains with mapped interfaces
      chain.1 =  length(setdiff(
        unique(ref_mpi$chain.1) ,
        unique(new_mpi$chain.1)
      )),
      # total number of ensembl prot ids with mapped interfaces
      ensembl.prot.id = length(setdiff(
        unique(ref_mpi$ensembl.prot.id) ,
        unique(new_mpi$ensembl.prot.id)
      )),
      # total number of different starting of each ensembl prot id mapped
      sstart = nrow(anti_join(
        unique(ref_mpi[, c("sstart", "ensembl.prot.id")]),
        unique(new_mpi[, c("sstart", "ensembl.prot.id")])
      )),
      # total number of different starting of each ensembl prot id mapped
      send = nrow(anti_join(
        unique(ref_mpi[, c("send", "ensembl.prot.id")]),
        unique(new_mpi[, c("send", "ensembl.prot.id")])
      )),
      # total number of mapped residues
      mapped.real.pos = nrow(ref_mpi) -  nrow(subset(mpi_compared, freq == 2)),
      # pdb id analysed
      id = ifelse( nrow(ref_mpi)> 0 , unique(ref_mpi$pdb.id.x), pdb_id),
      # identifier to plot venn diagrams afterwards
      data.id = "Reference(2017)"
    )
 
    # dataframe summarizing new results
    new_results <- data.frame(
      # total number of template chains with mapped interfaces
      chain =  length(setdiff(
        unique(new_mpi$chain),
        unique(ref_mpi$chain)
      )),
      # total number of interacting chains with mapped interfaces
      chain.1 = length(setdiff(
        unique(new_mpi$chain.1),
        unique(ref_mpi$chain.1)
      )),
      # total number of ensembl prot ids with mapped interfaces
      ensembl.prot.id = length(setdiff(
        unique(new_mpi$ensembl.prot.id),
        unique(ref_mpi$ensembl.prot.id)
      )),
      # total number of different starting of each ensembl prot id mapped
      sstart = nrow(anti_join(unique(new_mpi[, c("sstart", "ensembl.prot.id")]),
                              unique(ref_mpi[, c("sstart", "ensembl.prot.id")]))),
      # total number of different starting of each ensembl prot id mapped
      send = nrow(anti_join(unique(new_mpi[, c("send", "ensembl.prot.id")]),
                            unique(ref_mpi[, c("send", "ensembl.prot.id")]))),
      # total number of mapped residues
      mapped.real.pos =  nrow(new_mpi) -  nrow(subset(mpi_compared, freq == 2)),
      # pdb id analysed
      id = ifelse( nrow(new_mpi) > 0 , as.character(unique(new_mpi$pdb.id.x)), pdb_id),
      # identifier to plot venn diagrams afterwards
      data.id = "New(2019)"
    )
  
    # intersection between reference and new results
    intersect_results <- data.frame(
      # template chains
      chain = length(intersect(
        unique(ref_mpi$chain) ,
        unique(new_mpi$chain)
      )),
      # interacting chains
      chain.1 = length(intersect(
        unique(ref_mpi$chain.1) ,
        unique(new_mpi$chain.1)
      )),
      # ensembl protein ids
      ensembl.prot.id = length(intersect(
        unique(ref_mpi$ensembl.prot.id) ,
        unique(new_mpi$ensembl.prot.id)
      )),
      # start position
      sstart = nrow(inner_join(unique(new_mpi[, c("sstart", "ensembl.prot.id")]),
                               unique(ref_mpi[, c("sstart", "ensembl.prot.id")]))),
      # stop position
      send =  nrow(inner_join(unique(new_mpi[, c("send", "ensembl.prot.id")]),
                              unique(ref_mpi[, c("send", "ensembl.prot.id")]))),
      # residues that are exactly in the same location
      mapped.real.pos = length(which(mpi_compared$freq == 2)),
      # pdb id analysed
      id = as.character(unique(new_mpi$pdb.id.x)),
      # venn diagram id
      data.id = "Reference(2017)&New(2019)"
    )
    
    # store previous data frame into a single one (to be used to generate Venn diagrams afterwards)
    comparison_ref_new_table <- rbind(ref_results,
                                      new_results,
                                      intersect_results)
    
    return(comparison_ref_new_table)
  }