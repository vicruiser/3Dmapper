################################################################################################
#' Title  CALCULATION OF INTER-ATOMIC DISTANCES BETWEEN CHAINS OF A PDB FILE. 
#'
#' @param pdb_file 
#' @param comparison either 'protein', 'nucleic' or 'ligand'
#' @param atom_type_prot type of atom that we want to compare from 
#' template chain. By default all heavy atoms are selected
#' @param atom_type_other  type of atom that we want to compare from
#'  interacting chain. By default all heavy atoms are selected
#'
#' @return  list of distances matrices of all possible chains combinations
#' @export
#'
#' @examples
#' 

calc_PDB_dist <- function(pdb_file,
                          CompareTemplateChains,
                          CompareInteractionChains,
                          CompareTemplateChains.Vector,
                          CompareInteractionChains.Vector,
                          dist_threshold = 5,
                          type_of_interaction = "protein") {
 
  allDist <- dist.xyz(CompareTemplateChains.Vector[[1]],
                      CompareInteractionChains.Vector[[1]])
  
  filtered_pdb_dist <- which(allDist <= dist_threshold, arr.ind = T)
  
  if(nrow(filtered_pdb_dist) > 0){

    df_results <- data.frame(CompareTemplateChains[filtered_pdb_dist[,"row"],],
                             CompareInteractionChains[filtered_pdb_dist[,"col"],],
                             distance = allDist[filtered_pdb_dist],
                             interaction = type_of_interaction)

    return(df_results)
    
  } else {
    
    stop("No existing interaction between these two chains.")
    
  }

}
