#! /usr/bin/Rscript
################################################################################################
#' Title  CALCULATION OF INTER-ATOMIC DISTANCES BETWEEN CHAINS OF A PDB FILE.
#'
#' @param CompareTemplateChains
#' @param CompareInteractionChains
#' @param CompareTemplateChains.Vector
#' @param CompareInteractionChains.Vector
#' @param dist_threshold
#' @param type_of_interaction either 'protein', 'nucleic' or 'ligand'
#
#'
#' @return  list of distances matrices of all possible chains combinations
#' @export
#'
#' @examples
#'
options(echo = FALSE, verbose = F,warn = -1) 
calc_PDB_dist <- function(
                          CompareTemplateChains,
                          CompareInteractionChains,
                          CompareTemplateChains.Vector,
                          CompareInteractionChains.Vector,
                          dist_threshold,
                          type_of_interaction ) {

  allDist <- dist.xyz(CompareTemplateChains.Vector,#[[1]],
                      CompareInteractionChains.Vector) #[[1]])
  filtered_pdb_dist <- which(allDist <= dist_threshold, arr.ind = T)

  if (nrow(filtered_pdb_dist) > 0) {
    df_results <-
      data.frame(
        CompareTemplateChains[filtered_pdb_dist[, "row"], ],
        CompareInteractionChains[filtered_pdb_dist[, "col"], ],
        distance = allDist[filtered_pdb_dist],
        interaction = type_of_interaction
      )

    return(df_results)
    
  } else {
    df_results = data.frame()
    return(df_results)
    #stop("No existing interaction between these two chains.")
    
  }
}