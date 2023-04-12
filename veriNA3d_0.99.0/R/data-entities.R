#' List of NA containing PDB IDs and related data
#'
#' Number of each type of entity in a set of Nucleic Acid containing 
#' structures.
#'
#' @docType data
#'
#' @usage data(entities)
#'
#' @format An object of class data.frame
#' \describe{
#'     \item{pdbID:}{PDB ID.}
#'     \item{RNA:}{Number of different RNA entities.}
#'     \item{DNA:}{Number of different DNA entities.}
#'     \item{Hybrid:}{Number of different Hybrid DNA/RNA entities.}
#'     \item{PNA:}{Number of different PNA entities.}
#'     \item{Prot:}{Number of different Protein (L) entities.}
#'     \item{Dprot:}{Number of different Protein (D) entities.}
#'     \item{Ligands:}{Number of different ligands.}
#'     \item{Water:}{It's 1 when there's water in the structure 
#'                     and 0 when it is not.}
#'     \item{Other:}{Number of different "other" entities.}
#' }
#'
#' @return data.frame with a list and features of NA PDB IDs.
#'
"entities"
