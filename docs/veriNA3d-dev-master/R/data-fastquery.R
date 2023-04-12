#' List of NA containing PDB IDs and related data
#'
#' Presaved data to speed up some queries.
#'
#' @docType data
#'
#' @usage data(fastquery)
#'
#' @format An object of class data.frame with the following fields:
#' \describe{
#'     \item{pdbID:}{PDB ID.}
#'     \item{Technique:}{Experimental technique.}
#'     \item{Resol:}{Resolution. For NMR structures it contains an empty 
#'                             string.}
#'     \item{DNAclass:}{Output of classifyDNA function.}
#'     \item{RNAclassOver0:}{Output of classifyRNA with length=0 or length=1. If
#'                             the structure has one or two nucleotides, it is 
#'                             also considered an RNA containing structure.}
#'     \item{RNAclassOver2:}{Output of classifyRNA with length=3. RNA molecules
#'                             shorter than 3 are classified as NoRNA.}
#' }
#'
#' @return data.frame with a list and features of NA PDB IDs.
#'
"fastquery"
