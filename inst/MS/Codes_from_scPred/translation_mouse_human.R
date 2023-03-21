#' A translation matrix for human and mice Entrez gene IDs and symbols.
#' This list is primarily for the conversion of mouse to human gene annotation by mapping mouse genes to homologous human genes.
#' 
#' @format A matrix with 12825 rows and 4 variables:
#' \describe{
#'   \item{human_gene_symbol}{Human Entrez gene symbol}
#'   \item{human_entrez_ID}{Human Entrez gene ID}
#'   \item{mouse_gene_symbol}{Mouse Entrez gene symbol}
#'   \item{mouse_entrez_ID}{Mouse Entrez gene ID}
#'   ...
#' }
#' @source National Centre for Biotechnology Information, downloaded July 2017.
"translation_mouse_human"