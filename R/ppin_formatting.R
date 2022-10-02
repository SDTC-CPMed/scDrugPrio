#' Protein-protein interaction network (PPIN) formatting
#'
#' @param ppin matrix with 2 columns specifying all protein-protein interactions in the PPIN
#'
#' @return A 2 column matrix specifying the protein-protein interactions included in the PPINs largest connected component
#' @export
#'
#' @examples
#' load_example_data()
#' ppin_formatting(test_ppi)
#'
ppin_formatting <- function(ppin){

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \'igraph\' must be installed to use this function.",
      call. = FALSE
    )
  }

  library(igraph)

  # Self loops and incomplete interactions
  ppin <- ppin[!(ppin[,1]==ppin[,2]),] # no self loops
  ppin <- ppin[!is.na(ppin[,1]),] # no NA
  ppin <- ppin[!is.na(ppin[,2]),] # no NA

  # Extracting LCC
  mode(ppin) <- "character"
  ppin <- data.frame(node1 = ppin[,1], node2 = ppin[,2])
  ppin_graph <- graph_from_data_frame(ppin, directed = F) # igraph package
  lcc_genes <- extract_LCC(ppin_graph) # only ppin proteins that are connected to other ppin proteins
  print(paste("n of unique proteins/genes in PPIN: n = ", length(unique(lcc_genes))))

  ppin <- as.matrix(ppin)
  ppin <- ppin[ppin[,1] %in% lcc_genes,]
  ppin <- ppin[ppin[,2] %in% lcc_genes,]

  return(ppin)
}
