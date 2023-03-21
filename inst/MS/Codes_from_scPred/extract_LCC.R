#' Largest connected component
#'
#' Extraction of the largest connected component (LCC) within an igraph object
#'
#' @param network igraph object.
#'
#' @return Node labels of nodes included in LCC.
#' @export
#'
#' @examples
#' ppin <- as.data.frame(read.table(file = "data-raw/lit_ppi.txt", sep="\t", header = T, stringsAsFactors = F))
#' ppin_graph <- igraph::graph_from_data_frame(ppin, directed = F)
#' lcc_genes <- extract_LCC(ppin_graph) # only ppin proteins that are connected to other ppin proteins
#'
extract_LCC <- function(network){
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \'igraph\' must be installed to use this function.",
      call. = FALSE
    )
  }

  library(igraph)

  comp <- clusters(graph = network, mode = "weak")
  lcc <- which.max(comp$csize)
  vert_ids <- V(network)[comp$membership == lcc]
  subGr <- induced_subgraph(network, vids=vert_ids)

  return(V(subGr)$name)
}
