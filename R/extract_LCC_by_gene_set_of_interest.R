#' Extraction of the LCC formed by a gene set of interest in a PPIN
#'
#' @param network igraph object. Largest connceted component of all proteins in a protein-protein interaction network (PPIN).
#' @param gene_set_of_interest Vector specifying nodes in 'network'
#' @param plot Boolean, if TRUE creates plot of depicturing isolated components formed by 'gene_set_of_interest' in 'network'
#'
#' @return Vector specifying the node names of the LCC formed by the 'gene_set_of_interest'.
#' @export
#'
#' @examples
#' Checks for the LCC formed by 200 unique nodes in the network
#' ppin <- as.data.frame(read.table(file = "data-raw/lit_ppi.txt", sep="\t", header = T, stringsAsFactors = F))
#' goi <- unique(as.vector(ppin))[1:200]
#' ppin_graph <- igraph::graph_from_data_frame(ppin, directed = F)
#' lcc_genes <- extract_LCC_by_gene_set_of_interest(network = ppin_graph, gene_set_of_interest = goi) # only ppin proteins that are connected to other ppin proteins
#'
extract_LCC_by_gene_set_of_interest <- function(network, gene_set_of_interest, plot = F){
  if(!is.igraph(network)){
    stop(
      "Variable \'network\' not specified.",
      call. = FALSE
    )
  }
  if(!is.vector(gene_set_of_interest)){
    stop(
      "Variable \'gene_set_of_interest\' not specified.",
      call. = FALSE
    )
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \'igraph\' must be installed to use extract_LCC_by_gene_set_of_interest()",
      call. = FALSE
    )
  }
  library(igraph)
  if(plot){
    if (!requireNamespace("qgraph", quietly = TRUE)) {
      stop(
        "Package \'qgraph\' must be installed to use plot = TRUE in extract_LCC_by_gene_set_of_interest()",
        call. = FALSE
      )
    }
    library(qgraph)
  }

  V(network)$gene_set_of_interest <- 0
  V(network)$gene_set_of_interest[V(network)$name %in% gene_set_of_interest] <- 1
  gsoi <- V(network)[V(network)$gene_set_of_interest==1]
  subGr <- induced_subgraph(network, vids= gsoi)

  if(plot){
    e <- get.edgelist(subGr)
    e[,1] <- match(e[,1], V(subGr)$name)
    e[,2] <- match(e[,2], V(subGr)$name)
    mode(e) <- "numeric"
    l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subGr), area=8*(vcount(subGr)^2),repulse.rad=(vcount(subGr)^3.1))
    plot(subGr,layout=l,vertex.size=2,vertex.label=NA, vertex.color=c("grey","red")[1+(V(subGr)$gene_set_of_interest==1)])
  }

  subGr <- decompose(subGr, mode = "weak")

  if(length(subGr) > 0){
    out <- matrix(data = NA, nrow = length(V(network)), ncol = length(subGr))
    for(i in c(1:ncol(out))){
      temp <- V(subGr[[i]])$name
      out[1:length(temp),i] <- temp
    }
    if(all(colSums(!is.na(out)) == 1)){ # none of the elements in gene_set_of_interest form a connection with each other -> NO LCC
      out <- NA
    } else if (sum(colSums(!is.na(out)) > 1) == 1) { # only one column with all LCC genes
      out <- matrix(out[,colSums(!is.na(out)) > 1], ncol = 1)
    }
    if(!is.vector(out)){ # Decide which is the biggest connected component
      n_targets <- colSums(!is.na(out))
      out <- as.vector(out[,which.max(n_targets)])
      out <- out[!is.na(out)]
    }
  } else {
    warning(
      "No \'gene_set_of_interest\' found in 'network' and hence extract_LCC_by_gene_set_of_interest() will return NULL",
      call. = FALSE
    )
    out <- NULL
  }
  return(out)
}
