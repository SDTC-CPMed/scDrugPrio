#' NicheNet ligand activity analysis for all cell types
#' Will test cell type centrality based on NicheNet ligand activities. Results in an output file `Cell_type_centrality_summary.txt` that documents calculated centralities and properties of the intercellular disease model.
#'
#' @param all_ligand_activity Matrix. Outcome of `NicheNet_ligand_activity_analysis()`.
#' @param out_dir String. Indicates the full path to the folder in which results should be stored.
#'
#' @return Returns a matrix that summarizes cell type centralities in the multicellular disease module (MCDM)
#' @export
#'
NicheNet_cell_type_centrality <- function(all_ligand_activity, out_dir = ""){

  # Set working directory to out_dir
  set.seed(35)

  # CHECK IF BASIC REQUIREMENTS ARE FULFILLED
  ######################################################################

  if(out_dir != ""){
    if(substr(x = out_dir, start = nchar(out_dir), stop = nchar(out_dir)) != "/"){
      out_dir <- paste(out_dir, "/",sep="")
    }
  }

  if(is.null(all_ligand_activity)){
    stop(
      "Variable \'all_ligand_activity\' not specified.",
      call. = FALSE
    )
  }

  if (!requireNamespace("CINNA", quietly = TRUE)) {
    stop(
      "Package \'CINNA\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(CINNA)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \'igraph\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(igraph)


  # Perform centrality calculations
  ######################################################################
  all_ligand_activity <- as.matrix(all_ligand_activity)
  all_ligand_activity <- all_ligand_activity[as.numeric(all_ligand_activity[,4])>0,]

  # n cell types
  cell_types <- unique(c(all_ligand_activity[,ncol(all_ligand_activity)-1], all_ligand_activity[,ncol(all_ligand_activity)]))

  # Create graph
  g <- graph_from_edgelist(el = all_ligand_activity[,(ncol(all_ligand_activity)-1):ncol(all_ligand_activity)], directed = T)
  E(g)$weight <- as.numeric(all_ligand_activity[,4])

  out <- matrix(NA, nrow = 4, ncol = length(cell_types))


  rownames(out) <- c("node_degree_all", "node_degree_in", "node_degree_out", "closeness")
  colnames(out) <- as.character(V(g)$name)

  # Calculate centrality degree
  out[1,] <- centr_degree(g, mode = "all")$res # all
  out[2,] <- centr_degree(g, mode = "in")$res # in degree
  out[3,] <- centr_degree(g, mode = "out")$res # out degree

  # Calculate closeness
  out[4,] <- closeness(g, mode = "all")

  # Calculate centralities
  centrality_matrix <- calculate_centralities(g, include = c("eigenvector centralities", 
                                                             "K-core Decomposition", 
                                                             "Kleinberg's hub centrality scores",
                                                             "Laplacian Centrality",
                                                             "Leverage Centrality",
                                                             "Group Centrality",
                                                             "Local Bridging Centrality"),
                                              weights = E(g)$weight)

  centrality_matrix <- matrix(unlist(centrality_matrix), ncol = length(cell_types), byrow = T)
  colnames(centrality_matrix) <- as.character(V(g)$name)
  rownames(centrality_matrix) <- c("eigenvector centralities", "K-core Decomposition", "Kleinberg's hub centrality scores",
                                   "Laplacian Centrality", "Leverage Centrality", "Group Centrality", "Local Bridging Centrality")

  out <- rbind(out, centrality_matrix)

  write.table(out, file = paste(out_dir, "Cell_type_centrality_summary.txt", sep=""), sep="\t", col.names = NA, row.names = T)
  return(out)
}
