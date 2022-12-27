#' Intracellular drug traget celtrality calculation
#'
#' If no Largest Connected Component (LCC) can be found for a certain cell type, this cell type will be excluded completly from the calculations. As LCC we here defined a component of at least three proteins interacting with each other.
#'
#'
#' @param ppin Matrix with 2 columns, every row represents one protein-protein interaction in the protein-protein interaction network (PPIN) largest connected component (LCC). Output from ppin_formatting()
#' @param drug_target_matrix Matrix listing each drug's targets. One column per drug. Only drug targets that are included in 'ppin' are entered used. If drug targets that are not included in 'ppin' are included they will be removed prior to centrality calculations.
#' @param degs Matrix or vector. If matrix, DEGs for each cell type are stored in a column.
#' @param file_name Character string specifying the first part of the filename for the output.
#' @param cores number of cores used for computing.
#' @param centrality_alg string. Specifying the centrality algorithm that shall be used for calculation of intracellular centrality. Standard is 'eigenvector centralities'. The function 'calculate_centralities()' from the CINNA package is used for calculation of centralities and hence this parametar should be set according to what is available for this function. Other valid options include (but are not limited to): "K-core Decomposition", "Kleinberg's hub centrality scores", "Laplacian Centrality", "Leverage Centrality", "Group Centrality", "Local Bridging Centrality".
#' @param out_dir Character string specifying path to out directory where data will be saved.
#'
#' @return a matrix including the geometric mean drug target centrality in every cell types LCC (if existent) as well as the arithmetric mean of all geometric mean centrality scores. The later of which was by us used for prioritization among drugs.
#' @export
#'
#'
intracellular_drug_target_centrality <- function(ppin,
                                                 drug_target_matrix,
                                                 degs,
                                                 file_name = "",
                                                 cores = 1,
                                                 centrality_alg = "eigenvector centralities",
                                                 out_dir = "",
                                                 seed = 35){

  set.seed(seed)

  if(out_dir != ""){
    if(substr(x = out_dir, start = nchar(out_dir), stop = nchar(out_dir)) != "/"){
      out_dir <- paste(out_dir, "/",sep="")
    }
  }


  # CHECK IF BASIC REQUIREMENTS ARE FULFILLED
  ######################################################################

  if(is.null(ppin)){
    stop(
      "Variable \'ppin\' not specified.",
      call. = FALSE
    )
  }
  if(is.null(drug_target_matrix)){
    stop(
      "Variable \'drug_target_matrix\' not specified.",
      call. = FALSE
    )
  }
  if(is.null(degs)){
    stop(
      "Variable \'degs\' not specified.",
      call. = FALSE
    )
  }

  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop(
      "Package \'doParallel\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(doParallel)
  registerDoParallel(cores = cores)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \'igraph\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(igraph)

  if (!requireNamespace("CINNA", quietly = TRUE)) {
    stop(
      "Package \'CINNA\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(CINNA)


  # Check if input has correct format
  ######################################################################

  # Check 'ppin'
  if(is.matrix(ppin)){
    if(ncol(ppin) != 2){
      stop("\'ppin\' does not have two columns.", call. = FALSE)
    } else {
      all_genes <- unique(as.vector(ppin))
    }
  } else {
    stop("\'ppin\' is not a matrix.", call. = FALSE)
  }

  # Check 'drug_target_matrix'
  if(is.matrix(drug_target_matrix)){
    # check that all drug targets are included in PPIN
    all_drugt <- unique(as.vector(drug_target_matrix))
    all_drugt <- all_drugt[!is.na(all_drugt)]
    if(any(!all_drugt %in% all_genes)){
      warning("Not all targets specified in \'drug_target_matrix\' are found in \'ppin\' and hence they will be removed prior to the analysis.")
      for(i in 1:ncol(drug_target_matrix)){
        temp <- drug_target_matrix[,i]
        drug_target_matrix[,i] <- NA
        temp <- temp[temp %in% all_genes]
        if(length(temp)>0){
          drug_target_matrix[1:length(temp),i] <- temp
        }
      }
      if(any(colSums(!is.na(drug_target_matrix)) == 0)){
        warning("Not all drugs have targets in \'ppin\' and hence these drugs will be excluded from the analysis.")
        if(all(colSums(!is.na(drug_target_matrix)) == 0)){
          stop("No drug has targets in \'ppin\' and hence analysis can not be performed.", call. = FALSE)
        }
      }
      drug_target_matrix <- drug_target_matrix[rowSums(!is.na(drug_target_matrix)) > 0, colSums(!is.na(drug_target_matrix)) > 0]
    }
  } else {
    stop("\'drug_target_matrix\' is not a matrix.", call. = FALSE)
  }

  # Check 'degs'
  if(is.vector(degs)){
    if(any(degs %in% all_genes)){
      degs <- degs[degs %in% all_genes]
      degs <- matrix(degs[!is.na(degs)], ncol = 1)
    } else {
      stop("The vector supplied for the \'degs\' parameter does not include any values found in \'ppin\'.", call. = FALSE)
    }
  } else if (is.matrix(degs)) { # if degs = matrix
    for(i in 1:ncol(degs)){
      temp <- degs[,i]
      degs[,i] <- NA
      temp <- temp[temp %in% all_genes]
      if(length(temp)>0){
        degs[1:length(temp),i] <- temp
      }
    }
    if(all(colSums(!is.na(degs))==0)){
      stop("The matrix supplied for the \'degs\' parameter does not include any elements matchable with \'ppin\'.",call. = FALSE)
    }
    if(sum(colSums(!is.na(degs))>0)>1){
      if(sum(rowSums(!is.na(degs))>1)>1){
        degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
      } else {
        stop("The matrix supplied for the \'degs\' does only include one valid DEG per cell type which does not allow the LCC extraction'.",call. = FALSE)
      }
    } else {
      if(any(colSums(!is.na(degs)) == 1)){
        stop("The matrix supplied for the \'degs\' does only include one valid DEG per cell type which does not allow the LCC extraction'.",call. = FALSE)
      } else {
        warning(paste("\'degs\' does only includes elements matchable with \'ppin\' in the column: ", colnames(degs)[colSums(!is.na(degs))>0], sep = ""))
        degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
        degs <- matrix(degs, ncol = 1)
      }
    }
  } else {
    stop("\'degs\' is not a vector or matrix.", call. = FALSE)
  }

  if(length(centrality_alg)>1){
    stop("\'centrality_alg\' can not exceed length > 1.", call. = FALSE)
  }

  # FORMAT DATA
  ######################################################################

  ppin <- data.frame(node1 = as.character(ppin[,1]), node2 = as.character(ppin[,2]))
  ppi_graph <- graph_from_data_frame(ppin, directed = F) # igraph package


  # EXTRACT LCCs FOR EVERY CELL TYPE - disease neighbourhoods
  ######################################################################

  # Check if any cluster has less than 3 genes, which will not suffice for a LCC
  # Warn if any cell types does not have any LCC
  if(any(colSums(!is.na(degs)) < 3)){
    print(paste("The following clusters / cell types have too few genes to have a LCC:", paste(colnames(degs)[colSums(!is.na(degs)) < 3], collapse = ", "), sep =" "))
    degs <- degs[,colSums(!is.na(degs)) > 2]
  }

  LCCs <- foreach(i = c(1:ncol(degs)), .combine = "cbind") %dopar% {

    # extract LCC
    lcc_genes <- extract_LCC_by_gene_set_of_interest(network = ppi_graph, gene_set_of_interest = as.character(degs[!is.na(degs[,i]),i])) # only PPI proteins that are connected to other PPI proteins

    # check if two equal sized LCCs
    if(!is.null(ncol(lcc_genes))){ # several possible subgraphs
      if(ncol(lcc_genes)>1){ # several subgraphs
        if(nrow(lcc_genes)>1){
          lcc_genes <- lcc_genes[rowSums(!is.na(lcc_genes))>0,] # removes empty rows
          if(sum(colSums(!is.na(lcc_genes)) == max(colSums(!is.na(lcc_genes))))>1){ # two equal sized LCCs
            pos <- which(colSums(is.na(lcc_genes))==0)
            lcc_genes <- as.vector(lcc_genes[, pos[1]]) # 'lcc_genes' of the first discovered maximum subgraph will be used
          } else {
            lcc_genes <- as.vector(lcc_genes[, colSums(is.na(lcc_genes))==0])
          }
        } else { # no subgraphs with more than one element and a subgraph with only one element is not eligable for centrality calculations
          lcc_genes <- NA
        }
      } else { # only one subgraph
        lcc_genes <- lcc_genes[,1]
      }
    }

    # return
    return(c(lcc_genes, rep(NA, times = length(all_genes)-length(lcc_genes))))
  }
  colnames(LCCs) <- colnames(degs)


  # CALCULATE DRUG TARGET CENTRALITY IN LCCs
  ######################################################################

  # Define geometric mean function
  gm_mean = function(a){prod(a)^(1/length(a))}

  # Check if any valid LCC with at least 3 proteins can be extracted
  # Warn if any cell types does not have any LCC
  if(any(colSums(!is.na(LCCs)) < 3)){
    print(paste("The following clusters / cell types have no LCC:", paste(colnames(LCCs)[colSums(!is.na(LCCs)) < 3], collapse = ", "), sep =" "))
    if(sum(colSums(!is.na(LCCs)) > 3) > 1){
      LCCs <- LCCs[,colSums(!is.na(LCCs)) > 2]
    } else {
      temp <- colnames(LCCs)[colSums(!is.na(LCCs)) > 3]
      LCCs <- matrix(LCCs[,colSums(!is.na(LCCs)) > 2], ncol = 1)
      colnames(LCCs) <- temp
    }
  }

  # Find mean target centrality
  mean_target_centrality <- foreach(i = c(1:ncol(LCCs)), .combine = "cbind") %dopar% { # for every cell type specific LCC
    if(sum(!is.na(LCCs[,i]))>1){
      # Calculate centrality for proteins in disease neighborhood of this cluster
      ###########################################################################
      # make graph of disease neighboorhood LCC
      cand <- V(ppi_graph)[name %in% LCCs[,i]]
      g <- induced_subgraph(ppi_graph, vids = cand)

      # Calculate centrality
      centrality_matrix <- calculate_centralities(g, include = centrality_alg)
      centrality_matrix <- unlist(centrality_matrix)
      names(centrality_matrix) <- V(g)$name

      # Select only centralities of DEGs of this cluster
      cent <- vector()
      for(d in 1:ncol(drug_target_matrix)){ #for every drug candidate
        targets <- drug_target_matrix[,d]
        targets <- targets[!is.na(targets)]

        if(any(as.character(targets) %in% names(centrality_matrix))) {
          temp <- centrality_matrix[names(centrality_matrix)%in% targets]
        } else {
          temp <- 0
        }

        # Geometric mean of Eigenvector centrality of drug targets
        cent[d] <- gm_mean(temp)
      }
      return(cent)
    } else {
      return(rep(0, times = ncol(drug_target_matrix)))
    }
  }
  if(ncol(LCCs) == 1){
    mean_target_centrality <- matrix(mean_target_centrality, ncol = 1)
  }
  rownames(mean_target_centrality) <- colnames(drug_target_matrix)
  colnames(mean_target_centrality) <- colnames(LCCs)

  # Aggregate into mean centrality sum for all candidates in all LCCs
  mean_target_centrality <- cbind(mean_target_centrality, intracellular_centrality_mean = rowMeans(mean_target_centrality))

  # Save output
  write.table(mean_target_centrality, file = paste(out_dir, file_name, ".txt",sep=""),sep="\t", col.names = NA, row.names = T)

  return(mean_target_centrality)
}
