#' Closest Distance Network Drug Screening
#' Used for calculations of closest network proximity. Memory heavy.
#'
#' @param ppin Matrix with 2 columns, every row represents one protein-protein interaction in the protein-protein interaction network (PPIN) largest connected component (LCC). Output from ppin_formatting()
#' @param drug_target_matrix Matrix listing each drug's targets in a column.
#' @param disease_genes Vector specifying disease associated genes (DEGs, GWAS genes, OMIM genes, etc.)
#' @param file_name Character string specifying filename for the output.
#' @param disease_genes_lcc Boolean, if TRUE the LCC of disease genes is used for the calculations instead of disease genes.
#' @param min_bin_size Integer. The bin-size used for node degree adjustment of network calculations.
#' @param n_random_iterations number of permutations used to derive reference distribution of network distances.
#' @param cores number of cores used for computing.
#' @param prepare_recycle Boolean, if TRUE ppin_distances, bins and a cleaned version of drug_target_matrix are saved in  "/temp" folder for shortening computational time in coming iterations.
#' @param recycle Boolean, if TRUE ppin_distances and bins are loaded from "/temp". Decreases computational time. Can not be used if parallel sessions are accessing the same folder.
#' @param out_dir Character string specifying path to out directory where data will be saved.
#'
#' @return
#' @export
#'
#' @examples
#' load_example_data()
#' ppin <- ppin_formatting(lit_ppi)
#' drug_target_matrix <- prepare_drug_target_matrix_for_network_distance_calculation(ppin, drug_target_matrix, file_name = "in_lit_ppin", out_dir = "data/")
#' closest_distance_network_drug_screening(ppin = ppin, drug_target_matrix = drug_target_matrix, file_name = "test")
#'
average_closest_distance_network_drug_screening <- function(ppin,
                                               drug_target_matrix,
                                               disease_genes,
                                               file_name,
                                               disease_genes_lcc = F,
                                               min_bin_size = 100,
                                               n_random_iterations = 1000,
                                               cores = 1,
                                               prepare_recycle = F,
                                               recycle = F,
                                               out_dir = getwd()){

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
  if(is.null(disease_genes)){
    stop(
      "Variable \'disease_genes\' not specified.",
      call. = FALSE
    )
  }
  if(is.null(file_name)){
    stop(
      "Variable \'file_name\' not specified.",
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
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \'igraph\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(igraph)
  if(is.null(colnames(drug_target_matrix))){
    stop(
      "Variable \'drug_target_matrix\' must have valid colnames()",
      call. = FALSE
    )
  }
  if(ncol(ppin) != 2){
    stop(
      "Variable \'ppi\' has invalid dimensions",
      call. = FALSE
    )
  }
  if(!is.vector(disease_genes)){
    stop(
      "Variable \'disease_genes\' is not a vector",
      call. = FALSE
    )
  }

  # Set working directory to out_dir
  setwd(dir = out_dir)
  set.seed(35)


  # SET UP ENVIRONMENT
  ######################################################################

  # set up parallel running
  if(cores > 1){
    max_cores <- detectCores()
    if(cores <= max_cores){
      registerDoParallel(cores = cores)
    } else {
      stop(
        "Variable \'cores\' higher than available max available",
        call. = FALSE
      )
    }
  }

  # PREPARE INPUT OR LOAD PREPARED INPUT (IF RECYCLE = TRUE)
  ######################################################################

  if(recycle){
    print(list.files(path= "temp"))
    # ppin_dist
    ppin_dist <- read.table(file = "temp/ppin_distances.txt", sep="\t", header = F, stringsAsFactors = F)
    colnames(ppin_dist) <- as.character(ppin_dist[1,])
    ppin_dist <- ppin_dist[-1,]
    rownames(ppin_dist) <- as.character(ppin_dist[,1])
    ppin_dist <- as.matrix(ppin_dist[,-1])
    bins <- read.table(file = "temp/bins.txt", sep="\t", header = F, stringsAsFactors = F)
    # all_genes
    all_genes <- as.vector(read.table(file = "temp/all_genes.txt", sep="\t", header = F, stringsAsFactors = F)[,1])
    print("RECYCLING FILES")

  } else {

    dir.create(path = paste(getwd(), "/temp", sep=""))

    # Setting up ppin
    all_genes <- unique(c(unique(ppin[,1]), unique(ppin[,2])))
    ppin_graph <- graph_from_data_frame(data.frame(node1 = ppin[,1], node2 = ppin[,2]), directed = F) # igraph package
    if(prepare_recycle){  ################ SAVE ALL_GENES ################
      write.table(matrix(all_genes, ncol=1), file = "temp/all_genes.txt", sep="\t", col.names = F, row.names = F)
    }

    # Creating / Loading bins
    bins <- bin_creation_by_min_bin_size(ppin = ppin, min_bin_size = min_bin_size)
    if(prepare_recycle){  ################ SAVE BINS ################
      write.table(bins, file = "temp/bins.txt", sep="\t", col.names = F, row.names = F)
    }

    # Calculate distances between all proteins in LCC
    ppin_dist <- distances(ppin_graph, v = all_genes, to = all_genes)
    colnames(ppin_dist) <- as.character(colnames(ppin_dist))
    rownames(ppin_dist) <- as.character(rownames(ppin_dist))
    if(prepare_recycle){ ################ SAVE PPIN DISTANCES ################
      write.table(ppin_dist, file = "temp/ppin_distances.txt", sep="\t", col.names = NA, row.names = T)
    }

    if(!disease_genes_lcc){
      rm(ppin_graph)
    }
  }

  # PREPARE DISEASE GENES FOR CALCULATION
  ######################################################################

  disease_genes <- as.vector(disease_genes)
  disease_genes <- disease_genes[disease_genes %in% all_genes] # must be in PPIN

  if(disease_genes_lcc){ # if TRUE calculate LCC formed by disease_genes genes and use that instead
    print("Identifying LCC for disease genes in PPIN")
    disease_genes <- extract_LCC_by_gene_set_of_interest(ppin_graph, gene_set_of_interest = disease_genes, plot = F) # only ppin proteins that are connected to other ppin proteins
    if(is.null(disease_genes)){
      print("No LCC found for 'disease_genes'")
      return()
    }
  }

  # CALCULATION OF CLOSEST NETWORK DISTANCE
  ######################################################################

  print("RUNNING bin-adjusted reference distribution for average closest distances between random drug targets and random disease genes")
  # random_bin_drugs (1 column for every drug, 1 row for every n_random_iteration)
  random_bin_drugs <- foreach(i = c(1:n_random_iterations), .combine = rbind) %dopar% {
    # returns vector with n_random_iterations elements cotaining closest distances between disease and
    # (bin-adjusted) random drug tragets for a given drug
    random_drug_target_bin_adjusted_distances(bins, drug_target_matrix, disease_genes, ppin_dist, seed = i)
  }
  write.table(random_bin_drugs, file = paste("random_bin_drugs_",file_name,".txt", sep=""), sep = "\t", col.names = colnames(drug_target_matrix), row.names = F)


  # Calculation of average closest distances between
  print("RUNNING average closest distances between actual drug targets and actual disease genes")
  n_targets <- colSums(!is.na(drug_target_matrix))
  out <- foreach(i = 1:ncol(drug_target_matrix), .combine = rbind) %dopar% {
    temp <- rep(NA, times = 6)
    # n targets
    temp[1] <- n_targets[i]
    # Actual drug-disease distance
    drug_genes <- drug_target_matrix[,i]
    drug_genes <- drug_genes[!is.na(drug_genes)]
    if(length(drug_genes)>0){
      drug_genes <- as.character(unique(drug_genes))
      temp[2] <- average_closest_distance(ppin_dist, from = disease_genes, to = drug_genes)
    }
    # Mean random
    x <- as.numeric(random_bin_drugs[,i])
    temp[4] <- mean(x)
    # SD random
    temp[5] <- sd(x)
    # Z
    if(temp[2]!="NaN"){
      temp[3] <- (as.numeric(temp[2]) - as.numeric(temp[4]))/as.numeric(temp[5])
    } else {
      temp[3] <- NA
      temp[2] <- NA
    }
    # P
    temp[6] <- pnorm(as.numeric(temp[3]))
    # drug name
    temp <- c(colnames(drug_target_matrix)[i],temp)
    return(temp)
  }
  rm(n_targets, random_bin_drugs)
  colnames(out) <- c("Drug", "n_drug_targets", "dc", "zc", "mean(random dc)", "SD(random dc)", "P")
  write.table(out, file = paste("drug-disease_closest_distances_vs_random_bin_adjusted__",file_name,".txt", sep=""), sep = "\t", col.names = T, row.names = F)

  print("Network distance calculations = FINISHED")
  return(out)
}
