#' Calculation of closest distances between node degree adjusted, randomized drug target combinations and node degree adjusted, randomized disease genes for creation of a reference distribution
#'
#' @param bins bins used for network degree adjustment as previously obtained using bin_creation_by_min_bin_size()
#' @param drug_target_matrix Matrix listing each unique drug target combination in a column.
#' @param disease_genes Vector specifying disease associated genes in the PPIN.
#' @param ppin_dist Distances between all proteins included in the LCC of the PPIN.
#' @param seed Seed used for randomization of drug targets and disease genes. Needs to be altered for every iteration.
#'
#' @return Returns a vector containing average closest distances between every set of random (bin-adjusted) disease genes and random (bin-adjusted) drug targets.
#' @export
#'
random_drug_target_bin_adjusted_distances <- function(bins, drug_target_matrix, disease_genes, ppin_dist, seed){

  set.seed(seed)

  # which_bins_drugs -> columns = drugs, bins = rows
  which_bins_drugs <- foreach(k = c(1:ncol(drug_target_matrix)), .combine = cbind) %dopar% {
    # bins for drugs
    temp <- drug_target_matrix[,k]
    temp <- temp[!is.na(temp)]
    d.bins <- foreach(j = c(1:ncol(bins)), .combine = c) %do% {
      return(sum(temp%in%bins[,j]))
    }
    return(d.bins)
  }
  which_bins_disease <- vector()
  # bins for disease
  for(j in c(1:ncol(bins))){
    which_bins_disease <- c(which_bins_disease, sum(disease_genes %in% bins[!is.na(bins[,j]),j]))
  }

  # randomize disease model -> can be made more efficient by only randomizing disease_genes once for all drugs.
  random_disease <- vector()
  for(k in 1:ncol(bins)){
    if(which_bins_disease[k]>0){
      temp_bins <- bins[,k]
      temp_bins <- temp_bins[!is.na(temp_bins)]
      random_disease <-c(random_disease, sample(temp_bins, size = which_bins_disease[k], replace = F))
    }
  }
  if(length(random_disease)<length(disease_genes[!is.na(disease_genes)])){
    message("WARNING: Too few random disease genes!")
  }

  # randomize drug matrix & calculate average_closest_distance() between all 'random_drugs' and 'random_disease'
  if(length(disease_genes)>1){
    ppin_dist <- ppin_dist[rownames(ppin_dist)%in%random_disease,]
  }
  out <- foreach(k = c(1:ncol(which_bins_drugs)), .combine = c) %do% { # for every drug
    # create random drug
    random_drug <- vector()
    for(j in 1:ncol(bins)){
      if(which_bins_drugs[j,k]>0){
        temp_bins <- bins[,j]
        temp_bins <- temp_bins[!is.na(temp_bins)]
        random_drug <- c(random_drug, sample(temp_bins, size = which_bins_drugs[j,k], replace = F))
      }
    }
    if(length(random_drug)<sum(which_bins_drugs[,k])){
      message("WARNING: Too few random drug targets!")
    }
    return(average_closest_distance(ppin_dist, from = random_disease, to = random_drug))
  }
  return(out)
}
