#' Average closest distance
#'
#' Calculate average closest network distance (d(c)) between drug targets and the nearest disease protein
#' Proven useful for drug efficacy screening in Guney et al,  NatCom, 2016 (PMID: 26831545)
#'
#' @param ppin_dist Distance matrix between all proteins included in the largest connected component (LCC) of the protein-protein interaction network (PPIN).
#' @param from Character vector specifying genes between which distances should be calculated. E.g.: from = disease_genes
#' @param to Character vector specifying genes between which distances should be calculated. E.g.: to = drug_targets
#'
#' @return Double specifying the average closest network distance between elements in 'from' and elements in 'to'
#' @export
#'
average_closest_distance <- function(ppin_dist, from, to){
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop(
      "Package \'matrixStats\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(matrixStats)
  
  if(any(length(from)==0, length(to)==0)){
    warning("Empty vector supplied for from / to of function average_closest_distance(ppin_dist, from, to)", call. = FALSE)
    return(NA)
  }
  if(length(from)>1){
    ppin_dist <- ppin_dist[rownames(ppin_dist)%in%from,]

    if(length(to)>1){ # > 1 diseas protein & > 1 target protein
      ppin_dist <- ppin_dist[,colnames(ppin_dist)%in%to]
      return(mean(colMins(ppin_dist))) # average closest path

    } else { # only 1 target protein
      ppin_dist <- ppin_dist[,colnames(ppin_dist)%in%to]
      return(min(ppin_dist)) # closest path - no averaging needed
    }
  } else { # only 1 disease protein!
    if(length(to)>1){
      ppin_dist <- ppin_dist[,colnames(ppin_dist) %in% to]
      ppin_dist <- ppin_dist[rownames(ppin_dist)%in%from,]
      return(mean(ppin_dist)) # average closest path
    } else { # only one drug target & 1 disease protein
      return(ppin_dist[rownames(ppin_dist)%in% from, colnames(ppin_dist) %in% to])
    }
  }
}
