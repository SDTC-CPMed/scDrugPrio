#' Bin creation for network degree adjustments
#'
#' @param ppin Matrix with 2 columns, every row represents one protein-protein interaction in the protein-protein interaction network (PPIN) largest connected component (LCC). Output from ppin_formatting()
#' @param min_bin_size Integer. The bin-size used for node degree adjustment of network calculations.
#'
#' @return
#' @export
#'
bin_creation_by_min_bin_size <- function(ppin, min_bin_size){

  # calculate node degree for each protein/gene
  n_temp <- table(as.vector(ppin))
  n_temp <- sort(n_temp)

  # create bins
  bins <- matrix(NA, nrow = length(n_temp), ncol = length(n_temp)/min_bin_size)
  bins_col <- vector()
  degree <- unique(n_temp)
  degree <- degree[order(degree, decreasing = F)]
  temp <- vector()
  b <- 1
  for(i in c(1:length(degree))){
    temp <- c(temp, names(n_temp[n_temp == degree[i]]))
    if(length(temp)>=min_bin_size){
      bins[1:length(temp),b] <- temp
      bins_col <- c(bins_col, i)
      b <- b+1
      temp <- vector()
    }
  }
  if(length(temp)>0 & length(temp)< min_bin_size){ # in case last bin cant be filled up it is added to previous bin
    b <- b-1
    temp2 <- bins[,b]
    temp2 <- unique(temp2)
    temp2 <- temp2[!is.na(temp2)]
    temp2 <- c(temp2, temp)
    bins[1:length(temp2),b] <- temp2
  }
  rm(temp, temp2, b, degree, bins_col)
  bins <- bins[,colSums(!is.na(bins))>0]
  bins <- bins[rowSums(!is.na(bins))>0,]
  return(bins)
}
