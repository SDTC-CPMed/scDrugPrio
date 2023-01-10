#' Closest Distance Network Drug Screening
#'
#' Used for calculations of closest network proximity. Memory heavy.
#'
#' @param drugID Vector with IDs Matrix with 2 columns, every row represents one protein-protein interaction in the protein-protein interaction network (PPIN) largest connected component (LCC). Output from ppin_formatting()
#' @param target Vector with the same length as drugID that shall be sorted into different rows based on drugID. 
#' 
#' @return Matrix with length(unique(drugID)) columns. Each column specifies all targets for a certain drug. Column names specify which drugID the targets correspond to.
#' @export
#' 
create_drug_target_matrix <- function(drugID, target){
  library(doParallel)
  max <- max(table(drugID))
  unique_drugID <- unique(drugID)
  out <- foreach(i = c(1:length(unique_drugID)), .combine = "cbind") %do% {
    temp <- target[drugID %in% unique_drugID[i]]
    return(c(temp, rep(NA, times = max - length(temp))))
  }
  colnames(out) <- unique_drugID
  return(out)
}
