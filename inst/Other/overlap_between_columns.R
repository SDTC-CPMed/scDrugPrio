# Overlap between columns
###############################
# Variable input:
###############################
# input = matrix; lists all DEGs for a cluster in a column.
#
# out = matrix; lists all DEGs found in 1, 2, 3, ... all clusters in a column.
overlapping_genes <- function(input){
  # make sure no duplicates are found in a column.
  for(i in c(1:ncol(input))){
    if(any(duplicated(input[,i]))){
      temp <- unique(input[,i])
      input[,i] <- NA
      input[1:length(temp),i] <- temp
    }
  }
  # calculate overlap
  out <- matrix(data = NA, nrow = nrow(input)*ncol(input), ncol = ncol(input))
  input <- table(input)
  for(i in c(1:ncol(out))){
    temp <- names(input[input >= i])
    if(length(temp)>0){
      out[1:length(temp),i] <- temp
    }
  }
  out <- out[rowSums(is.na(out))<ncol(out),]
  colnames(out) <- c("non-overlapping", paste("found_in_", 2:(ncol(out)-1),"_clusters", sep=""), "found_in_all")
  out <- out[,colSums(!is.na(out))>0]
  return(out)
}
