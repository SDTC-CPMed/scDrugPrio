#' Final drug prioritization 
#'
#' Ranking of candidates passing selection criteria. Candidates are primarily ranked by the sum of the intercellular centrality scores of clusters in which they are considered candidates. Secondarily, candidates are ranked by their cluster-wide average intracellular centrality.
#'
#' @param fc_evaluation Matrix with 2 columns, every row represents one protein-protein interaction in the protein-protein interaction network (PPIN) largest connected component (LCC). Output from ppin_formatting()
#' @param pos_clusterID Integer specifying the column in which the cluster IDs / cell type labels can be found that are to be matched with names in inter_cent.
#' @param pos_DrugID Integer specifying the column in which the cluster IDs / cell type labels can be found that are to be matched with names in inter_cent.  
#' @param keep Vector specifying which columns (as indicated by column number) to keep in the final drug ranking file. 
#' @param inter_cent Vector specifying the centrality score of each cluster. The name of each element indicates which cluster/cell type the centrality score belongs to. Names must be matchable with cluster identifiers in 'fc_evaluation'.  
#' @param intra_cent Vector specifying the cluster-wide average of the intracellular centralities. Names indicate which drug the centrality score belongs to. Names must be matchable with column 1 in 'fc_evaluation'.
#' @param file_name Character string specyfing the file name by which the final ranking will be saved.
#' @param out_dir Character string specifying path to out directory where the final ranking will be saved.
#'
#' @return a matrix
#' @export
#'
final_drug_prioritization <- function(fc_evaluation,
                                      pos_DrugID = 1,
                                      pos_clusterID = ncol(fc_evaluation),
                                      keep = c(pos_DrugID, pos_clusterID),
                                      inter_cent,
                                      intra_cent,
                                      file_name,
                                      out_dir = ""){
  
  
  if(out_dir != ""){
    if(substr(x = out_dir, start = nchar(out_dir), stop = nchar(out_dir)) != "/"){
      out_dir <- paste(out_dir, "/",sep="")
    }
  }
  
  # CHECK IF BASIC REQUIREMENTS ARE FULFILLED
  ######################################################################
  
  if(is.null(fc_evaluation)){
    stop(
      "Variable \'fc_evaluation\' not specified.",
      call. = FALSE
    )
  }
  if(is.null(pos_clusterID)){
    stop(
      "Variable \'pos_clusterID\' not specified.",
      call. = FALSE
    )
  } else {
    if(!is.numeric(pos_clusterID)){
      stop(
        "Variable \'pos_clusterID\' is not an integer",
        call. = FALSE
      )
    }
  }
  if(is.null(pos_DrugID)){
    stop(
      "Variable \'pos_DrugID\' not specified.",
      call. = FALSE
    )
  } else {
    if(!is.numeric(pos_DrugID)){
      stop(
        "Variable \'pos_DrugID\' is not an integer",
        call. = FALSE
      )
    }
  }
  if(is.null(keep)){
    stop(
      "Variable \'keep\' not specified.",
      call. = FALSE
    )
  } else {
    if(!is.vector(keep)){
      if(!is.numeric(keep)){
        stop(
          "Variable \'keep\' is not numeric",
          call. = FALSE
        )
      }
    } else {
      temp <- vector()
      for(i in keep){
        temp <- c(temp, !is.numeric(i))
      }
      if(any(temp)){
        stop(
          "Variable \'keep\' is not numeric",
          call. = FALSE
        )
      }
    }
  }
  if(is.null(inter_cent)){
    stop(
      "Variable \'inter_cent\' not specified.",
      call. = FALSE
    )
  }
  if(is.null(intra_cent)){
    stop(
      "Variable \'intra_cent\' not specified.",
      call. = FALSE
    )
  }
  if(is.null(file_name)){
    stop(
      "Variable \'file_name\' not specified.",
      call. = FALSE
    )
  }
  
  # MATCH DATA
  ######################################################################
  
  # Make sure to put drug ID in first column
  keep <- c(pos_DrugID, keep)
  
  # Calculate inter_cent for each drug
  temp <- fc_evaluation[, c(pos_DrugID, pos_clusterID)]
  unique_drugs <- unique(temp[,1])
  temp2 <- vector()
  for(i in 1:length(unique_drugs)){
    clust <- temp[temp[,1] %in% unique_drugs[i],2]
    temp2[i] <- sum(inter_cent[names(inter_cent)%in% clust])
  }
  inter_cent <- cbind(unique_drugs, temp2)
  
  # Select only unique drugs from FC evaluation file
  out <- fc_evaluation[,keep]
  out <- out[!duplicated(out[,1]),]
  
  # Match with inter_cent
  out <- cbind(out, Inter_cent = inter_cent[match(out[,1], inter_cent[,1]),2])
  
  # Match with intra_cent 
  out <- cbind(out, Intra_cent = intra_cent[match(out[,1], names(intra_cent))])
  
  # Final score
  out <- cbind(out, combined_centrality_score = as.numeric(out[,(ncol(out)-1)])+10^-1*as.numeric(out[,ncol(out)]))
  out <- out[order(as.numeric(out[,ncol(out)]), decreasing = T),]
  
  # Rank
  out <- cbind(out, rank = rank(1000 - as.numeric(out[,ncol(out)]), ties.method = "average"))
  out <- out[,-1]
  
  write.table(out, file = paste(out_dir, file_name, ".txt",sep=""), sep="\t", col.names = T, row.names = F)
  
  return(out)
}