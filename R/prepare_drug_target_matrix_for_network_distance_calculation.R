#' Preparation of drug target matrix for network distance analysis
#' Removes all drug targets not found in PPIN and thereafter identifies unique drug target combinations.
#' Unique drug target combinations are then used for calculation of network distances to save computational time and
#' can thereafter be translated back to individual drugs.
#'
#' @param ppin Matrix with 2 columns, every row represents one protein-protein interaction in the protein-protein interaction network (PPIN) largest connected component (LCC). Output from ppin_formatting()
#' @param drug_target_matrix Matrix listing each drug's targets in a column.
#' @param file_name Character string specifying filename for the output.
#' @param out_dir Character string specifying path to out directory where data will be saved.
#'
#' @return Returns drug_target_matrix including only unique drug target combinations. The matching between unique drug target combinations and
#'  their respective drugs are documented in paste(out_dir,"/SAME_DRUG_TARGETS_",file_name,".txt",sep="")
#' @export
#'
prepare_drug_target_matrix_for_network_distance_calculation <- function(ppin, drug_target_matrix, file_name, out_dir){
  
  if(out_dir != ""){
    if(substr(x = out_dir, start = nchar(out_dir), stop = nchar(out_dir)) != "/"){
      out_dir <- paste(out_dir, "/",sep="")
    }
  }
  
  ppin <- ppin_formatting(ppin) # extract LCC of PPIN
  all_genes <- unique(as.vector(ppin))

  # formatting drug matrix to only include drugs with at least one target found in ppin
  drug_target_matrix <- as.matrix(drug_target_matrix)
  for(i in c(1:ncol(drug_target_matrix))){
    temp <- as.character(drug_target_matrix[,i])
    drug_target_matrix[,i] <- NA
    temp <- unique(temp[temp%in%all_genes])
    if(length(temp)>0){
      drug_target_matrix[1:length(temp),i] <- temp
    }
  }
  rm(temp, i)
  drug_target_matrix <- drug_target_matrix[rowSums(!is.na(drug_target_matrix))>0,colSums(!is.na(drug_target_matrix))>0]

  # finding drugs with unique drug target combinations
  same_drugs <- foreach(i = 1:ncol(drug_target_matrix), .combine = rbind) %do% {
    temp <- unique(drug_target_matrix[!is.na(drug_target_matrix[,i]),i])
    out <- vector()
    for(j in 1:ncol(drug_target_matrix)){
      temp2 <- unique(drug_target_matrix[!is.na(drug_target_matrix[,j]),j])
      if(length(temp) == length(temp2)){
        if(all(temp %in% temp2)){
          out <- rbind(out, c(colnames(drug_target_matrix)[i], colnames(drug_target_matrix)[j]))
        }
      }
    }
    return(out)
  }
  colnames(same_drugs) <- c("Drug1", "Drug2")
  
  # select only unique drug target combinations
  remove <- vector()
  for(i in 1:(ncol(drug_target_matrix)-1)){
    temp <- drug_target_matrix[!is.na(drug_target_matrix[,i]),i]
    for(j in (i+1):ncol(drug_target_matrix)){
      temp2 <- drug_target_matrix[!is.na(drug_target_matrix[,j]),j]
      if(length(temp)==length(temp2)){
        if(all(temp %in% temp2)){
          remove <- c(remove, j)
        }
      }
    }
  }
  print(paste("n unique drugs included in analysis: n = ", ncol(drug_target_matrix),sep="")) # n = 1840
  unique_drug_target_matrix <- drug_target_matrix[,-unique(remove)]
  print(paste("n drugs with unique drug target combinations: n = ",ncol(unique_drug_target_matrix),sep="")) # n = 1204
  print(paste("n unique drugs that can be extracted from UNIQUE_DRUG_TARGET_COMBINATION_", file_name, ".txt using SAME_DRUG_TARGETS_",
              file_name,".txt: n = ", length(unique(c(same_drugs[same_drugs[,1]%in%colnames(unique_drug_target_matrix),2],
                                                      same_drugs[same_drugs[,2]%in%colnames(unique_drug_target_matrix),1]))),sep="")) # n = 1840
  
  # sorting drugs by n drug targets
  n_targets <- colSums(!is.na(unique_drug_target_matrix))
  unique_drug_target_matrix <- unique_drug_target_matrix[,order(n_targets, decreasing = T)]
  rm(n_targets)
  print("Drug matrix = FORMATTED")

  same_drugs <- same_drugs[same_drugs[,1] %in% colnames(drug_target_matrix),]
  write.table(same_drugs, file = paste(out_dir, "SAME_DRUG_TARGETS_", file_name, ".txt",sep=""), sep="\t", col.names = T, row.names = F)
  
  # Save file and return unique drug target combinations
  write.table(unique_drug_target_matrix, file = paste(out_dir, "UNIQUE_DRUG_TARGET_COMBINATION_", file_name, ".txt",sep=""), sep="\t", col.names = T, row.names = F)
  return(unique_drug_target_matrix)
}

