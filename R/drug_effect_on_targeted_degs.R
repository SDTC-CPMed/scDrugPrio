# Matching pharmacological action and fold change of targeted DEGs to potential candidates drug targets
#' Title
#'
#' @param input matrix with 8 columns in ```mode(degs) <- "character"```, listing individual drug's network distances. Column 1 = DrugBank ID (character string), column 2 = known drug for this disease (boolean), column 3 = n of drug targets (integer), column 4 = dc (double), column 5 = zc (double), column 6 = mean random dc (double), column 7 = SD of random dc (double), column 8 = P value for significance of dc (double). Only column 1 is used for merging information in this function, though the matrix is reshuffled before return and therefore requires all 8 columns.
#' @param degs matrix with 6 columns in ```mode(degs) <- "character"```, listing DEGs associated with the input file. Prserves output format of Seurat based DEG calculation. Column 1 = gene name in annotation that matches 'drug_effect', column 3 = average log fold change (avg_logFC), column 6 = adjusted P value.
#' @param drug_effect matrix with 3 columns in ```mode(degs) <- "character"```, listing drug target information derived from DrugBank. Column 1 = DrugBankID used for matching with 'input', column 2 = gene name in annotation that matches with 'degs', column 3 = pharmacological action on the drug target.
#'
#' @return Returns matrix for manual evaluation of drug effect on targeted DEGs. Matches all matrices and determins fold change of targeted DEGs.
#' @export
#'
drug_effect_on_targeted_degs <- function(input, degs, drug_effect){
  if(!is.null(input)){
    # Format drug effect matrix
    out <- foreach(i = 1:nrow(input), .combine = rbind) %do% {
      temp <- input[i,1] # DrugBankID
      if(sum(drug_effect[,1]%in%temp)>1){
        temp <- drug_effect[drug_effect[,1]%in%temp,] # all drug targets
      } else {
        if(sum(drug_effect[,1]%in%temp)==0){
          warning(paste("No drug targets found for ", temp,"! Check for abnormality in gene annotation translation!",sep=""))
          return(c(NA, NA, NA))
        }
        temp <- matrix(drug_effect[drug_effect[,1] == temp,], nrow = 1) # all drug targets (n = 1)
      }

      # Remove duplicates:
      if(nrow(temp)>1){
        temp <- temp[!duplicated(temp[,2]),]
        if(is.null(nrow(temp))){
          if(length(temp)>0){
            temp <- matrix(temp, nrow = 1)
          } else {
            temp <- rep(NA, times = ncol(drug_effect))
          }
        }
      }

      # Effect on drug targets?
      t <- paste(as.character(temp[,2]), " (", as.character(temp[,3]),")", sep="")
      t <- paste(t,  collapse = ", ")

      # DEG being counteracted?
      t2 <- cbind(as.character(temp[,2]), NA, NA)
      t2[,2] <- degs[match(t2[,1],degs[,1]),3] # match deg gene names to drug target gene names and get avg_logFC values
      t2[as.numeric(t2[,2])>0,3] <- "up"
      t2[as.numeric(t2[,2])<0,3] <- "down"
      t3 <- sum(is.na(t2[,3])) # sum of drug targets that are not targeting DEGs
      if(sum(!is.na(t2[,3]))< 2){
        if(sum(!is.na(t2[,3])) == 1){
          t2 <- matrix(t2[!is.na(t2[,3]),], nrow = 1)
        } else {
          t2 <- matrix(NA, ncol = 3, nrow = 1)
        }
      } else {
        t2 <- t2[!is.na(t2[,3]),] # only include DEGs
      }
      t2 <- paste(as.character(t2[,1]), " (", t2[,3], ")", sep="")
      t2 <- paste(t2, collapse = ", ")
      return(c(t3, t2, t))
    }
    #colnames(out) <- c("n_drug_targets_that_are_non_DEGs", "Fold_change_of_drug_targets", "Pharmacological_effect_on_drug_targets")
    if(nrow(input) > 1){
      input <- cbind(input[,c(1,4,6,7,5,8,2,3)],out[,1],NA,NA,NA,NA,NA, out[,2:3])
      colnames(input) <- c("Drug", "dc", "mean_dc_random", "SD_random_dc", "zc", "P", "known_disease_specific_drug", "n_drug_targets_total",
                           "n_drug_targets_that_are_non_DEGs","n_drug_targets_counteracting_disease", "n_drug_targets_mimmicking_disease",
                           "Drug_targets_counteracting_fold_change_DEGs_%","Drug_targets_mimmicking_fold_change_DEGs_%", "Drug_targets_targeting_non_DEGs_%",
                           "Fold_change_of_drug_targets", "Pharmacological_effect_on_drug_targets")

    } else if (is.vector(input) || nrow(input) == 1){
      input <- c(input[c(1,4,6,7,5,8,2,3)],out[1],NA,NA,NA,NA,NA, out[2:3])
      input <- matrix(input, nrow = 1)
      colnames(input) <- c("Drug", "dc", "mean_dc_random", "SD_random_dc", "zc", "P", "known_disease_specific_drug", "n_drug_targets_total",
                           "n_drug_targets_that_are_non_DEGs","n_drug_targets_counteracting_disease", "n_drug_targets_mimmicking_disease",
                           "Drug_targets_counteracting_fold_change_DEGs_%","Drug_targets_mimmicking_fold_change_DEGs_%", "Drug_targets_targeting_non_DEGs_%",
                           "Fold_change_of_drug_targets", "Pharmacological_effect_on_drug_targets")
    }
    return(input)
  } else {
    warning("'input' is NULL")
    return(NULL)
  }
}
