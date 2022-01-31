#' Data formatting of basic files used for network calculations
#'
#' @return Creates a drug target matrix based on DrugBank input
#' @export
#'
#' @examples data_analysis_and_formatting()
#'
data_analysis_and_formatting <- function(){

  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop(
      "Package \'doParallel\' must be installed to use this function.",
      call. = FALSE
    )
  }

  library(doParallel)


  # DrugBank filtering and creation of drug_target_matrix (DrugBank_drug_target_matrix.txt)
  #########################################################

  if(!file.exists("data/DrugBank_drug_target_matrix.txt")){
    transl <- as.matrix(read.table(file = "data-raw/HGNC_transl.txt", sep="\t", header = T, stringsAsFactors = F))
    transl <- transl[,c(2,6)]
    transl <- transl[!is.na(transl[,1]),]
    transl <- transl[!is.na(transl[,2]),]
    transl[,2] <- as.numeric(transl[,2])
    drugs <- as.matrix(read.table(file = "data-raw/all_drug_targets_drug_bank.txt", sep="\t", header = T, stringsAsFactors = F)) # n drugs = 13339
    drugs <- drugs[grepl(pattern = "approved", drugs[,3]),] # only approved drugs; n drugs = 4021
    drugs <- drugs[grepl(pattern = "Humans", drugs[,8]),] # approved for humans; n drugs = 1964
    drugs <- drugs[!is.na(drugs[,6]),] # needs known drug targets; n drugs = 1864
    drugs[,5] <- NA
    drugs[,8] <- NA
    drugs[,5] <- transl[match(drugs[,6], transl[,1]),2]
    drugs <- drugs[!is.na(drugs[,5]),] # only drugs with drug targets that can be translated to Entrez ID; n drugs = 1844

    # create drug target matrix
    unique_drugs <- unique(drugs[,1])
    drug_target_matrix <- foreach(i = 1:length(unique(drugs[,1])), .combine = cbind) %do% {
      temp <- unique(drugs[drugs[,1]%in% unique_drugs[i],5]) # all drug targets for a given drug
      temp <- matrix(c(temp, rep(NA, times = 1000-length(temp))), ncol = 1)
      colnames(temp) <- unique_drugs[i]
      return(temp)
    }
    drug_target_matrix <- drug_target_matrix[rowSums(!is.na(drug_target_matrix))>0,]
    write.table(drug_target_matrix, file = "data/DrugBank_drug_target_matrix.txt", sep="\t",col.names = T, row.names = F)
    rm(unique_drugs, temp, drugs, transl, i, drug_target_matrix)
  }

  #
  #########################################################
}
