#
# COMBINE ALL INDIVIDUAL CD PATIENTS DRUG SUMMARY FOR FC CRITERIA CHECKING
#
#################################################################################

fp <- getwd()
setwd(paste(fp, "/..",sep=""))

library(doParallel)

# INPUT
lit_data <- as.matrix(read.table(file = "Output/Literature_search_drugs/Filtered_drug_effects_220704.txt",sep="\t", header = T, stringsAsFactors = F, quote = ""))
lit_data <- lit_data[,1:3]

# GET SUMMARY FOR EVERY PATIENT (except 1 & 10 which are already done)
#################################################################################
out <- foreach(i = c(2:9,11), .combine = "rbind") %do% {
  
  #temp <- as.matrix(read.table(file = paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_",i,"/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion_patient_",i,".txt", sep=""),
  #                   sep= "\t", stringsAsFactors = F, header = T, quote = ""))
  
  if(file.exists(paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_",i,"/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt", sep=""))){
    temp <- as.matrix(read.table(file = paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_",i,"/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt", sep=""),
                                 sep= "\t", stringsAsFactors = F, header = T))
    temp_col_name <- unlist(strsplit(colnames(temp), split = "X."))
    temp_col_name <- unlist(strsplit(temp_col_name, split = "\\.$"))
    colnames(temp) <- temp_col_name
    
    # add additional columns + headers  
    temp <- cbind(temp[,1:9], 
                  n_drug_targets_mimmicking_disease = NA, 
                  n_drug_targets_mimmicking_disease = NA,
                  counteracting_perc = NA,
                  mimicking_perc = NA,
                  outside_model_perc = NA,
                  temp[,10:13],
                  Additional_information_targets = NA)
    
    # now fill in information from previous literature search
    pos <- match(temp[,1], lit_data[,1])
    pos <- pos[!is.na(pos)]
    temp[temp[,1]%in% lit_data[,1],c(16,19)] <- lit_data[pos,2:3]
    
    # save
    write.table(x = temp, file = paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_",i,"/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion_patient_",i,".txt", sep=""),
                sep= "\t", col.names = T, row.names = F)
    
    # return
    return(cbind(temp, patient = paste("Patient_",i,sep="")))
  } else {
    print(paste("DOES NOT EXIST: i = ",i,sep=""))
    return()
  }
}

# ADDITIONAL LITERATURE SEARCH NEEDED TO CHECK FC CRITERIA?
#################################################################################
out <- cbind(out, calculated_n_of_direct_DEGs = NA, additional_lit_search_needed = "FALSE")
for(i in 1:nrow(out)){
  temp <- out[i,15]
  temp <- unlist(strsplit(temp, split = "), "))
  temp <- !grepl(pattern = "\\(NA", temp) # TRUE = a direct target DEG
  out[i,21] <- sum(temp)
  
  temp2 <- out[i,16]
  temp2 <- unlist(strsplit(temp2, split = "), "))
  temp2 <- !grepl(pattern = "\\(NA", temp2) # TRUE = we do know the drug action
  
  if(length(temp) == length(temp2)){
    temp <- temp - temp2 
    # 1 = DEG for which pharmacological action is not known
    # 0 = not a DEG or a DEG for which we know the drug action
    # -1 = not a DEG but we know the drug action 
    
    if(any(temp == 1)){
      out[i,22] <- "TRUE"
    }
  } else {
    out[i,22] <- "CAN NOT BE EVALUATED AUTOMATICALLY"
  }
}

# SAVE
write.table(out, file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/COMBINED_EVALUATION_FC_for_all_but_pat1_and_10.txt",sep="\t", col.names = T, row.names = F)

rm(list = ls())