#
# COMBINE CD PATIENTS DRUG WITH AVAILABLE LITERATURE SEARCH FOR FC CRITERIA CHECKING
#
#################################################################################

fp <- getwd()
setwd(paste(fp, "/..",sep=""))

library(doParallel)
library(readxl)

# INPUT
#################################################################################
lit_data <- as.matrix(read.table(file = "Output/Literature_search_drugs/Filtered_drug_effects_220704.txt",sep="\t", header = T, stringsAsFactors = F, quote = ""))
lit_data <- lit_data[,1:3]

individual_cd_data <- as.matrix(read_xlsx(path = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/COMBINED_EVALUATION_FC_for_all_but_pat1_and_10_evaluated.xlsx", sheet = 1))


# Combine input files into new file
#################################################################################

individual_cd_data <- individual_cd_data[,c(1,16,19)]
individual_cd_data <- individual_cd_data[!is.na(individual_cd_data[,3]),]
individual_cd_data <- individual_cd_data[nrow(individual_cd_data):1,]
individual_cd_data <- individual_cd_data[!duplicated(individual_cd_data[,1]),]
colnames(individual_cd_data) <- NULL

lit_data <- rbind(lit_data[!(lit_data[,1] %in% individual_cd_data[,1]),], individual_cd_data)

write.table(lit_data, file ="Output/Literature_search_drugs/Filtered_drug_effects_221013.txt", sep="\t", col.names = T)

# Match with FC summary file
#################################################################################


temp <- as.matrix(read.table(file = "CD_batch_corrected/Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt", sep= "\t", stringsAsFactors = F, header = T))

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
write.table(x = temp, file = "CD_batch_corrected/Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion_evaluated.txt",sep= "\t", col.names = T, row.names = F)
out <- temp    

# ADDITIONAL LITERATURE SEARCH NEEDED TO CHECK FC CRITERIA?
out <- cbind(out, additional_lit_search_needed = "FALSE")
for(i in 1:nrow(out)){
  temp <- out[i,15]
  temp <- unlist(strsplit(temp, split = "), "))
  temp <- !grepl(pattern = "\\(NA", temp) # TRUE = a direct target DEG
  
  temp2 <- out[i,16]
  temp2 <- unlist(strsplit(temp2, split = "), "))
  temp2 <- !grepl(pattern = "\\(NA", temp2) # TRUE = we do know the drug action
  
  if(length(temp) == length(temp2)){
    temp <- temp - temp2 
    # 1 = DEG for which pharmacological action is not known
    # 0 = not a DEG or a DEG for which we know the drug action
    # -1 = not a DEG but we know the drug action 
    
    if(any(temp == 1)){
      out[i,20] <- "TRUE"
    }
  } else {
    out[i,20] <- "CAN NOT BE EVALUATED AUTOMATICALLY"
  }
}

# SAVE
write.table(x = out, file = "CD_batch_corrected/Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion_evaluated.txt",sep= "\t", col.names = T, row.names = F)

rm(list = ls())