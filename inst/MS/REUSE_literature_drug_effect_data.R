#
# COMBINE ALL INDIVIDUAL CD PATIENTS DRUG SUMMARY FOR FC CRITERIA CHECKING
#
#################################################################################

library(doParallel)

# INPUT
lit_data <- as.matrix(read.table(file = "../Input/Literature_search_drugs/Filtered_drug_effects_220704.txt",sep="\t", header = T, stringsAsFactors = F, quote = ""))
lit_data <- lit_data[,1:3]

# GET FC SUMMARY FOR MS DATA
temp <- as.matrix(read.table(file = "../Output/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria.txt",sep= "\t", stringsAsFactors = F, header = T))

# add additional columns + headers  
temp <- cbind(temp[,1:11], 
              n_drug_targets_counteracting_disease = NA, 
              n_drug_targets_mimmicking_disease = NA,
              counteracting_perc = NA,
              mimicking_perc = NA,
              outside_model_perc = NA,
              temp[,12:14],
              Additional_information_targets = NA)

# now fill in information from previous literature search
pos <- match(temp[,1], lit_data[,1])
pos <- pos[!is.na(pos)]
temp[temp[,1]%in% lit_data[,1],c(18,20)] <- lit_data[pos,2:3]

# save
write.table(x = temp, file = "../Output/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria_checked.txt",sep= "\t", col.names = T, row.names = F)

rm(list = ls())