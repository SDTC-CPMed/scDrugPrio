#
# COMBINE CD PATIENTS DRUG WITH AVAILABLE LITERATURE SEARCH FOR FC CRITERIA CHECKING
#
#################################################################################

library(doParallel)
library(readxl)

# INPUT
#################################################################################
lit_data <- as.matrix(read.table(file = "../Input//Literature_search_drugs/Filtered_drug_effects_220704.txt",sep="\t", header = T, stringsAsFactors = F, quote = ""))
lit_data <- lit_data[,1:3]

# Responders
#################################################################################
temp <- as.matrix(read_xlsx(path = "../Output/FC_criteria_checking/Responder/SUMMARY_only_drugs_passing_network_criteria.xlsx", sheet = 1))

# now fill in information from previous literature search
pos <- match(temp[,1], lit_data[,1])
pos <- pos[!is.na(pos)]
temp[temp[,1]%in% lit_data[,1],c(18,19)] <- lit_data[pos,2:3]

# save
write.table(x = temp, file = "../Output/FC_criteria_checking/Responder/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt",sep= "\t", col.names = T, row.names = F)
out <- temp    


# Non-responders
#################################################################################
temp <- as.matrix(read_xlsx(path = "../Output/FC_criteria_checking/Non-responder/SUMMARY_only_drugs_passing_network_criteria.xlsx", sheet = 1))

# now fill in information from previous literature search
pos <- match(temp[,1], lit_data[,1])
pos <- pos[!is.na(pos)]
temp[temp[,1]%in% lit_data[,1],c(18,19)] <- lit_data[pos,2:3]

# save
write.table(x = temp, file = "../Output/FC_criteria_checking/Non-responder/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt",sep= "\t", col.names = T, row.names = F)

rm(list = ls())