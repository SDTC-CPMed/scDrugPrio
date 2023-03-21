#
# Check which MS drugs in PPIN 
#
# 221030
# SAMUEL SCHAEFER
#################################################

setwd("C:/Users/samsc76/OneDrive - Linköpings universitet/Old RA mouse project - Barabasi/SIGMA_13012020/MS_batch_corrected/drug_prediction_R")

# INPUT
###############################

drugs <- read.table(file = "../Input/Drugs/MS_drugs_from_DrugBank.txt",sep="\t", header = T, stringsAsFactors = F)
drugs_in_ppin <- read.table(file ="../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt",sep="\t", header = T, stringsAsFactors = F)
unique_drug_targets <- read.table(file="../Input/Drugs/drug_targets_unique_literature_ppi.txt",sep="\t", header = T, stringsAsFactors = F)

drugs_in_ppin <- drugs_in_ppin[drugs_in_ppin[,1]%in% colnames(unique_drug_targets),]

included <- unique(drugs[drugs[,1] %in% drugs_in_ppin[,2],2])