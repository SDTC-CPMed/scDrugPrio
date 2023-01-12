#
#
# Prepare files for establishment of AIA MCDM ligand importance 
#
#

RA_drugs <- as.matrix(read.table(file = "../Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", stringsAsFactors = F, header = T, quote = ""))
drugs  <- as.matrix(read.table(file ="../Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", header = T, stringsAsFactors = F))
drugs <- drugs[drugs[,1] %in% RA_drugs[,1],6]
drugs <- unique(drugs)

dir.create(path = "../Output/Establishment_of_importance_of_MDCM_ligands")
write.table(matrix(drugs,ncol = 1), file = "../Output/Establishment_of_importance_of_MDCM_ligands/RA_drug_tragets.txt",sep="\t", col.names = "Drug_targets")


ligands <- read.table(file = "../Output/NicheNet/all_ligand_activity.txt", sep="\t", header = T, stringsAsFactors = F)
ligands <- ligands[as.numeric(ligands[,4])>0,]
write.table(ligands, file = "../Output/Establishment_of_importance_of_MDCM_ligands/all_ligand_activity.txt", sep="\t", col.names = T, row.names = F)