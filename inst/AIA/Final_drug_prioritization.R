
# FINAL DRUG RANKING
# Run on R 4.0.4
####################################################################################

#setwd("/data/samsc76/Doctis_data/R")
source("Codes_from_scPred/final_drug_prioritization.R")

dir.create("../Output/Final_ranking", showWarnings = F)

# Analysis
########################################

fc_eval <- as.matrix(read.table(file = "../Output/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt", sep = "\t", header = T, stringsAsFactors = F))
fc_eval[fc_eval == ""] <- NA
fc_eval <- fc_eval[rowSums(!is.na(fc_eval))>0, colSums(!is.na(fc_eval))>0]
fc_eval <- fc_eval[as.numeric(fc_eval[,12]) > 0,]
temp <- fc_eval[,19]
temp <- unlist(strsplit(temp, split = "FC_INDIVIDUAL_DRUGS_"))
temp <- unlist(strsplit(temp, split = "_res=0.6_dims=32_k=15"))
fc_eval[,19] <- temp

inter_cent <- read.table(file = "../Output/NicheNet/Cell_type_centrality_summary.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(inter_cent) <- inter_cent[,1]
inter_cent <- inter_cent[,-1]
inter_cents <- inter_cent[5,]
names(inter_cents) <- colnames(inter_cent)

intra_cents <- as.matrix(read.table(file = "../Output/Intracellular_centrality/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", stringsAsFactors = F, header = T))
rownames(intra_cents) <- intra_cents[,1]
intra_cents <- intra_cents[,-1]
intra_cents <- intra_cents[,ncol(intra_cents)]

drug_rank <- final_drug_prioritization(fc_evaluation = fc_eval, 
                                       pos_DrugID = 2,
                                       pos_clusterID = 19, 
                                       keep = c(1:4,18), 
                                       inter_cent = inter_cents, 
                                       intra_cent = intra_cents, 
                                       file_name = "FINAL_drug_ranking", 
                                       out_dir = "../Output/Final_ranking")

head(drug_rank)
