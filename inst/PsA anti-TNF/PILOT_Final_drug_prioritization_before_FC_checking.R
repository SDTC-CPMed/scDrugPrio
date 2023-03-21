
# PILOT FINAL DRUG RANKING
# Run on R 4.0.4
####################################################################################

#setwd("/data/samsc76/Doctis_data/R")
source("Codes_from_scPred/final_drug_prioritization.R")

dir.create("../Output/Final_ranking", showWarnings = F)

# Responder
########################################

fc_eval <- as.matrix(read.table(file = "../Output/FC_criteria_checking/Responder/SUMMARY_only_drugs_passing_network_criteria.txt", sep = "\t", header = T, stringsAsFactors = F))
# ADDITION:
fc_eval <- cbind(fc_eval[,1:11], NA, NA, NA, NA, NA, NA, fc_eval[,12:14])
#fc_eval <- fc_eval[fc_eval[,13] > 0,]

temp <- fc_eval[,20]
temp <- unlist(strsplit(temp, split = "FC_INDIVIDUAL_DRUGS_"))
temp <- unlist(strsplit(temp, split = "_R_res=0.8_k=15"))
fc_eval[,20] <- temp

inter_cent <- read.table(file = "../Output/NicheNet/Responder/Cell_type_centrality_summary.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(inter_cent) <- inter_cent[,1]
inter_cent <- inter_cent[,-1]
inter_cents <- inter_cent[5,]
names(inter_cents) <- colnames(inter_cent)

intra_cents <- as.matrix(read.table(file = "../Output/Intracellular_centrality/Responder/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", stringsAsFactors = F, header = T))
rownames(intra_cents) <- intra_cents[,1]
intra_cents <- intra_cents[,-1]
intra_cents <- intra_cents[,ncol(intra_cents)]

dir.create("../Output/Final_ranking/Responder", showWarnings = F)

drug_rank <- final_drug_prioritization(fc_evaluation = fc_eval, 
                                       pos_DrugID = 2,
                                       pos_clusterID = 20, 
                                       keep = c(1:4,19), 
                                       inter_cent = inter_cents, 
                                       intra_cent = intra_cents, 
                                       file_name = "FINAL_drug_ranking_Responder", 
                                       out_dir = "../Output/Final_ranking/Responder")

head(drug_rank)

temp <- drug_rank[drug_rank[,3] == T,]


# Non-responder
########################################

fc_eval <- as.matrix(read.table(file = "../Output/FC_criteria_checking/Non-responder/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt", sep = "\t", header = T, stringsAsFactors = F))
fc_eval <- fc_eval[fc_eval[,13] > 0,]
temp <- fc_eval[,21]
temp <- unlist(strsplit(temp, split = "FC_INDIVIDUAL_DRUGS_"))
temp <- unlist(strsplit(temp, split = "_NR_res=0.8_k=15"))
fc_eval[,21] <- temp

inter_cent <- read.table(file = "../Output/NicheNet/Non_responder/Cell_type_centrality_summary.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(inter_cent) <- inter_cent[,1]
inter_cent <- inter_cent[,-1]
inter_cents <- inter_cent[5,]
names(inter_cents) <- colnames(inter_cent)

intra_cents <- as.matrix(read.table(file = "../Output/Intracellular_centrality/Non_responder/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", stringsAsFactors = F, header = T))
rownames(intra_cents) <- intra_cents[,1]
intra_cents <- intra_cents[,-1]
intra_cents <- intra_cents[,ncol(intra_cents)]

dir.create("../Output/Final_ranking/Non_responder", showWarnings = F)

drug_rank <- final_drug_prioritization(fc_evaluation = fc_eval, 
                                       pos_DrugID = 2,
                                       pos_clusterID = 21, 
                                       keep = c(1:4,19), 
                                       inter_cent = inter_cents, 
                                       intra_cent = intra_cents, 
                                       file_name = "FINAL_drug_ranking_Non_responder", 
                                       out_dir = "../Output/Final_ranking/Non_responder")

head(drug_rank)
