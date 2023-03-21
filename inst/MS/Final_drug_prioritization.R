
# FINAL DRUG RANKING
# Run on R 4.0.4
####################################################################################

#setwd("/data/samsc76/MS_batch_corrected/R")
source("Codes_from_scPred/final_drug_prioritization.R")

dir.create("../Output/Final_ranking", showWarnings = F)

# Analysis
########################################

fc_eval <- as.matrix(read.table(file = "../Output/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt", sep = "\t", header = T, stringsAsFactors = F, quote = ""))
fc_eval <- fc_eval[as.numeric(fc_eval[,12]) > 0,]
temp <- fc_eval[,19]
temp <- unlist(strsplit(temp, split = "FC_INDIVIDUAL_DRUGS_"))
temp <- unlist(strsplit(temp, split = "_res=0.35_dims=32_k=15"))
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

dir.create("../Output/Final_ranking", showWarnings = F)

n_clusters <- table(fc_eval[,1]) 
fc_eval <- cbind(fc_eval, n_clusters = n_clusters[match(fc_eval[,1], names(n_clusters))])

drug_rank <- final_drug_prioritization(fc_evaluation = fc_eval, 
                                       pos_DrugID = 2,
                                       pos_clusterID = 19, 
                                       keep = c(1:4,18,21), 
                                       inter_cent = inter_cents, 
                                       intra_cent = intra_cents, 
                                       file_name = "FINAL_drug_ranking", 
                                       out_dir = "../Output/Final_ranking")

head(drug_rank)
print("Precision in top 100 ranked drugs:")
print(sum(drug_rank[as.numeric(drug_rank[,9]) <= 100,3] == "TRUE") * 100 / nrow(drug_rank[as.numeric(drug_rank[,9]) <= 100,]))
print("n captured MS drugs overall:")
print(sum(drug_rank[,3] == "TRUE"))

# combine with previous literature knowledge
########################################

prev_prediction <- read.table(file = "../../Output/MS/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_MS.txt", sep = "\t", header = T, stringsAsFactors = F)
print(paste("Approved MS drugs n = ", sum(prev_prediction[,4] == "MS_drug"),sep=""))
print(paste("Total n candidates = ", nrow(prev_prediction),sep=""))

prev_prediction <- as.matrix(prev_prediction[, c(1,3,2,44:47)])
colnames(prev_prediction)[3] <- "non_batch_corrected_rank"
prev_prediction[prev_prediction == ""] <- NA

drug_rank <- cbind(drug_rank, prev_prediction[match(drug_rank[,1], prev_prediction[,1]) ,3:ncol(prev_prediction)])

print("Precision in top 100 ranked drugs:")
print(sum(drug_rank[as.numeric(drug_rank[,9]) <= 100 & !is.na(drug_rank[,12]),12] == "Yes") * 100 / nrow(drug_rank[as.numeric(drug_rank[,9]) <= 100,]))

write.table(drug_rank, file = "../Output/Final_ranking/FINAL_drug_ranking_with_previous_lit_evidence.txt", sep="\t", col.names = T, row.names = F)
