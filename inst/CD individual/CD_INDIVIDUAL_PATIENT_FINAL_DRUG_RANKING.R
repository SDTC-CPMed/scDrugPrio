
#######################################################################################
# Drug ranking script for processing after all selection criteria has been run
# FOR CD DATA SET
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/../", sep=""))
library(doParallel)

##############################################################################################################################################################################
# PATIENT 10
##############################################################################################################################################################################

#######################################################################################
# LOAD files
#######################################################################################
# File with all drugs passing mean d(c) and P value criteria used for manually controlling FC-criteria
fc_criteria <- read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion_patient_10.txt",
                          sep="\t", quote = "\\", stringsAsFactors = F, fill =F, header = T)
fc_criteria <- as.matrix(fc_criteria)
fc_criteria <- fc_criteria[,-c((ncol(fc_criteria)-1):ncol(fc_criteria))]
fc_criteria[fc_criteria[,10]=="",10] <- NA
fc_criteria <- fc_criteria[fc_criteria[,1]!="",] # removes empty rows at bottom of excel files

# DrugBank
drugbank <- as.matrix(read.table(file = "Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", stringsAsFactors = F, header = T))

# Drug ID matching with unique drug target combination used in prediction
drugID_matching <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F))

# Unique drug target combination matrix
unique_target_combinations <- as.matrix(read.table(file = "Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F))
unique_target_combinations <- colnames(unique_target_combinations)

# Centrality
centrality_scores <- as.matrix(read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/NicheNet/Cell-cell_centrality_summary_CD_P10.txt", sep="\t", header = T, stringsAsFactors = F))
centrality_scores <- centrality_scores[5,-1] #choosing only eigenvector centrality
temp_names <- names(centrality_scores)
temp_names <- unlist(strsplit(temp_names, split = "_CD_P10")) # perhaps not needed for other data sets
centrality_scores <- as.numeric(centrality_scores)
names(centrality_scores) <- temp_names
rm(temp_names)

#######################################################################################
# Process drugs that pass FC criteria
#######################################################################################
temp <- fc_criteria[is.na(fc_criteria[,10]),]
print("FC-criteria could not be checked in all cluster for:")
print(unique(temp[,1])) # 0 drugs, 0 clusters
final_drugs <- fc_criteria[!is.na(fc_criteria[,10]),] # remove drugs for which FC-criteria could not be checked
final_drugs <- final_drugs[as.numeric(final_drugs[,10])>0,]
print("n of drugs for which FC-criteria could be checked in some clusters but not all:")
print(sum(unique(temp) %in% unique(final_drugs[,1]))) # 0 drugs
print(intersect(unique(temp), unique(final_drugs[,1])))

#######################################################################################
# Creation of output file
#######################################################################################

# following columns:
# drug name + n drug targets + known CD drug?
#######################################################################

unique_drugs <- final_drugs[!duplicated(final_drugs[,1]),c(1,7,8)]
# DrugBank ID
#######################################################################

drug_IDs <- drugbank[match(unique_drugs[,1],drugbank[,2]),1:2]
#print(cbind(unique_drugs[,1], drug_IDs))
print(paste("Complete matching?", any(!is.na(drug_IDs[,1])), sep=" "))
drug_IDs <- drug_IDs[,1]
unique_drugs <- cbind(unique_drugs[,1], drug_IDs, unique_drugs[,2:3])
colnames(unique_drugs)[1:2] <- c("Drug", "DrugBank_ID")

# mean closest distance to each cluster, one column for every cluster
#######################################################################

# detect which clusters had DEGs
file_pattern <- "drug-disease_closest_distances_vs_random_bin_adjusted__CD_DCA_MAST_DEGs_literature_PPI_"
file_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/literature_PPI/"
lf <- list.files(path = file_path, pattern = file_pattern)
lf <- lf[!grepl(lf, pattern = "_top_")]
clusters <- unlist(strsplit(lf, split = file_pattern))
clusters <- unlist(strsplit(clusters, split= ".txt"))

# create matrix for matching drug ID with drug ID used for the unique drug target combination in the network prediction
match_IDs <- cbind(drug_IDs, NA)
colnames(match_IDs) <- c("Actual_drugID", "DrugID_for_unique_target_combination")
for(i in 1:nrow(match_IDs)){
  temp <- unique(c(drugID_matching[drugID_matching[,2]==match_IDs[i,1],1], drugID_matching[drugID_matching[,1]==match_IDs[i,1],2]))
  temp <- temp[temp %in% unique_target_combinations]
  if(length(temp)>1){
    print(match_IDs[i,1])
    print("Too many matched drugs!!")
  } else {
    match_IDs[i,2] <- temp
  }
}

# create output matrix
mean_dist <- matrix(NA, nrow = nrow(unique_drugs), ncol = length(clusters))
rownames(mean_dist) <- unique_drugs[,1]
colnames(mean_dist) <- clusters

for(i in 1:length(clusters)){ # once for every cluster/column
  # given that we use DEG thresholds the correct file has to be chosen for mean network distances
  lf <- list.files(path = file_path)
  lf <- lf[grepl(lf, pattern = file_pattern)]
  lf1 <- lf[grepl(lf, pattern = paste(clusters[i],".txt",sep=""))] # protects for example from Cluster_3 being matched with files of Cluster_31
  lf2 <- lf[grepl(lf, pattern = paste(clusters[i],"_",sep=""))]
  lf <- unique(c(lf1,lf2))
  if(any(grepl(lf,pattern = "_top_3500"))){
    lf <- lf[grepl(lf, pattern = "top_3500")]
  } else {
    lf <- lf[!grepl(lf, pattern = "_top_")]
  }
  pred_mat <- as.matrix(read.table(file = paste(file_path, lf, sep=""), sep="\t", header = T, stringsAsFactors = F))  
  pred_mat <- pred_mat[pred_mat[,1]%in%match_IDs[,2],]
  mean_dist[,i] <- pred_mat[match(match_IDs[,2],pred_mat[,1]),3]
}
colnames(mean_dist) <- paste("mean_dc_", colnames(mean_dist),sep="")

rm(pred_mat, i, match_IDs, file_pattern)

# mean of mean closest distance
#######################################################################

mean_dist <- as.matrix(mean_dist)
mode(mean_dist) <- "numeric"
mean_dist <- cbind(mean_dist, rowMeans(mean_dist))
colnames(mean_dist)[ncol(mean_dist)] <- "mean_of_mean_dc"

# n of clusters in which drug fulfills selection criteria
#######################################################################

n_clusters_fulfilled <- table(final_drugs[,1])
n_clusters_fulfilled <- n_clusters_fulfilled[match(unique_drugs[,1],names(n_clusters_fulfilled))]

# fullfills selection criteria in which clusters? one column for every cluster, indicated by TRUE or FALSE if criteria is passed
# + combined Eigenvector centrality score for all clusters in which selection criteria is passed
#######################################################################

fulfilled_criteria <- matrix(F, nrow = nrow(unique_drugs), ncol = length(clusters))
rownames(fulfilled_criteria) <- unique_drugs[,1]
colnames(fulfilled_criteria) <- clusters

combined_centrality <- vector()
centrality_scores <- centrality_scores[match(colnames(fulfilled_criteria), names(centrality_scores))]

for(i in 1:nrow(fulfilled_criteria)){
  temp <- final_drugs[final_drugs[,1] %in% rownames(fulfilled_criteria)[i],17] # extract in which clusters drug fulfilled criteria
  temp <- as.numeric(temp)
  temp <- paste("Cluster_", temp, sep="")
  fulfilled_criteria[i,colnames(fulfilled_criteria)%in% temp] <- T
  
  combined_centrality[i] <- sum(centrality_scores[names(centrality_scores)%in% temp])
}
mode(fulfilled_criteria) <- "character"
colnames(fulfilled_criteria) <- paste("Passes_selection_criteria_",colnames(fulfilled_criteria),sep="")
fulfilled_criteria <- cbind(fulfilled_criteria, n_clusters_fulfilled, combined_centrality)
colnames(fulfilled_criteria)[(ncol(fulfilled_criteria)-1):ncol(fulfilled_criteria)] <- c("Fulfilled_criteria_in_n_clusters", "Sum_of_targeted_clusters_Eigenvector_centrality")
rm(combined_centrality, centrality_scores, temp, i)

# COMBINE ALL MATRICES INTO 1
#######################################################################
out <- cbind(unique_drugs, mean_dist, fulfilled_criteria)
out <- out[order(as.numeric(out[,ncol(out)]), as.numeric(out[,(5+length(clusters))]), decreasing = T),] # ordered by 1) combined Eigenvector centrality and 2) cluster-wide mean of mean d(c)
write.table(out, file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_INDIVIDUAL_CD_PATIENT_10.txt", sep="\t", col.names = T, row.names = F)

setwd(fp)
rm(list = ls())


##############################################################################################################################################################################
# PATIENT 1
##############################################################################################################################################################################

fp <- getwd()
setwd(paste(fp, "/../", sep=""))
library(doParallel)

#######################################################################################
# LOAD files
#######################################################################################
# File with all drugs passing mean d(c) and P value criteria used for manually controlling FC-criteria
fc_criteria <- read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/literature_PPI/SUMMARY/COMBINED_EVALUATION_FILES_for_checking_FC_criterion_patient_1.txt",
                          sep="\t", quote = "\\", stringsAsFactors = F, fill =T, header = T)
fc_criteria <- as.matrix(fc_criteria)
fc_criteria <- fc_criteria[,-c((ncol(fc_criteria)-1):ncol(fc_criteria))]
fc_criteria[fc_criteria[,10]=="",10] <- NA
fc_criteria <- fc_criteria[fc_criteria[,1]!="",] # removes empty rows at bottom of excel files

# DrugBank
drugbank <- as.matrix(read.table(file = "Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", stringsAsFactors = F, header = T))

# Drug ID matching with unique drug target combination used in prediction
drugID_matching <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F))

# Unique drug target combination matrix
unique_target_combinations <- as.matrix(read.table(file = "Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F))
unique_target_combinations <- colnames(unique_target_combinations)

# Centrality
centrality_scores <- as.matrix(read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/NicheNet/Cell-cell_centrality_summary_CD_P1.txt", sep="\t", header = T, stringsAsFactors = F))
centrality_scores <- centrality_scores[5,-1] #choosing only eigenvector centrality
temp_names <- names(centrality_scores)
temp_names <- unlist(strsplit(temp_names, split = "_CD_P1")) # perhaps not needed for other data sets
centrality_scores <- as.numeric(centrality_scores)
names(centrality_scores) <- temp_names
rm(temp_names)

#######################################################################################
# Process drugs that pass FC criteria
#######################################################################################
temp <- fc_criteria[is.na(fc_criteria[,10]),]
print("FC-criteria could not be checked in all cluster for:")
print(unique(temp[,1])) # 0
final_drugs <- fc_criteria[!is.na(fc_criteria[,10]),] # remove drugs for which FC-criteria could not be checked
final_drugs <- final_drugs[as.numeric(final_drugs[,10])>0,]
print("n of drugs for which FC-criteria could be checked in some clusters but not all:")
print(sum(unique(temp) %in% unique(final_drugs[,1]))) # 0 drugs
print(intersect(unique(temp), unique(final_drugs[,1])))

#######################################################################################
# Creation of output file
#######################################################################################

# following columns:
# drug name + n drug targets + known CD drug?
#######################################################################

unique_drugs <- final_drugs[!duplicated(final_drugs[,1]),c(1,7,8)]
# DrugBank ID
#######################################################################

drug_IDs <- drugbank[match(unique_drugs[,1],drugbank[,2]),1:2]
#print(cbind(unique_drugs[,1], drug_IDs))
print(paste("Complete matching?", any(!is.na(drug_IDs[,1])), sep=" "))
drug_IDs <- drug_IDs[,1]
unique_drugs <- cbind(unique_drugs[,1], drug_IDs, unique_drugs[,2:3])
colnames(unique_drugs)[1:2] <- c("Drug", "DrugBank_ID")

# mean closest distance to each cluster, one column for every cluster
#######################################################################

# detect which clusters had DEGs
file_pattern <- "drug-disease_closest_distances_vs_random_bin_adjusted__CD_DCA_MAST_DEGs_literature_PPI_"
file_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/literature_PPI/"
lf <- list.files(path = file_path, pattern = file_pattern)
lf <- lf[!grepl(lf, pattern = "_top_")]
clusters <- unlist(strsplit(lf, split = file_pattern))
clusters <- unlist(strsplit(clusters, split= ".txt"))

# create matrix for matching drug ID with drug ID used for the unique drug target combination in the network prediction
match_IDs <- cbind(drug_IDs, NA)
colnames(match_IDs) <- c("Actual_drugID", "DrugID_for_unique_target_combination")
for(i in 1:nrow(match_IDs)){
  temp <- unique(c(drugID_matching[drugID_matching[,2]==match_IDs[i,1],1], drugID_matching[drugID_matching[,1]==match_IDs[i,1],2]))
  temp <- temp[temp %in% unique_target_combinations]
  if(length(temp)>1){
    print(match_IDs[i,1])
    print("Too many matched drugs!!")
  } else {
    match_IDs[i,2] <- temp
  }
}

# create output matrix
mean_dist <- matrix(NA, nrow = nrow(unique_drugs), ncol = length(clusters))
rownames(mean_dist) <- unique_drugs[,1]
colnames(mean_dist) <- clusters

for(i in 1:length(clusters)){ # once for every cluster/column
  # given that we use DEG thresholds the correct file has to be chosen for mean network distances
  lf <- list.files(path = file_path)
  lf <- lf[grepl(lf, pattern = file_pattern)]
  lf1 <- lf[grepl(lf, pattern = paste(clusters[i],".txt",sep=""))] # protects for example from Cluster_3 being matched with files of Cluster_31
  lf2 <- lf[grepl(lf, pattern = paste(clusters[i],"_",sep=""))]
  lf <- unique(c(lf1,lf2))
  if(any(grepl(lf,pattern = "_top_3500"))){
    lf <- lf[grepl(lf, pattern = "top_3500")]
  } else {
    lf <- lf[!grepl(lf, pattern = "_top_")]
  }
  pred_mat <- as.matrix(read.table(file = paste(file_path, lf, sep=""), sep="\t", header = T, stringsAsFactors = F))  
  pred_mat <- pred_mat[pred_mat[,1]%in%match_IDs[,2],]
  mean_dist[,i] <- pred_mat[match(match_IDs[,2],pred_mat[,1]),3]
}
colnames(mean_dist) <- paste("mean_dc_", colnames(mean_dist),sep="")

rm(pred_mat, i, match_IDs, file_pattern)

# mean of mean closest distance
#######################################################################

mean_dist <- as.matrix(mean_dist)
mode(mean_dist) <- "numeric"
mean_dist <- cbind(mean_dist, rowMeans(mean_dist))
colnames(mean_dist)[ncol(mean_dist)] <- "mean_of_mean_dc"

# n of clusters in which drug fulfills selection criteria
#######################################################################

n_clusters_fulfilled <- table(final_drugs[,1])
n_clusters_fulfilled <- n_clusters_fulfilled[match(unique_drugs[,1],names(n_clusters_fulfilled))]

# fullfills selection criteria in which clusters? one column for every cluster, indicated by TRUE or FALSE if criteria is passed
# + combined Eigenvector centrality score for all clusters in which selection criteria is passed
#######################################################################

fulfilled_criteria <- matrix(F, nrow = nrow(unique_drugs), ncol = length(clusters))
rownames(fulfilled_criteria) <- unique_drugs[,1]
colnames(fulfilled_criteria) <- clusters

combined_centrality <- vector()
centrality_scores <- centrality_scores[match(colnames(fulfilled_criteria), names(centrality_scores))]

for(i in 1:nrow(fulfilled_criteria)){
  temp <- final_drugs[final_drugs[,1] %in% rownames(fulfilled_criteria)[i],17] # extract in which clusters drug fulfilled criteria
  temp <- as.numeric(temp)
  temp <- paste("Cluster_", temp, sep="")
  fulfilled_criteria[i,colnames(fulfilled_criteria)%in% temp] <- T
  
  combined_centrality[i] <- sum(centrality_scores[names(centrality_scores)%in% temp])
}
mode(fulfilled_criteria) <- "character"
colnames(fulfilled_criteria) <- paste("Passes_selection_criteria_",colnames(fulfilled_criteria),sep="")
fulfilled_criteria <- cbind(fulfilled_criteria, n_clusters_fulfilled, combined_centrality)
colnames(fulfilled_criteria)[(ncol(fulfilled_criteria)-1):ncol(fulfilled_criteria)] <- c("Fulfilled_criteria_in_n_clusters", "Sum_of_targeted_clusters_Eigenvector_centrality")
rm(combined_centrality, centrality_scores, temp, i)

# COMBINE ALL MATRICES INTO 1
#######################################################################
out <- cbind(unique_drugs, mean_dist, fulfilled_criteria)
out <- out[order(as.numeric(out[,ncol(out)]), as.numeric(out[,(5+length(clusters))]), decreasing = T),] # ordered by 1) combined Eigenvector centrality and 2) cluster-wide mean of mean d(c)
write.table(out, file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_INDIVIDUAL_CD_PATIENT_1.txt", sep="\t", col.names = T, row.names = F)

rm(list = ls())


