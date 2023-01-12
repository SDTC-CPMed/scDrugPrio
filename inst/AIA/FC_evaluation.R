
# FC EVALUATION
# Run on R 4.0.4
####################################################################################

source("Codes_from_scPred/separate_unique_drug_target_combinations_into_individual_drugs.R")
source("Codes_from_scPred/pharma_effect_on_fold_change.R")
source("Codes_from_scPred/combine_evaluation_files.R")
dir.create("../Output/FC_criteria_checking", showWarnings = F)

# INPUT
same_drugs <- as.matrix(read.table(file ="../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F))
drugs <- read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F)
ra_drugs <- unique(read.table(file = "../Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T, stringsAsFactors = F, quote = "")[,1])
# lets translate the raw drug bank matrix back to human symbols in order to make the interpretation a little bit easier
drug_bank <- read.table(file = "../Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", header = T, stringsAsFactors = F)

# TRANSLATE DEG FILES TO HUMAN ENTREZ SYMBOLS
dir.create(path = "../Output/DCA_MAST_DEGs/TRANSLATED")
transl <- read.table(file = "../Input/Human-mouse_homologs/transl.txt", sep="\t", header = T)
lf <- list.files(path = "../Output/DCA_MAST_DEGs//", pattern = "Cluster_", full.names = T)
lf_name <- unlist(strsplit(lf, split = "../Output/DCA_MAST_DEGs///"))
lf_name <- lf_name[lf_name != ""]
for(i in 1:length(lf)){
  temp <- read.table(file = lf[i], sep="\t", header = T)
  temp <- cbind(transl[match(temp[,1],transl[,3]),1], temp)
  if(sum(!is.na(temp[,1]))>1){
    temp <- temp[!is.na(temp[,1]),-2]
  } else if(sum(!is.na(temp[,1])) == 1){
    temp_col <- colnames(temp)
    temp <- matrix(temp[!is.na(temp[,1]),-2], nrow = 1)
    colnames(temp) <- temp_col
  } else {
    next
  }
  write.table(temp, file = paste("../Output/DCA_MAST_DEGs/TRANSLATED/", lf_name[i], sep=""), sep="\t", col.names = T, row.names = F)
}


# Literature PPIN
#########################################

# seperate unique drug target combinations into individual drugs again
out_dir <- "../Output/Network_distances/"
lf  <- list.files(path = out_dir, pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__")
separate_unique_drug_target_combinations_into_individual_drugs(files = lf, 
                                                               remove_part_of_file_name_for_output = "drug-disease_closest_distances_vs_random_bin_adjusted__",
                                                               same_drugs = same_drugs, 
                                                               disease_specific_drugs = ra_drugs, 
                                                               output_file_name_add_on = "INDIVIDUAL_DRUGS_", 
                                                               in_dir = out_dir, 
                                                               out_dir = out_dir)

# Let's apply FC criteria
#############################

# translate DEG files
lf <- list.files(path = "../Output/DCA_MAST_DEGs/TRANSLATED", pattern = "Cluster_", full.names = T)
# get paths to average closest netork distances
drug_dists <- list.files("../Output/Network_distances/", pattern = "INDIVIDUAL_DRUGS_", full.names = T)
# remove deg files that are empty
temp <- colnames(read.table(file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_ID.txt", header = T))
temp <- paste(temp, "_", sep = "")
temp2 <- rep(F, times = length(lf))
for(i in 1:length(temp)){
  temp2[grepl(lf, pattern = temp[i])] <- T
}
lf <- lf[temp2]
rm(temp, temp2)

# save names
save_name <- unlist(strsplit(lf, split = ".txt"))
save_name <- unlist(strsplit(save_name, split = "../Output/DCA_MAST_DEGs/TRANSLATED/"))
save_name <- save_name[save_name != ""]
save_name <- paste("FC_INDIVIDUAL_DRUGS_", save_name, sep="")
# make output folder
out_dir <- "../Output/FC_criteria_checking/"
dir.create(out_dir, showWarnings = F)

fc_evaluation <- pharma_effect_on_fold_change(drug_dist_files = drug_dists, 
                                              deg_files = lf, 
                                              pharma_effect = drug_bank[,c(1,6:7)], 
                                              save_names = save_name, 
                                              out_dir = out_dir)

# Apply drug selection network criteria
fc_evaluation <- fc_evaluation[as.numeric(fc_evaluation[,4]) < 1,]
fc_evaluation <- fc_evaluation[as.numeric(fc_evaluation[,8]) < 0.05,]
fc_evaluation <- fc_evaluation[order(fc_evaluation[,1]),]
fc_evaluation <- cbind(Drug_name = drug_bank[match(fc_evaluation[,1], drug_bank[,1]),2], fc_evaluation)
write.table(fc_evaluation, file = paste(out_dir, "SUMMARY_only_drugs_passing_network_criteria.txt",sep=""), sep="\t", col.names = T, row.names = F)


# HuRI
#########################################
