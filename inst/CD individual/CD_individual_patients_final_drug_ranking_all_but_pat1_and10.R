
#######################################################################################
# Drug ranking script for processing after all selection criteria has been run
# FOR CD DATA SET - all individuals patients but patient 1 and 10 that were processed earlier
# BY SAMUEL SCHAEFER
#######################################################################################

source("From scDrugPrio 220925/final_drug_prioritization.R")

fp <- getwd()
setwd(paste(fp, "/../", sep=""))
library(doParallel)

#######################################################################################
# General input data
#######################################################################################

# DrugBank
drugbank <- as.matrix(read.table(file = "Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", stringsAsFactors = F, header = T))

# Drug ID matching with unique drug target combination used in prediction
drugID_matching <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F))

# Unique drug target combination matrix
unique_target_combinations <- as.matrix(read.table(file = "Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F))
unique_target_combinations <- colnames(unique_target_combinations)

#######################################################################################
# Seperate FC criteria files by patient
#######################################################################################

out_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"

all_file <- readxl::read_xlsx(path = paste(out_path, "COMBINED_EVALUATION_FC_for_all_but_pat1_and_10_evaluated.xlsx", sep = ""), sheet = 1, col_names = T, col_types =  "text")
patient <- unique(all_file$patient)

all_file <- all_file[,colSums(is.na(all_file)) != nrow(all_file)]
all_file$n_drug_targets_counteracting_disease <- as.numeric(all_file$n_drug_targets_counteracting_disease)
all_file$n_drug_targets_mimmicking_disease <- as.numeric(all_file$n_drug_targets_mimmicking_disease)

for(pat in patient){
    temp <- all_file[all_file$patient %in% pat, ]
    temp <- as.matrix(temp)
    temp <- cbind(temp[,1], drugbank[match(as.character(temp[,1]), drugbank[,2]),1], temp[,2:ncol(temp)])
    write.table(temp, file = paste(out_path, pat, "/", pat, "_FC_criteria_evaluated.txt", sep =""), sep ="\t", col.names = T, row.names = F)
}

for(pat in patient){
  
  # LOAD FC file
  fc_file <- as.matrix(read.table(file = paste(out_path, pat, "/", pat, "_FC_criteria_evaluated.txt", sep =""), sep ="\t", header =  T, stringsAsFactors = F))
  fc_file[,18] <- paste("Cluster_",fc_file[,18],sep="")
  
  # LOAD intercellular centrality
  inter_cents <- as.matrix(read.table(file = paste(out_path, pat, "/NicheNet/Cell_type_centrality_summary.txt", sep = ""), sep="\t", header = T, stringsAsFactors = F))
  inter_cents <- inter_cents[5,-1] #choosing only eigenvector centrality
  temp_names <- names(inter_cents)
  inter_cents <- as.numeric(inter_cents)
  names(inter_cents) <- temp_names
  rm(temp_names)
  
  # LOAD intracellular centrality
  intra_cents <- as.matrix(read.table(file = paste(out_path, pat, "/Intracellular_centrality/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = ""), sep="\t", header = T, stringsAsFactors = F))
  temp_names <- intra_cents[,1]
  intra_cents <- as.numeric(intra_cents[,ncol(intra_cents)])
  names(intra_cents) <- temp_names
  rm(temp_names)
  
  # FINAL DRUG RANKING
  dir.create(path = paste(out_path, pat, "/Drug_ranking", sep =""))
  final <- final_drug_prioritization(fc_evaluation = fc_file,
                            pos_DrugID = 2,
                            pos_clusterID = 18,
                            keep = c(1,2,8,17,19),
                            inter_cent = inter_cents, 
                            intra_cent = intra_cents, 
                            file_name = "FINAL_drug_ranking", 
                            out_dir = paste(out_path, pat, "/Drug_ranking/", sep=""))
}

patient <- paste("Patient_",c(1,10),sep="")

for(pat in patient){
  
  print(pat)
  
  lf <- list.files(path = paste(out_path, pat,"/literature_PPI/SUMMARY/", sep = ""), pattern = ".xlsx", full.names = T)
  lf <- lf[grepl(lf, pattern = "COMBINED_EVALUATION_FILES_for_checking_FC_criterion_patient_")]
  # LOAD FC file
  temp <- readxl::read_xlsx(path = lf[1], sheet = 1, col_names = T, col_types =  "text")
  temp <- as.matrix(temp)
  temp <- cbind(temp[,1], drugbank[match(as.character(temp[,1]), drugbank[,2]),1], temp[,2:ncol(temp)])
  write.table(temp, file = paste(out_path, pat, "/", pat, "_FC_criteria_evaluated.txt", sep =""), sep ="\t", col.names = T, row.names = F)
  
  # LOAD STANDARD FC FILE
  fc_file <- as.matrix(read.table(file = paste(out_path, pat, "/", pat, "_FC_criteria_evaluated.txt", sep =""), sep ="\t", header =  T, stringsAsFactors = F))
  fc_file[,18] <- paste("Cluster_", fc_file[,18],sep="")
  
  # LOAD intercellular centrality
  inter_cents <- as.matrix(read.table(file = paste(out_path, pat, "/NicheNet/Cell_type_centrality_summary.txt", sep = ""), sep="\t", header = T, stringsAsFactors = F))
  inter_cents <- inter_cents[5,-1] #choosing only eigenvector centrality
  temp_names <- names(inter_cents)
  inter_cents <- as.numeric(inter_cents)
  names(inter_cents) <- temp_names
  rm(temp_names)
  
  # LOAD intracellular centrality
  intra_cents <- as.matrix(read.table(file = paste(out_path, pat, "/Intracellular_centrality/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = ""), sep="\t", header = T, stringsAsFactors = F))
  temp_names <- intra_cents[,1]
  intra_cents <- as.numeric(intra_cents[,ncol(intra_cents)])
  names(intra_cents) <- temp_names
  rm(temp_names)
  
  # FINAL DRUG RANKING
  dir.create(path = paste(out_path, pat, "/Drug_ranking", sep =""))
  final <- final_drug_prioritization(fc_evaluation = fc_file,
                                     pos_DrugID = 2,
                                     pos_clusterID = 18,
                                     keep = c(1,2,8,17,19),
                                     inter_cent = inter_cents, 
                                     intra_cent = intra_cents, 
                                     file_name = "FINAL_drug_ranking", 
                                     out_dir = paste(out_path, pat, "/Drug_ranking/", sep=""))
}




