# Drug prediction based on DEGs calculated based on DCA normalized data with MAST
library(doParallel)

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
set.seed(35)

# Get human-mouse gene translation  
transl <- read.table(file = "Input/Human-mouse_homologs/transl.txt", sep="\t", header = T)
transl <- transl[!is.na(transl[,3]),]

# Load DCA + MAST DEGs
degs <- read.table("Output/DCA_MAST_DEGs/Summary_sig_adj_MAST_DEGs_log(1.5)_all_clusters_res=0.6_dims=32.txt", sep="\t", header = T)
degs <- as.matrix(degs)
mode(degs) <- "character"
for(i in 1:ncol(degs)){
  temp <- degs[,i]
  degs[,i] <- NA
  temp <- transl[transl[,3]%in%temp,2]
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
rm(temp)

# Set up output directories
dir.create("Output/old_RA_mouse")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis of RA DEGs")
out.dir <- paste(fp, "/../Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
setwd(paste(fp, "/..", sep=""))

# LOAD PPI
ppi <- as.matrix(read.table(file = "Input/literature_PPI/ppi.txt", sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
ppi <- ppi[ppi[,1]%in%transl[,2],]
ppi <- ppi[ppi[,2]%in%transl[,2],]
print("PPI LOADED")

# Load drug_matrix
drugs <- as.matrix(read.table("Input/DrugBank/all_drug_targets_drug_bank.txt", sep="\t", header = T))
drugs <- drugs[drugs[,8] == "Humans",]
drugs <- drugs[grepl(drugs[,3],pattern = "approved"),]
drugs[,6] <- transl[match(drugs[,6], transl[,1]),2]
drugs <- drugs[!is.na(drugs[,6]),]
unique_drugs <- unique(drugs[,1])
drugs <- foreach(i = c(1:length(unique_drugs)), .combine = "cbind") %do% {
  temp <- drugs[drugs[,1] == unique_drugs[i],6]
  temp <- c(unique_drugs[i],temp[order(as.numeric(temp), decreasing = F)], rep(NA, times = 500-length(temp)))
}
colnames(drugs) <- drugs[1,]
drugs <- drugs[-1,]
drugs <- drugs[rowSums(!is.na(drugs))>0,]
# keep only unique drug target combinations
combinations <- vector()
for(i in 1:ncol(drugs)){
  temp <- drugs[!is.na(drugs[,i]),i]
  if(length(temp) == 1){
    combinations <- c(combinations, temp)
  } else {
    combinations <- c(combinations, paste(temp, collapse = "_"))
  }
}
names(combinations) <- colnames(drugs)
unique_combinations <- table(combinations)
same_drugs <- foreach(i = c(1:length(unique_combinations)), .combine = "rbind") %do% {
  temp <- names(combinations)[combinations == names(unique_combinations)[i]]
  return(cbind(rep(temp[1], times = length(temp)), temp))
}
write.table(x = same_drugs, file = "Input/Same_drugs_for_drug_targets_translated_by_human_mouse_homologs.txt",sep="\t", col.names = F, row.names = F)
drugs <- drugs[,colnames(drugs) %in% same_drugs[,1]]
write.table(x = drugs, file = "Input/Drugs_translated_by_human_mouse_homologs.txt",sep="\t", col.names = T, row.names = F)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
print("DRUGS LOADED")
rm(temp, i, combinations, unique_combinations, transl, unique_drugs)

# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table(file = "Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 16 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)
print("RA DRUGS LOADED")

# Analysis
save_name <- "RA_DCA_MAST_DEGs_literature_PPI_MOUSE_GENES_ONLY"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes")
  }
  rm(out)
}


'
# start drug prediction HuRI PPI
##################################################################################################################
print("START literature PPI network analysis of RA DEGs")
out.dir <- paste(fp, "/../Output/old_RA_mouse/DCA_MAST_DEGs_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
setwd(paste(fp, "/..", sep=""))

# LOAD PPI
ppi <- as.matrix(read.table(file = "Input/HuRI_PPI/ppi.txt", sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Load drug_matrix 
drugs <- read.table(paste("Input/Drugs/drug_targets_unique_HuRI_ppi.txt" sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table(file = "Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_HuRI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 10 unique drug combinations found in HuRI PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "RA_DCA_MAST_DEGs_HuRI_PPI"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes")
  }
  rm(out)
}
'

# remove everything when done
rm(list = ls())
