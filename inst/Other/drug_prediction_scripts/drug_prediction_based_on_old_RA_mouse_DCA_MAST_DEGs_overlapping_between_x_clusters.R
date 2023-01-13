# Drug prediction based on DEGs (calculated with DCA + MAST) that overlap between x clusters + LCCs of these overlapping DEGs

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
set.seed(35)

degs <- read.table("Output/DCA_MAST_DEGs/TRANSLATED_Summary_sig_adj_MAST_DEGs_all_clusters_res=0.6_dims=32.txt", sep="\t", header = T)
degs <- as.matrix(degs)
mode(degs) <- "character"
source(paste(fp, "/overlap_between_columns.R",sep="")) # contains overlapping_genes() function
degs <- overlapping_genes(degs)
degs <- degs[,-1]
degs <- degs[rowSums(!is.na(degs))>0,]

# Set up output directories
dir.create("Output/old_RA_mouse")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions/HuRI_PPI")


# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis of overlapping RA DEGs")
out.dir <- paste(fp, "/../Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
setwd(paste(fp, "/..", sep=""))

# LOAD PPI
ppi <- as.matrix(read.table(file = "Input/literature_PPI/ppi.txt", sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
print("PPI LOADED")
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
print("DRUGS LOADED")
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table(file = "Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 16 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)
print("RA DRUGS LOADED")

# Analysis
save_name <- "OVERLAPPING_RA_DCA_MAST_DEGs_literature_PPI"
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

save_name <- "OVERLAPPING_RA_DCA_MAST_DEGs_LCCs_literature_PPI"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = T, out.dir = out.dir,
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



# start drug prediction literature PPI
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
drugs <- read.table("Input/Drugs/drug_targets_unique_HuRI_ppi.txt", sep="\t", header = T)
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
save_name <- "OVERLAPPING_RA_DCA_MAST_DEGs_HuRI_PPI"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: OVERLAPPING DEGs from ", colnames(degs)[i], sep=""))
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

save_name <- "OVERLAPPING_RA_DCA_MAST_DEGs_LCCs_HuRI_PPI"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: OVERLAPPING DEGs LCCs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = T, out.dir = out.dir,
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

# remove everything when done
rm(list = ls())
