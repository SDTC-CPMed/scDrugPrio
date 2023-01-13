#
#
# Drug prediction based on DRESS skin DEGs calculated based on DCA denoised data with MAST
#
##################################################################################################################
fp <- getwd()
setwd(paste(fp, "/..", sep=""))
set.seed(35)

# Load DCA + MAST DEGs
degs <- read.table("Input/DRESS GSE132802/Skin cells/DCA_MAST_DEGs_logFC_1.5_DRESS_SKIN/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt", sep="\t", header = T)
degs <- as.matrix(degs)
mode(degs) <- "character"
# Set up output directories
dir.create("Output/DRESS_skin")
dir.create("Output/DRESS_skin/DCA_MAST_DEGs_predictions")
dir.create("Output/DRESS_skin/DCA_MAST_DEGs_predictions/literature_PPI")
dir.create("Output/DRESS_skin/DCA_MAST_DEGs_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis for DRESS_skin DEGs")
out.dir <- paste(getwd(), "/Output/DRESS_skin/DCA_MAST_DEGs_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()

# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # removes "temp" folder if exists
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/DRESS_drugs_from_articles.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #25 drugs, 18 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "DRESS_skin_DCA_MAST_DEGs_literature_PPI"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 10, recycle = rec, n_random_iterations = 1000)
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



# start drug prediction HuRI PPI
##################################################################################################################
print("START HuRI PPI network analysis of DRESS_skin DEGs")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/DRESS_skin/DCA_MAST_DEGs_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()

source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/HuRI_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_HuRI_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/DRESS_drugs_from_articles.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_HuRI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #15 drugs, 10 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "DCA_MAST_DEGs_HuRI_PPI"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 10, recycle = rec, n_random_iterations = 1000)
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
