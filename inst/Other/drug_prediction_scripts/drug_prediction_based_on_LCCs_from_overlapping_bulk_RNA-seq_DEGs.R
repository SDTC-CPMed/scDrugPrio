# Drug prediction based on DEGs overlapping between data sets (blood and synovium)

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

set.seed(35)

# Load bulk RNA-seq data
bulk_synovium <- read.table("Input/Bulk_data/synovium_DEGs_GSE55235.txt", sep="\t", header = T)
bulk_synovium <- bulk_synovium[order(bulk_synovium[,6], decreasing = F),]
bulk_synovium <- bulk_synovium[bulk_synovium[,6]<0.05,1]

bulk_blood <- read.table("Input/Bulk_data/whole_blood_DEGs_GSE93272.txt", sep="\t", header = T)
bulk_blood <- bulk_blood[order(bulk_blood[,6], decreasing = F),]
bulk_blood <- bulk_blood[bulk_blood[,6]<0.05,1]

# Combine gene lists
degs <- cbind(c(bulk_synovium, rep(NA, times = 10000-length(bulk_synovium))), c(bulk_blood, rep(NA, times = 10000-length(bulk_blood))))
colnames(degs) <- c("Synovial_DEGs_bulk_RNA-seq", "Whole_blood_DEGs_bulk_RNA-seq")
degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
degs <- as.matrix(degs)
mode(degs) <- "character"
rm(bulk_synovium, bulk_blood)
# Overlapping DEGs
degs <- intersect(degs[,1], degs[,2])
degs <- matrix(degs, ncol = 1)
colnames(degs) <- "bulk_DEGs_overlapping"

# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/predefined_RA_drugs.txt", sep="\t", header = F)[,1])

# Set up output directories
dir.create("Output/Other_predictions/literature_PPI")
dir.create("Output/Other_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/Other_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])

# Analysis
save_name <- "LCCs_formed_by_overlapping_bulk_RNA-seq_DEGs_between_data_sets_in_literature_PPI_1000"
print(save_name)
rec <- F

for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = T, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes that can be mapped to PPI")
  }
  rm(out)
}



# start drug prediction HuRI PPI
##################################################################################################################
print("START HuRI PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/Other_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/HuRI_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Analysis
save_name <- "LCCs_formed_by_overlapping_bulk_RNA-seq_DEGs_between_data_sets_in_HuRI_PPI_1000"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = T, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes that can be mapped to PPI")
  }
  rm(out)
}


# remove everything when done
rm(list = ls())
