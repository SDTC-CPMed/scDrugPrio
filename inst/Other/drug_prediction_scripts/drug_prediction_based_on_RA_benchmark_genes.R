# Drug prediction based on GWAS, OMIM, GWAS & OMIM, bulk RNA-seq data (as benchmark)

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

set.seed(35)

# Load GWAS
gwas <- read.table(file = "Input/GWAScat/TRANSLATED_GWAS_RA_P%3C1e-8.txt", sep="\t", header = T)
gwas <- as.vector(gwas[,1])
# Load OMIM
omim <- as.matrix(read.table(file = "Input/OMIM/OMIM.txt", sep="\t", header = T))
omim <- as.vector(omim)
# Load microarray data
bulk_synovium <- read.table("Input/Bulk_data/synovium_DEGs_GSE55235.txt", sep="\t", header = T)
bulk_synovium <- bulk_synovium[order(bulk_synovium[,6], decreasing = F),]
bulk_synovium <- bulk_synovium[bulk_synovium[,6]<0.05,1]

bulk_blood <- read.table("Input/Bulk_data/whole_blood_DEGs_GSE93272.txt", sep="\t", header = T)
bulk_blood <- bulk_blood[order(bulk_blood[,6], decreasing = F),]
bulk_blood <- bulk_blood[bulk_blood[,6]<0.05,1]

# Combine gene lists
degs <- cbind(c(gwas, rep(NA, times = 10000-length(gwas))), c(omim, rep(NA, times = 10000-length(omim))), c(gwas, omim, rep(NA, times = 10000-length(gwas)-length(omim))), 
              c(bulk_synovium, rep(NA, times = 10000-length(bulk_synovium))), c(bulk_blood, rep(NA, times = 10000-length(bulk_blood))))
colnames(degs) <- c("GWAS", "OMIM", "GWAS_&_OMIM", "Synovial_DEGs", "Whole_blood_DEGs")
degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
degs <- as.matrix(degs)
mode(degs) <- "character"
rm(gwas,omim, bulk_synovium, bulk_blood)

# Set up output directories
dir.create("Output/RA_benchmark_predictions/literature_PPI")
dir.create("Output/RA_benchmark_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis for RA benchmark DEGs")
out.dir <- paste(getwd(), "/Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
setwd(paste(fp, "/..", sep=""))
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 16 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "RA_benchmark_genes_literature_PPI"
print(save_name)
rec <- F

for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 10, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes mapped to PPI")
  }
  rm(out)
}



# start drug prediction HuRI PPI
##################################################################################################################
print("START HuRI PPI network analysis for RA benchmark DEGs")
out.dir <- paste(getwd(), "/Output/RA_benchmark_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
setwd(paste(fp, "/..", sep=""))
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
predefined_list <- as.vector(read.table("Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_HuRI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 10 unique drug combinations found in HuRI PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "RA_benchmark_genes_HuRI_PPI"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 15, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes mapped to PPI")
  }
  rm(out)
}


# remove everything when done
rm(list = ls())
