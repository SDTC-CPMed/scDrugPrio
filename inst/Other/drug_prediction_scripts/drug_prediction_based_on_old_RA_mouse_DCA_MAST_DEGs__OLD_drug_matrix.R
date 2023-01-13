# Drug prediction based on DEGs calculated based on DCA normalized data with MAST

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

set.seed(35)

# Load DCA + MAST DEGs
degs <- read.table("Output/DCA_MAST_DEGs/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt", sep="\t", header = T, stringsAsFactors = F)
degs <- as.matrix(degs)
mode(degs) <- "character"
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_old.txt", sep="\t", header = T, stringsAsFactors = F)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- read.table("Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T, stringsAsFactors = F, quote = "")
predefined_list <- as.vector(predefined_list[,1])

# Set up output directories
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI")
dir.create("Output/old_RA_mouse/DCA_MAST_DEGs_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis with OLD DRUG MATRIX")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/old_RA_mouse/DCA_MAST_DEGs_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# make sure drugs are in PPI
unique_proteins <- unique(c(ppi[,1],ppi[,2]))
for(i in 1:ncol(drugs)){
  temp <- drugs[,i]
  temp <- temp[!is.na(temp)]
  temp <- temp[temp %in% unique_proteins]
  drugs[,i] <- NA
  if(length(temp)>0){
    drugs[1:length(temp),i] <- temp
  }
}
drugs <- drugs[rowSums(!is.na(drugs))>0, colSums(!is.na(drugs))>0]
rm(temp, unique_proteins)

# Analysis
save_name <- "DCA_MAST_DEGs_in_literature_PPI_OldDrugMatrix35"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 14, recycle = rec, n_random_iterations = 1000)
  rec <- F
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
print("START HuRI PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/DCA_MAST_DEGs_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/HuRI_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Analysis
save_name <- "DCA_MAST_DEGs_in_HuRI_PPI_1000"
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs from ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  rec <- T
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
