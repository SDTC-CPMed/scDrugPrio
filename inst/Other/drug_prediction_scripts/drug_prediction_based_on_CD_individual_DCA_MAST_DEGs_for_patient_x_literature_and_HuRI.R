#
#
# Drug prediction based on CD DEGs calculated based on DCA denoised data with MAST
#
##################################################################################################################
patient <- "Patient_"

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
set.seed(35)

# Load DCA + MAST DEGs
DEG <- read.table(file = paste("Input/CD GSE134809/Individual_patients/",patient,"/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt", sep=""), sep="\t", header = T)
DEG <- as.matrix(DEG)
mode(DEG) <- "character"

# Set up output directories
dir.create("Output/CD")
dir.create("Output/CD/Individual_patients")
dir.create("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions")
dir.create(paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/", patient,sep=""))
dir.create(paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/", patient, "/literature_PPI",sep=""))
dir.create(paste("Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/", patient, "/HuRI_PPI",sep=""))

# start drug prediction literature PPI
##################################################################################################################
print(paste("START literature PPI network analysis for CD DEGs of ",patient, sep=""))
out.dir <- paste(getwd(), "/Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/", patient, "/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()

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
predefined_list <- as.vector(read.table("Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", header = T, quote = "")[,1]) # 28 drugs
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) # 67 drugs with same drug target combination as CD drugs
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)] # 15 unique drug combinations found in literature PPI
rm(same_drugs)
# limit analysis to 4000 most significant DEGs found in PPI
degs <- DEG
for(i in 1:ncol(degs)){
  temp <- degs[,i]
  temp <- temp[temp%in%unique(c(ppi[,1], ppi[,2]))]
  degs[,i] <- NA
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
if(all(colSums(!is.na(degs))<2) | ((ncol(degs)) >= (sum((colSums(!is.na(degs))==0))+1))){ # removing row/column could set matrix object to vector() if only 1 row filled with DEGs or only one column filled with DEGs
  if(all(colSums(!is.na(degs))<2)){ # having only one DEG in every cluster
    if((ncol(degs)) == (sum((colSums(!is.na(degs))==0))+1)){
      if(any(colSums(!is.na(degs)))>1){
        degs <- degs[rowSums(!is.na(degs))>0,]
      } else {
        degs <- rbind(degs[rowSums(!is.na(degs))>0,],NA)
      }
    }
  }
  if(sum(colSums(!is.na(degs))>0)==1){ # having only one column with DEGs - handled later in the in silico drug prediction script
    degs <- rbind(degs[rowSums(!is.na(degs))>0,],NA)
  } else {
    print("NO DEGs IN PPI NETWORK - ANALYSIS ABORTED")
  }
} else {
  degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
}
if(nrow(degs)>4000){
  degs <- degs[1:4000,]
}
rm(temp)

# Analysis
save_name <- "CD_DCA_MAST_DEGs_literature_PPI"
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
print(paste("START HuRI PPI network analysis of CD DEGs of ",patient, sep=""))
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/", patient, "/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()

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
predefined_list <- as.vector(read.table("Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", header = T,quote = "")[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_HuRI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #15 drugs, 10 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)
# limit analysis to 4000 most significant DEGs found in PPI
degs <- DEG
for(i in 1:ncol(degs)){
  temp <- degs[,i]
  temp <- temp[temp%in%unique(c(ppi[,1], ppi[,2]))]
  degs[,i] <- NA
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
if(all(colSums(!is.na(degs))<2) | ((ncol(degs)) >= (sum((colSums(!is.na(degs))==0))+1))){ # removing row/column could set matrix object to vector() if only 1 row filled with DEGs or only one column filled with DEGs
  if(all(colSums(!is.na(degs))<2)){ # having only one DEG in every cluster
    if((ncol(degs)) == (sum((colSums(!is.na(degs))==0))+1)){
      if(any(colSums(!is.na(degs)))>1){
        degs <- degs[rowSums(!is.na(degs))>0,]
      } else {
        degs <- rbind(degs[rowSums(!is.na(degs))>0,],NA)
      }
    }
  }
  if(sum(colSums(!is.na(degs))>0)==1){ # having only one column with DEGs - handled later in the in silico drug prediction script
    degs <- rbind(degs[rowSums(!is.na(degs))>0,],NA)
  } else {
    print("NO DEGs IN PPI NETWORK - ANALYSIS ABORTED")
  }
} else {
  degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
}
if(nrow(degs)>4000){
  degs <- degs[1:4000,]
}
rm(temp)


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
