
#######################################################################################
# FOR CD DATA SET - all individuals patients 
# BY SAMUEL SCHAEFER
# Run in R 4.0.4
#######################################################################################

source("From scDrugPrio 220925/intracellular_drug_centrality.R")
source("From scDrugPrio 220925/extract_LCC_by_gene_set_of_interest.R")
source("From scDrugPrio 220925/ppin_formatting.R")
source("From scDrugPrio 220925/extract_LCC.R")

fp <- getwd()
setwd(paste(fp, "/../", sep=""))
library(doParallel)

#######################################################################################
# INPUT
#######################################################################################

# Patients
patient <- paste("Patient_", 1:11, sep = "")

#PPIN
lit_ppi <- as.matrix(read.table(file = "Input/literature_PPI/ppi.txt", sep ="\t", header = T, stringsAsFactors = F))
ppin <- ppin_formatting(lit_ppi)

# transl
transl <- read.table(file = "Input/HGNC translation matrix 201108/transl.txt", header = T, stringsAsFactors = F, sep="\t")
transl <- transl[,c(2,6)]

# CD drugs
cd_drug_info <- read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", header = T, sep="\t", stringsAsFactors = F, quote = "")
cd_drugs <- unique(cd_drug_info[,1])

# Drug target matrix
drugs <- read.table(file = "Input/Drugs/drug_targets_unique_literature_ppi.txt", header = T, stringsAsFactors = F, sep="\t")
same_drugs <- read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
#same_drugs <- same_drugs[same_drugs[,2] %in% cd_drugs,]
same_drugs <- same_drugs[same_drugs[,1] %in% colnames(drugs),]
if(ncol(cd_drug_info) < 4){
  cd_drug_info <- cbind(cd_drug_info, unique_drug_target_combination = NA)
}
cd_drug_info[cd_drug_info[,1]%in% same_drugs[,2],4] <- "x"
write.table(cd_drug_info, file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", col.names = T, row.names = F)

drugs <- drugs[,colnames(drugs)%in% same_drugs[,1]]
drugs <- as.matrix(drugs)

# Check if drugs have valid target in PPIN
for(i in 1:ncol(drugs)){
  print(paste(colnames(drugs)[i], "has valid targets:", all(drugs[!is.na(drugs[,i]),i] %in% as.vector(ppin)), sep=" "))
}
#write.table(same_drugs, "../Input/Drugs/SAME_DRUGS_cd.txt", sep="\t", col.names = T, row.names = F)
#write.table(drugs, "../Input/Drugs/DRUG_TARGET_MATRIX_cd.txt", sep="\t", col.names = T, row.names = F)

#######################################################################################
# Intracellular centrality for all patients
#######################################################################################

out_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
in_path <- "Input/CD GSE134809/Individual_patients/"

for(pat in patient){
  
  print(paste("PROCESSING: ", pat, sep=""))
  
  out_dir <- paste(out_path,pat,"/Intracellular_centrality",sep="")
  dir.create(out_dir, showWarnings = F)
  
  # DEGs
  degs <- as.matrix(read.table(file = paste(in_path, pat, "/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt", sep=""),
                                            header = T, stringsAsFactors = F, sep="\t"))
  mode(degs) <- "numeric"
  
  # Select only top 3500 most significant DEGs
  if(nrow(degs)>3500){
    degs <- degs[1:3500,]
  }
  
  # run intracellular centrality analysis
  intra_cent <- intracellular_drug_target_centrality(ppin = ppin,
                                                     drug_target_matrix = drugs,
                                                     degs = degs,
                                                     file_name = "CD_lit_PPIN_intracellular_centralities",
                                                     centrality_alg = "eigenvector centralities",
                                                     out_dir = out_dir)
  
  #print(head(intra_cent))
  
  # seperate into individual drugs
  intra_cent <- cbind(same_drugs[,2], intra_cent[match(same_drugs[,1], rownames(intra_cent)),])
  rownames(intra_cent) <- same_drugs[,2]
  write.table(intra_cent, file = paste(out_dir, "/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep=""), sep = "\t", col.names = NA, row.names = T)
}




