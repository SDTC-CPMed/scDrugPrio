
# INTRACELLULAR CENTRALITY
# Run in R 4.0.4
######################################################################

#setwd("/data/samsc76/SIGMA_13012020/CD_batch_corrected/drug_prediction_R")
dir.create("../Output/Intracellular_centrality", showWarnings = F)
source("Codes_from_scPred/intracellular_drug_centrality.R")
source("Codes_from_scPred/extract_LCC_by_gene_set_of_interest.R")
source("Codes_from_scPred/ppin_formatting.R")
source("Codes_from_scPred/extract_LCC.R")


#PPIN
lit_ppi <- as.matrix(read.table(file = "../Input/literature_PPI/ppi.txt", sep ="\t", header = T, stringsAsFactors = F))
ppin <- ppin_formatting(lit_ppi)

# transl
transl <- read.table(file = "../Input/HGNC translation matrix 201108/transl.txt", header = T, stringsAsFactors = F, sep="\t")
transl <- transl[,c(2,6)]

# MS drugs
cd_drug_info <- read.table(file = "../Input/Drugs/CD_drugs_from_DrugBank_201015.txt", header = T, sep="\t", stringsAsFactors = F, quote = "")
cd_drugs <- unique(cd_drug_info[,1])

# Drug target matrix
drugs <- read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt", header = T, stringsAsFactors = F, sep="\t")
same_drugs <- read.table(file = "../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
same_drugs <- same_drugs[same_drugs[,2] %in% cd_drugs,]
same_drugs <- same_drugs[same_drugs[,1] %in% colnames(drugs),]
cd_drug_info[cd_drug_info[,1]%in% same_drugs[,2],4] <- "x"
colnames(cd_drug_info)[4] <- "targets_in_lit_PPI"
write.table(cd_drug_info, file = "../Input/Drugs/CD_drugs_from_DrugBank.txt", sep="\t", col.names = T, row.names = F)
same_drugs <- read.table(file = "../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
same_drugs <- same_drugs[same_drugs[,1] %in% colnames(drugs),]


#drugs <- drugs[,colnames(drugs)%in% same_drugs[,1]]
drugs <- as.matrix(drugs)

# Intracellular centrality for all drugs
########################################

# out_dir
out_dir <- "../Output/Intracellular_centrality"

# DEGs
degs <- as.matrix(read.table(file = "../Input/CD/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt", header = T, stringsAsFactors = F, sep="\t"))

# Select only top 1800 most significant DEGs
if(nrow(degs)>1800){
  degs <- degs[1:1800,]
}

# run intracellular centrality analysis
intra_cent <- intracellular_drug_target_centrality(ppin = ppin,
                                                   drug_target_matrix = drugs,
                                                   degs = degs,
                                                   file_name = "CD_drugs_lit_PPIN_intracellular_centralities",
                                                   centrality_alg = "eigenvector centralities",
                                                   out_dir = out_dir)

print(head(intra_cent))

# seperate into individual drugs
intra_cent <- cbind(same_drugs[,2], intra_cent[match(same_drugs[,1], rownames(intra_cent)),])
rownames(intra_cent) <- same_drugs[,2]
write.table(intra_cent, file = "../Output/Intracellular_centrality/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", col.names = NA, row.names = T)

