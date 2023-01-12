
# INTRACELLULAR CENTRALITY
# Run in R 4.0.4
######################################################################

#setwd("/data/samsc76/Doctis_data_v3/R")
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

# PsA drugs
ra_drug_info <- read.table(file = "../Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", header = T, sep="\t", stringsAsFactors = F, quote = "")
ra_drugs <- unique(ra_drug_info[,1])

# Drug target matrix
drugs <- read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt", header = T, stringsAsFactors = F, sep="\t")
same_drugs <- read.table(file = "../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
ra_drug_info <- cbind(ra_drug_info, in_DrugBank = NA)
ra_drug_info[ra_drug_info[,1]%in% same_drugs[,2],4] <- "x"
write.table(ra_drug_info, file = "../Input/Drugs/ra_drugs_from_DrugBank.txt", sep="\t", col.names = T, row.names = F)

drugs <- drugs[,colnames(drugs)%in% same_drugs[,1]]
drugs <- as.matrix(drugs)

# Check if drugs have valid target in PPIN
for(i in 1:ncol(drugs)){
  print(paste(colnames(drugs)[i], "has valid targets:", all(drugs[!is.na(drugs[,i]),i] %in% as.vector(ppin)), sep=" "))
}

# Responder
########################################

# DEGs
degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_ID.txt", header = T, stringsAsFactors = F))
mode(degs) <- "numeric"

# run intracellular centrality analysis
dir.create("../Output/Intracellular_centrality", showWarnings = F)
intra_cent <- intracellular_drug_target_centrality(ppin = ppin,
                                                   drug_target_matrix = drugs,
                                                   degs = degs,
                                                   file_name = "RA_drugs_lit_PPIN_intracellular_centralities",
                                                   centrality_alg = "eigenvector centralities",
                                                   out_dir = "../Output/Intracellular_centrality/")

print(head(intra_cent))

# seperate into individual drugs
intra_cent <- cbind(same_drugs[,2], intra_cent[match(same_drugs[,1], rownames(intra_cent)),])
rownames(intra_cent) <- same_drugs[,2]
write.table(intra_cent, file = "../Output/Intracellular_centrality/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", col.names = NA, row.names = T)

