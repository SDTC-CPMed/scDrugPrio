
# INTRACELLULAR CENTRALITY
# Run in R 4.0.4
######################################################################

#setwd("/home/samsc76/scDrugPrio/SIGMA_201213/Doctis_anti_IL17_v1/R")
dir.create("../Output/Intracellular_centrality", showWarnings = F)
source("Codes_from_scPred/intracellular_drug_centrality.R")
source("Codes_from_scPred/extract_LCC_by_gene_set_of_interest.R")
source("Codes_from_scPred/ppin_formatting.R")
source("Codes_from_scPred/extract_LCC.R")

max_degs <- 3000

#PPIN
lit_ppi <- as.matrix(read.table(file = "../Input/literature_PPI/ppi.txt", sep ="\t", header = T, stringsAsFactors = F))
ppin <- ppin_formatting(lit_ppi)

# transl
transl <- read.table(file = "../Input/HGNC translation matrix 201108/transl.txt", header = T, stringsAsFactors = F, sep="\t")
transl <- transl[,c(2,6)]

# PsA drugs
psa_drug_info <- read.table(file = "../Input/Drugs/PsA_drugs_from_DrugBank.txt", header = T, sep="\t", stringsAsFactors = F)
psa_drugs <- unique(psa_drug_info[,1])

# Drug target matrix
drugs <- as.matrix(read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt", header = T, stringsAsFactors = F, sep="\t"))
same_drugs <- read.table(file = "../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
same_drugs <- same_drugs[same_drugs[,1] %in% colnames(drugs),]
psa_drug_info[psa_drug_info[,1]%in% same_drugs[,2],5] <- "x"
write.table(psa_drug_info, file = "../Input/Drugs/PsA_drugs_from_DrugBank.txt", sep="\t", col.names = T, row.names = F)

# Check if drugs have valid target in PPIN
for(i in 1:ncol(drugs)){
  print(paste(colnames(drugs)[i], "has valid targets:", all(drugs[!is.na(drugs[,i]),i] %in% as.vector(ppin)), sep=" "))
}


# Responder
########################################

# out_dir
out_dir <- "../Output/Intracellular_centrality/Responder"
dir.create(out_dir, showWarnings = F)

# DEGs
degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_R.txt", header = T, stringsAsFactors = F, sep="\t"))

# make sure ppin and degs have the same gene annotation system so they can be matched
for(i in 1:ncol(degs)){
  temp <- transl[match(x = degs[,i], table =  transl[,1]),2]
  temp <- temp[!is.na(temp)]
  degs[,i] <- NA
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
degs <- degs[rowSums(!is.na(degs))>0, colSums(!is.na(degs))>0]
mode(degs) <- "numeric"
write.table(degs, file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_R_TRANSLATED.txt", col.names = T, row.names = F, sep = "\t")

# Select only top max_degs most significant DEGs
if(nrow(degs)>max_degs){
  degs <- degs[1:max_degs,]
}

# run intracellular centrality analysis
dir.create("../Output/Intracellular_centrality", showWarnings = F)
intra_cent <- intracellular_drug_target_centrality(ppin = ppin,
                                                   drug_target_matrix = drugs,
                                                   degs = degs,
                                                   file_name = "PsA_drugs_lit_PPIN_intracellular_centralities",
                                                   centrality_alg = "eigenvector centralities",
                                                   out_dir = out_dir)

print(head(intra_cent))

# seperate into individual drugs
intra_cent <- cbind(same_drugs[,2], intra_cent[match(same_drugs[,1], rownames(intra_cent)),])
rownames(intra_cent) <- same_drugs[,2]
write.table(intra_cent, file = "../Output/Intracellular_centrality/Responder/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", col.names = NA, row.names = T)

# Non-Responder
########################################

# out_dir
out_dir <- "../Output/Intracellular_centrality/Non_responder"
dir.create(out_dir, showWarnings = F)

# DEGs
degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_NR.txt", header = T, stringsAsFactors = F, sep="\t"))

# make sure ppin and degs have the same gene annotation system so they can be matched
for(i in 1:ncol(degs)){
  temp <- transl[match(x = degs[,i], table =  transl[,1]),2]
  temp <- temp[!is.na(temp)]
  degs[,i] <- NA
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
degs <- degs[rowSums(!is.na(degs))>0, colSums(!is.na(degs))>0]
mode(degs) <- "numeric"
write.table(degs, file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_NR_TRANSLATED.txt", col.names = T, row.names = F, sep = "\t")

# Select only top max_degs most significant DEGs
if(nrow(degs)>max_degs){
  degs <- degs[1:max_degs,]
}

# run intra-cellular centrality analysis
intra_cent <- intracellular_drug_target_centrality(ppin = ppin,
                                                   drug_target_matrix = drugs,
                                                   degs = degs,
                                                   file_name = "PsA_drugs_lit_PPIN_intracellular_centralities",
                                                   centrality_alg = "eigenvector centralities",
                                                   out_dir = out_dir)
print(head(intra_cent))

# seperate into individual drugs
intra_cent <- cbind(same_drugs[,2], intra_cent[match(same_drugs[,1], rownames(intra_cent)),])
rownames(intra_cent) <- same_drugs[,2]
write.table(intra_cent, file = "../Output/Intracellular_centrality/Non_responder/INDIVIDUAL_DRUGS_intracellular_centrality_lit_PPIN.txt", sep = "\t", col.names = NA, row.names = T)


