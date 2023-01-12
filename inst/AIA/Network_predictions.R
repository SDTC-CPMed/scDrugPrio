
# NETWORK-BASED DRUG SCREENING
# Run in R 4.0.4
######################################################################

dir.create("../Output/Network_distances", showWarnings = F)
source("Codes_from_scPred/average_closest_distance_network_drug_screening.R")
source("Codes_from_scPred/average_closest_distance.R")
source("Codes_from_scPred/bin_creation_by_min_bin_size.R")
source("Codes_from_scPred/random_drug_target_bin_adjusted_distances.R")
source("Codes_from_scPred/extract_LCC_by_gene_set_of_interest.R")
source("Codes_from_scPred/ppin_formatting.R")
source("Codes_from_scPred/extract_LCC.R")
source("Codes_from_scPred/combine_evaluation_files.R")

# INPUT
lit_ppi <- as.matrix(read.table(file = "../Input/literature_PPI/ppi.txt", sep ="\t", header = T, stringsAsFactors = F))
ppin <- ppin_formatting(lit_ppi)
same_drugs <- read.table(file ="../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)
drugs <- read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F)

set.seed(23)

# LITERATURE PPIN
########################################

out_dir <- "../Output/Network_distances"
dir.create(out_dir, showWarnings = F)

degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_ID.txt", header = T, stringsAsFactors = F))
# if(nrow(degs)>3000){
#   degs <- degs[1:3000,]
# }

for(i in 1:ncol(degs)){ # For every cell types DEGs
  print(colnames(degs)[i])
  # calculate average closest distances between all drugs and the degs
  # if possible on your machine this should be parallelized (cores > 1), though the function will be memory heavy
  average_closest_distance_network_drug_screening(ppin = ppin,
                                                  drug_target_matrix = drugs,
                                                  disease_genes = as.numeric(degs[,i]),
                                                  file_name = colnames(degs)[i],
                                                  cores = 30,
                                                  out_dir = out_dir)
}

# select output files based on file names
lf <- list.files(path = out_dir, pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__Cluster_", full.names = T)
# labels
label <- unlist(strsplit(x = lf, split = "drug-disease_closest_distances_vs_random_bin_adjusted__"))
label <- label[seq(from = 2, to = length(label), by = 2)]
label <- unlist(strsplit(x = label, split = ".txt"))

out <- combine_evaluation_files(files = lf, label_for_files = label, output_file_name = "SUMMARY_drug_dist.txt", cores = 1)
print(head(out))


# HuRI
########################################
# 
# out_dir <- "../Output/Network_distances/HuRI"
# dir.create(out_dir, showWarnings = F)
# 
# degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_NR_TRANSLATED.txt", header = T, stringsAsFactors = F, sep = "\t"))
# # if(nrow(degs)>3000){
# #   degs <- degs[1:3000,]
# # }
# 
# for(i in 1:ncol(degs)){ # For every cell types DEGs
#   print(colnames(degs)[i])
#   # calculate average closest distances between all drugs and the degs
#   # if possible on your machine this should be parallelized (cores > 1), though the function will be memory heavy
#   average_closest_distance_network_drug_screening(ppin = ppin,
#                                                   drug_target_matrix = drugs,
#                                                   disease_genes = as.numeric(degs[,i]),
#                                                   file_name = colnames(degs)[i],
#                                                   cores = 14,
#                                                   out_dir = out_dir)
# }
# 
# # select output files based on file names
# lf <- list.files(path = out_dir, pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__Cluster_", full.names = T)
# # labels
# label <- unlist(strsplit(x = lf, split = "drug-disease_closest_distances_vs_random_bin_adjusted__"))
# label <- label[seq(from = 2, to = length(label), by = 2)]
# label <- unlist(strsplit(x = label, split = ".txt"))
# 
# out <- combine_evaluation_files(files = lf, label_for_files = label, output_file_name = paste(out_dir,"/SUMMARY_drug_dist.txt",sep=""), cores = 1)
# print(head(out))
# 
# 
# 
