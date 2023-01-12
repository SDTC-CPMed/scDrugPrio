#
#
# Validation of centrality by GWAS enrichment
#
# BY: SAMUEL SCHAEFER
###################################################

library(Seurat)

# Load data
degs <- read.table(file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_ID.txt", header = T)

drug_targets <- read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F)
known_RA_drugs <- read.delim(file = "../Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T, stringsAsFactors = F)[,1]
same_drugs <- read.table(file = "../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", stringsAsFactors = F)

known_RA_drugs <- unique(same_drugs[same_drugs[,1] %in% known_RA_drugs, 2])
drug_targets <- drug_targets[,colnames(drug_targets) %in% known_RA_drugs]
drug_targets <- unlist(drug_targets)
drug_targets <- unique(drug_targets[!is.na(drug_targets)])

all_genes <- readRDS(file = "../Input/DCA_adjusted_matrix/joint_RA.rds")
all_genes <- rownames(all_genes@assays$RNA@counts)

transl <- read.table(file = "../Input/Human-mouse_homologs/transl.txt", sep="\t", header = T)
all_genes <- transl[transl[,3] %in% all_genes,2]


# Fisher exact tests for checking enrichement of RA drug targets among B cell DEGs
t <- matrix(0, nrow = 2, ncol = 2)
colnames(t) <- c("RA", "Other")
rownames(t) <- c("in", "out")

# Cluster 10
cl <- which(colnames(degs) == "Cluster_10")
t[1,1] <- sum(degs[,cl] %in% drug_targets) 
t[1,2] <- sum(!is.na(degs[,cl]))-t[1,1]
t[2,1] <- length(drug_targets) - t[1,1]
t[2,2] <- length(all_genes) - sum(t[1,]) -t[2,1]
fisher.test(t, alternative = "greater", conf.int = T) # P = 0.000275

# Cluster 12
cl <- which(colnames(degs) == "Cluster_12")
t[1,1] <- sum(degs[,cl] %in% drug_targets) 
t[1,2] <- sum(!is.na(degs[,cl]))-t[1,1]
t[2,1] <- length(drug_targets) - t[1,1]
t[2,2] <- length(all_genes) - sum(t[1,]) -t[2,1]
fisher.test(t, alternative = "greater", conf.int = T) # P = 0.004733

# Cluster 1
cl <- which(colnames(degs) == "Cluster_1")
t[1,1] <- sum(degs[,cl] %in% drug_targets) 
t[1,2] <- sum(!is.na(degs[,cl]))-t[1,1]
t[2,1] <- length(drug_targets) - t[1,1]
t[2,2] <- length(all_genes) - sum(t[1,]) -t[2,1]
fisher.test(t, alternative = "greater", conf.int = T) # P = 0.9996

# Cluster 15
cl <- which(colnames(degs) == "Cluster_15")
t[1,1] <- sum(degs[,cl] %in% drug_targets) 
t[1,2] <- sum(!is.na(degs[,cl]))-t[1,1]
t[2,1] <- length(drug_targets) - t[1,1]
t[2,2] <- length(all_genes) - sum(t[1,]) -t[2,1]
fisher.test(t, alternative = "greater", conf.int = T) # P = 1
