
# NICHENET ANALYSIS
# Run on R 4.0.4
######################################################################

# Load clusteres cd_all data
#setwd("/data/samsc76/SIGMA_13012020/CD_batch_corrected/drug_prediction_R")
cd_all <- readRDS(file = "../Input/CD/CD_batch_corrected_clustered.rds")
rm(list = ls()[ls() != "cd_all"])

source("Codes_from_scPred/background_genes_NicheNet.R")
source("Codes_from_scPred/NicheNet_cell_type_centrality.R")
source("Codes_from_scPred/NicheNet_ligand_activity_analysis.R")
library(doParallel)
library(Seurat)

id <- Idents(cd_all)
id <- cbind(cell_ID = as.character(names(id)), cluster_id = as.numeric(as.character(id)))


# format cell IDs to fit in function
cell_IDs <- cd_all$seurat_clusters
unique_cl <- unique(as.numeric(as.character(cell_IDs)))
cell_IDs <- foreach(i = c(1:length(unique_cl))) %do% {
  temp <- cell_IDs[as.character(cell_IDs) == unique_cl[i]]
  return(names(temp))
}
names(cell_IDs) <- paste("Cluster_", unique_cl, sep="")

# Set out_dir
out_dir <- "../Output/NicheNet"
dir.create(path = out_dir, showWarnings = F)

# background gene calculation
bg <- background_genes_NicheNet(data = as.matrix(cd_all@assays$RNA@counts), cell_IDs = cell_IDs)

# translate bg
transl <- read.table(file = "../Input/HGNC translation matrix 201108/transl.txt", header = T, sep = "\t", stringsAsFactors = F)
transl <- transl[,c(7,2)]
transl <- transl[!is.na(transl[,1]),]
transl <- transl[!is.na(transl[,2]),]
for(i in 1:ncol(bg)){
  temp <- bg[!is.na(bg[,i]),i]
  temp <- transl[transl[,1]%in% temp,2]
  bg[,i] <- NA
  if(length(temp)>0){
    bg[1:length(temp),i] <- temp
  }
  rm(temp)
}
write.table(bg, file = paste(out_dir, "/background_genes.txt", sep=""), sep="\t", col.names = T, row.names = F)

# DEGs
degs <- as.matrix(read.table(file = "../Input/CD/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_HGNC_SYMBOL_sig_MAST_DEGs_log(1.5)_all_clusters.txt", header = T, stringsAsFactors = F, sep="\t"))
degs <- degs[,colSums(!is.na(degs))>0]

# NicheNet ligand activity
all_ligand_activity <- NicheNet_ligand_activity_analysis(degs = degs, background_genes = bg, out_dir = out_dir, cores = 30)
print(head(all_ligand_activity))

#Centrality analysis
intercellular_centrality <- NicheNet_cell_type_centrality(all_ligand_activity, out_dir = out_dir)
print(intercellular_centrality)

