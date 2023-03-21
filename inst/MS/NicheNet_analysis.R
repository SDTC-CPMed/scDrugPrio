
# NICHENET ANALYSIS
# Run on R 4.0.4
######################################################################

# Load clusteres ms_all data
#setwd("/data/samsc76/SIGMA_13012020/MS_batch_corrected/drug_prediction_R")
ms_all <- readRDS(file = "../Input/MS/MS_batch_corrected_clustered.rds")
rm(list = ls()[ls() != "ms_all"])

source("Codes_from_scPred/background_genes_NicheNet.R")
source("Codes_from_scPred/NicheNet_cell_type_centrality.R")
source("Codes_from_scPred/NicheNet_ligand_activity_analysis.R")
library(doParallel)
library(Seurat)

id <- Idents(ms_all)
id <- cbind(cell_ID = as.character(names(id)), cluster_id = as.numeric(as.character(id)))


# format cell IDs to fit in function
cell_IDs <- ms_all$seurat_clusters
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
bg <- background_genes_NicheNet(data = as.matrix(ms_all@assays$RNA@counts), cell_IDs = cell_IDs)

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
degs <- as.matrix(read.table(file = "../Input/MS/DCA_MAST_DEGs_logFC_1.5_MS/TRANSLATED_HGNC_SYMBOL_sig_MAST_DEGs_log(1.5)_all_clusters.txt", header = T, stringsAsFactors = F, sep="\t"))
degs <- degs[,colSums(!is.na(degs))>0]
#degs <- degs[1:1800,] # limiting the number of DEGs

# NicheNet ligand activity
all_ligand_activity <- NicheNet_ligand_activity_analysis(degs = degs, background_genes = bg, out_dir = out_dir, cores = 30)
print(head(all_ligand_activity))

#Centrality analysis
intercellular_centrality <- NicheNet_cell_type_centrality(all_ligand_activity = all_ligand_activity, out_dir = out_dir)
print(intercellular_centrality)

###############################################################################
# Output for Cytoscape
###############################################################################
all_ligand <- read.table(file = "../Output/NicheNet/all_ligand_activity.txt", sep = "\t", header = T, stringsAsFactors = F)
all_ligand <- all_ligand[as.numeric(all_ligand[,4]) > 0,]
color <- read.table(file = "../Output/MS/Cluster_colors.txt", sep="\t", header = T, stringsAsFactors = F)


cl1 <- unique(all_ligand[,5])
cl2 <- unique(all_ligand[,6])
cl <- unique(c(cl1, cl2))

# summed pearsson coefficients for plotting
out_ligand <- foreach(s = c(1:length(cl1)), .combine = "rbind") %do% {
  if(any(all_ligand[,5] %in% cl1[s])){
    out <- foreach(t = c(1:length(cl2)), .combine = "rbind") %do% {
      if(any((all_ligand[,6] %in% cl2[t]) & (all_ligand[,5] %in% cl1[s]))){
        if(sum((all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t])) > 1){
          temp <- all_ligand[(all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t]),]
          temp <- c(NA, NA, nrow(temp), sum(as.numeric(temp[,4])), cl1[s], cl2[t])
          names(temp) <- NULL
          return(matrix(temp, nrow = 1))
        } else if(sum((all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t])) == 0){
          return(rep(NA, times = ncol(all_ligand)))
        } else {
          temp <- as.vector(all_ligand[(all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t]),])
          temp[,3] <- 1 # n of ligand interactions
          names(temp) <- NULL
          return(matrix(temp, nrow = 1))
        }
      } else {
        return(matrix(rep(NA, times = ncol(all_ligand)), nrow = 1))
      }
    }
    colnames(out) <- NULL
    return(out)
  } else {
    return(matrix(rep(NA, times = ncol(all_ligand)), nrow = 1))
  }
}
colnames(out_ligand) <- colnames(all_ligand)
colnames(out_ligand)[3] <- "n_ligand_interactions"
out_ligand <- out_ligand[!is.na(out_ligand[,5]),]

# add node colors through extra rows
temp <-cbind(rep("node_color",times = length(cl)), NA, NA, NA, cl, NA)
colnames(temp) <- colnames(out_ligand)
out_ligand <- rbind(out_ligand, temp)
out_ligand <- cbind(out_ligand, color[match(out_ligand[,5], paste("Cluster_",as.numeric(color[,1]),sep="")),4])
colnames(out_ligand)[7] <- "color"

out_ligand <- as.data.frame(out_ligand[,-2])
out_ligand$test_ligand <- as.character(out_ligand$test_ligand)
out_ligand$n_ligand_interactions <- as.numeric(out_ligand$n_ligand_interactions)
out_ligand$pearson <- as.numeric(out_ligand$pearson)
out_ligand$Sender <- as.character(out_ligand$Sender)
out_ligand$Target <- as.character(out_ligand$Target)
out_ligand$color <- as.character(out_ligand$color)
out_ligand$color2 <- as.character(out_ligand$color)

write.table(x = out_ligand, file = "../Output/NicheNet/MCDM_MS_own_colors_221027_pearson_summed.txt", sep = "\t", col.names = T, row.names = F)
