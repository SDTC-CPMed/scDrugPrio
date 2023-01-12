#
#
# Test heterogeneity between patients and RA mice
#
# BY SAMUEL SCHAEFER
#
############################################################################################################################################################
library(dplyr)
library(Seurat)
library(R.filesets)
set.seed(54)
fp <- getwd()
setwd(paste(fp, "/..", sep=""))
dir.create(path = "Output/Comparison_heterogenity")

############################################################################################################################################################
# Check seperation of sick RA mice in latent space
############################################################################################################################################################
dir.data <- 'Input/DCA_adjusted_matrix/'
lat <- read.table(paste(dir.data, "/latent.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)
lat <- lat[,grepl(colnames(lat), pattern = "Sick")]

id <- vector()
for(i in 1:6){
  id[grepl(colnames(lat),pattern =paste("_mouse_",i,"_",sep=""))] <- i
}
names(id) <- colnames(lat)
id <- factor(x = id, levels = c(1:6), ordered = T)

joint_RA <- CreateSeuratObject(counts = lat, project = "Drug_prediction_denoised", min.cells = 1, min.features = 1)
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
joint_RA@reductions$pca <- pca
Idents(joint_RA) <- id

#joint_RA <- FindNeighbors(joint_RA, dims = 1:nrow(lat), k = 20)
#joint_RA <- FindClusters(joint_RA, resolution = 0.6)
joint_RA <- RunTSNE(joint_RA,dims.use = 1:nrow(lat))

pdf(file = "Output/Comparison_heterogenity/Clustering_of_sick_RA_cells_colored_by_mouse_ID.pdf", width = 5, height = 4)
TSNEPlot(joint_RA, label = F)
dev.off()

