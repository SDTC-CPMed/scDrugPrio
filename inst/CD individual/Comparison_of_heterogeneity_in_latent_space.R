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

############################################################################################################################################################
# Check seperation of sick samples from CD patients in latent space
############################################################################################################################################################
dir.data <- 'Input/CD GSE134809/CD_all/'
lat <- read.table(paste(dir.data, "/latent.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)
lat <- lat[,grepl(colnames(lat), pattern = "Involved")]

# identify patients
temp <- cbind(c("GSM3972009","GSM3972010","GSM3972011","GSM3972012",
                "GSM3972013","GSM3972014","GSM3972016","GSM3972015",
                "GSM3972017","GSM3972018","GSM3972020","GSM3972019",
                "GSM3972022","GSM3972021","GSM3972024","GSM3972023",
                "GSM3972026","GSM3972025","GSM3972028","GSM3972027",
                "GSM3972030","GSM3972029"),c("Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl"))
temp <- cbind(temp, c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11))
temp <- temp[temp[,2]%in% "Infl",]

id <- vector()
for(i in 1:nrow(temp)){
  id[grepl(colnames(lat),pattern = temp[i,1])] <- as.numeric(temp[i,3])
}
names(id) <- colnames(lat)
id <- factor(x = id, levels = c(1:11), ordered = T)

joint_CD <- CreateSeuratObject(counts = lat, project = "Drug_prediction_denoised", min.cells = 1, min.features = 1)
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
joint_CD@reductions$pca <- pca
Idents(joint_CD) <- id

#joint_RA <- FindNeighbors(joint_RA, dims = 1:nrow(lat), k = 20)
#joint_RA <- FindClusters(joint_RA, resolution = 0.6)
joint_CD <- RunTSNE(joint_CD,dims.use = 1:nrow(lat))

pdf(file = "Output/Comparison_heterogenity/Clustering_of_cells_in_sick_CD_patient_samples_colored_by_patient_ID.pdf", width = 5, height = 4)
TSNEPlot(joint_CD, label = F)
dev.off()

############################################################################################################################################################
# Check seperation of sick samples from CD patients in latent space AFTER BATCH CORRECTION
############################################################################################################################################################
rm(joint_CD)
dir.data <- 'CD_batch_corrected/Input/CD/Results_CD_all_integrated/'
lat <- read.table(paste(dir.data, "/latent.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)
lat <- lat[,grepl(colnames(lat), pattern = "Involved")]

# identify patients
temp <- cbind(c("GSM3972009","GSM3972010","GSM3972011","GSM3972012",
                "GSM3972013","GSM3972014","GSM3972016","GSM3972015",
                "GSM3972017","GSM3972018","GSM3972020","GSM3972019",
                "GSM3972022","GSM3972021","GSM3972024","GSM3972023",
                "GSM3972026","GSM3972025","GSM3972028","GSM3972027",
                "GSM3972030","GSM3972029"),c("Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl","Infl","Uninfl","Infl","Uninfl",
                                             "Infl","Uninfl"))
temp <- cbind(temp, c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11))
temp <- temp[temp[,2]%in% "Infl",]

id <- vector()
for(i in 1:nrow(temp)){
  id[grepl(colnames(lat),pattern = temp[i,1])] <- as.numeric(temp[i,3])
}
names(id) <- colnames(lat)
id <- factor(x = id, levels = c(1:11), ordered = T)

joint_CD <- CreateSeuratObject(counts = lat, project = "Drug_prediction_denoised", min.cells = 1, min.features = 1)
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
joint_CD@reductions$pca <- pca
Idents(joint_CD) <- id

#joint_RA <- FindNeighbors(joint_RA, dims = 1:nrow(lat), k = 20)
#joint_RA <- FindClusters(joint_RA, resolution = 0.6)
joint_CD <- RunTSNE(joint_CD,dims.use = 1:nrow(lat))

pdf(file = "Output/Comparison_heterogenity/Clustering_of_cells_in_sick_CD_patient_samples_colored_by_patient_ID_BATCH_CORRECTED.pdf", width = 5, height = 4)
TSNEPlot(joint_CD, label = F)
dev.off()


############################################################################################################################################################
# Check seperation of sick samples from MS patients in latent space
############################################################################################################################################################
dir.data <- 'Input/MS GSE138266/MS_Martin/MS_all/'
lat <- read.table(paste(dir.data, "/latent.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)
lat <- lat[,grepl(colnames(lat), pattern = "_MS_")]

# identify patients
temp <- t(matrix(c("GSM4104122", "MS",          
          "GSM4104123", "MS",         
          "GSM4104124", "MS",        
          "GSM4104125", "MS",       
          "GSM4104126", "MS",      
          "GSM4104127", "MS",     
          "GSM4104128", "IIH",    
          "GSM4104129", "IIH",   
          "GSM4104130", "IIH",  
          "GSM4104131", "IIH", 
          "GSM4104132", "IIH",
          "GSM4104133", "IIH"), nrow = 2))
temp <- cbind(temp, c(1:nrow(temp)))
temp <- temp[temp[,2]%in% "MS",]

id <- vector()
for(i in 1:nrow(temp)){
  id[grepl(colnames(lat),pattern = temp[i,1])] <- as.numeric(temp[i,3])
}
names(id) <- colnames(lat)
id <- factor(x = id, levels = c(1:nrow(temp)), ordered = T)

joint_MS <- CreateSeuratObject(counts = lat, project = "Drug_prediction_denoised", min.cells = 1, min.features = 1)
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
joint_MS@reductions$pca <- pca
Idents(joint_MS) <- id

#joint_RA <- FindNeighbors(joint_RA, dims = 1:nrow(lat), k = 20)
#joint_RA <- FindClusters(joint_RA, resolution = 0.6)
joint_MS <- RunTSNE(joint_MS,dims.use = 1:nrow(lat))

pdf(file = "Output/Comparison_heterogenity/Clustering_of_cells_in_sick_MS_patient_samples_colored_by_patient_ID.pdf", width = 5, height = 4)
TSNEPlot(joint_MS, label = F)
dev.off()

# Recreate real cluster plot
#####################################################################
setwd("Input/MS GSE138266/")

ms_all <- readRDS('MS_Martin/MS_all_v2.rds') # counts = DCA adjusted data

# New parameters
ms_all <- FindNeighbors(ms_all, k.param = 10)
ms_all <- FindClusters(ms_all, resolution = 0.20)

temp <- t(matrix(c("GSM4104122", "MS",          
                   "GSM4104123", "MS",         
                   "GSM4104124", "MS",        
                   "GSM4104125", "MS",       
                   "GSM4104126", "MS",      
                   "GSM4104127", "MS",     
                   "GSM4104128", "IIH",    
                   "GSM4104129", "IIH",   
                   "GSM4104130", "IIH",  
                   "GSM4104131", "IIH", 
                   "GSM4104132", "IIH",
                   "GSM4104133", "IIH"), nrow = 2))
temp <- cbind(temp, c(1:nrow(temp)))
id <- vector()
for(i in 1:nrow(temp)){
  id[grepl(colnames(lat),pattern = temp[i,1])] <- temp[i,1]
}
names(id) <- colnames(lat)
id <- factor(x = id, levels = temp[,1], ordered = T)


# Create cellID / cluster ID matrix
id2 <- cbind(as.character(names(Idents(ms_all))), as.character(Idents(ms_all)))
colnames(id2) <- c("Cell_ID", "Cluster_ID")

temp <- unlist(strsplit(id2[,1], split="_"))
temp <- temp[seq(from = 1,to = length(temp), by = 3)]
id2 <- cbind(id2, temp)

library(scales)
ms_all@active.ident <- factor(x = id2[,3], levels = levels(id), ordered = T)
#ms_all <- RunTSNE(ms_all,dims.use = 1:nrow(lat))
pdf(file = "../../Output/Comparison_heterogenity/Patient_colored_in_MS_clustering.pdf", width = 6, height = 4)
TSNEPlot(ms_all, cols = hue_pal()(12))
dev.off()
#FeaturePlot(ms_all, reduction = "tsne")

############################################################################################################################################################
# Check seperation of sick samples from MS patients in latent space AFTER BATCH CORRECTION
############################################################################################################################################################

dir.data <- 'MS_batch_corrected/Input/MS/Results_MS_all_integrated/'
lat <- read.table(paste(dir.data, "/latent.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)
lat <- lat[,grepl(colnames(lat), pattern = "_MS_")]

# identify patients
temp <- t(matrix(c("GSM4104122", "MS",          
                   "GSM4104123", "MS",         
                   "GSM4104124", "MS",        
                   "GSM4104125", "MS",       
                   "GSM4104126", "MS",      
                   "GSM4104127", "MS",     
                   "GSM4104128", "IIH",    
                   "GSM4104129", "IIH",   
                   "GSM4104130", "IIH",  
                   "GSM4104131", "IIH", 
                   "GSM4104132", "IIH",
                   "GSM4104133", "IIH"), nrow = 2))
temp <- cbind(temp, c(1:nrow(temp)))
temp <- temp[temp[,2]%in% "MS",]

id <- vector()
for(i in 1:nrow(temp)){
  id[grepl(colnames(lat),pattern = temp[i,1])] <- as.numeric(temp[i,3])
}
names(id) <- colnames(lat)
id <- factor(x = id, levels = c(1:nrow(temp)), ordered = T)

joint_MS <- CreateSeuratObject(counts = lat, project = "Drug_prediction_denoised", min.cells = 1, min.features = 1)
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
joint_MS@reductions$pca <- pca
Idents(joint_MS) <- id

#joint_RA <- FindNeighbors(joint_RA, dims = 1:nrow(lat), k = 20)
#joint_RA <- FindClusters(joint_RA, resolution = 0.6)
joint_MS <- RunTSNE(joint_MS,dims.use = 1:nrow(lat))

pdf(file = "Output/Comparison_heterogenity/Clustering_of_cells_in_sick_MS_patient_samples_colored_by_patient_ID_BATCH_CORRECTED.pdf", width = 5, height = 4)
TSNEPlot(joint_MS, label = F)
dev.off()



