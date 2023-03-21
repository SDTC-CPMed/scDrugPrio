
# Run in R 4.0.4 due to cluster limitations.

#setwd("/data/samsc76/Doctis_data_v3/R")
#install.packages("Seurat")
library(Seurat)
library(doParallel)
registerDoParallel(cores = 15)
set.seed(21)

# get DCA output
###############################################################################

psa <- readRDS("../Input/PsA_data/PsA_all_integrated_singlet_patients.rds")
rownames(psa@assays$RNA@counts)[1:5]

meta <- psa@meta.data
print("Healthy_controls:")
print(unique(paste(meta[,5], "_", meta[,12],sep="")))
print("PsA patients: response_patientID")
print(unique(paste(meta[,26], "_", meta[,13],sep=""))[order(unique(paste(meta[,26], "_", meta[,13],sep="")))])
meta <- cbind(rownames(meta), meta$Patients, meta$reference, meta$response, meta$seurat_clusters)


#psa <- DietSeurat(object = psa, graphs = NULL, )
psa <- RunUMAP(psa, dims = 1:32)

dir.create("../Output/Clustering", showWarnings = F)
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.4.pdf")
TSNEPlot(psa, group.by = "RNA_snn_res.0.4")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.4_healthy.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "healthy", 1], group.by = "RNA_snn_res.0.4")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.4_sick.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "sick", 1], group.by = "RNA_snn_res.0.4")
dev.off()

# get number of cells / cluster
###############################################################################

get_n_cells <- function(){
  id <- Idents(psa)
  meta[,5] <- as.character(id[match(meta[,1], names(id))])
  
  out <- foreach(i =c(1:length(unique(id))), .combine = "cbind" )%dopar%{
    t <- vector()
    temp <- meta[meta[,5] %in% (as.character(unique(id))[i]),3:4]
    t[1] <- sum(temp[,1] %in% "healthy")
    temp <- temp[temp[,1] %in% "sick",2]
    t[2] <- sum(temp %in% "NR")
    t[3] <- sum(temp %in% "R")
    return(t)
  }
  colnames(out) <- paste("Cluster_",as.character(unique(id)), sep="")
  rownames(out) <- c("Healthy", "NR", "R")
  return(out)
}
out <- get_n_cells()
out

dir.create("../Output/Clustering", showWarnings = F)
write.table(out, file = "../Output/Clustering/Martin_initial_clustering.txt", sep="\t", col.names = NA, row.names = T)

# change cluster parameters
###############################################################################

# k = 15 and resolution = 0.4 in the above loaded RDS
psa <- FindNeighbors(psa, k.param = 15, nn.method = "annoy")
psa <- FindClusters(psa, resolution = 0.45)
out <- get_n_cells()
out
write.table(out, file = "../Output/Clustering/Clustering_k=15_res=0.45.txt", sep="\t", col.names = NA, row.names = T)
id <- Idents(psa)
id <- cbind(cell_id = names(id), cluster_id = as.character(id))
write.table(id, file = "Cell_identity_k=15_res=0.45.txt",sep="\t", col.names = T, row.names = F)

# remake plots
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45.pdf")
TSNEPlot(psa)
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_healthy.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "healthy", 1])
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_sick.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "sick", 1])
dev.off()

pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_UMAP.pdf")
UMAPPlot(psa)
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_healthy_UMAP.pdf")
UMAPPlot(psa, cells = meta[meta[,3] %in% "healthy", 1])
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_sick_UMAP.pdf")
UMAPPlot(psa, cells = meta[meta[,3] %in% "sick", 1])
dev.off()

saveRDS(object=psa, file="../Input/PsA_data/PsA_all_clustered.rds")
#psa <- readRDS(file = "../Input/PsA_data/PsA_all_clustered.rds")

# Cell typing
###############################################################################

library(matrixStats)
library(reshape2)
library(ggplot2)
library(viridis)
library(textshape)
library(doParallel)

# load marker genes
marker_genes <- read.table(file = "../Input/Marker_genes/CD_cell_type_markers_original_article.txt",sep="\t", header = T, stringsAsFactors = F)[,2]
marker_genes <- c("TRAC", marker_genes)

# extract count-adjusted expression values of marker_genes
X <- psa@assays$RNA@counts
X <- as.matrix(X)
temp <- colSums(X)
X <- t(t(X)/temp)
X <- X[rownames(X)%in% marker_genes,]
X <- as.matrix(X)

# cluster order
id <- Idents(psa)
id <- cbind(cell_ID = names(id), cluster_ID = as.numeric(as.character(id)))
id[,2] <- floor(as.numeric(as.character(id[,2])))

temp_id <- as.factor(id[,2])
names(temp_id) <- id[,1]
psa@active.ident <- temp_id

unique_cl <- unique(as.numeric(as.character(psa@active.ident)))
 

# Summarize expression scores
out <- foreach(cl = c(1:length(cluster_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%marker_genes, colnames(X)%in%id[id[,2]==cluster_order[cl],1]]
  exp <- t(exp[match(marker_genes, rownames(exp)),])
  exp <- cluster_matrix(exp, dim = "row")
  return(rbind(exp, matrix(NA, nrow = 75, ncol = ncol(exp))))
}
colnames(out) <- marker_genes
out <- out[,colSums(!is.na(out))>0]
rownames(out) <- nrow(out):1

# Create heatmap for all markers
mode(out) <- "numeric" 
print(floor(min(out[!is.na(out)])*10^6)/10^6) # 1e-6
print(floor(max(out[!is.na(out)])*1000)/1000) # 0.007

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- nrow(out)
for(i in (nrow(out)-1):2){
  if(i%in%temp){
    if(!((i-1) %in% temp)){
      pos <- c(pos, i)
    }
    if(!((i+1)%in%temp)){
      pos <- c(pos,i)
    }
  }
}
y_breaks <- vector()
for(i in 1:length(cluster_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cluster_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = 1:max(out$Var1), ordered = T)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")
out$value[out$value < 1e-5] <- 1e-5
out$value[out$value > 1e-3] <- 1e-3

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(10)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-5,1e-3)) +
  scale_y_discrete(breaks = floor(y_breaks), labels = y_labels) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_blank(), # bg of the panel
    plot.background = element_blank(), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "../Output/Clustering/Expression_scores_for_all_CD_cell_type_markers.pdf", width = 25, height = 12)
heatmap
dev.off()

png(file = "../Output/Clustering/Expression_scores_for_all_CD_cell_type_markers.png", width = 35, height = 17, units = "cm", res = 300)
heatmap
dev.off()


# Mean z-score expression heatmap
###############################################################################

# calculate mean expression per cluster
out <- foreach(cl = c(1:length(unique_cl)), .combine = "rbind") %do% {
  exp <- X[, colnames(X)%in%id[id[,2] == unique_cl[cl],1]]
  exp <- exp[match(marker_genes, rownames(exp)),]
  return(rowMeans(exp))
}
colnames(out) <- marker_genes
rownames(out) <- unique_cl
out <- out[,colSums(!is.na(out))>0]

# Z - score
mode(out) <- "numeric"
temp <-colSds(out)
out <- t(t(out) - colMeans(out)) # center around 0
out <- t(t(out)/temp) # make z-scores

out <- melt(out)
out$Var1 <- as.character(out$Var1)
out$Var1 <- factor(x = paste("Cluster_",as.character(out$Var1),sep=""), levels = paste("Cluster_",unique_cl[match(rev(cluster_order), unique_cl)],sep=""), ordered = T)

thresh <-1.5

my_breaks <- seq(from = -thresh, to = thresh, by = 0.5)
my_labels <- seq(from = -thresh, to = thresh, by = 0.5)

out[out[,3] < -thresh,3] <- -thresh
out[out[,3] > thresh,3] <- thresh

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) +
  geom_raster(na.rm = T) + scale_fill_gradientn(colours = viridis(12), name = "count", breaks = my_breaks, labels = my_labels, limits = c(-thresh,thresh)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) #, axis.text.y=as.character(), axis.ticks.y=element_blank())

pdf(file = "../Output/Clustering/Mean_expression_scores_for_all_CD_cell_type_markers.pdf", width = 25, height = 12)
heatmap
dev.off()

png(file = "../Output/Clustering/Mean_expression_scores_for_all_CD_cell_type_markers.png", width = 35, height = 12, units = "cm", res = 300)
heatmap
dev.off()


cluster_order <- c(9,24,5,10,20,17,0,4,11,6,35,22,23,15,21,31, #CD4
                   7,1,18,26,19, # CD8
                   32, # ILC
                   12,8,30,29, # B
                   3,16,27,13,14,2,33,28,25, #MNP
                   34) #DC
cluster <- cbind(cluster_order, 
                 c(rep("CD4+", times = 16), rep("CD8+", times = 5), "ILC", rep("B", times = 4), rep("MNP", times = 9), "DC"),
                 c(rep("Naive CD4+ T+", times = 10), rep("Early activated CD4+ T", times = 6), 
                   rep("Naive CD8+", times = 3), rep("Early activated CD8+ T", times = 2),
                   "ILC",
                   rep("Naive B", times = 4),
                   rep("Neutrophil", times = 6), rep("Monocyte", times = 2), "NK",
                   "DC")
                 )
colnames(cluster) <- c("Cluster_ID", "Major_cell_type", "Minor_cell_type")

write.table(cluster, file = "../Output/Clustering/Cell_typing_Doctis_v5_anti_TNF.txt",sep="\t", col.names = T, row.names = F)

# DEG calculations
###############################################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MAST")

id <- Idents(psa)
id <- cbind(cell_ID = as.character(names(id)), cluster_id = floor(as.numeric(as.character(id))), degs_id = floor(as.numeric(as.character(id))))
unique_cl <- sort(as.numeric(unique(id[,2])))

meta <- psa@meta.data
meta <- cbind(rownames(meta), meta$Patients, meta$reference, meta$response, meta$seurat_clusters)
id[id[,1] %in% meta[meta[,3] %in% "sick", 1],3] <- as.numeric(id[id[,1] %in% meta[meta[,3] %in% "sick", 1],3]) + 0.6
id[id[,1] %in% meta[meta[,4] %in% "NR", 1],3] <- as.numeric(id[id[,1] %in% meta[meta[,4] %in% "NR", 1],3]) - 0.5 # NR = X.1
#id[id[,1] %in% meta[meta[,4] %in% "R", 1],3] <- as.numeric(id[id[,1] %in% meta[meta[,4] %in% "R", 1],3]) + 0.5 #R = X.6 

# set new ID in psa
temp_id <- as.factor(id[,3])
names(temp_id) <- id[,1]
psa@active.ident <- temp_id

# test VlnPlots
pdf(file = "../Output/Clustering/VlnPlot_TNF.pdf", width = 25, height = 12)
VlnPlot(psa, features = "TNF", slot = "counts", pt.size = 0)
dev.off()
pdf(file = "../Output/Clustering/VlnPlot_NR3C1.pdf", width = 25, height = 12)
VlnPlot(psa, features = "NR3C1", slot = "counts", pt.size = 0)
dev.off()
pdf(file = "../Output/Clustering/VlnPlot_NR3C1.pdf", width = 25, height = 12)
VlnPlot(psa, features = "NR3C1", slot = "counts", pt.size = 0)
dev.off()
pdf(file = "../Output/Clustering/VlnPlot_NR3C1.pdf", width = 25, height = 12)
VlnPlot(psa, features = "NR3C1", slot = "counts", pt.size = 0)
dev.off()
pdf(file = "../Output/Clustering/VlnPlot_NR3C1.pdf", width = 25, height = 12)
VlnPlot(psa, features = "NR3C1", slot = "counts", pt.size = 0)
dev.off()

# unique deg_ids
temp_id_unique <- sort(as.numeric(unique(as.character(temp_id))))
print(temp_id_unique)

# Exclude clusters with less than 3 healthy or sick cells from calculation
# => excludes no cluster
print(table(Idents(psa)))

# Output directory
dir.create(path="../Output/DCA_MAST_DEGs", showWarnings = F)

for(i in 1:(length(unique_cl))){
  print(paste("Checking cluster:", unique_cl[i],sep=" "))
  if(any(temp_id %in% unique_cl[i])){
    if(any(temp_id %in% temp_id[temp_id %in% (unique_cl[i]+0.1)])){
      # Non-responders
      ########################################
      #if(!unique_cl%in%c()){ # exclude clusters as they have no genes with logFC.threshold > log(1.5)
      temp <- FindMarkers(psa,
                          slot = "counts",
                          test.use = "MAST", 
                          ident.1 = temp_id[temp_id %in% (unique_cl[i]+0.1)],
                          ident.2 = temp_id[temp_id %in% unique_cl[i]], # important that this is a factor vector with cell labels as names
                          only.pos = F,
                          random.seed = 3,
                          pseudocount.use = 0,
                          logfc.threshold = log(1.5),
                          min.pct = 0.1,
                          min.cells.group = 3)
    }
    
    if(exists("temp")){
      print("Saving ...")
      write.table(temp, file = paste("../Output/DCA_MAST_DEGs/Cluster_", unique_cl[i], "_NR_res=0.45_k=15.txt",sep=""),sep="\t",col.names=NA,row.names=T)
      rm(temp)
    }
    
    print("Responder:")
    # Responders
    ########################################
    if(any(temp_id %in% temp_id[temp_id %in% (unique_cl[i]+0.6)])){ # exclude clusters as they have no genes with logFC.threshold > log(1.5)
      if(unique_cl[i] != 26){ # no DEGs with logFC > log(1.5) which causes an error
        temp <- FindMarkers(psa,
                            slot = "counts",
                            test.use = "MAST", 
                            ident.1 = temp_id[temp_id %in% (unique_cl[i]+0.6)], # important that this is a factor vector with cell labels as names
                            ident.2 = temp_id[temp_id %in% unique_cl[i]],
                            only.pos = F,
                            random.seed = 3,
                            pseudocount.use = 0,
                            logfc.threshold = log(1.5),
                            min.pct = 0.1,
                            min.cells.group = 3)
        
      }
    }
    
    if(exists("temp")){
      print("Saving ...")
      write.table(temp, file = paste("../Output/DCA_MAST_DEGs/Cluster_", unique_cl[i], "_R_res=0.45_k=15.txt",sep=""),sep="\t",col.names=NA,row.names=T)
      rm(temp)
    }
  }
}


# DEG summary - non-responder
###############################################################################

outdir = "../Output/DCA_MAST_DEGs"
files <- list.files(path = outdir, pattern = "_NR_res=0.45_k=15.txt")

degs.clusters.h.vs.s <- matrix(NA, nrow = 35000, ncol = length(files))

for(i in 1:ncol(degs.clusters.h.vs.s)){
  temp <- read.table(file = paste(outdir, "/", files[i], sep=""), sep="\t", header = T)
  temp <- temp[temp$p_val_adj < 0.05,]
  temp <- as.character(temp[order(as.numeric(temp$p_val_adj), decreasing = F),1]) # order by significance
  if(length(temp)>0){
    degs.clusters.h.vs.s[1:length(temp),i] <- temp
  }
  rm(temp)
}
colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = "_NR_res=0.45_k=15.txt"))

# remove empty rows and columns
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,]
head(degs.clusters.h.vs.s)

# save
write.table(degs.clusters.h.vs.s, file = paste(outdir, "/SUMMARY_MAST_DEGs_all_clusters_NR.txt", sep=""),sep="\t", col.names = T, row.names = F)


# DEG summary - responder
###############################################################################

outdir = "../Output/DCA_MAST_DEGs"
files <- list.files(path = outdir, pattern = "_R_res=0.45_k=15.txt")

degs.clusters.h.vs.s <- matrix(NA, nrow = 35000, ncol = length(files))

for(i in 1:ncol(degs.clusters.h.vs.s)){
  temp <- read.table(file = paste(outdir, "/", files[i], sep=""), sep="\t", header = T)
  temp <- temp[temp$p_val_adj < 0.05,]
  temp <- as.character(temp[order(as.numeric(temp$p_val_adj), decreasing = F),1]) # order by significance
  if(length(temp)>0){
    degs.clusters.h.vs.s[1:length(temp),i] <- temp
  }
  rm(temp)
}
colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = "_R_res=0.45_k=15.txt"))

# remove empty rows and columns
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,]
head(degs.clusters.h.vs.s)

# save
write.table(degs.clusters.h.vs.s, file = paste(outdir, "/SUMMARY_MAST_DEGs_all_clusters_R.txt", sep=""),sep="\t", col.names = T, row.names = F)

# Find alternative gene targets
genes <- rownames(x = psa@assays$RNA@counts)
genes <- genes[grepl(genes, pattern = "TNF")]

# Check DCA denoised expression of "genes" (TNF related genes) in clusters seperated by non-responder / responder / healthy
###############################################################################
library(textshape)
library(reshape2)
library(ggplot2)
library(viridis)

# extract count-adjusted expression values of marker_genes
X <- psa@assays$RNA@counts
X <- as.matrix(X)
temp <- colSums(X)
X <- t(t(X)/temp)
X <- X[rownames(X)%in% genes,]
X <- as.matrix(X)

# get meta information 
meta <- psa@meta.data
meta <- meta[,c(4,12,13,26)]

# divide cells by group (NR / R / healthy)
healthy <- rownames(meta[meta[,2] %in% "healthy",])
nr <- meta[meta[,2] %in% "sick",]
nr <- rownames(meta[meta[,4] %in% "NR",])
r <- meta[meta[,2] %in% "sick",]
r <- rownames(meta[meta[,4] %in% "R",])

# cluster order
id <- Idents(psa)
id <- cbind(cell_ID = names(id), cluster_ID = floor(as.numeric(as.character(id))))
id <- cbind(id, group = NA, cluster = id[,2])

# For dividing clusters by NR / R / healthy
# id[id[,1] %in% healthy,3] <- "Ctrl"
# id[id[,1] %in% nr,3] <- "NR"
# id[id[,1] %in% r,3] <- "R"
# 
# id <- cbind(id, paste(id[,2],"_", id[,3],sep=""))

# Cluster order
unique_cl <- unique(id[,4])
cluster_order <- unique_cl[order(unique_cl)]

# Summarize expression scores
out <- foreach(cl = c(1:length(cluster_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%genes, colnames(X)%in%id[id[,4]==cluster_order[cl],1]]
  exp <- t(exp[match(genes, rownames(exp)),])
  exp <- cluster_matrix(exp, dim = "row")
  return(rbind(exp, matrix(NA, nrow = 75, ncol = ncol(exp))))
}
colnames(out) <- genes
out <- out[,colSums(!is.na(out))>0]
rownames(out) <- nrow(out):1

# Create heatmap for all markers
mode(out) <- "numeric" 
print(floor(min(out[!is.na(out)])*10^7)/10^7) # 7e-7
print(floor(max(out[!is.na(out)])*1000)/1000) # 0.001

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- nrow(out)
for(i in (nrow(out)-1):2){
  if(i%in%temp){
    if(!((i-1) %in% temp)){
      pos <- c(pos, i)
    }
    if(!((i+1)%in%temp)){
      pos <- c(pos,i)
    }
  }
}
y_breaks <- vector()
for(i in 1:length(cluster_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cluster_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = 1:max(out$Var1), ordered = T)
my_breaks <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
my_labels <- c("1e-8", "1e-7", "1e-6", "1e-5", "1e-4", "1e-3")
out$value[out$value < 1e-7] <- 1e-7
out$value[out$value > 1e-3] <- 1e-3

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(10)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-7,1e-3)) +
  scale_y_discrete(breaks = floor(y_breaks), labels = y_labels) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_blank(), # bg of the panel
    plot.background = element_blank(), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "../Output/Clustering/Expression_scores_for_TNF_related_genes.pdf", width = 25, height = 12)
heatmap
dev.off()

png(file = "../Output/Clustering/Expression_scores_for_TNF_related_genes.png", width = 35, height = 17, units = "cm", res = 300)
heatmap
dev.off()

