
# Run in R 4.0.4 due to cluster limitations.

#setwd("/local/data1/user/samsc76/Doctis_anti_IL17_v2_artificial_target/R")
#install.packages("Seurat")
library(Seurat)
library(ggplot2)
library(doParallel)
registerDoParallel(cores = 15)
set.seed(21)

# get DCA output
###############################################################################

psa <- readRDS("../Input/PsA_data/PsA_all_integrated_IL17_singlet_patients.rds")
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
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45.pdf")
TSNEPlot(psa, group.by = "RNA_snn_res.0.45")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_healthy.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "healthy", 1], group.by = "RNA_snn_res.0.45")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_sick.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "sick", 1], group.by = "RNA_snn_res.0.45")
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

# k = 15 and resolution = 0.45 in the above loaded RDS
psa <- FindNeighbors(psa, k.param = 15, nn.method = "annoy")
psa <- FindClusters(psa, resolution = 0.45)
out <- get_n_cells()
out
write.table(out, file = "../Output/Clustering/Clustering_k=15_res=0.45.txt", sep="\t", col.names = NA, row.names = T)
id <- Idents(psa)
id <- cbind(cell_id = names(id), cluster_id = as.character(id))
write.table(id, file = "Cell_identity_k=15_res=0.45.txt",sep="\t", col.names = T, row.names = F)

# remake plots
'pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45.pdf")
TSNEPlot(psa, group.by = "RNA_snn_res.0.45")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_healthy.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "healthy", 1], group.by = "RNA_snn_res.0.45")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_sick.pdf")
TSNEPlot(psa, cells = meta[meta[,3] %in% "sick", 1], group.by = "RNA_snn_res.0.45")
dev.off()'

pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_UMAP.pdf")
UMAPPlot(psa, group.by = "RNA_snn_res.0.45")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_healthy_UMAP.pdf")
UMAPPlot(psa, cells = meta[meta[,3] %in% "healthy", 1], group.by = "RNA_snn_res.0.45")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.45_sick_UMAP.pdf")
UMAPPlot(psa, cells = meta[meta[,3] %in% "sick", 1], group.by = "RNA_snn_res.0.45")
dev.off()

#saveRDS(object=psa, file="../Input/PsA_data/PsA_all_clustered.rds")
psa <- readRDS(file = "../Input/PsA_data/PsA_all_clustered.rds")

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
unique_cl <- unique(as.numeric(as.character(psa@active.ident)))
cluster_order <- c(7,15,0,1,14,12,13,3,5,8,18,10,6,16,9,11,20,2,4,19,17)

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
print(floor(max(out[!is.na(out)])*100)/100) # 0.1

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

pdf(file = "../Output/Clustering/Expression_scores_for_all_CD_cell_type_markers_smaller.pdf", width = 12, height = 6.5)
heatmap
dev.off()

jpeg(filename = "../Output/Clustering/Expression_scores_for_all_CD_cell_type_markers.jpeg", width = 30, height = 15, units = "cm", res = 300, pointsize = 7)
plot(heatmap)
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

thresh <-2

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




#######################################################################################################################################################
# INPUT FINAL CELL TYPING
#######################################################################################################################################################

cluster_order <- c(7,15,0,1,14,12, # naive CD4
                   13,3,5, # activated CD4+
                   8, # naive CD8+
                   18,10,6, # activated CD8+
                   16,9, # B
                   11,20,2,4, # granulocytes
                   19, # Monocyte
                   17) #DC

cellt_clusters <- cluster_order
names(cellt_clusters) <- c(rep("T",times = 13),  rep("B", times = 2), rep("Neutrophil", times = 4), "Monocyte", "DC")

cellt_clusters <- cbind(cellt_clusters, names(cellt_clusters), 
                        c(paste("Naive CD4+ T ", 1:6, sep=""), 
                          paste("Activated CD4+ T ", 1:3, sep=""),
                          "Naive CD8+ T",
                          paste("Activated CD8+ T ", 1:3, sep=""),
                          paste("B cell ", 1:2, sep=""),
                          paste("Neutrophil ", 1:4, sep=""),
                          "Monocyte",
                          "Dentritic cell")) 

write.table(cellt_clusters, file = "../Output/Clustering/Cell_typing_by_marker_expression.txt", sep="\t", col.names = c("Cluster", "Major_cell_type", "Cell_subtype"), row.names = F)

#######################################################################################################################################################
# Barplots comparing n of cells between inflammed and uninflammed
#######################################################################################################################################################
library(scales)
library(reshape2)
print(nrow(cellt_clusters)) # = 21

col_palette2 <- c("#660000","#800000", "#FF0000", "#e7772b", "#f49728", "#FF7F50", # Reds
                  "#e8c124", "#f6ed31", "#F0E68C", # Yellows
                  
                  "#000000", "#808080", "#A9A9A9", "#D3D3D3", # Black & greys 4
                  
                  "#000080", "#2005A5", "#0000FF", "#2b58a3", "#1F65CC", "#3686D3", "#4197D9","#4AA4DE",  # Blues 8
                  "#4169E1", "#6495ED", "#1E90FF", "#00BFFF", "#87CEFA", "#ADD8E6", "#2f9da7", "#008080", # Blues 8
                  "#553c93","#82318f","#882865","#4B0082", # Purples
                  
                  "#006400", "#228B22", "#556B2F","#4F7942","#6B8E23","#9ACD32","#BFE890","#AFE1AF", # Greens 8
                  "#32CD32", "#2E8B57", "#3CB371", "#66CDAA", "#7CFC00", "#00FF00", "#8FBC8F", "#ADFF2F", "#9EF38F") # Greens 8

show_col(col_palette2)

col_to_cluster <- c(
  # naive CD4
  7, col_palette2[26],
  15, col_palette2[25],
  0, col_palette2[24],
  1, col_palette2[22],
  14, col_palette2[20],
  12, col_palette2[18],    
  # activated CD4+
  13, col_palette2[16],
  3, col_palette2[15],
  5, col_palette2[14],  
  # naive CD8+
  8, col_palette2[33],    
  # activated CD8+
  18, col_palette2[30],
  10, col_palette2[31],
  6, col_palette2[32],
  # B
  16, col_palette2[43],
  9, col_palette2[42],
  # granulocytes
  11, col_palette2[8],
  20, col_palette2[7],
  2, col_palette2[5],
  4, col_palette2[4],
  # Monocyte
  19, col_palette2[3],    
  #DC
  17, col_palette2[2])   
col_to_cluster <- t(matrix(t(col_to_cluster), nrow = 2))
show_col(col_to_cluster[,2])

cellt_clusters <- cbind(cellt_clusters, col_to_cluster[match(as.numeric(cellt_clusters[,1]), as.numeric((col_to_cluster[,1]))),2])
write.table(cellt_clusters, file = "../Output/Clustering/Cluster_colors.txt",sep="\t",col.names =  T, row.names = F)

# ALL clusters
psa <- RunUMAP(psa, dims = 1:ncol(psa@reductions$pca@cell.embeddings))

pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_all_cells.pdf", width = 12, height = 6)
TSNEPlot(psa, label.size = 0) + scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                                                      labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_all_cells_UMAP.pdf", width = 12, height = 6)
UMAPPlot(psa, label.size = 0) + scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                                                      labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

# Sick cells only
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_PsA_sick_cells.pdf", width = 12, height = 6)
TSNEPlot(psa, label.size = 0, cells = meta[meta[,3] %in% "sick", 1]) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_PsA_sick_cells_UMAP.pdf", width = 12, height = 6)
UMAPPlot(psa, label.size = 0, cells = meta[meta[,3] %in% "sick", 1]) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

# Healthy cells only
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_PsA_healthy_cells.pdf", width = 12, height = 6)
TSNEPlot(psa, label.size = 0, cells = meta[meta[,3] %in% "healthy", 1]) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_PsA_healthy_cells_UMAP.pdf", width = 12, height = 6)
UMAPPlot(psa, label.size = 0, cells = meta[meta[,3] %in% "healthy", 1]) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()


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
id[id[,1] %in% meta[meta[,4] %in% "NR", 1],3] <- as.numeric(id[id[,1] %in% meta[meta[,4] %in% "NR", 1],3]) - 0.5 # NR = X.2
#id[id[,1] %in% meta[meta[,4] %in% "R", 1],3] <- as.numeric(id[id[,1] %in% meta[meta[,4] %in% "R", 1],3]) + 0.5 #R = X.6 

# set new ID in psa
temp_id <- as.factor(id[,3])
names(temp_id) <- id[,1]
psa@active.ident <- temp_id

# test VlnPlots
#VlnPlot(psa, features = "TNF", slot = "counts", pt.size = 0)
#VlnPlot(psa, features = "NR3C1", slot = "counts", pt.size = 0)
#dev.off()

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
      if(unique_cl[i] != 17){ # exclude clusters as they have no genes with logFC.threshold > log(1.5)
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


# get artificial targets
il17_genes <- rownames(psa@assays$RNA@counts)
il17_genes <- il17_genes[grepl(il17_genes, pattern = "IL17")]
write.table(il17_genes, file = "../Input/Drugs/articicial_IL17_targets.txt",sep="\t",col.names = F, row.names = F)

# Check DCA denoised expression of il17 genes
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
X <- X[rownames(X)%in% il17_genes,]
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
  exp <- X[rownames(X)%in%il17_genes, colnames(X)%in%id[id[,4]==cluster_order[cl],1]]
  exp <- t(exp[match(il17_genes, rownames(exp)),])
  exp <- cluster_matrix(exp, dim = "row")
  return(rbind(exp, matrix(NA, nrow = 75, ncol = ncol(exp))))
}
colnames(out) <- il17_genes
out <- out[,colSums(!is.na(out))>0]
rownames(out) <- nrow(out):1

# Create heatmap for all markers
mode(out) <- "numeric" 
print(floor(min(out[!is.na(out)])*10^6)/10^6) # 2e-6
print(floor(max(out[!is.na(out)])*10^4)/10^4) # 2e-4

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
my_breaks <- c(1e-6, 1e-5, 1e-4)
my_labels <- c("1e-6", "1e-5", "1e-4")
out$value[out$value < 1e-6] <- 1e-6
out$value[out$value > 1e-4] <- 1e-4

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(10)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-4)) +
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

pdf(file = "../Output/Clustering/Expression_scores_for_IL17_related_genes.pdf", width = 25, height = 12)
heatmap
dev.off()

png(file = "../Output/Clustering/Expression_scores_for_IL17_related_genes.png", width = 35, height = 17, units = "cm", res = 300)
heatmap
dev.off()

