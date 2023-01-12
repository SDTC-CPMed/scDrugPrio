
# Run in R 4.0.4 due to cluster limitations.

library(dplyr)
library(Seurat)
library(R.filesets)
library(doParallel)
library(scales)

set.seed(54)
fp <- getwd()
dir.create("../Output/Clustering")
dir.create("../Output/DCA_MAST_DEGs")


#####################################################################################################################
# Load data
#####################################################################################################################
dir.data <- '../Input/DCA_adjusted_matrix/'
denoised_DCA <- read.table(paste(dir.data, "mean.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(denoised_DCA) <- denoised_DCA[,1]
denoised_DCA <- denoised_DCA[,-1]
colnames(denoised_DCA) <- denoised_DCA[1,]
denoised_DCA <- denoised_DCA[-1,]

lat <- read.table(paste(dir.data, "latent.tsv", sep = ""), sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)

######################################################################################################################
# Pre-processing
######################################################################################################################
joint_RA <- CreateSeuratObject(counts = denoised_DCA, project = "Drug_prediction_denoised", min.cells = 1, min.features = 1)

#######################################################################################################################
# Dimensional reduction
# no linear PCA used as DCA has a built in non-linear PCA => latent features => PCs
# instead of PCs = latent features from DCA are entered into the Seurat object
######################################################################################################################
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
joint_RA@reductions$pca <- pca

# Cluster cells
joint_RA <- FindNeighbors(joint_RA, dims = 1:nrow(lat), k = 15)
joint_RA <- FindClusters(joint_RA, resolution = 0.6)

# Non-linear dimensional reduction
joint_RA <- RunTSNE(joint_RA, dims = 1:nrow(lat))
joint_RA <- RunUMAP(joint_RA, dims = 1:nrow(lat))

# make plots
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.6.pdf")
TSNEPlot(joint_RA)
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.6_Healthy.pdf")
TSNEPlot(joint_RA, cells = grep("Healthy", Cells(joint_RA), value =T), group.by = "RNA_snn_res.0.8")
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.6_Sick.pdf")
TSNEPlot(joint_RA, cells = grep("Sick", Cells(joint_RA), value =T))
dev.off()

pdf(file = "../Output/Clustering/Clustering_k=15_res=0.6_UMAP.pdf")
UMAPPlot(joint_RA)
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.6_Healthy_UMAP.pdf")
UMAPPlot(joint_RA, cells = grep("Healthy", Cells(joint_RA), value =T))
dev.off()
pdf(file = "../Output/Clustering/Clustering_k=15_res=0.6_Sick_UMAP.pdf")
UMAPPlot(joint_RA, cells = grep("Sick", Cells(joint_RA), value =T))
dev.off()

# Save clustering outcomes
saveRDS(joint_RA, file = "../Input/AIA/AIA_clustered_k=15_res=0.6.rds", version = 2)
#joint_RA <- loadRDS(file = "../Input/AIA/AIA_clustered_k=15_res=0.6.rds")

id <- Idents(joint_RA)
id <- cbind(as.character(names(id)), as.character(id))
colnames(id) <- c("Cell_ID", "Cluster_ID")
write.table(id, file = "../Input/AIA/Cell_identity.txt", sep="\t", col.names = T, row.names = F)

temp <- unlist(strsplit(id[,1], split = "Joint_"))
temp <- strsplit(temp, split = "_")
temp <- temp[seq(from = 2, to = length(temp), by = 2)]
temp <- foreach(i = c(1:length(temp)), .combine = "c") %do% {
  return(paste(temp[[i]][1:3], collapse = "_"))
}
temp <- cbind(id, temp)
mice <- unique(temp[,3])
cl <- unique(temp[,2])
temp <- foreach(i = c(1:length(cl)), .combine = "rbind")%do%{
  temp2 <- temp[temp[,2] %in% cl[i],3]
  temp2 <- table(temp2)
  out <- rep(0, times = length(mice))
  names(out) <- mice
  out[names(out)%in% names(temp2)] <- temp2[order(names(out)[names(out) %in% names(temp2)], names(temp2))]
  return(out)
}
rownames(temp) <- paste("Cluster_",cl,sep="")
write.table(temp, file = "../Input/AIA/Cells_by_sample.txt", sep="\t", col.names = NA, row.names = T)
rm(mice, cl, out, temp)


######################################################################################################################
# Cluster information
######################################################################################################################

id <- Idents(joint_RA)
clust.info <- matrix(NA, nrow=length(unique(id)), ncol = 5)
colnames(clust.info) <- c("Cluster", "n_Sick_cells", "n_Healthy_cells", "Sick/Healthy", "total_n_cells")
for(i in 0:(length(unique(id))-1)){
  temp <- id[id==i]
  clust.info[i+1,1] <- i
  clust.info[i+1,2] <- sum(grepl("Sick", names(temp)))
  clust.info[i+1,3] <- sum(grepl("Healthy", names(temp)))
  clust.info[i+1,4] <- clust.info[i+1,2]/clust.info[i+1,3]
  clust.info[i+1,5] <- clust.info[i+1,2] + clust.info[i+1,3]
}
write.table(clust.info, file = "../Input/AIA/Frequency_healthy_vs_sick_clusters.txt", sep="\t", col.names = T, row.names = F)
rm(temp, clust.info, i, id)

######################################################################################################################
# DEG calculation - logfc.threshold = log(1.5)
######################################################################################################################
id <- as.matrix(read.table(file = "../Input/AIA/Cell_identity.txt", sep="\t", header = T))
temp_id <- as.numeric(id[,2])
names(temp_id) <- id[,1]
temp_id[grepl("Healthy", id[,1])] <- temp_id[grepl("Healthy", id[,1])] + 0.1 # mark cells from control joints (uninflamed)  
temp_id <- as.factor(temp_id)
joint_RA@active.ident <- temp_id
# unique clusters
temp_id_unique <- sort(as.numeric(unique(as.character(temp_id))))

# Exclude clusters with only healthy/sick cells from calculation
# => excludes no clusters
for(i in 0:(length(unique(id[,2]))-1)){
  if(any(!(i%in%temp_id_unique),!((i+0.1)%in%temp_id_unique))){
    print(paste("Exclude cluster ", i, " from DEG calculation", sep=""))
    temp_id_unique <- temp_id_unique[!temp_id_unique%in% c(i,i+0.1)]
  }
}
# Exclude clusters with less than 3 healthy or sick cells from calculation
# => excludes no cluster
temp <- table(Idents(joint_RA))
if(any(temp<3)){
  print(paste("Exclude clusters:", names(temp[temp<3]), sep=" "))
  temp <- as.numeric(names(temp[temp<3]))
  temp_id_unique <- temp_id_unique[!temp_id_unique %in% c(temp, temp+0.1, temp-0.1)]
}

# DEG function
calculate_DEGs <- function(temp_id_unique, outdir, joint_RA, pseudo.count = 0, logfc = log(1.5)){
  unique_cl <- unique(floor(temp_id_unique))
  
  # Running MAST DEGs analysis
  for(i in 1:length(unique_cl)){
    print("Checking cluster:")
    print(unique_cl[i])
    if((unique_cl[i] + 0.1) %in% temp_id_unique){
      temp <- FindMarkers(joint_RA,
                          slot = "counts",
                          test.use = "MAST", 
                          ident.1 = unique_cl[i], # important that this is a factor vector with cell labels as names
                          ident.2 = (unique_cl[i]+0.1),
                          only.pos = F,
                          random.seed = 3,
                          pseudocount.use = pseudo.count,
                          logfc.threshold = logfc,
                          min.pct = 0.1,
                          min.cells.group = 3)
      
      if(exists("temp")){
        print("Saving ...")
        write.table(temp, file = paste(outdir, "/Cluster_", unique_cl[i], "_res=0.6_dims=32_k=15.txt",sep=""),sep="\t",col.names=NA,row.names=T)
        rm(temp)
      }
    }
  }
  
  #######################################
  # Running MAST DEGs summary
  #######################################
  files <- list.files(path = outdir)
  files <- files[grepl(files, pattern = "_res=0.6_dims=32_k=15")]
  degs.clusters.h.vs.s <- matrix(NA, nrow = 30000, ncol = length(files))
  for(i in 1:ncol(degs.clusters.h.vs.s)){
    temp <- read.table(file = paste(outdir, "/", files[i], sep=""), sep="\t", header = T)
    temp <- temp[temp$p_val_adj < 0.05,]
    temp <- as.character(temp[order(as.numeric(temp$p_val_adj), decreasing = F),1]) # order by significance
    if(length(temp)>0){
      degs.clusters.h.vs.s[1:length(temp),i] <- temp
    }
    rm(temp)
  }
  colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = "_res=0.6_dims=32_k=15.txt"))
  degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,colSums(!is.na(degs.clusters.h.vs.s))>0]
  write.table(degs.clusters.h.vs.s, file = paste(outdir, "/SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters.txt", sep=""),sep="\t", col.names = T, row.names = F)
  return(degs.clusters.h.vs.s)
}

# run on DCA denoised, log-tranformed data - excludes clusters as they have no genes with logFC.threshold > log(1.5)
degs.clusters.h.vs.s <- calculate_DEGs(temp_id_unique, outdir = "../Output/DCA_MAST_DEGs", 
                                       joint_RA, pseudo.count = 0, logfc = log(1.5))

# Translate DEGs to human Entrez ID
for(i in 1:ncol(degs.clusters.h.vs.s)){
  temp <- degs.clusters.h.vs.s[,i]
  degs.clusters.h.vs.s[,i] <- NA
  temp <- transl[match(temp, transl[,3]),2]
  temp <- temp[!is.na(temp)]
  
  if(length(temp)>0){
    degs.clusters.h.vs.s[1:length(temp),i] <- temp
  }
}
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,colSums(!is.na(degs.clusters.h.vs.s))>0]
write.table(degs.clusters.h.vs.s, file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_ID.txt")


# Cell typing
###############################################################################

library(matrixStats)
library(reshape2)
library(ggplot2)
library(viridis)
library(textshape)

transl <- read.table(file ="../Input/Human-mouse_homologs/transl.txt", sep="\t", header = T)

# extract count-adjusted expression values of marker_genes
X <- joint_RA@assays$RNA@counts
X <- as.matrix(X)
temp <- colSums(X)
X <- t(t(X)/temp)
X <- X[rownames(X)%in% marker_genes,]
X <- as.matrix(X)

# load marker genes
marker_genes <- read.table(file = "../Input/Marker_genes/Markers_for_cell_typing_heatmap_RA.txt",sep="\t", header = T, stringsAsFactors = F)
marker_genes <- marker_genes[order(marker_genes[,3]),]
marker_genes <- marker_genes[!is.na(marker_genes[,3]),1]

# cluster order
id <- Idents(joint_RA)
id <- cbind(cell_ID = names(id), cluster_ID = as.numeric(as.character(id)))
id[,2] <- floor(as.numeric(as.character(id[,2])))

temp_id <- as.factor(id[,2])
names(temp_id) <- id[,1]
joint_RA@active.ident <- temp_id

unique_cl <- unique(as.numeric(as.character(joint_RA@active.ident)))
#cluster_order <- unique_cl
cluster_order <- c(18,16,8, # Erythrocytes
                   3,13,7,10, # Stem cells + progenitors
                   9,19,1,20, # B cells
                   0,15,5,11,12, # T cells
                   4,2,6,14,17 # Myeloid clusters 1 to 5
)

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
print(floor(min(out[!is.na(out)])*10^7)/10^7) # 5e-7
print(floor(max(out[!is.na(out)])*100)/100) # 0.41

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
my_breaks <- c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
my_labels <- c("1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")
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

pdf(file = "../Output/Clustering/Expression_scores_for_all_RA_cell_type_markers.pdf", width = 25, height = 12)
heatmap
dev.off()

png(file = "../Output/Clustering/Expression_scores_for_all_RA_cell_type_markers.png", width = 35, height = 17, units = "cm", res = 300)
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

thresh <-1

my_breaks <- seq(from = -thresh, to = thresh, by = 0.5)
my_labels <- seq(from = -thresh, to = thresh, by = 0.5)

out[out[,3] < -thresh,3] <- -thresh
out[out[,3] > thresh,3] <- thresh

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) +
  geom_raster(na.rm = T) + scale_fill_gradientn(colours = viridis(12), name = "count", breaks = my_breaks, labels = my_labels, limits = c(-thresh,thresh)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) #, axis.text.y=as.character(), axis.ticks.y=element_blank())

pdf(file = "../Output/Clustering/Mean_expression_scores_for_all_RA_cell_type_markers.pdf", width = 25, height = 12)
heatmap
dev.off()

png(file = "../Output/Clustering/Mean_expression_scores_for_all_RA_cell_type_markers.png", width = 35, height = 12, units = "cm", res = 300)
heatmap
dev.off()


# Colorblind friendly colors
###############################################################################

cluster_order <- c(18,16,8, # Erythrocytes
                   3,13,7,10, # Stem cells + progenitors
                   9,19,1,20, # B cells
                   0,15,5,11,12, # T cells
                   4,2,6,14,17 # Myeloid clusters 1 to 5
)
names(cluster_order) <- c(rep("Erythrocyte",times = 3),rep("Stem cells",times=4), rep("B", times = 4), rep("T", times = 5), rep("Myeloid cluster",times=5)) 
cluster_order <- cbind(cluster_order, names(cluster_order), c("Erythrocyte, stadium 1", "Erythrocyte, stadium 2", "Erythrocyte, stadium 3",
                                                                 "Hematopoetic stem cell","Multipotent progenitor","Lymphoid progenitor","Myeloid progenitor",
                                                                 "B cell, immature","B cell, early activated", "B cell, activated","Plasma cell",
                                                                 "Cd4+ T cell, naive", "Cd4+ Il17f+ T cell", "Cd4+ Il2+ T cell", "Cd8+, inactive", "Cd8+, active",
                                                                 "Myeloid cluster, stadium 1","Myeloid cluster, stadium 2", "Myeloid cluster, stadium 3",
                                                                 "Myeloid cluster, stadium 4","Myeloid cluster, stadium 5"))
col_palette <- c("#660000","#e92925", "#e7772b", "#f49728","#e8c124",
                 "#f6ed31","#AFE1AF","#BFE890","#4F7942","#2f9da7",
                 "#2b58a3","#553c93","#82318f","#882865","#4B0082", 
                 "#2005A5","#000080", "#000000", "#A9A9A9", "#D3D3D3",
                 "#2b58a3", "#3686D3","#4AA4DE","#87CEFA","#2F9DA7",
                 "#009999","#6666FF","#4545E2","#000080","#553c93",
                 "#82318f","#882865","#4B0082","#330066", "#2E8B57")
show_col(col_palette)

col_palette <- col_palette[c(18,19,20, # 3 erythrocyte
              30,26,27,1, # 4 stem
              7,8,35,9, # 4 B cells
              23,21,17, # 3 CD4+ T cells
              13,15, # 2 CD8+ T cells
              2:6)]
cluster_order <- cbind(cluster_order, col_palette)


#col_palette2 <- c("#556B2F","#006400","#7CFC00","#00FF00","#4F7942","#6B8E23","#9ACD32","#ADFF2F","#228B22","#BFE890","#AFE1AF","#2E8B57","#3CB371","#66CDAA","#32CD32")
#show_col(col_palette2)

temp_id <- floor(as.numeric(as.character(Idents(joint_RA))))
temp_id <- as.factor(temp_id)
names(temp_id) <- names(Idents(joint_RA))
joint_RA@active.ident <- temp_id

# ALL clusters
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly.pdf", width = 8, height = 4)
TSNEPlot(joint_RA, label.size = 0) +  scale_color_manual(values = as.character(cluster_order[order(as.numeric(cluster_order[,1])),4]), labels = as.character(cluster_order[order(as.numeric(cluster_order[,1])),3]))
dev.off()
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_UMAP.pdf", width = 8, height = 4)
UMAPPlot(joint_RA, label.size = 0) + scale_color_manual(values = as.character(cluster_order[order(as.numeric(cluster_order[,1])),4]), labels = as.character(cluster_order[order(as.numeric(cluster_order[,1])),3]))
dev.off()

# Sick cells only
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_sick.pdf", width = 8, height = 4)
TSNEPlot(joint_RA, label.size = 0, cells = grep("Sick", Cells(joint_RA), value =T)) +  
  scale_color_manual(values = as.character(cluster_order[order(as.numeric(cluster_order[,1])),4]), labels = as.character(cluster_order[order(as.numeric(cluster_order[,1])),3]))
dev.off()
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_sick_UMAP.pdf", width = 8, height = 4)
UMAPPlot(joint_RA, label.size = 0, cells = grep("Sick", Cells(joint_RA), value =T)) + 
  scale_color_manual(values = as.character(cluster_order[order(as.numeric(cluster_order[,1])),4]), labels = as.character(cluster_order[order(as.numeric(cluster_order[,1])),3]))
dev.off()

# Healthy cells only
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_healthy.pdf", width = 8, height = 4)
TSNEPlot(joint_RA, label.size = 0, cells = grep("Healthy", Cells(joint_RA), value =T)) +  
  scale_color_manual(values = as.character(cluster_order[order(as.numeric(cluster_order[,1])),4]), labels = as.character(cluster_order[order(as.numeric(cluster_order[,1])),3]))
dev.off()
pdf(file = "../Output/Clustering/Cluster_plot_own_color_palette_colorblind_friendly_healthy_UMAP.pdf", width = 8, height = 4)
UMAPPlot(joint_RA, label.size = 0, cells = grep("Healthy", Cells(joint_RA), value =T)) + 
  scale_color_manual(values = as.character(cluster_order[order(as.numeric(cluster_order[,1])),4]), labels = as.character(cluster_order[order(as.numeric(cluster_order[,1])),3]))
dev.off()

write.table(cluster_order, file = "../Output/Clustering//Cluster_colors.txt",sep="\t",col.names =  T, row.names = F)

#######################################################################################################################################################
# Overview figure % of cells per cluster for individual mice 
#######################################################################################################################################################

# Create cellID / cluster ID matrix
id <- cbind(as.character(names(Idents(joint_RA))), as.character(Idents(joint_RA)))
colnames(id) <- c("Cell_ID", "Cluster_ID")
#id <- as.matrix(read.table(file = "CD_all/Cell_ID_to_cluster_ID_k=15_res=0.8.txt", sep="\t", header = T, stringsAsFactors = F))
# split by patients
temp <- strsplit(id[,1], split = "_")
temp2 <- vector()
temp3 <- vector()
for(i in 1:nrow(id)){
  temp2[i] <- temp[[i]][2]
  temp3[i] <- temp[[i]][4]
}
id <- cbind(id[,1], temp2, temp3, id[,2])
colnames(id) <- c("Cell", "Sample", "Patient", "Cluster")
rm(temp, temp2, temp3, i)
id[,3] <- paste(id[,2], "_mouse_", id[,3], sep ="")

# Make table showing n cells per patient and cluster
unique_pat <- unique(id[,3])
out <- matrix(0, nrow = length(unique_pat), ncol = length(unique(id[,4])))
temp <- unique(id[,4])
temp <- sort(as.numeric(temp))
colnames(out) <- temp

for(i in 1:length(unique_pat)){
  temp <- table(id[id[,3]==unique_pat[i],4])
  temp <- temp[order(as.numeric(names(temp)))]
  out[i,colnames(out)%in%as.character(as.numeric(names(temp)))] <- temp
}
colnames(out) <- paste("Cluster_",colnames(out),sep="")
rownames(out) <- unique_pat
write.table(out, file = "../Output/Clustering/AIA_n_cells_of_individual_samples_per_cluster.txt", sep="\t", row.names = T, col.names = NA)

# in percent
out <- 100*out/rowSums(out)
write.table(out, file = "../Output/Clustering/AIA_%_of_cells_of_individual_samples_per_cluster.txt", sep="\t", row.names = T, col.names = NA)

out <- cbind(rownames(out),out) 
colnames(out)[1] <- "Patient"
rownames(out) <- NULL
out <- as.data.frame(out)
out$Patient <- factor(x = out$Patient, levels = unique_pat, ordered = T)
for(i in 2:ncol(out)){
  out[,i] <- as.numeric(as.character(out[,i]))
}
temp <- unlist(strsplit(as.character(out$Patient), split = "_"))
out$disease <- temp[seq(from = 1, to = length(temp), by = 3)]
out <- melt(out, id.vars = c("Patient", "disease"))

# order by cluster ID
out$variable <- factor(x = out$variable, levels = paste("Cluster_",0:(length(unique(Idents(joint_RA)))-1),sep=""), ordered = T)

# order by major cell types
out$variable <- factor(x = out$variable, levels = paste("Cluster_",cluster_order[,1],sep=""), ordered = T)

pdf(file = "../Output/Clustering/AIA_percentage_of_total_cells_per_cluster_individual_mice_ordered_by_major_cell_type.pdf", width = 6, height = 6)
ggplot(out, aes(x = Patient, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  theme(axis.text.x = element_text(angle = 90, colour = "black"),
        panel.background = element_blank()) + #+ facet_grid(~ disease)
  scale_fill_manual(values = as.character(cluster_order[,4]))
dev.off()


