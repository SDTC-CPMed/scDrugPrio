

#####################################################################################################################
# Clustering of CD (Martin et al, cell) lamina propria ileal cells
# BY MARTIN SMELIK
#####################################################################################################################

#setwd("/data/samsc76/SIGMA_13012020/CD_batch_corrected/drug_prediction_R")
dir.create("../Input/CD")
dir.create("../Input/CD/CD_all")
dir.create("../Input/CD/Individual_patients")
dir.create("../Output/CD/")

library(dplyr)
library(Seurat)
#install_version("R.filesets", version = "2.12.1")
library(R.filesets)
set.seed(54)

#*INSET MARTIN CODE*

#####################################################################################################################
# DEGs and cell typing of CD (Martin et al, cell) lamina propria ileal cells
# BY SAMUEL SCHAEFER
#####################################################################################################################
library(doParallel)
library(matrixStats)
#fp <- "/data/sharedData/scPred_validation_data_sets_SamueL_drug_prediction_project/CD/"
fp <- getwd()
setwd(paste(fp, "/../Input/CD",sep=""))


cd_all <- loadRDS('CD_all_integrated.rds') # counts = DCA adjusted data

# Martins old clustering parameters
#cd_all <- FindNeighbors(cd_all, k.param = 15)
#cd_all <- FindClusters(cd_all, resolution = 0.25)

# New parameters
cd_all <- FindNeighbors(cd_all, k.param = 15)
cd_all <- FindClusters(cd_all, resolution = 0.8)

saveRDS(cd_all, file = "CD_batch_corrected_clustered.rds")

# Create cellID / cluster ID matrix
id <- cbind(as.character(names(Idents(cd_all))), as.character(Idents(cd_all)))
colnames(id) <- c("Cell_ID", "Cluster_ID")
write.table(id, file = "CD_all/Cell_ID_to_cluster_ID_k=15_res=0.8.txt", sep="\t", col.names = T, row.names = F)

#######################################
# Load data
#######################################
TSNEPlot(cd_all, label = T, label.size = 3)
TSNEPlot(cd_all, cells = grep("Involved", Cells(cd_all), value =T)) # only cells from sick
TSNEPlot(cd_all, cells = grep("Uninvolved", Cells(cd_all), value =T)) # only cells from healthy volunteers
dev.off()

temp_id <- Idents(cd_all)
for(i in 1:length(unique(temp_id))){
  print(paste("Investigating Cluster ", unique(temp_id)[i], sep=""))
  temp <- temp_id[temp_id == unique(temp_id)[i]]
  print(paste("n of cells ", length(temp), " n of healthy cells ", sum(grepl(pattern = "Involved", x = names(temp))), 
              " n of sick cells ", sum(grepl(pattern = "Uninvolved", x = names(temp))), sep=""))
}

#######################################
# DEG calculation - logfc.threshold = log(1.5)
#######################################
temp_id <- as.numeric(id[,2])
names(temp_id) <- id[,1]
temp_id[grepl("Uninvolved", id[,1])] <- temp_id[grepl("Uninvolved", id[,1])] + 0.1 # mark cells from univolved areas  
temp_id <- as.factor(temp_id)
cd_all@active.ident <- temp_id
# unique clusters
temp_id_unique <- sort(as.numeric(unique(as.character(temp_id))))

# Exclude clusters with only healthy/sick cells from calculation
for(i in 0:(length(unique(id[,2]))-1)){
  if(any(!(i%in%temp_id_unique),!((i+0.1)%in%temp_id_unique))){
    print(paste("Exclude cluster ", i, " from DEG calculation", sep=""))
    temp_id_unique <- temp_id_unique[!temp_id_unique%in% c(i,i+0.1)]
  }
}
# Exclude clusters with less than 3 healthy or sick cells from calculation
# => excludes no cluster
temp <- table(Idents(cd_all))
if(any(temp<3)){
  print(paste("Exclude clusters:", names(temp[temp<3]), sep=" "))
  temp <- as.numeric(names(temp[temp<3]))
  temp_id_unique <- temp_id_unique[!temp_id_unique %in% c(temp, temp+0.1, temp-0.1)]
}

# Output directory
dir.create(path="DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells")

# Running MAST DEGs analysis
for(i in 1:(length(temp_id_unique)/2)){
  print("Checking cluster:")
  pos <- i*2-1
  print(temp_id_unique[pos:(pos+1)])
  if(!temp_id_unique[pos]%in%c(1,2,6,11,13,15,30,32)){ # exclude clusters as they have no genes with logFC.threshold > log(1.5)
    temp <- FindMarkers(cd_all,
                        slot = "counts",
                        test.use = "MAST", 
                        ident.1 = temp_id_unique[pos], # important that this is a factor vector with cell labels as names
                        ident.2 = temp_id_unique[pos+1],
                        only.pos = F,
                        random.seed = 3,
                        pseudocount.use = 0,
                        logfc.threshold = log(1.5),
                        min.pct = 0.1,
                        min.cells.group = 3)
  }
  
  if(exists("temp")){
    print("Saving ...")
    write.table(temp, file = paste("DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/Cluster_", temp_id_unique[pos], "_res=0.8_dims=32_k=15.txt",sep=""),sep="\t",col.names=NA,row.names=T)
    rm(temp)
  }
}

#######################################
# Running MAST DEGs summary
#######################################
files <- list.files(path = "DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/")
files <- files[grepl(files, pattern = "_res=0.8_dims=32_k=15")]
degs.clusters.h.vs.s <- matrix(NA, nrow = 15000, ncol = length(files))
for(i in 1:ncol(degs.clusters.h.vs.s)){
  temp <- read.table(file = paste("DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/", files[i], sep=""), sep="\t", header = T)
  temp <- temp[as.numeric(temp$p_val_adj) < 0.05,]
  temp <- as.character(temp[order(as.numeric(temp$p_val_adj), decreasing = F),1]) # order by significance
  print(length(temp))
  if(length(temp)>0){
    degs.clusters.h.vs.s[1:length(temp),i] <- temp
  }
  rm(temp)
}
colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = "_res=0.8_dims=32_k=15.txt"))
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,colSums(!is.na(degs.clusters.h.vs.s))>0]
write.table(degs.clusters.h.vs.s, file = "DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters.txt",sep="\t", col.names = T, row.names = F)

#######################################
# Translate MAST DEGs for further analysis - translate to Entrez Gene IDs
#######################################
transl <- as.matrix(read.table(file = "../HGNC translation matrix 201108/transl.txt", sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,6]),] # Entrez Gene ID
transl <- transl[!is.na(transl[,7]),] # Ensemble Gene ID
print(head(transl))

out <- matrix(NA, ncol = ncol(degs.clusters.h.vs.s), nrow = nrow(degs.clusters.h.vs.s))
colnames(out) <- colnames(degs.clusters.h.vs.s)
for(i in 1:ncol(degs.clusters.h.vs.s)){
  print(paste("Working on", colnames(degs.clusters.h.vs.s)[i], sep = " "))
  temp <- degs.clusters.h.vs.s[,i]
  temp <- temp[!is.na(temp)]
  print(paste("n of DEGs: ", length(temp), sep=""))
  temp <- transl[match(temp, transl[,7]),6] # keeps the order
  temp <- unique(temp[!is.na(temp)])
  print(paste("translated to n of DEGs: ", length(temp), sep=""))
  out[1:length(temp),i] <- temp
}
out <- out[rowSums(!is.na(out))>0,colSums(!is.na(out))>0]
write.table(out, file = "DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt",sep="\t", col.names = T, row.names = F)

#######################################################################################################################################################
# Cell typing
#######################################################################################################################################################

cd_all <- readRDS(file = "CD_batch_corrected_clustered.rds")

# Adjust for sequencing depth
X <- cd_all@assays$RNA@counts
X <- as.matrix(X)
temp <- colSums(X)
X <- t(t(X)/temp)
cd_all@assays$RNA@counts <- X

#######################################
# CREATE HEATMAPS FOR CELL TYPE MARKERS FROM MARTIN et al - mean expression per cluster for major cell type markers
#######################################
# Load translation matrix
transl <- as.matrix(read.table(file = "../HGNC translation matrix 201108/transl.txt", sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,2]),] # Official approved HGNC symbol
transl <- transl[!is.na(transl[,7]),] # Ensemble Gene ID
# load markers
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,1]),1] # n = 29
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 
# create cell ID matrix
id <- Idents(cd_all)
id <- cbind(names(id),as.numeric(as.character(id)))
unique_cl <- unique(id[,2])
unique_cl <- unique_cl[order(as.numeric(unique_cl), decreasing = F)]
# extract count-adjusted expression values
X <- cd_all@assays$RNA@counts
X <- X[rownames(X)%in% transl[,7],]
X <- as.matrix(X)
# calculate mean expression per cluster
out <- foreach(cl = c(1:length(unique_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==unique_cl[cl],1]]
  exp <- exp[match(markers[,2], rownames(exp)),]
  return(rowMeans(exp))
}
colnames(out) <- markers[,1]
rownames(out) <- unique_cl
out <- out[,colSums(!is.na(out))>0]

# Create heatmap
library(viridis)
library(gplots)
library(tidyverse)
mode(out) <- "numeric" 

out <- log10(out)
max_exp <- ceiling(max(out))
min_exp <- floor(min(out))

color.palette <- c(colorRampPalette(viridis(20)[1:7])(n = 10), colorRampPalette(viridis(20)[8:20])(n = 10))
col_breaks <- seq(from=min_exp, to = max_exp, by = 0.2)

print(length(col_breaks))
print(length(color.palette))

pdf(file = "Marker_genes/Expression_scores_for_major_cell_type_markers.pdf", paper = "a4") # save as horisontal A3
heatmap.2(out, trace = "none", 
          Colv = F,
          dendrogram = "row",
          col = color.palette, 
          symbreaks = F,
          breaks = col_breaks,
          key = T, 
          density.info = "none",
          keysize = 1,
          symkey = F,
          na.rm = F,
          margins = c(6,4),
          cexCol = 0.8,
          colsep = NA,
          rowsep = NA,
          sepwidth = c(0.001,0.001), 
          sepcolor = NULL) # try to remove white lines between blocks - not possible since heatmap.2 is reliant on image()
dev.off()

#######################################
# ALL CELL TYPE MARKERS FROM ORIGINAL ARTICLE - cell specific expression + dendrogram of all cell type markers
#######################################

# load markers again
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,2]),2] # n = 29 / 187
markers <- c("TRAC", markers)
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 
# calculate mean expression per cluster
out <- foreach(cl = c(1:length(unique_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==unique_cl[cl],1]]
  exp <- exp[match(markers[,2], rownames(exp)),]
  return(rowMeans(exp))
}
colnames(out) <- markers[,1]
rownames(out) <- unique_cl
out <- out[,colSums(!is.na(out))>0]
# Create heatmap for all markers
library(ggplot2)
library(reshape2)
library(ggdendro)
mode(out) <- "numeric" 

# dendrogram
out.dendro <- as.dendrogram(hclust(d = dist(x = log10(out))))
dendro.plot <- ggdendrogram(data = out.dendro, rotate = T)
pdf(file = "Marker_genes/Dendro_for_expression_scores_for_all_cell_type_markers.pdf") # save as horisontal A3
print(dendro.plot)
dev.off()
# get order of clusters
out.order <- order.dendrogram(out.dendro)

# heatmap
min_exp <- floor(min(out)*10^6)/10^6
out <- melt(out)
out[,1] <- as.character(out[,1])
out$Var1 <- factor(x = as.character(out$Var1), levels = unique_cl[out.order], ordered = T)

my_breaks <- c(min_exp, min_exp*10, min_exp*10^2, min_exp*10^3, min_exp*10^4)
my_labels <- c(as.character(min_exp), "6e-05", "6e-04", "6e-03", "6e-02")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = T) + scale_fill_gradientn(colours = c("black",viridis(10)), name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))#, axis.text.y=as.character(), axis.ticks.y=element_blank())  

pdf(file = "Marker_genes/Mean_expression_scores_of_clusters_for_all_cell_type_markers.pdf", width = 23, height = 10) # save as horisontal A3
heatmap
dev.off()

#######################################
# CREATE HEATMAP INCLUDING EVERY CELL TYPE - COUNT ADJUSTED DCA DENOISED EXPRESSION SCORES
#######################################
library(textshape)
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,2]),2] # n = 187
markers <- c("TRAC", markers)
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# 34-2 = T, 33 = ILC, 23-32 = plasma, 16-40 = B, 18-9 = MNP + pDC, 36-37 = stroma / glial
cluster_order <- c(33,30,27,26,7,20,21,5,0,4,3,2,11,1,32,22,9,12,37,28,16,23,10,29,31,15,14,34,13,6,39,17,8,35,25,18,38,19,24,36)

# collect all cells for cluster
out <- foreach(cl = c(1:length(cluster_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cluster_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- cluster_matrix(exp, dim = "row")
  return(rbind(exp, matrix(NA, nrow = 75, ncol = ncol(exp))))
}
colnames(out) <- markers[,1]
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

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(10)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  scale_y_discrete(breaks = floor(y_breaks), labels = y_labels) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Expression_scores_for_all_cell_type_markers.pdf", width = 25, height = 12)
heatmap
dev.off()

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F, show.legend = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(10)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  scale_y_discrete(breaks = floor(y_breaks), labels = y_labels) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
  )

png(file = "Marker_genes/Expression_scores_for_all_cell_type_markers.png", width = 52, height = 24, units = "cm", res = 600)
heatmap
dev.off()


#######################################
# INPUT MAJOR CELL TYPES FOR CLUSTERS TO DO CELL TYPE SPECIFIC HEATMAPS
#######################################
cellt_clusters <- c(33,30,27,26,7,20,21,5,0,4,3,2,11,1,32,22,9,12,37,28,16,23,10,29,31,15,14,34,13,6,39,17,8,35,25,18,38,19,24,36)
names(cellt_clusters) <- c(rep("T",times = 14), "ILC", rep("Plasma", times = 10), rep("B", times = 6), "MNP", "MNP", rep("Stroma_glia",times=7)) 

#######################################
# CREATE HEATMAP FOR MNP (mononuclear phagocytes) specific cell type markers (Figure 2A in original article) - Z every single cell + Z mean
#######################################
library(RColorBrewer)
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,3]),3] # n = 41
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "MNP"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colSums(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
cl_order <- rownames(out)
out <- cluster_matrix(out, dim = "row", method = "ward.D")

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
min_exp <- floor(min(out[!is.na(out)])) # -0.93
max_exp <- floor(max(out[!is.na(out)])) # 0.56

out[out < -0.5] <- -0.5
out[out>0.5] <- 0.5

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_col_breaks <- seq(from = -0.5, to = 0.5, by = 0.1)
my_col_labels <- as.character(seq(from = -0.5, to = 0.5, by = 0.1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_col_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_col_labels) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_of_mean_expression_scores_for_MNP_cell_type_markers.pdf", width = 7, height = 2)
heatmap
dev.off()

# collect all cells for clusters
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "mcquitty")
  exp <- rbind(exp, matrix(NA, nrow = 25, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
min_exp <- floor(min(out[!is.na(out)])) # -6
max_exp <- floor(max(out[!is.na(out)])) # 3

out[out< -2] <- -2
out[out>2] <- 2

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = order_rownames, ordered = T)
my_col_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_y_discrete(breaks = floor(y_breaks), labels = as.character(y_labels)) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_expression_scores_for_MNP_cell_type_markers.pdf", width = 10, height = 6)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR T cells specific cell type markers (Figure 2F in orginal article) - Z every single cell + Z mean
#######################################
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,7]),7] # n = 17
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "T"  | names(cellt_clusters) == "ILC"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colSums(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
cl_order <- rownames(out)
out <- cluster_matrix(out, dim = "row", method = "ward.D")

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -7
print(max(out[!is.na(out)])) # 4

out[out < -2] <- -2
out[out>2] <- 2

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_breaks, labels = my_labels) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_of_mean_expression_scores_for_T_cell_type_markers.pdf", width = 5, height = 4)
heatmap
dev.off()

# collect expression for all cells
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "mcquitty")
  exp <- rbind(exp, matrix(NA, nrow = 100, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -7
print(max(out[!is.na(out)])) # 4

out[out< -2] <- -2
out[out>2] <- 2

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = order_rownames, ordered = T)
my_col_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_y_discrete(breaks = floor(y_breaks), labels = as.character(y_labels)) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_expression_scores_for_T_cell_type_markers.pdf", width = 5, height = 4)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR T cells specific cell type markers - Z every single cell + Z mean
#######################################
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,12]),12] # n = 85
markers <- c("TRAC", markers)
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "T"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colSums(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
out <- cluster_matrix(out, dim = "row", method = "ward.D")
cl_order <- rownames(out)

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -4.7
print(max(out[!is.na(out)])) # 1.7

out[out < -1.5] <- -1.5
out[out>1.5] <- 1.5

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- seq(from = -1.5, to = 1.5, by = 0.5)
my_labels <- as.character(my_breaks)

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_breaks, labels = my_labels) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_mean_expression_scores_for_T_cell_sub_type_markers.pdf", width = 15, height = 4)
heatmap
dev.off()

# collect expression for all cells
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "mcquitty")
  exp <- rbind(exp, matrix(NA, nrow = 100, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -7.7
print(max(out[!is.na(out)])) # 5.7

out[out< -2] <- -2
out[out>2] <- 2

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = order_rownames, ordered = T)
my_col_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_y_discrete(breaks = floor(y_breaks), labels = as.character(y_labels)) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_expression_scores_for_T_cell_sub_type_markers.pdf", width = 15, height = 10)
heatmap
dev.off()

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F, show.legend = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_y_discrete(breaks = floor(y_breaks), labels = as.character(y_labels)) +
  theme(
    #axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

png(file = "Marker_genes/Z_expression_scores_for_T_cell_sub_type_markers.png", width = 35, height = 20, units = "cm", res = 600)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR T cells specific cell type markers (Figure S3E in orginal article) - COUNT ADJUSTED DCA DENOISED EXPRESSION SCORES
#######################################
# collect all cells for cluster
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "mcquitty")
  exp <- rbind(exp, matrix(NA, nrow = 100, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(floor(min(out[!is.na(out)])*10^6)/10^6) # 4e-6
print(floor(max(out[!is.na(out)])*100)/100) # 0.05

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2,1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Expression_scores_for_all_T_cell_sub_type_markers.pdf", width = 15, height = 10)
heatmap
dev.off()

# collect mean expression for every cluster
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  return(colMeans(exp))
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- cl_order

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # 1.2e-5
print(max(out[!is.na(out)])) # 4e-3

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2,1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Mean_expression_scores_for_all_T_cell_sub_type_markers.pdf", width = 15, height = 4)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR MNP cells specific cell type markers 2 (Figure 2D in original article) - every single cell 
#######################################
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,6]),6] # n = 29
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "MNP"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colSums(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
cl_order <- rownames(out)
out <- cluster_matrix(out, dim = "row", method = "ward.D")

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -1.5
print(max(out[!is.na(out)])) # 0.7

out[out < -0.5] <- -0.5
out[out>0.5] <- 0.5

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- seq(from = -0.5, to = 0.5, by = 0.1)
my_labels <- as.character(seq(from = -0.5, to = 0.5, by = 0.1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "transparent", name = "count", breaks = my_breaks, labels = my_labels) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_of_mean_expression_scores_for_MNP_cell_type_markers_2.pdf", width = 7, height = 2)
heatmap
dev.off()

# collect expression for all cells
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "mcquitty")
  exp <- rbind(exp, matrix(NA, nrow = 25, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -7
print(max(out[!is.na(out)])) # 4

out[out< -2] <- -2
out[out>2] <- 2

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = order_rownames, ordered = T)
my_col_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_y_discrete(breaks = floor(y_breaks), labels = as.character(y_labels)) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_expression_scores_for_MNP_cell_type_markers_2.pdf", width = 5, height = 4)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR STROMA & GLIA cells specific cell type markers (Figure 2H in original article) - every single cell 
#######################################
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,9]),9] # n = 24
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "Stroma_glia"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colSums(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
cl_order <- rownames(out)
out <- cluster_matrix(out, dim = "row", method = "ward.D")

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -5.5
print(max(out[!is.na(out)])) # 2.1

out[out < -2] <- -2
out[out>2] <- 2

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_breaks, labels = my_labels) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_of_mean_expression_scores_for_STROMA_GLIA_cell_type_markers.pdf", width = 5, height = 3)
heatmap
dev.off()

# collect expression for all cells
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "mcquitty")
  exp <- rbind(exp, matrix(NA, nrow = 50, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# expression per mean
means <- colMeans(out[!is.na(out[,1]),])
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -7.4
print(max(out[!is.na(out)])) # 4.6

out[out< -2] <- -2
out[out>2] <- 2

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = order_rownames, ordered = T)
my_col_breaks <- seq(from = -2, to = 2, by = 1)
my_labels <- as.character(seq(from = -2, to = 2, by = 1))

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #scale_fill_gradientn(colours = viridis(20,option = "B"), na.value = "transparent", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), na.value = "black", name = "count", breaks = my_col_breaks, labels = my_labels) +
  scale_y_discrete(breaks = floor(y_breaks), labels = as.character(y_labels)) +
  theme(
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Z_expression_scores_for_STROMA_GLIA_cell_type_markers.pdf", width = 5, height = 4)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR B cell specific cell type markers (Figure S3A in original article) - COUNT ADJUSTED DCA DENOISED EXPRESSION SCORES
#######################################
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,11]),11] # n = 71
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 
markers[1,2] <- "ENSG00000132465" # JCHAIN was not translated but crucial to cell typing

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "Plasma" | names(cellt_clusters) == "B"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colSums(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
out <- cluster_matrix(out, dim = "row", method = "mcquitty")
cl_order <- rownames(out)
cl_order <- cl_order[length(cl_order):1]

# collect all cells for cluster
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- cluster_matrix(exp, dim = "row", method = "ward.D")
  exp <- rbind(exp, matrix(NA, nrow = 50, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # 2.7e-6
print(max(out[!is.na(out)])) # 0.12

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2,1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Expression_scores_for_all_B_cell_sub_type_markers.pdf", width = 15, height = 10)
heatmap
dev.off()

# collect mean expression for every cluster
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  return(colMeans(exp))
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- cl_order

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # 1.2e-5
print(max(out[!is.na(out)])) # 4e-3

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2,1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Mean_expression_scores_for_all_B_cell_sub_type_markers.pdf", width = 15, height = 4)
heatmap
dev.off()

#######################################
# CREATE HEATMAP FOR STROMA_GLIA cell specific cell type markers (Figure 2G in original article) - COUNT ADJUSTED DCA DENOISED EXPRESSION SCORES
#######################################
markers <- read.table(file = "Marker_genes/Cell_type_markers_original_article.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,8]),8] # n = 82
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 

# collect mean expression for clusters
# considered_cl <- cellt_clusters[names(cellt_clusters) == "Stroma_glia"]
# out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
#   exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
#   exp <- t(exp[match(markers[,2], rownames(exp)),])
#   exp <- exp[,colSums(!is.na(exp))>0]
#   exp <- colSums(exp)
#   return(exp)
# }
# colnames(out) <- markers[match(colnames(out), markers[,2]),1]
# rownames(out) <- as.character(considered_cl)
# out <- cluster_matrix(out, dim = "row", method = "ward.D2")
cl_order <- c("36","26","19","39","20","37","25")
cl_order <- cl_order[length(cl_order):1]

# collect all cells for cluster
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  #exp <- cluster_matrix(exp, dim = "row", method = "ward.D2")
  exp <- rbind(exp, matrix(NA, nrow = 20, ncol = ncol(exp)))
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- 1:nrow(out)
order_rownames <- rownames(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # 3.7e-6
print(max(out[!is.na(out)])) # 0.12

# set breaks
temp <- rownames(out)[which(!is.na(out[,1]))] # finds valid cells position
pos <- vector()
pos[1] <- 1
for(i in 2:(nrow(out)-1)){
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
for(i in 1:length(cl_order)){
  y_breaks[i] <- mean(c(pos[(i*2-1)],pos[(i*2)]))
}
y_labels <- cl_order

# format for heatmap
out <- melt(out)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2,1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Expression_scores_for_all_STROMA_GLIA_cell_sub_type_markers.pdf", width = 15, height = 10)
heatmap
dev.off()

# collect mean expression for every cluster
out <- foreach(cl = c(1:length(cl_order)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==cl_order[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  return(colMeans(exp))
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- cl_order

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # 1.2e-5
print(max(out[!is.na(out)])) # 4e-3

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2,1e-1)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-6,1e-1)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Mean_expression_scores_for_all_STROMA_GLIA_cell_sub_type_markers.pdf", width = 15, height = 4)
heatmap
dev.off()




#######################################################################################################################################################
# INPUT FINAL CELL TYPING
#######################################################################################################################################################
cellt_clusters <- cbind(cellt_clusters, names(cellt_clusters), c(
  "Highly activated T cells","Tregs","Naive/CM T cells","Naive/CM T cells","Tissue-resident T cells","Highly activated T cells","Naive/CM T cells", #22
  "CD8+ Cytotoxic T cells","CD8+ Cytotoxic T cells","Tissue-resident CD4+ T cells","CD8+ Tissue-resident T cells", "CD8+ Tissue-resident T cells", #3
  "CD8+ Tissue-resident T cells", "Type 1/3 cytokine Tissue-resident T cells", "Innate lymphoid cells", # 33
  "Plasmablasts","IgA plasma cells","IgG plasma cells","IgM plasma cells","IgM plasma cells","IgA plasma cells","IgA plasma cells","IgA plasma cells", #11
  "IgA plasma cells","IgG plasma cells","CD27- Naive B cells", "CD27- Naive B cells", "Memory B cell", "Memory B cell", "Memory B cell", "CD27+ Memory B cell", #40
  "Dentritic cells", "Activated Macrophages", "ACKR1+ Endothelial cells","Fibroblasts","Fibroblasts","Activated Fibroblasts","Smooth muscle cells","Glial cells","Smooth muscle cells"))

write.table(cellt_clusters, file = "Marker_genes/Cell_typing_by_marker_expression.txt", sep="\t", col.names = c("Cluster", "Major_cell_type", "Cell_subtype"), row.names = F)

# Activated Fibroblasts
temp <- id[id[,2]==39,1]
temp <- t(matrix(unlist(strsplit(temp, split = "_")), nrow = 3))
print(table(temp[,1]))
# Activated T cells
temp <- id[id[,2] %in% c(21,34),1]
temp <- temp[grepl(x = temp, pattern = "Involved")]
temp <- t(matrix(unlist(strsplit(temp, split = "_")), nrow = 3))
print(table(temp[,1]))
# Activated Macrophages
temp <- id[id[,2]==9,1]
temp <- temp[grepl(x = temp, pattern = "Involved")]
temp <- t(matrix(unlist(strsplit(temp, split = "_")), nrow = 3))
print(table(temp[,1]))
# IgG plasma cells
temp <- id[id[,2] %in% c(32, 13),1]
temp <- temp[grepl(x = temp, pattern = "Involved")]
temp <- t(matrix(unlist(strsplit(temp, split = "_")), nrow = 3))
print(table(temp[,1]))

#######################################################################################################################################################
# Barplots comparing n of cells between inflammed and uninflammed
#######################################################################################################################################################
library(scales)
library(reshape2)
cellt_clusters <- c(33,30,27,26,7,20,21,5,0,4,3,2,11,1,32,22,9,12,37,28,16,23,10,29,31,15,14,34,13,6,39,17,8,35,25,18,38,19,24,36)
names(cellt_clusters) <- c(rep("T",times = 14), "ILC", rep("Plasma", times = 10), rep("B", times = 6), "MNP", "MNP", rep("Stroma_glia",times=7)) 
cellt_clusters <- cbind(cellt_clusters, names(cellt_clusters), c(
  "Highly activated T cells","Tregs","Naive/CM T cells","Naive/CM T cells","Tissue-resident CD4+ T cells","Highly activated T cells","Naive/CM T cells", #22
  "CD8+ Cytotoxic T cells","CD8+ Cytotoxic T cells","Tissue-resident CD4+ T cells","CD8+ Tissue-resident T cells", "CD8+ Tissue-resident T cells", #3
  "CD8+ Tissue-resident T cells", "Type 1/3 cytokine Tissue-resident T cells", "Innate lymphoid cells", # 33
  "Plasmablasts","IgA plasma cells","IgG plasma cells","IgM plasma cells","IgM plasma cells","IgA plasma cells","IgA plasma cells","IgA plasma cells", #11
  "IgA plasma cells","IgG plasma cells","CD27- Naive B cells", "CD27- Naive B cells", "Memory B cell", "Memory B cell", "Memory B cell", "CD27+ Memory B cell", #40
  "Dentritic cells", "Activated Macrophages", "ACKR1+ Endothelial cells","Fibroblasts","Fibroblasts","Activated Fibroblasts","Smooth muscle cells","Glial cells","Smooth muscle cells"))


print(nrow(cellt_clusters)) # = 40

col_palette2 <- c("#660000","#800000", "#FF0000", "#e7772b", "#f49728", "#FF7F50", # Reds
                  "#e8c124", "#f6ed31", "#F0E68C", # Yellows
                  
                  "#000000", "#808080", "#A9A9A9", "#D3D3D3", # Black & greys 4
                  
                  "#000080", "#2005A5", "#0000FF", "#2b58a3", "#1F65CC", "#3686D3", "#4197D9","#4AA4DE",  # Blues 8
                  "#4169E1", "#6495ED", "#1E90FF", "#00BFFF", "#87CEFA", "#ADD8E6", "#2f9da7", "#008080", # Blues 8
                  "#553c93","#82318f","#882865","#4B0082", # Purples
                  
                  "#006400", "#228B22", "#556B2F","#4F7942","#6B8E23","#9ACD32","#BFE890","#AFE1AF", # Greens 8
                  "#32CD32", "#2E8B57", "#3CB371", "#66CDAA", "#7CFC00", "#00FF00", "#8FBC8F", "#ADFF2F") # Greens 8

show_col(col_palette2)

# need 6 B, 10 Plasma, -> Green 16
# 2 MNP, 1 ILC,  -> 2 grey +1 Black
# 7 stroma, -> Red + yellow
# 14 T - > Blues + Purple 14

col_to_cluster <- c(
  "5", col_palette2[30], "CD8+ T cytotoxic 1",
  "0", col_palette2[31], "CD8+ T cytotoxic 2",
  "3", col_palette2[32], "CD8+ tissue resident 1",
  "2", col_palette2[33], "CD8+ tissue resident 2",
  "11", col_palette2[17], "CD8+ tissue resident 3",
  "33", col_palette2[14], "Highly activated T cells 1",
  "20", col_palette2[15], "Highly activated T cells 2",
  "30", col_palette2[16], "T regs",
  "27", col_palette2[19], "Naive/CM T cells 1",
  "26", col_palette2[21], "Naive/CM T cells 2",
  "21", col_palette2[22], "Naive/CM T cells 3",
  "4", col_palette2[24], "CD4+ Tissue resident memory 1",
  "7", col_palette2[25], "CD4+ Tissue resident memory 2",
  "1", col_palette2[26], "Type 1/3 cytokine Tissue-resident T cells",
  
  "15", col_palette2[40], "CD27- Naive B cells 1",
  "14", col_palette2[41], "CD27- Naive B cells 2",
  "12", col_palette2[34], "IgG plasma cells 1",
  "31", col_palette2[35], "IgG plasma cells 2",
  "37", col_palette2[46], "IgM plasma cells 1",
  "28", col_palette2[47], "IgM plasma cells 2",
  "9", col_palette2[36],  "IgA plasma cells 1",
  "16", col_palette2[37], "IgA plasma cells 2",
  "23", col_palette2[38], "IgA plasma cells 3",
  "10", col_palette2[39], "IgA plasma cells 4",
  "29", col_palette2[49], "IgA plasma cells 5",
  "39", col_palette2[42], "CD27+ Memory B cell",
  "34", col_palette2[43], "Memory B cell 1",
  "13", col_palette2[44], "Memory B cell 2",
  "6", col_palette2[45], "Memory B cell 3",
  "22", col_palette2[48], "Plasmablasts",
  
  "32", col_palette2[10], "Innate lymphoid cells",
  "17", col_palette2[12], "Dentritic cells",
  "8", col_palette2[11], "Activated Macrophages",
  
  "25", col_palette2[1], "Fibroblasts 1",
  "18", col_palette2[2], "Fibroblasts 2",
  "38", col_palette2[3], "Activated Fibroblasts",
  "19", col_palette2[6], "Smooth muscle cells 1",
  "36", col_palette2[5], "Smooth muscle cells 2",
  "24", col_palette2[7], "Glial cells",
  "35", col_palette2[9], "ACKR1+ Endothelial cells"
)                  
col_to_cluster <- t(matrix(t(col_to_cluster), nrow = 3))

cellt_clusters <- cbind(cellt_clusters, col_to_cluster[match(as.numeric(cellt_clusters[,1]), as.numeric((col_to_cluster[,1]))),2:3])

# ALL clusters
pdf(file = "../../Output/CD/Cluster_plot_own_color_palette_colorblind_friendly_all_cells.pdf", width = 12, height = 6)
TSNEPlot(cd_all, label.size = 0) + scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                                                      labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),5]))
dev.off()

# Sick cells only
cl <- Idents(cd_all)
cl <- cl[names(cl) %in% grep("Involved",names(cl),value = T)]
cl <- as.numeric(as.character(unique(cl)))
temp <- cellt_clusters[cellt_clusters[,1]%in%cl,]
pdf(file = "../../Output/CD/Cluster_plot_own_color_palette_colorblind_friendly_CD_inflammed_cells.pdf", width = 12, height = 6)
TSNEPlot(cd_all, label.size = 0, cells = grep("Involved", Cells(cd_all), value =T)) + 
  scale_color_manual(values = as.character(temp[order(as.numeric(temp[,1])),4]),
                     labels = as.character(temp[order(as.numeric(temp[,1])),5]))
dev.off()

# Healthy cells only
cl <- Idents(cd_all)
cl <- cl[names(cl) %in% grep("Involved",names(cl),value = T)]
cl <- as.numeric(as.character(unique(cl)))
temp <- cellt_clusters[cellt_clusters[,1]%in%cl,]
pdf(file = "../../Output/CD/Cluster_plot_own_color_palette_colorblind_friendly_CD_unaffected_cells.pdf", width = 12, height = 6)
TSNEPlot(cd_all, label.size = 0, cells = grep("Uninvolved", Cells(cd_all), value =T)) + 
  scale_color_manual(values = as.character(temp[order(as.numeric(temp[,1])),4]),
                     labels = as.character(temp[order(as.numeric(temp[,1])),5]))
dev.off()

write.table(cellt_clusters, file = "../../Output/CD/Cluster_colors.txt",sep="\t",col.names =  T, row.names = F)

# Show cell types as percentage of total cells from sample split by patients
###############################################################################################################################################
id <- Idents(cd_all)
id <- cbind(names(id), as.numeric(as.character(id)))
temp <- strsplit(id[,1], split = "_")
temp2 <- vector()
temp3 <- vector()
for(i in 1:nrow(id)){
  temp2[i] <- temp[[i]][2]
  temp3[i] <- temp[[i]][1]
}
id <- cbind(id[,1], temp2, temp3, id[,2])
colnames(id) <- c("Cell", "Sample", "Patient", "Cluster")
rm(temp, temp2, temp3, i)

# input pat_data
pat_data <- cbind(c("GSM3972009","GSM3972010","GSM3972011","GSM3972012",
                    "GSM3972013","GSM3972014","GSM3972016","GSM3972015",
                    "GSM3972017","GSM3972018","GSM3972020","GSM3972019",
                    "GSM3972022","GSM3972021","GSM3972024","GSM3972023",
                    "GSM3972026","GSM3972025","GSM3972028","GSM3972027",
                    "GSM3972030","GSM3972029"),
                  c("Inflammed","Unaffected","Inflammed","Unaffected","Inflammed","Unaffected","Inflammed",
                    "Unaffected","Inflammed","Unaffected","Inflammed","Unaffected","Inflammed","Unaffected",
                    "Inflammed","Unaffected","Inflammed","Unaffected","Inflammed","Unaffected",
                    "Inflammed","Unaffected"),
                  c(paste("Patient", c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11),sep=" ")))

# add pat ID
id <- cbind(id, pat_data[match(id[,3],pat_data[,1]),3])
# add cell type a cell belongs to
id <- cbind(id, cellt_clusters[match(id[,4], cellt_clusters[,1]),c(2,5,4)])
# new colnames
colnames(id) <- c("Cell_ID", "Disease_status", "Sample_ID", "Cluster", "Patient_ID", "Major_cell_type", "Cell_type", "clust_color")

# create matrix summarizing n cells per cell type for each sample
unique_sample <- unique(id[,3])
out <- matrix(0, nrow = length(unique_sample), ncol = length(unique(id[,7])))
temp <- unique(id[,7])
temp <- sort(temp)
colnames(out) <- temp

for(i in 1:length(unique_sample)){
  temp <- table(id[id[,3]==unique_sample[i],7])
  temp <- temp[sort(names(temp))]
  out[i,colnames(out)%in%as.character(names(temp))] <- temp
}
rownames(out) <- unique_sample
write.table(out, file = "Individual_patients/n_cells_of_individual_CD_samples_per_cell_type.txt", sep="\t", row.names = T, col.names = NA)

# Percentage of total cells
temp <- out
temp <- 100*temp/rowSums(temp) # expressed as % of total cells
temp <- melt(temp)
temp <- cbind(temp, pat_data[match(rownames(out), pat_data[,1]),2:3])
temp <- cbind(temp, col_to_cluster[match(as.character(temp[,2]), col_to_cluster[,3]),2])
colnames(temp)[(ncol(temp)-2):ncol(temp)] <- c("Disease", "Patient", "Color")
temp <- as.data.frame(temp)
temp$Var2 <- factor(temp$Var2, levels = col_to_cluster[,3], ordered = T)
temp$Patient <- factor(temp$Patient, levels = c(paste("Patient ", 1:11,sep="")), ordered = T)

pdf(file = "Percentage_of_total_cells_per_cluster_individual_CD_patients.pdf", width = 14, height = 6)

ggplot(temp, aes(x = Disease, y = value, fill = Var2)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  facet_grid(~ Patient, margins = F) + 
  scale_fill_manual("Cell types", values = col_to_cluster[,2]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, margin = margin(-13,0,0,0), colour = "black", family = "Helvetica"),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", family = "Helvetica", face = "bold"),
        axis.title = element_blank())

dev.off()
rm(temp)

# Combine cell types for sub cell type plots
#####################################################################################################################################################################################
id <- Idents(cd_all)
id <- cbind(names(id), as.numeric(as.character(id)))
temp <- strsplit(id[,1], split = "_")
temp2 <- vector()
temp3 <- vector()
for(i in 1:nrow(id)){
  temp2[i] <- temp[[i]][2]
  temp3[i] <- temp[[i]][1]
}
id <- cbind(id[,1], temp2, temp3, id[,2])
colnames(id) <- c("Cell", "Sample", "Patient", "Cluster")
rm(temp, temp2, temp3, i)

# input pat_data
pat_data <- cbind(c("GSM3972009","GSM3972010","GSM3972011","GSM3972012",
                    "GSM3972013","GSM3972014","GSM3972016","GSM3972015",
                    "GSM3972017","GSM3972018","GSM3972020","GSM3972019",
                    "GSM3972022","GSM3972021","GSM3972024","GSM3972023",
                    "GSM3972026","GSM3972025","GSM3972028","GSM3972027",
                    "GSM3972030","GSM3972029"),
                  c("Inflammed","Unaffected","Inflammed","Unaffected","Inflammed","Unaffected","Inflammed",
                    "Unaffected","Inflammed","Unaffected","Inflammed","Unaffected","Inflammed","Unaffected",
                    "Inflammed","Unaffected","Inflammed","Unaffected","Inflammed","Unaffected",
                    "Inflammed","Unaffected"),
                  c(paste("Patient", c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11),sep=" ")))

# add pat ID
id <- cbind(id, pat_data[match(id[,3],pat_data[,1]),3])
# add cell type a cell belongs to
id <- cbind(id, cellt_clusters[match(id[,4], cellt_clusters[,1]),c(2,3,4)])
# new colnames
colnames(id) <- c("Cell_ID", "Disease_status", "Sample_ID", "Cluster", "Patient_ID", "Major_cell_type", "Cell_type", "clust_color")

# create matrix summarizing n cells per cell type for each sample
unique_sample <- unique(id[,3])
out <- matrix(0, nrow = length(unique_sample), ncol = length(unique(id[,7])))
temp <- unique(id[,7])
temp <- sort(temp)
colnames(out) <- temp

for(i in 1:length(unique_sample)){
  temp <- table(id[id[,3]==unique_sample[i],7])
  temp <- temp[sort(names(temp))]
  out[i,colnames(out)%in%as.character(names(temp))] <- temp
}
rownames(out) <- unique_sample
write.table(out, file = "Individual_patients/n_cells_of_individual_CD_samples_per_combined_cell_type.txt", sep="\t", row.names = T, col.names = NA)

# function for bar plots
bar_plots_cell_type_proportions_for_individual_patients <- function(in_matrix, name, label_order, pat_order, width, height){
  mode(in_matrix) <- "numeric"
  out <- as.data.frame(in_matrix)
  out$Patient <- rownames(out)
  out <- melt(out, id.vars = "Patient")
  out$variable <- factor(x = out$variable, levels = label_order, ordered = T)
  out$Patient <- factor(x = out$Patient, levels = pat_order, ordered = T)
  
  pdf(file = paste("Individual_patients/Cell_type_proportions_for_individual_patients__",name,".pdf",sep=""), width = width, height = height)
  temp <- ggplot(out, aes(x = Patient, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'stack') + scale_fill_viridis_d() +
    theme(
      axis.text.x = element_text(angle=90, hjust=1, color = "black", face = "bold"),
      axis.text.y = element_text(hjust=1.5, color = "black", face = "bold"),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
      axis.line.y = element_line(colour = "black")
    )
  plot(temp)
  dev.off()
} 

#cellt_clusters[,3] <- cellt_clusters[,5]

# prepare matrix for plotting stacked bar graphs for plasma cells
temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="Plasma",3])]
temp <- temp/rowSums(temp) # get percentages
# seperate by inflammed/uninflammed
temp_infl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Inflammed",],]
rownames(temp_infl) <- pat_data[match(rownames(temp_infl), pat_data[,1]),3]
temp_uninfl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Unaffected",],]
rownames(temp_uninfl) <- pat_data[match(rownames(temp_uninfl), pat_data[,1]),3]
rm(temp)
#
# GENERAL PATIENT ORDER - based on similarity
#
pat_order <- paste("Patient ", c(1,4,3,7,2,6,9,10,5,8,11), sep="")
# plot for inflammed plasma cells
label_order <- c("IgG plasma cells", "IgM plasma cells", "IgA plasma cells", "Plasmablasts")
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_infl,
                                                        label_order = label_order,
                                                        name = "plasma_cells_infl",
                                                        pat_order = pat_order,
                                                        width = 5, height = 4)
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_uninfl,
                                                        label_order = label_order,
                                                        name = "plasma_cells_uninfl",
                                                        pat_order = pat_order,
                                                        width = 5, height = 4)

# prepare matrix for plotting stacked bar graphs for T cells
temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="T",3])]
temp <- temp/rowSums(temp) # get percentages
# seperate by inflammed/uninflammed
temp_infl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Inflammed",],]
rownames(temp_infl) <- pat_data[match(rownames(temp_infl), pat_data[,1]),3]
temp_uninfl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Unaffected",],]
rownames(temp_uninfl) <- pat_data[match(rownames(temp_uninfl), pat_data[,1]),3]
rm(temp)
# plot for inflammed plasma cells
label_order <- c("Type 1/3 cytokine Tissue-resident T cells", "Tissue-resident CD4+ T cells", "CD8+ Tissue-resident T cells", "Tregs",
                 "Naive/CM T cells","CD8+ Cytotoxic T cells", "Highly activated T cells")
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_infl,
                                                        label_order = label_order,
                                                        name = "T_cells_infl",
                                                        pat_order = pat_order,
                                                        width = 7, height = 4)
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_uninfl,
                                                        label_order = label_order,
                                                        name = "T_cells_uninfl",
                                                        pat_order = pat_order,
                                                        width = 7, height = 4)

# prepare matrix for plotting stacked bar graphs for Stroma/Glia cells
temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="Stroma_glia",3])]
temp <- temp/rowSums(temp) # get percentages
# seperate by inflammed/uninflammed
temp_infl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Inflammed",],]
rownames(temp_infl) <- pat_data[match(rownames(temp_infl), pat_data[,1]),3]
temp_uninfl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Unaffected",],]
rownames(temp_uninfl) <- pat_data[match(rownames(temp_uninfl), pat_data[,1]),3]
rm(temp)
# plot for inflammed plasma cells
label_order <- c("Smooth muscle cells", "Glial cells", "Fibroblasts","ACKR1+ Endothelial cells", "Activated Fibroblasts")
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_infl,
                                                        label_order = label_order,
                                                        name = "Stroma_Glia_cells_infl",
                                                        pat_order = pat_order,
                                                        width = 5, height = 4)
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_uninfl,
                                                        label_order = label_order,
                                                        name = "Stroma_Glia_cells_uninfl",
                                                        pat_order = pat_order,
                                                        width = 5, height = 4)

# prepare matrix for plotting stacked bar graphs for myeloid cells
temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="MNP",3])]
temp <- temp/rowSums(temp) # get percentages
# seperate by inflammed/uninflammed
temp_infl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Inflammed",],]
rownames(temp_infl) <- pat_data[match(rownames(temp_infl), pat_data[,1]),3]
temp_uninfl <- temp[rownames(temp)%in% pat_data[pat_data[,2]=="Unaffected",],]
rownames(temp_uninfl) <- pat_data[match(rownames(temp_uninfl), pat_data[,1]),3]
rm(temp)
# plot for inflammed plasma cells
label_order <- c("Dentritic cells", "Activated Macrophages")
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_infl,
                                                        label_order = label_order,
                                                        name = "myeloid_cells_infl",
                                                        pat_order = pat_order,
                                                        width = 5, height = 4)
bar_plots_cell_type_proportions_for_individual_patients(in_matrix = temp_uninfl,
                                                        label_order = label_order,
                                                        name = "myeloid_cells_uninfl",
                                                        pat_order = pat_order,
                                                        width = 5, height = 4)


#######################################################################################################################################################
# Get geometric mean for GIMATs
#######################################################################################################################################################
geo_mean <- vector()
geo_cols <- c("IgG plasma cells", "Highly activated T cells", "ACKR1+ Endothelial cells","Activated Fibroblasts", "Activated Macrophages")
# Plasma cells
temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="Plasma",3])]
temp <- temp/rowSums(temp) # get percentages
geo_mean <- cbind(geo_mean, temp[,colnames(temp) %in% geo_cols])

temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="T",3])]
temp <- temp/rowSums(temp) # get percentages
geo_mean <- cbind(geo_mean, temp[,colnames(temp) %in% geo_cols])

temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="Stroma_glia",3])]
temp <- temp/rowSums(temp) # get percentages
geo_mean <- cbind(geo_mean, temp[,colnames(temp) %in% geo_cols])

temp <- out[,colnames(out) %in% unique(cellt_clusters[cellt_clusters[,2]=="MNP",3])]
temp <- temp/rowSums(temp) # get percentages
geo_mean <- cbind(geo_mean, temp[,colnames(temp) %in% geo_cols])

temp <- vector()
gm_mean = function(a){prod(a)^(1/length(a))}
for(i in 1:nrow(geo_mean)){
  temp[i] <- gm_mean(as.vector(geo_mean[i,]))
}
temp2 <- rowMeans(geo_mean)
geo_mean <- cbind(geo_mean, temp, temp2)
colnames(geo_mean) <- c(geo_cols, "geo_mean", "arithmetric_mean")

geo_mean_infl <- geo_mean[rownames(geo_mean)%in% pat_data[pat_data[,2]=="Inflammed",],6:7]
rownames(geo_mean_infl) <- pat_data[match(rownames(geo_mean_infl), pat_data[,1]),3]
geo_mean_infl <- geo_mean_infl[match(pat_order, rownames(geo_mean_infl)),]
geo_mean_uninfl <- geo_mean[rownames(geo_mean)%in% pat_data[pat_data[,2]=="Unaffected",],6:7]
rownames(geo_mean_uninfl) <- pat_data[match(rownames(geo_mean_uninfl), pat_data[,1]),3]
geo_mean_uninfl <- geo_mean_uninfl[match(pat_order, rownames(geo_mean_uninfl)),]

pdf(file = "Individual_patients/GIMATs_geo_metric_means_for_patients_inflamed.pdf", width = 4, height = 3)
barplot(geo_mean_infl[,1], xlab = "", las = 2, main = "", sub = "", ylab = "Geometric mean (frequencies)", ylim = c(0,0.4))
dev.off()

pdf(file = "Individual_patients/GIMATs_arithmetric_metric_means_for_patients_inflamed.pdf", width = 4, height = 3)
barplot(geo_mean_infl[,2], xlab = "", las = 2, main = "", sub = "", ylab = "Mean (frequencies)", ylim = c(0,0.5))
dev.off()

pdf(file = "Individual_patients/GIMATs_geo_metric_means_for_patients_uninflamed.pdf", width = 4, height = 3)
barplot(geo_mean_uninfl[,1], xlab = "", las = 2, main = "", sub = "", ylab = "Geometric mean (frequencies)", ylim = c(0,0.4))
dev.off()

pdf(file = "Individual_patients/GIMATs_arithmetric_metric_means_for_patients_uninflamed.pdf", width = 4, height = 3)
barplot(geo_mean_uninfl[,2], xlab = "", las = 2, main = "", sub = "", ylab = "Mean (frequencies)", ylim = c(0,0.5))
dev.off()

