
#####################################################################################################################
# Clustering of MS (Schafflick et al, Nat Com) leukocytes #GSE138266
# BY SAMUEL SCHAEFER
#####################################################################################################################

#setwd("/data/samsc76/SIGMA_13012020/MS_batch_corrected/drug_prediction_R")
library(matrixStats)
library(reshape2)
library(ggplot2)


fp <- getwd()
setwd(paste(fp, "/../Input/MS",sep=""))

#######################################
# Load data
#######################################

ms_all <- readRDS(file = "MS_batch_corrected.rds") # counts = batch corrected DCA adjusted data, PCs = latent features

# New parameters
ms_all <- FindNeighbors(ms_all, k.param = 15)
ms_all <- FindClusters(ms_all, resolution = 0.35)
print(length(unique(Idents(ms_all))))

# Create cellID / cluster ID matrix
id <- cbind(as.character(names(Idents(ms_all))), as.character(Idents(ms_all)))
colnames(id) <- c("Cell_ID", "Cluster_ID")
write.table(id, file = "Cell_ID_to_cluster_ID_k=15_res=0.35.txt", sep="\t", col.names = T, row.names = F)

# split by patients
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

# Make table showing n cells per sample and cluster
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
write.table(out, file = "MS_n_cells_of_individual_samples_per_cluster.txt", sep="\t", row.names = T, col.names = NA)

# in percent
out <- 100*out/rowSums(out)
write.table(out, file = "MS_%_of_cells_of_individual_samples_per_cluster.txt", sep="\t", row.names = T, col.names = NA)

# make stacked bar graph
out <- cbind(rownames(out),out) 
colnames(out)[1] <- "GEO_ID"
# rename values in column GEO_ID based on "MS" or "IIH" sample
temp <- cbind(c("GSM4104122","GSM4104123","GSM4104124","GSM4104125",
                "GSM4104126","GSM4104127","GSM4104128","GSM4104129",
                "GSM4104130","GSM4104131","GSM4104132","GSM4104133"),
              c(rep("MS",times = 6), rep("IIH",times = 6)))
out[,1] <- temp[match(temp[,1], out[,1]),2]
out <- cbind(out, NA)
out[,ncol(out)] <- c(paste("Sample", c(1,2,3,4,5,6),sep=" "), paste("Sample",c(1,2,3,4,5,6),sep=" "))
colnames(out)[ncol(out)] <- "Patient"
rownames(out) <- NULL
out <- as.data.frame(out)
out$Patient <- factor(x = out$Patient, levels = paste("Sample", c(1,2,3,4,5,6),sep=" "), ordered = T)
for(i in 2:18){
  out[,i] <- as.numeric(as.character(out[,i]))
}

# P value for cellular heterogenity between sick MS patients
temp <- out
temp <- temp[temp[,1] == "MS",-1]
rownames(temp) <- temp[,ncol(temp)]
temp <- as.matrix(temp[,-ncol(temp)])
mode(temp) <- "numeric"
chisq <- chisq.test(temp)
print(chisq)

out <- melt(out, id.vars = c("Patient", "GEO_ID"))

pdf(file = "Percentage_of_total_cells_per_cluster_individual_MS_patients.pdf", width = 6, height = 5)

ggplot(out, aes(x = Patient, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ GEO_ID, margins = F) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, margin = margin(-13,0,0,0), colour = "black", family = "Helvetica"),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", family = "Helvetica", face = "bold"),
        axis.title = element_blank())

dev.off()

# Reload ID
id <- cbind(as.character(names(Idents(ms_all))), as.character(Idents(ms_all)))
colnames(id) <- c("Cell_ID", "Cluster_ID")

#######################################
# Clusters
#######################################
dir.create(path = "../../Output/Clustering")

pdf(file = "../../Output/Clustering/Clustering_k=15_res=0.35.pdf")
TSNEPlot(ms_all)
dev.off()
pdf(file = "../../Output/Clustering/Clustering_k=15_res=0.35_healthy.pdf")
TSNEPlot(ms_all, cells = grep("_IIH_", Cells(ms_all), value =T)) # only cells from control
dev.off()
pdf(file = "../../Output/Clustering/Clustering_k=15_res=0.35_sick.pdf")
TSNEPlot(ms_all, cells = grep("_MS_", Cells(ms_all), value =T)) # only cells from sick
dev.off()

temp_id <- Idents(ms_all)
for(i in 1:length(unique(temp_id))){
  temp <- temp_id[temp_id == unique(temp_id)[i]]
  print(paste("Investigating Cluster ", unique(temp_id)[i], "       n of cells ", length(temp), " n of healthy cells ", sum(grepl(pattern = "_IIH_", x = names(temp))), 
              " n of sick cells ", sum(grepl(pattern = "_MS_", x = names(temp))), sep=""))
}

#######################################################################################################################################################
# CELL TYPE MARKERS - cell specific expression + dendrogram of all cell type markers
#######################################################################################################################################################
library(ggplot2)
library(reshape2)
library(viridis)
library(textshape)

# Adjust for sequencing depth
X <- ms_all@assays$RNA@counts
X <- as.matrix(X)
temp <- colSums(X)
X <- t(t(X)/temp)

# Check marker genes
dir.create(path = "Marker_genes", showWarnings = F)

unique_cl <- unique(id[,2])
unique_cl <- unique_cl[order(as.numeric(unique_cl), decreasing = F)]

# Load translation matrix
transl <- as.matrix(read.table(file = "../HGNC translation matrix 201108/transl.txt", sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,2]),] # Official approved HGNC symbol
transl <- transl[!is.na(transl[,7]),] # Ensemble Gene ID

# load markers again
markers <- read.table(file = "Marker_genes/Marker genes collection.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,1]),1] # n = 93
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID
markers <- markers[!is.na(markers[,2]),] # n = 92

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
mode(out) <- "numeric" 

# set cluster order to form groups
out.order <- match(as.character(rev(c(10,5,3,1,6,15,13,7,2,0,9,12,4,14,19,18,17,8,11,16,20))),unique_cl)

# heatmap
print(max(out)) # 0.002
min_exp <- floor(min(out)*10^7)/10^7 
print(min_exp) # 5e-7
out <- melt(out)
out[,1] <- as.character(out[,1])
out$Var1 <- factor(x = as.character(out$Var1), levels = unique_cl[out.order], ordered = T)

my_breaks <- c(1e-6,1e-5,1e-4,1e-3)#, min_exp*10^4)#, min_exp*10^5)
my_labels <- c("1e-06","1e-05", "1e-04", "1e-03")#, "3e-02")#, "3e-01")
out$value[out$value>1e-3] <- 1e-3
out$value[out$value<1e-5] <- 1e-5

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = T) + scale_fill_gradientn(colours = c("black",viridis(4)), 
                                                name = "count", 
                                                trans = "log", 
                                                breaks = my_breaks, 
                                                labels = my_labels, 
                                                limits = c(1e-5,1e-3)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))#, axis.text.y=as.character(), axis.ticks.y=element_blank())  

pdf(file = "Marker_genes/Mean_expression_scores_of_clusters_for_all_cell_type_markers.pdf", width = 17, height = 6) # save as horisontal A3
heatmap
dev.off()

#######################################################################################################################################################
# CREATE HEATMAP INCLUDING EVERY CELL TYPE - COUNT ADJUSTED DCA DENOISED EXPRESSION SCORES
#######################################################################################################################################################
cluster_order <- c(10,5,3,1,6,15,13,7,2,0,9,12,4,14,19,18,17,8,11,16,20)

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
print(floor(min(out[!is.na(out)])*10^8)/10^8) # 1e-8
print(floor(max(out[!is.na(out)])*1000)/1000) # 0.006

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
my_breaks <- c(1e-8,1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
my_labels <- c("1e-8","1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2")

out$value[out$value>1e-3] <- 1e-3
out$value[out$value<1e-5] <- 1e-5

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", 
                       breaks = my_breaks, 
                       labels = my_labels, 
                       limits = c(1e-5,1e-3)) +
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

#######################################################################################################################################################
# INPUT MAJOR CELL TYPES FOR CLUSTERS TO DO CELL TYPE SPECIFIC HEATMAPS
#######################################################################################################################################################

cellt_clusters <- c(10,5,3,1,6,15,13,7,2, # CD4
                    0,9,12,4,14, #CD8
                    19, #NK
                    18,17, # B
                    8,11,16, #MNP
                    20) #DC
names(cellt_clusters) <- c(rep("T",times = 14), "NK", "Plasma","B","MNP","MNP","MNP","DC")

#######################################################################################################################################################
# CREATE HEATMAP FOR T cells specific cell type markers - Z every single cell + Z mean
#######################################################################################################################################################
markers <- read.table(file = "Marker_genes/T_sub_set_markers_combined.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,1]),1] 
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 
markers <- markers[!is.na(markers[,2]),]

# collect mean expression for clusters
considered_cl <- cellt_clusters[names(cellt_clusters) == "T"] #| names(cellt_clusters) == "ILC"]
out <- foreach(cl = c(1:length(considered_cl)), .combine = "rbind") %do% {
  exp <- X[rownames(X)%in%markers[,2], colnames(X)%in%id[id[,2]==considered_cl[cl],1]]
  exp <- t(exp[match(markers[,2], rownames(exp)),])
  exp <- exp[,colSums(!is.na(exp))>0]
  exp <- colMeans(exp)
  return(exp)
}
colnames(out) <- markers[match(colnames(out), markers[,2]),1]
rownames(out) <- as.character(considered_cl)
cl_order <- rev(rownames(out))
out <- cluster_matrix(out, dim = "row", method = "ward.D")

# expression per mean
means <- colMeans(out)
out <- t(t(out)/means)
out <- log2(out)

# Create heatmap for all markers
mode(out) <- "numeric" 
print(min(out[!is.na(out)])) # -3.08
print(max(out[!is.na(out)])) # 1.72

out[out < -1.5] <- -1.5
out[out>1.5] <- 1.5

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- seq(from = -1.5, to = 1.5, by = 0.5)
my_labels <- as.character(seq(from = -1.5, to = 1.5, by = 0.5))

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

pdf(file = "Marker_genes/Z_of_mean_expression_scores_for_T_cell_type_markers.pdf", width = 8, height = 3.5)
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
print(min(out[!is.na(out)])) # -8.05
print(max(out[!is.na(out)])) # 4.56

out[out< -1.5] <- -1.5
out[out>1.5] <- 1.5

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
my_col_breaks <- seq(from = -1.5, to = 1.5, by = 0.5)
my_labels <- as.character(seq(from = -1.5, to = 1.5, by = 0.5))

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

pdf(file = "Marker_genes/Z_expression_scores_for_T_cell_type_markers.pdf", width = 8, height = 4)
heatmap
dev.off()

#######################################################################################################################################################
# CREATE HEATMAP FOR T cells specific cell type markers - DCA expression every single cell + DCA mean expression
#######################################################################################################################################################
markers <- read.table(file = "Marker_genes/T_sub_set_markers_combined.txt", sep = "\t", stringsAsFactors = F, header = T) # Corresponding to locally saved supplementary excel file by Martin et al
markers[markers == ""] <- NA
markers <- markers[!is.na(markers[,1]),1]
markers <- cbind(markers, transl[match(markers, transl[,2]),7]) # translate from NCBI symbols to ensembl gene ID 
markers <- markers[!is.na(markers[,2]),]

# collect all cells for cluster
cl_order <- rev(c(10,5,3,1,6,15,13,7,2,0,9,12,4,14))

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
print(floor(min(out[!is.na(out)])*10^6)/10^6) # 1e-6
print(floor(max(out[!is.na(out)])*1000)/1000) # 0.007

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
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3,1e-2)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2")

out$value[out$value > 1e-3] <- 1e-3
out$value[out$value<1e-5] <- 1e-5

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)), na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-5,1e-3)) +
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

pdf(file = "Marker_genes/Expression_scores_for_T_cell_sub_type_markers_2.pdf", width = 8, height = 4)
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
print(min(out[!is.na(out)])) # 1.39e-5
print(max(out[!is.na(out)])) # 3e-3

# format for heatmap
out <- melt(out)
out$Var1 <- factor(x = out$Var1, levels = cl_order, ordered = T)
my_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3)
my_labels <- c("1e-6", "1e-5", "1e-4", "1e-3")

out$value[out$value>1e-3] <- 1e-3
out$value[out$value<1e-5] <- 1e-5

heatmap <- ggplot(out, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster(na.rm = F) + 
  #geom_tile(aes(fill = value), colour = "transparent") +
  scale_fill_gradientn(colours = c("black",viridis(4)),
                       na.value = "transparent", name = "count", trans = "log", breaks = my_breaks, labels = my_labels, limits = c(1e-5,1e-3)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, color = "black"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )      

pdf(file = "Marker_genes/Mean_expression_scores_for_T_cell_sub_type_markers_2.pdf", width = 8, height = 3.5)
heatmap
dev.off()


#######################################################################################################################################################
# DEGs calculation  
#######################################################################################################################################################

temp_id <- as.numeric(id[,2])
names(temp_id) <- id[,1]
temp_id[grepl("_IIH_", id[,1])] <- temp_id[grepl("_IIH_", id[,1])] + 0.1 # mark cells from patients with IIH  
temp_id <- as.factor(temp_id)
ms_all@active.ident <- temp_id
# unique clusters
temp_id_unique <- sort(as.numeric(unique(as.character(temp_id))))

# Exclude clusters with only healthy/sick cells from calculation
for(i in 0:(length(unique(id[,2]))-1)){
  if(any(!(i%in%temp_id_unique),!((i+0.1)%in%temp_id_unique))){
    print(paste("Exclude cluster ", i, " from DEG calculation", sep=""))
    temp_id_unique <- temp_id_unique[!temp_id_unique%in% c(i,(i+0.1))]
  }
}
# Exclude clusters with less than 3 healthy or sick cells from calculation
temp <- table(Idents(ms_all))
if(any(temp<2)){ # excludes no cluster 
  temp <- as.numeric(names(temp[temp<=3]))
  temp_id_unique <- temp_id_unique[!temp_id_unique %in% c(temp, temp+0.1, temp-0.1)]
}

# Output directory
dir.create(path="DCA_MAST_DEGs_logFC_1.5_MS", showWarnings = F)

# Running MAST DEGs analysis
for(i in 1:(length(temp_id_unique)/2)){
  print("Checking cluster:")
  pos <- i*2-1
  print(temp_id_unique[pos:(pos+1)])
  temp <- FindMarkers(ms_all,
                      slot = "counts",
                      test.use = "MAST", 
                      ident.1 = temp_id_unique[pos], # important that this is a factor vector with cluster labels as names
                      ident.2 = temp_id_unique[pos+1],
                      only.pos = F,
                      random.seed = 3,
                      pseudocount.use = 0,
                      logfc.threshold = log(1.5),
                      min.pct = 0.1,
                      min.cells.group = 3
  )
  if(exists("temp")){
    print("Saving ...")
    write.table(temp, file = paste("DCA_MAST_DEGs_logFC_1.5_MS/Cluster_", temp_id_unique[pos], "_res=0.35_dims=32_k=15.txt",sep=""),sep="\t",col.names=NA,row.names=T)
    rm(temp)
  }
}

#######################################
# Running MAST DEGs summary
#######################################
files <- list.files(path = "DCA_MAST_DEGs_logFC_1.5_MS")
files <- files[grepl(files, pattern = "_res=0.35_dims=32_k=15")]
degs.clusters.h.vs.s <- matrix(NA, nrow = 15000, ncol = length(files))
for(i in 1:ncol(degs.clusters.h.vs.s)){
  temp <- read.table(file = paste("DCA_MAST_DEGs_logFC_1.5_MS/", files[i], sep=""), sep="\t", header = T, stringsAsFactors = F)
  temp <- temp[as.numeric(temp$p_val_adj) < 0.05,]
  temp <- as.character(temp[order(as.numeric(temp$p_val_adj), decreasing = F),1]) # order by significance
  print(length(temp))
  if(length(temp)>0){
    degs.clusters.h.vs.s[1:length(temp),i] <- temp
  }
  rm(temp)
}
colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = "_res=0.35_dims=32_k=15.txt"))
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,colSums(!is.na(degs.clusters.h.vs.s))>0]
write.table(degs.clusters.h.vs.s, file = "DCA_MAST_DEGs_logFC_1.5_MS/SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters.txt",sep="\t", col.names = T, row.names = F)

#######################################
# Translate MAST DEGs for further analysis - translate to Entrez Gene IDs
#######################################
transl <- as.matrix(read.table(file = paste(fp, "/../Input/HGNC translation matrix 201108/transl.txt", sep=""), sep="\t", header = T, stringsAsFactors = F))
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
write.table(out, file = "DCA_MAST_DEGs_logFC_1.5_MS/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt",sep="\t", col.names = T, row.names = F)

#######################################
# Translate MAST DEGs & background genes for further analysis - translate to official approved HGNC symbols
#######################################
transl <- as.matrix(read.table(file = paste(fp, "/../Input/HGNC translation matrix 201108/transl.txt", sep=""), sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,2]),] # Official approved HGNC symbol
transl <- transl[!is.na(transl[,7]),] # Ensemble Gene ID
print(head(transl))

# Translate DEGs
###########################
out <- matrix(NA, ncol = ncol(degs.clusters.h.vs.s), nrow = nrow(degs.clusters.h.vs.s))
colnames(out) <- colnames(degs.clusters.h.vs.s)
for(i in 1:ncol(degs.clusters.h.vs.s)){
  print(paste("Working on", colnames(degs.clusters.h.vs.s)[i], sep = " "))
  temp <- degs.clusters.h.vs.s[,i]
  temp <- temp[!is.na(temp)]
  print(paste("n of DEGs: ", length(temp), sep=""))
  temp <- transl[match(temp, transl[,7]),2] # keeps order
  temp <- unique(temp[!is.na(temp)])
  print(paste("translated to n of DEGs: ", length(temp), sep=""))
  out[1:length(temp),i] <- temp
}
out <- out[rowSums(!is.na(out))>0,colSums(!is.na(out))>0]
write.table(out, file = "DCA_MAST_DEGs_logFC_1.5_MS/TRANSLATED_HGNC_SYMBOL_sig_MAST_DEGs_log(1.5)_all_clusters.txt",sep="\t", col.names = T, row.names = F)

#######################################################################################################################################################
# INPUT FINAL CELL TYPING
#######################################################################################################################################################
library(scales)
cellt_clusters <- c(10,5,3,1,6,15,13,7,2, # CD4
                    0,9,12,4,14, #CD8
                    19, #NK
                    18,17, # B
                    8,11,16, #MNP
                    20) #DC
names(cellt_clusters) <- c(rep("CD4 T",times=9), rep("CD8 T", times = 5), "NK", rep("B",times = 2), rep("MNP",times=3), "DC") 


cellt_clusters <- cbind(cellt_clusters, names(cellt_clusters), 
                        c("Resting CD4+ Tm 1",
                          "Resting CD4+ Tm 2", 
                          "Resting CD4+ Tm 3",
                          "Resting CD4+ Tm 4",
                          "Naive CD4+ T 1",
                          "Naive CD4+ T 2",
                          "Early activated CD4+ memory T cell",
                          "Activated CD4+ memory T cell 1",
                          "Activated CD4+ memory T cell 2",
                          
                          "Naive CD8+ T cell",
                          "Early activated CD8+",
                          "Activated CD8+ effector memory 1",
                          "Activated CD8+ effector memory 2",
                          "Activated terminally differentiated effector CD8+ T cell",
                          
                          "Natural killer cell",
                          
                          "Plasma B cell", 
                          "Naive B cell",
                          
                          "Mononuclear phagocyte 1", 
                          "Mononuclear phagocyte 2",
                          "Mononuclear phagocyte 3",
                          
                          "Dentritic cell")
                        )

col_palette <- c( "#2b58a3","#3686D3","#4AA4DE","#87CEFA","#2F9DA7","#009999",
                  "#6666FF","#4545E2","#000080", # CD4+
                  "#553c93","#82318f","#882865","#4B0082","#330066", # CD8+
                  "#000000", # NK
                  "#4F7942", "#BFE890", # plasma + naive B
                  "#e8c124","#e7772b","#e92925", # MNP
                  "#A9A9A9") # DC
show_col(col_palette)

cellt_clusters <- cbind(cellt_clusters, col_palette)
write.table(cellt_clusters, file = "../../Output/MS/Cluster_colors.txt",sep="\t",col.names =  T, row.names = F)

temp_id <- floor(as.numeric(as.character(Idents(ms_all))))
temp_id <- as.factor(temp_id)
names(temp_id) <- names(Idents(ms_all))
ms_all@active.ident <- temp_id

# ALL clusters
pdf(file = "../../Output/MS/Cluster_plot_own_color_palette_colorblind_friendly_all_cells.pdf", width = 13, height = 6)
TSNEPlot(ms_all, label.size = 0) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

# Sick cells only
pdf(file = "../../Output/MS/Cluster_plot_own_color_palette_colorblind_friendly_MS_cells.pdf", width = 13, height = 6)
TSNEPlot(ms_all, label.size = 0, cells = grep("MS", Cells(ms_all), value =T)) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

# Healthy cells only
pdf(file = "../../Output/MS/Cluster_plot_own_color_palette_colorblind_friendly_IIH_cells.pdf", width = 13, height = 6)
TSNEPlot(ms_all, label.size = 0, cells = grep("IIH", Cells(ms_all), value =T)) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()


ms_all <- RunUMAP(ms_all, dims = 1:ncol(ms_all@reductions$pca@cell.embeddings))

# ALL clusters
pdf(file = "../../Output/MS/Cluster_plot_own_color_palette_colorblind_friendly_all_cells_UMAP.pdf", width = 13, height = 6)
UMAPPlot(ms_all, label.size = 0) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

# Sick cells only
pdf(file = "../../Output/MS/Cluster_plot_own_color_palette_colorblind_friendly_MS_cells_UMAP.pdf", width = 13, height = 6)
UMAPPlot(ms_all, label.size = 0, cells = grep("MS", Cells(ms_all), value =T)) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

# Healthy cells only
pdf(file = "../../Output/MS/Cluster_plot_own_color_palette_colorblind_friendly_IIH_cells_UMAP.pdf", width = 13, height = 6)
UMAPPlot(ms_all, label.size = 0, cells = grep("IIH", Cells(ms_all), value =T)) + 
  scale_color_manual(values = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),4]),
                     labels = as.character(cellt_clusters[order(as.numeric(cellt_clusters[,1])),3]))
dev.off()

#######################################################################################################################################################
# Cell proportion plot in color friendly colors
#######################################################################################################################################################

# Create cellID / cluster ID matrix
id <- cbind(as.character(names(Idents(ms_all))), as.character(Idents(ms_all)))
colnames(id) <- c("Cell_ID", "Cluster_ID")

# split by patients
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

# Make table showing n cells per sample and cluster
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

# in percent
out <- 100*out/rowSums(out)
write.table(rbind(colnames(out),cellt_clusters[match(colnames(out),paste("Cluster_",cellt_clusters[,1],sep="")),3],out), file = "MS_%_of_cells_of_individual_samples_per_cluster.txt", sep="\t", row.names = T, col.names = NA)

# make stacked bar graph
out <- cbind(rownames(out),out) 
# alternative to above: out <- read.table(file = "MS_%_of_cells_of_individual_samples_per_cluster.txt", sep="\t", header = T)
colnames(out)[1] <- "GEO_ID"
# rename values in column GEO_ID based on "MS" or "IIH" sample
temp <- cbind(c("GSM4104122","GSM4104123","GSM4104124","GSM4104125",
                "GSM4104126","GSM4104127","GSM4104128","GSM4104129",
                "GSM4104130","GSM4104131","GSM4104132","GSM4104133"),
              c(rep("MS",times = 6), rep("IIH",times = 6)))
out[,1] <- temp[match(temp[,1], out[,1]),2]
out <- cbind(out, NA)
out[,ncol(out)] <- c(paste("Sample", c(1,2,3,4,5,6),sep=" "), paste("Sample",c(1,2,3,4,5,6),sep=" "))
colnames(out)[ncol(out)] <- "Patient"
rownames(out) <- NULL
out <- as.data.frame(out)
out$Patient <- factor(x = out$Patient, levels = paste("Sample", c(1,2,3,4,5,6),sep=" "), ordered = T)
for(i in 2:(ncol(out)-1)){
  out[,i] <- as.numeric(as.character(out[,i]))
}
out <- melt(out, id.vars = c("Patient", "GEO_ID"))
out$variable <- as.character(cellt_clusters[match(out$variable, paste("Cluster_",cellt_clusters[,1],sep="")), 3])
out$variable <- factor(x = out$variable, levels = cellt_clusters[,3], ordered = T)

pdf(file = "Percentage_of_total_cells_per_cluster_individual_MS_patients_color_blind_friendly.pdf", width = 12, height = 5)

ggplot(out, aes(x = Patient, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ GEO_ID, margins = F) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, margin = margin(-13,0,0,0), colour = "black", family = "Helvetica"),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", family = "Helvetica", face = "bold"),
        axis.title = element_blank()) +
  scale_fill_manual(values = as.character(cellt_clusters[,4]))

dev.off()




saveRDS(ms_all, file = "MS_batch_corrected_clustered.rds")
rm(list = ls())

