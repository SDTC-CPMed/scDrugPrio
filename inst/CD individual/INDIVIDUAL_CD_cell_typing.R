#
#
# Cell type clusters of individual patients
#
################################################################################################

library(doParallel)
library(matrixStats)



# Derive clustered cell ID data for each patient
pat <- paste("patient_",1:11,sep="")
Pat <- paste("Patient_",1:11,sep="")

dir <- "../Input/CD GSE134809/Individual_patients/"

cluster_out <- foreach(i = c(1:length(pat)), .combine = "rbind") %do% {
  temp <- cbind(read.table(file = paste(dir, "Cell_identities_",pat[i],".txt",sep=""),sep="\t", header = T), Pat[i], stringsAsFactors = F)
  return(temp)
}
cluster_out[,1] <- gsub(pattern = "-", replacement = "\\.", x = cluster_out[,1])
head(cluster_out)

# Derive pooled cell typing
pooled_clust <- read.table(file = "../CD_batch_corrected/Input/CD/CD_all/Cell_ID_to_cluster_ID_k=15_res=0.8.txt", sep="\t", header = T, stringsAsFactors = F)
pooled_cell_types <- read.table(file = "../CD_batch_corrected/Output/CD/Cluster_colors.txt",sep="\t", header = T, stringsAsFactors = F)

# Divide pooled clustering into clusters
cl <- as.numeric(pooled_cell_types[,1])
pooled_id <- foreach(i = c(1:length(cl)), .combine = "cbind") %do% {
  temp <- pooled_clust[pooled_clust[,2] %in% cl[i],1]
  temp <- c(temp, rep(NA, times = 10000-length(temp)))
  return(temp)
}
pooled_id <- pooled_id[rowSums(!is.na(pooled_id))>0,]
colnames(pooled_id) <- cl

# For every patient, check overlap between pooled clusters and individual clusters
out <- foreach(i = c(1:length(Pat)), .combine = "rbind") %do% {
  temp <- cluster_out[cluster_out[,3] %in% Pat[i],]
  ind_cl <- sort(unique(temp[,2])) 
  temp_out <- foreach(x = c(1:length(ind_cl)), .combine = "rbind") %do% {
    overlap <- vector() 
    for(z in 1:ncol(pooled_id)){
      overlap[z] <- sum(pooled_id[!is.na(pooled_id[,z]),z] %in% temp[temp[,2] %in% ind_cl[x],1])
    }
    return(overlap)
  }
  rownames(temp_out) <- paste(Pat[i], "__", ind_cl,sep="")
  return(temp_out)
}
colnames(out) <- pooled_cell_types[,5]

write.table(out, file ="../Input/CD GSE134809/Individual_patients/Overlap_with_batch_corrected_pooled_cell_typing.txt", col.names = T, row.names = T, sep="\t")

# Divide pooled clustering into major cell types
pooled_cell_types <- cbind(pooled_cell_types, 
                           c("Treg", rep("CD4+ T", times = 9), rep("CD8+ T", times = 4), "IgM plasma", 
                             rep("IgG plasma", times = 5), rep("IgA plasma", times = 3), rep("Memory B", times = 3), 
                             rep("Naive B", times = 4), "Memory B",rep("MNP", times = 4), rep("Stroma cell", times = 5))
                           )

cl <- unique(pooled_cell_types[,6])
pooled_id <- foreach(i = c(1:length(cl)), .combine = "cbind") %do% {
  temp <- pooled_clust[pooled_clust[,2] %in% pooled_cell_types[pooled_cell_types[,6]%in% cl[i],1], 1]
  temp <- c(temp, rep(NA, times = 30000-length(temp)))
  return(temp)
}
pooled_id <- pooled_id[rowSums(!is.na(pooled_id))>0,]
colnames(pooled_id) <- cl

# For every patient, check overlap between pooled clusters and individual clusters
out <- foreach(i = c(1:length(Pat)), .combine = "rbind") %do% {
  temp <- cluster_out[cluster_out[,3] %in% Pat[i],]
  ind_cl <- sort(unique(temp[,2])) 
  temp_out <- foreach(x = c(1:length(ind_cl)), .combine = "rbind") %do% {
    overlap <- vector() 
    for(z in 1:ncol(pooled_id)){
      overlap[z] <- sum(pooled_id[!is.na(pooled_id[,z]),z] %in% temp[temp[,2] %in% ind_cl[x],1])
    }
    return(overlap)
  }
  rownames(temp_out) <- paste(Pat[i], "__", ind_cl,sep="")
  return(temp_out)
}
colnames(out) <- cl

# Set color according to above cell typing
colors <- cbind(Treg = c("#0000FF",NA,NA,NA,NA),
                CD4 = c("#2005A5","#1E90FF","#87CEFA","#3686D3", "#00BFFF"),
                CD8 = c("#553c93","#82318f","#882865", "#4B0082", NA),
                IgM = c("#228B22", NA, NA, NA, NA),
                IgG = c("#006400", "#556B2F","#6B8E23", "#9ACD32", NA),
                IgA = c("#7CFC00","#00FF00","#9EF38F", NA, NA),
                MemB = c("#32CD32","#2E8B57","#3CB371","#66CDAA", NA),
                NaiveB = c("#8FBC8F","#BFE890", "#AFE1AF","#ADFF2F", NA),
                MNP = c("#D3D3D3","#808080","#000000", NA, NA),
                Stroma = c("#FF0000", "#660000","#f49728","#F0E68C", NA))

out2 <- foreach(i = c(1:length(Pat)), .combine = "rbind") %do% {
  
  temp <- out[which(grepl(pattern = paste(Pat[i],"__",sep=""),x = rownames(out))),]
  xy <- which(temp == rowMaxs(temp), arr.ind = T)
  ct <- unique(xy[,2])
  
  labels <- matrix(NA, ncol = 2, nrow = nrow(temp))
  for(x in 1:length(ct)){
    n <- sum(xy[,2] == ct[x])
    if(n > 5){
      warning(paste("Too many similar cell types at i = ", i, " and x = ", x,sep=""))
    }
    labels[xy[xy[,2] %in% ct[x],1],1] <- colors[1:n,ct[x]]
    if(n == 1){
      labels[xy[xy[,2] %in% ct[x],1],2] <- colnames(out)[ct[x]]
    } else {
      labels[xy[xy[,2] %in% ct[x],1],2] <- paste(colnames(out)[ct[x]], " ", 1:n,sep="")
    }
  }
  colnames(labels) <- c("Color", "Cell_type")
  write.table(cbind(Cluster = rownames(temp),labels), file = paste("../Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/", Pat[i], "/Automatic_cluster_colors.txt", sep=""),
              sep="\t", col.names = T, row.names = F)
  
  return(labels)
}
colnames(out2) <- c("Color", "Cell_type")
out <- cbind(out, out2)

write.table(out, file ="../Input/CD GSE134809/Individual_patients/Overlap_with_batch_corrected_pooled_major_cell_typing.txt", col.names = T, row.names = T, sep="\t")

