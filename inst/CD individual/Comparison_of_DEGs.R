#
#
#
#
#

library(doParallel)

set.seed(54)

#setwd("/data/samsc76/SIGMA_13012020/drug_prediction_R")
#setwd("~/Library/CloudStorage/OneDrive-Linköpingsuniversitet/Old RA mouse project - Barabasi/SIGMA_13012020/drug_prediction_R")

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

# ACTUAL DEGS BUT RELOADED

lf <- list.files(path = "Output/DCA_MAST_DEGs", pattern = "^MAST_DEGs_log", full.names = T)
actual_DEGs <- foreach(i = c(1:length(lf)), .combine = "cbind") %do% {
  temp <- as.matrix(read.table(file = lf[i], sep="\t", header = T, stringsAsFactors = F))
  temp <- temp[order(as.numeric(temp[,6])),]
  if(sum(temp[,6] < 0.05)>2){
    temp <- temp[as.numeric(temp[,6]) < 0.05,]
  } else if(sum(temp[,6] < 0.05) == 1){
    temp <- matrix(temp[as.numeric(temp[,6]) < 0.05,], nrow = 1)
  } else {
    
  }
  temp <- temp[,1]
  if(length(temp) < 20000){
    temp <- c(temp, rep(NA, times = 20000 - length(temp)))
  }
  return(temp)
}
temp <- list.files(path = "Output/DCA_MAST_DEGs", pattern = "^MAST_DEGs_log")
temp <- unlist(strsplit(x = temp, split = "MAST_DEGs_log\\(1\\.5\\)_threshold_cluster_"))
temp <- unlist(strsplit(x = temp, split = "_res=0.6_dims=32.txt"))
colnames(actual_DEGs) <- paste("Cluster_", temp, sep ="")
actual_DEGs <- actual_DEGs[rowSums(!is.na(actual_DEGs))>0,]
write.table(actual_DEGs, file = "Output/DCA_MAST_DEGs/Ordered_summary_DEGs_for_comparisons_221018.txt",sep="\t", col.names = T, row.names = F)

transl <- read.table(file = "Input/NCBI_annotation/transl.txt", sep="\t", header = T, stringsAsFactors = F)

transl_actual_DEGs <- actual_DEGs
for(i in 1:ncol(transl_actual_DEGs)){
  temp <- transl_actual_DEGs[,i]
  temp <- temp[!is.na(temp)]
  temp <- transl[match(temp, transl[,3]),1]
  temp <- temp[!is.na(temp)]
  temp <- c(temp, rep(NA, times = nrow(transl_actual_DEGs)-length(temp)))
  transl_actual_DEGs[,i] <- temp
}
transl_actual_DEGs <- transl_actual_DEGs[rowSums(!is.na(transl_actual_DEGs))>0, colSums(!is.na(transl_actual_DEGs))>0]

# FILE USED FOR ALL DRUG PREDICTIONS

used_degs <- read.table("Output/DCA_MAST_DEGs/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt", sep="\t", header = T, stringsAsFactors = F)
used_degs <- as.matrix(used_degs)
mode(used_degs) <- "numeric"

ordered_used_degs <- read.table("Output/DCA_MAST_DEGs/TRANSLATED_Summary_sig_adj_MAST_DEGs_all_clusters_res=0.6_dims=32.txt", sep="\t", header = T)
ordered_used_degs <- as.matrix(ordered_used_degs)
mode(ordered_used_degs) <- "numeric"

clustering <- read.table(file = "Output/Seurat_clusters/Cell_identity.txt",sep="\t", header = T)



