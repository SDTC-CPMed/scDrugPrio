#
# Combine DEG lists CD for Supplementary data
#
# SAMUEL SCHAEFER
# 2022-01-21
#################################################


library(doParallel)


# load transl
transl <- as.matrix(read.table(file = "../Input/HGNC translation matrix 201108/transl.txt", sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,7]),] # Ensembl Gene ID
print(head(transl))

setwd("../Output/DCA_MAST_DEGs/")

# Identify DEG files
lf <- list.files(pattern = "Cluster_")

degs <- foreach(i = c(1:length(lf)), .combine = "rbind") %do% {
  temp <- read.delim(file = lf[i], stringsAsFactors = F, header = T)
  if(nrow(temp)>0){
    cluster <- unlist(strsplit(unlist(strsplit(lf[i], split = "Cluster_")), split = "_res=0.45_k=15.txt"))
    temp <- cbind(temp, paste("Cluster ",cluster, sep=""))
    return(temp)
  } else {
    return(rep(NA, times = 4))
  }
}

# Translate to Entrez ID and human gene symbol
degs <- cbind(degs, transl[match(degs[,1],transl[,2]),6])
degs <- degs[,c(1,8,2,3,6,7)]
colnames(degs) <- c("Entrez_gene_ID", "Entrez_symbol", "p_val", "avg_logFC", "p_val_FDR_adj", "Cluster")

write.table(degs, file = "ALL_DEGs_AND_TRANSLATIONS_FOR_SUPPLEMENT.txt",sep="\t", col.names = T, row.names = F)

