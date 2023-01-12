#
# Combine DEG lists RA for Supplementary data
#
# SAMUEL SCHAEFER
#
#################################################

library(doParallel)

# load transl
transl <- as.matrix(read.table(file = "../Input/Human-mouse_homologs/transl.txt", sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,2]),] # Entrez Gene ID
transl <- transl[!is.na(transl[,3]),] # Mouse gene symbol
print(head(transl))

setwd("../Output/DCA_MAST_DEGs/")

# Identify DEG files  
lf <- list.files(pattern = ".txt")
lf <- lf[grepl(lf, pattern = "Cluster_")]

degs <- foreach(i = c(1:length(lf)), .combine = "rbind") %do% {
  temp <- read.table(file = lf[i], stringsAsFactors = F, header = T)
  cluster <- unlist(strsplit(unlist(strsplit(lf[i], split = "Cluster_")), split = "_res=0.6_dims=32_k=15.txt"))
  temp <- cbind(temp, paste("Cluster ",cluster, sep=""))
  return(temp)
}

# Translate to mouse & human Entrez ID and human gene symbol
degs <- cbind(degs, transl[match(degs[,1],transl[,3]),c(1,2,4)])
degs <- degs[,c(1,8,10,9,2,3,6,7)]
colnames(degs) <- c("mouse_gene_symbol", "human_gene_symbol", "mouse_entrez_ID", "human_entrez_ID", "p_val", "avg_logFC", "p_val_FDR_adj", "Cluster")

write.table(degs, file = "ALL_DEGs_AND_TRANSLATIONS_FOR_SUPPLEMENT.txt",sep="\t", col.names = T, row.names = F)

