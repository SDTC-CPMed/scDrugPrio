

# NR vs R DEGs for Doctis presentation
# Run in R 4.0.4 due to cluster limitations.

#install.packages("Seurat")
library(Seurat)
library(doParallel)
registerDoParallel(cores = 15)
set.seed(21)

# get clustered PSA data
###############################################################################

# k = 15 and resolution = 0.4 in the above loaded RDS
#psa <- FindNeighbors(psa, k.param = 15, nn.method = "annoy")
#psa <- FindClusters(psa, resolution = 0.8)
#saveRDS(object=psa, file="../Input/PsA_data/PsA_all_clustered.rds")
psa <- readRDS(file = "../Input/PsA_data/PsA_all_clustered.rds")


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

# unique deg_ids
temp_id_unique <- sort(as.numeric(unique(as.character(temp_id))))
print(temp_id_unique)

# Exclude clusters with less than 3 healthy or sick cells from calculation
# => excludes no cluster
print(table(Idents(psa)))

# Output directory
dir.create(path="../Output/DCA_MAST_DEGs_R_vs_NR", showWarnings = F)

for(i in 1:(length(unique_cl))){
  print(paste("Checking cluster:", unique_cl[i],sep=" "))
  if(any(temp_id %in% (unique_cl[i]+0.6))){
    if(any(temp_id %in% temp_id[temp_id %in% (unique_cl[i]+0.1)])){
      
      temp <- FindMarkers(psa,
                          slot = "counts",
                          test.use = "MAST", 
                          ident.1 = temp_id[temp_id %in% (unique_cl[i]+0.1)],
                          ident.2 = temp_id[temp_id %in% (unique_cl[i]+0.6)], # important that this is a factor vector with cell labels as names
                          only.pos = F,
                          random.seed = 3,
                          pseudocount.use = 0,
                          logfc.threshold = log(1.5),
                          min.pct = 0.1,
                          min.cells.group = 3)
    }
    
    if(exists("temp")){
      print("Saving ...")
      write.table(temp, file = paste("../Output/DCA_MAST_DEGs_R_vs_NR/Cluster_", unique_cl[i], "_R_vs_NR_res=0.8_k=15.txt",sep=""),sep="\t",col.names=NA,row.names=T)
      rm(temp)
    }
  }
}


# DEG summary - responder vs non-responder
###############################################################################

outdir = "../Output/DCA_MAST_DEGs_R_vs_NR"
files <- list.files(path = outdir, pattern = "_R_vs_NR_res=0.8_k=15.txt")

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
colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = "_R_vs_NR_res=0.8_k=15.txt"))

# remove empty rows and columns
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,]
head(degs.clusters.h.vs.s)

# save
write.table(degs.clusters.h.vs.s, file = paste(outdir, "/SUMMARY_MAST_DEGs_all_clusters_R_vs_NR.txt", sep=""),sep="\t", col.names = T, row.names = F)

