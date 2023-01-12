
# NICHENET ANALYSIS
# Run on R 4.0.4
######################################################################

# Load clusteres joint_RA data
#joint_RA <- readRDS(file = "../Input/joint_RA_data/joint_RA_all_clustered.rds")
rm(list = ls()[ls() != "joint_RA"])

source("Codes_from_scPred/background_genes_NicheNet.R")
source("Codes_from_scPred/NicheNet_cell_type_centrality.R")
source("Codes_from_scPred/NicheNet_ligand_activity_analysis.R")
library(doParallel)
library(Seurat)

id <- Idents(joint_RA)
id <- cbind(cell_ID = as.character(names(id)), cluster_id = as.numeric(as.character(id)), degs_id = "Sick")
id[id[,1] %in% grep("Healthy", Cells(joint_RA), value =T),3] <- "Healthy"
print(table(id[,3]))

###############################################################################
# Analysis
###############################################################################

# format cell IDs to fit in function
cell_IDs <- joint_RA$seurat_clusters
unique_cl <- unique(as.numeric(as.character(cell_IDs)))
cell_IDs <- foreach(i = c(1:length(unique_cl))) %do% {
  temp <- cell_IDs[as.character(cell_IDs) == unique_cl[i]]
  temp <- names(temp)
  temp <- temp[temp %in% id[id[,3] %in% c("Healthy", "Sick"),1]]
  return(temp)
}
names(cell_IDs) <- paste("Cluster_", unique_cl, sep="")

# Set out_dir
out_dir <- "../Output/NicheNet"
dir.create(path = out_dir, showWarnings = F)

# transl
transl <- read.table(file = "../Input/Human-mouse_homologs/transl.txt",sep="\t", header = T)

# background gene calculation
bg <- background_genes_NicheNet(data = as.matrix(joint_RA@assays$RNA@counts), cell_IDs = cell_IDs)
for(i in 1:ncol(bg)){
  temp <- bg[,i]
  bg[,i] <- NA
  temp <- transl[transl[,3] %in% temp,1]
  if(length(temp)>0){
    bg[1:length(temp),i] <- temp
  }
}
bg <- bg[rowSums(!is.na(bg))>0,colSums(!is.na(bg))>0]
write.table(bg, file = paste(out_dir, "/background_genes.txt", sep=""), sep="\t", col.names = T, row.names = F)

# DEGs
degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters.txt", header = T, stringsAsFactors = F, sep="\t"))
for(i in 1:ncol(degs)){
  temp <- degs[,i]
  degs[,i] <- NA
  temp <- transl[transl[,3] %in% temp,1]
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
write.table(degs, file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_symbol.txt", sep="\t", col.names = T, row.names = F)


# NicheNet ligand activity
all_ligand_activity <- NicheNet_ligand_activity_analysis(degs = degs, background_genes = bg, out_dir = out_dir, cores = 30)
print(head(all_ligand_activity))

#Centrality analysis
intercellular_centrality <- NicheNet_cell_type_centrality(all_ligand_activity, out_dir = out_dir)
print(intercellular_centrality)


# Output for Cytoscape
########################################

all_ligand <- read.table(file = "../Output/NicheNet/all_ligand_activity.txt", sep = "\t", header = T, stringsAsFactors = F)
all_ligand <- all_ligand[as.numeric(all_ligand[,4]) > 0,]
color <- read.table(file = "../Output/Clustering/Cluster_colors.txt", sep="\t", header = T, stringsAsFactors = F)


cl1 <- unique(all_ligand[,5])
cl2 <- unique(all_ligand[,6])
cl <- unique(c(cl1, cl2))

# summed pearsson coefficients for plotting
out_ligand <- foreach(s = c(1:length(cl1)), .combine = "rbind") %do% {
  
  if(any(all_ligand[,5] %in% cl1[s])){
    out <- foreach(t = c(1:length(cl2)), .combine = "rbind") %do% {
      
      if(any((all_ligand[,6] %in% cl2[t]) & (all_ligand[,5] %in% cl1[s]))){
        if(sum((all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t])) > 1){
          temp <- all_ligand[(all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t]),]
          temp <- c(NA, NA, nrow(temp), sum(as.numeric(temp[,4])), cl1[s], cl2[t])
          names(temp) <- NULL
          return(matrix(temp, nrow = 1))
        } else if(sum((all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t])) == 0){
          return(rep(NA, times = ncol(all_ligand)))
        } else {
          temp <- as.vector(all_ligand[(all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t]),])
          temp[3] <- 1 # n of ligand interactions
          names(temp) <- NULL
          return(matrix(temp, nrow = 1))
        }
      } else {
        return(matrix(rep(NA, times = ncol(all_ligand)), nrow = 1))
      }
      
    }
    
    colnames(out) <- NULL
    return(out)
  } else {
    return(matrix(rep(NA, times = ncol(all_ligand)), nrow = 1))
  }
}
colnames(out_ligand) <- colnames(all_ligand)
colnames(out_ligand)[3] <- "n_ligand_interactions"
out_ligand <- out_ligand[!is.na(out_ligand[,5]),]

# add node colors through extra rows
temp <-cbind(rep("node_color",times = length(cl)), NA, NA, NA, cl, NA)
colnames(temp) <- colnames(out_ligand)
out_ligand <- rbind(out_ligand, temp)
out_ligand <- cbind(out_ligand, color[match(out_ligand[,5], paste("Cluster_",as.numeric(color[,1]),sep="")),4])
colnames(out_ligand)[7] <- "color"

out_ligand <- as.data.frame(out_ligand[,-2])
out_ligand$test_ligand <- as.character(out_ligand$test_ligand)
out_ligand$n_ligand_interactions <- as.numeric(out_ligand$n_ligand_interactions)
out_ligand$pearson <- as.numeric(out_ligand$pearson)
out_ligand$Sender <- as.character(out_ligand$Sender)
out_ligand$Target <- as.character(out_ligand$Target)
out_ligand$color <- as.character(out_ligand$color)
out_ligand$color2 <- as.character(out_ligand$color)

write.table(x = out_ligand, file = "../Output/NicheNet/MCDM_joint_AIA_own_colors_pearson_summed.txt", sep = "\t", col.names = T, row.names = F)

rm(list = ls())



