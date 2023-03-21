
#######################################################################################
# FOR CD DATA SET - all individuals patients 
# BY SAMUEL SCHAEFER
#######################################################################################

source("From scDrugPrio 220925/background_genes_NicheNet.R")
source("From scDrugPrio 220925/NicheNet_ligand_activity_analysis.R")
source("From scDrugPrio 221026/NicheNet_cell_type_centrality.R")

fp <- getwd()
setwd(paste(fp, "/../", sep=""))
library(doParallel)


#######################################################################################
# NicheNet for all patients
#######################################################################################
patient <- paste("Patient_", 1:11, sep = "")
out_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
in_path <- "Input/CD GSE134809/Individual_patients/"

# Translate MAST DEGs for further analysis - translate to gene symbols
transl <- as.matrix(read.table(file = "Input/HGNC translation matrix 201108/transl.txt", sep="\t", header = T, stringsAsFactors = F))
transl <- transl[!is.na(transl[,2]),] # Entrez Gene Symbol
transl <- transl[!is.na(transl[,7]),] # Ensemble Gene ID
print(head(transl))

# Load Seurat object
denoised_DCA_clustered <- readRDS(file = "Input/CD GSE134809/CD_all_v2.rds")

for(pat in patient){
  
  print(paste("Processing: ", pat, sep=""))
  set.seed(4)
  
  # format cell IDs to fit in function
  id <- read.table(file = paste(in_path,"Cell_identities_", gsub("P", "p", pat),".txt", sep=""), sep = "\t", header = T, stringsAsFactors = F)
  unique_cl <- unique(id[,2])
  cell_IDs <- foreach(i = c(1:length(unique_cl))) %do% {
    temp <- id[as.character(id[,2]) == unique_cl[i],1]
    return(temp)
  }
  names(cell_IDs) <- paste("Cluster_", unique_cl, sep="")
  
  # get DEGs
  degs_transl <- as.matrix(read.table(file = paste(in_path, pat, "/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_HGNC_SYMBOL_sig_MAST_DEGs_log(1.5)_all_clusters.txt", sep =""), sep="\t", header = T, stringsAsFactors = F))
  
  # background gene calculation
  bg <- background_genes_NicheNet(data = as.matrix(denoised_DCA_clustered@assays$RNA@counts), cell_IDs = cell_IDs)
  
  # Translate background genes
  for(i in 1:ncol(bg)){
    temp <- transl[match(x = bg[,i], table =  transl[,7]),2]
    temp <- temp[!is.na(temp)]
    bg[,i] <- NA
    if(length(temp)>0){
      bg[1:length(temp),i] <- temp
    }
  }
  bg <- bg[rowSums(!is.na(bg))>0, colSums(!is.na(bg))>0]
  
  # NicheNet ligand activity
  dir.create(path = paste(out_path,pat,"/NicheNet",sep=""), showWarnings = F)
  all_ligand_activity <- NicheNet_ligand_activity_analysis(degs = degs_transl, background_genes = bg, out_dir = paste(out_path,pat,"/NicheNet",sep=""), cores = 5)
  #print(head(all_ligand_activity))
  
  # Calculate centrality for cell types based on ligand activity outcomes
  intercellular_centrality <- NicheNet_cell_type_centrality(all_ligand_activity = all_ligand_activity, out_dir = paste(out_path,pat,"/NicheNet",sep=""))
  #print(intercellular_centrality)
}



###############################################################################
# Output for Cytoscape
###############################################################################

patient <- paste("Patient_", 1:11, sep = "")
out_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
in_path <- "Input/CD GSE134809/Individual_patients/"

for(pat in patient){
  all_ligand <- read.table(file = paste(out_path,pat,"/NicheNet/all_ligand_activity.txt",sep=""), sep = "\t", header = T, stringsAsFactors = F)
  if(pat %in% c("Patient_1", "Patient_10")){
    color <- read.table(file = paste(out_path,pat,"/Cluster_colors.txt",sep=""), sep="\t", header = T, stringsAsFactors = F)[, - c(2:3)]
  } else {
    color <- read.table(file = paste(out_path,pat,"/Automatic_cluster_colors.txt",sep=""), sep="\t", header = T, stringsAsFactors = F)
    temp <- unlist(strsplit(color[,1], split = paste(pat, "__",sep="")))
    temp <- temp[temp != ""]
    color[,1] <- temp
  }

  cl1 <- unique(all_ligand[,6])
  cl2 <- unique(all_ligand[,7])
  cl <- unique(c(cl1, cl2))
  all_ligand <- all_ligand[as.numeric(all_ligand[,4]) > 0,-5]
  
  # summed pearsson coefficients for plotting
  out_ligand <- foreach(s = c(1:length(cl1)), .combine = "rbind") %do% {
    
   if(any(all_ligand[,5] %in% cl1[s])){
      out <- foreach(t = c(1:length(cl2)), .combine = "rbind") %do% {
        if(any((all_ligand[,6] %in% cl2[t]) & (all_ligand[,5] %in% cl1[s]))){
          if(sum((all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t])) >= 1){
            temp <- all_ligand[(all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t]),]
            temp <- c(NA, NA, nrow(temp), sum(as.numeric(temp[,4])), cl1[s], cl2[t])
            names(temp) <- NULL
            return(matrix(temp, nrow = 1))
          } else if(sum((all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t])) == 0){
            return(rep(NA, times = ncol(all_ligand)))
          } else {
            temp <- as.vector(all_ligand[(all_ligand[,5] %in% cl1[s]) & (all_ligand[,6] %in% cl2[t]),])
            temp[,3] <- 1 # n of ligand interactions
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
  out_ligand <- cbind(out_ligand, color[match(out_ligand[,5], paste("Cluster_",as.numeric(color[,1]),sep="")),2])
  colnames(out_ligand)[7] <- "color"
  
  out_ligand <- as.data.frame(out_ligand[,-2])
  out_ligand$test_ligand <- as.character(out_ligand$test_ligand)
  out_ligand$n_ligand_interactions <- as.numeric(out_ligand$n_ligand_interactions)
  out_ligand$pearson <- as.numeric(out_ligand$pearson)
  out_ligand$Sender <- as.character(out_ligand$Sender)
  out_ligand$Target <- as.character(out_ligand$Target)
  out_ligand$color <- as.character(out_ligand$color)
  out_ligand$color2 <- as.character(out_ligand$color)
  
  write.table(x = out_ligand, file = paste(out_path,pat,"/NicheNet/",pat,"_CD_MCDM_own_colors_pearson_summed.txt",sep=""), sep = "\t", col.names = T, row.names = F)
}




