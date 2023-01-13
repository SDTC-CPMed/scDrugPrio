# Drug prediction based on pathways (sig enriched with DCA MAST DEGs) that overlap between x clusters

library(doParallel)
fp <- getwd()
setwd(paste(fp, "/..", sep=""))

set.seed(35)

# Load KEGG pathways
kegg_pathways <- as.matrix(read.table(file = "Input/KEGG_pathways/TRANSLATED_kegg_genes_by_pathway.txt", sep="\t", header = T))

# Load background genes (all genes measured with scRNAseq) - depending on that 
bg <- read.table("Input/Human-mouse_homologs/mouse_genes_from_scRNAseq.txt", sep="\t", header = F)
bg <- as.vector(bg[,1])

# Load DCA + MAST DEGs
degs <- read.table("Output/DCA_MAST_DEGs/TRANSLATED_Summary_sig_adj_MAST_DEGs_all_clusters_res=0.6_dims=32.txt", sep="\t", header = T)
degs <- as.matrix(degs)
mode(degs) <- "character"
# CALCULATE: Enriched pathways for all cluster specific DEGs
source(paste(fp, "/pathway_enrichment.R", sep=""))
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  # pathway enrichment
  sig_pathw <- pathway_enrichment_function(pathways = kegg_pathways, degs = degs[,i], bg = bg, cores = 16, return_p_all_pathways = T, return_genes_sig_pathways = F)
  if(i == 1){
    pathways <- sig_pathw
  } else {
    pathways <- rbind(pathways, sig_pathw)
  }
}
rm(sig_pathw,i)
# CALCULATE: Overlapping pathways
pathways <- pathways[pathways[,2]<0.05,]
pathways <- table(pathways[,1])
write.table(cbind(names(pathways), pathways), file = "Output/DCA_MAST_DEGs_predictions/KEGG_pathways_enriched_by_multiple_cell_types.txt", sep="\t",
            col.names = c("Pathway", "enriched_in_x_cell_types"), row.names = F)
pathways <- pathways[pathways > 1]
registerDoParallel(cores = 16)
pathways <- foreach(i = c(1:(max(pathways)-1)), .combine = cbind) %dopar% {
  pathways <- names(pathways[pathways == i+1])
  temp <- kegg_pathways[,pathways]
  temp <- as.vector(temp)
  temp <- unique(temp[!is.na(temp)])
  
  temp <- c(temp, rep(NA, times = 8000-length(temp)))
  return(temp)
}
pathways <- pathways[rowSums(!is.na(pathways))>0,]
colnames(pathways) <- paste("Pathways_enriched_by_",2:(ncol(pathways)+1), "_cell_types",sep="")
degs <- pathways
rm(pathways, kegg_pathways)

# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/predefined_RA_drugs.txt", sep="\t", header = F)[,1])

# Set up output directories
dir.create("Output/DCA_MAST_DEGs_predictions/literature_PPI")
dir.create("Output/DCA_MAST_DEGs_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/DCA_MAST_DEGs_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])

# Analysis
save_name <- "Enriched_pathways_overlapping_between_clusters_in_literature_PPI_1000"
print(save_name)
rec <- F

for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes that can be mapped to PPI")
  }
  rm(out)
}



# start drug prediction HuRI PPI
##################################################################################################################
print("START HuRI PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/DCA_MAST_DEGs_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/HuRI_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Analysis
save_name <- "Enriched_pathways_overlapping_between_clusters_in_HuRI_PPI_1000"
print(save_name)
rec <- F

for(i in 1:ncol(degs)){
  print(paste("CALCULATING: DEGs ", colnames(degs)[i], sep=""))
  out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[,i]), disease_module_lcc = F, out.dir = out.dir,
                                           disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
  #rec <- T
  if(!(length(out)==1) & !all(is.na(out))){
    out <- precision_and_recall_from_z_score(predefined_list, out)
    print(out)
    write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    print("No disease genes that can be mapped to PPI")
  }
  rm(out)
}


# remove everything when done
rm(list = ls())
