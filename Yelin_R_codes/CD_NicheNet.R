#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(tidyverse)
library(stringr)
library(readxl)
library(plyr)

#########################
#1. Single dataset####
#########################
#change it when changing source
datasource <- "CD"

inputdir = paste0(getwd(),"data/CD NicheNet")
outputdir = paste(getwd(),"/Results/CD",sep="")
if (dir.exists(outputdir)==FALSE){
  dir.create(outputdir, recursive = T)
  print('outputdir was created')
  }
  
#### Step 0:  NicheNet’s ligand-target prior model  ####
#This model denotes the prior potential that a particular ligand might regulate the expression of a specific target gene.
ligand_target_matrix <- readRDS("ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]  # target genes in rows, ligands in columns

all_expressed_genes <- read.table(file = "./CD NicheNet/TRANSLATED_HGNC_SYMBOLS_background_genes_all_clusters_DCA_adj_data.txt", 
                                  header= T, sep="\t",stringsAsFactor = T) 

all_expressed_genes[1:5,1:5]
colnames(all_expressed_genes)
dim(all_expressed_genes)

#add prefix to each cluster
colnames(all_expressed_genes) <- paste(datasource, colnames(all_expressed_genes), sep = "_")
colnames(all_expressed_genes)
dim(all_expressed_genes)

all_expressed_genes_CD <- all_expressed_genes

DEGs_human <- read.table(file="./CD NicheNet/TRANSLATED_HGNC_SYMBOL_sig_MAST_DEGs_log(1.5)_all_clusters.txt",
                         header=T, sep="\t", stringsAsFactors = T)

DEGs_human[1:5,1:5]
colnames(DEGs_human) <- paste(datasource, colnames(DEGs_human), sep = "_")
colnames(DEGs_human)
dim(DEGs_human)

DEGs_human_CD <- DEGs_human


#### looping for ligands activity for single dataset ####
all_ligand_activity_1 <- data.frame(matrix(ncol=6,nrow=0))
colnames(all_ligand_activity_1) <- c("test_ligand","auroc","aupr","pearson", "Sender", "Target")

for(i in 1:dim(DEGs_human)[2]){ 
  
  for(j in 1:dim(DEGs_human)[2]){
    sender_cluster = colnames(DEGs_human)[i]
    receiver_cluster = colnames(DEGs_human)[j]   
    
    expressed_genes_sender = as.vector(DEGs_human[,i])
    expressed_genes_receiver = as.vector(all_expressed_genes[,receiver_cluster])
    
    # gene set of interest - setDEG as gene set of interest ####
    geneset_oi <- as.vector(DEGs_human[, j] %>% .[. %in% rownames(ligand_target_matrix)] )
    
    # set background expressed genes  ####
    background_expressed_genes = as.vector(expressed_genes_receiver %>%
                                             .[. %in% rownames(ligand_target_matrix)])
    
    # Define a set of potential ligands ####
    # Putative ligand-receptor links were gathered from NicheNet’s ligand-receptor data sources.
    lr_network = readRDS("lr_network.rds")
    
    # If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions 
    # and only keep ligand-receptor interactions that are described in curated databases. 
    # To do this: uncomment following line of code:
    # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
    
    ligands = lr_network %>% pull(from) %>% unique()        
    expressed_ligands = intersect(ligands,expressed_genes_sender)  
    
    if (length(expressed_ligands)==0) {
      next
    }
    
    receptors = lr_network %>% pull(to) %>% unique()
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    if (length(expressed_receptors)==0) {
      next
    }
    
    lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)      #选取ligand 和receptor都要有表达的
    head(lr_network_expressed)
    
    potential_ligands = lr_network_expressed %>% pull(from) %>% unique()  #只有lr_network在当前数据中有表达的ligands才考虑作为potentially active ligands for the NicheNet analysis
    head(potential_ligands)
    
    ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                  background_expressed_genes = background_expressed_genes, 
                                                  ligand_target_matrix = ligand_target_matrix, 
                                                  potential_ligands = potential_ligands)
    
    ligand_activities %>% arrange(-pearson)
    ligand_activities <- cbind(ligand_activities,Sender=sender_cluster, Target=receiver_cluster)
    all_ligand_activity_1 <- rbind(all_ligand_activity_1,ligand_activities)
  }
}

View(all_ligand_activity_1)

#### replace a whole column with parts of their characters ####
a<-all_ligand_activity_1
all_ligand_activity_1<-a
Sender <- strsplit(all_ligand_activity_1$Sender,"_")
all_ligand_activity_1$Sender <- lapply(Sender, function(Sender) Sender [[4]])
all_ligand_activity_1$Sender <- paste( all_ligand_activity_1$Sender, datasource,sep = "_")
#all_ligand_activity_1$Sender <- as.numeric(unlist(all_ligand_activity_1$Sender)) # conver a list to numeric

Target <- strsplit(all_ligand_activity_1$Target,"_")
all_ligand_activity_1$Target <- lapply(Target, function(Target) Target [[4]])
all_ligand_activity_1$Target <- paste(all_ligand_activity_1$Target, datasource,sep = "_")
#all_ligand_activity_1$Target <- as.numeric(unlist(all_ligand_activity_1$Target))
 
write.table(all_ligand_activity_1, 
            file=paste(outputdir,"/all_ligand_activity_",datasource,".txt",sep=""),
            sep="\t", quote=F, col.names=T,row.names=F)

dim(table(all_ligand_activity_1$Sender,all_ligand_activity_1$Target))
 
cluster_intera_ye <- data.frame(table(all_ligand_activity_1$Sender,all_ligand_activity_1$Target))
colnames(cluster_intera_ye) = c("Sender", "Target", "No.Interactions")

write.table(cluster_intera_ye, 
            file=paste(outputdir,"/cluster_intera_",datasource,".txt",sep=""), 
            sep="\t", quote=F, col.names=T,row.names=F)
dim(cluster_intera_ye)
View(cluster_intera_ye)

#Centrality####
source("2.1 NicheNet_centrality_single_dataset.R")


