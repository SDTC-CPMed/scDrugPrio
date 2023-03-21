#
#
# Make plots establishing the value of rank aggregation parameters. - CD
########################################################################################################

library(doParallel)
library(ggplot2)
library(reshape2)
fp <- getwd()
setwd(paste(fp, "/..", sep=""))

##############################################################################################################################################
# CD
##############################################################################################################################################

CD_drugs <- as.matrix(read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", stringsAsFactors = F, header = T, quote = ""))
dir <- "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/"
lf <- list.files(path = dir)
lf <- lf[grepl(pattern = "INDIVIDUAL_DRUGS_prediction_based_on_CD_DCA_MAST_DEGs_literature_PPI_Cluster_", x = lf)]

# Read in silico drug efficacy output
out <- foreach(i = (1:length(lf)), .combine = "rbind") %do% {
  infile <- as.matrix(read.table(file = paste(dir, lf[i], sep=""), sep="\t",header = T, stringsAsFactors = F))
  temp <- unlist(strsplit(strsplit(lf[i], split = "_")[[1]][13], split = ".txt")) # get cluster number
  #infile <- cbind(infile, paste("Cluster_",temp, sep=""))
  infile <- cbind(infile, temp)
  return(infile)
}
colnames(out)[9] <- "Cluster"
cl <- unique(out[,9]) # clusters

fc_criterium_file <- read.table(file = paste(dir, "COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt",sep=""),sep="\t",stringsAsFactors = F,header = T, quote = "", fill = T)
z_d_fc_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  if(any(as.numeric(fc_criterium_file[,17])==cl[i])){
    temp <- fc_criterium_file[as.numeric(fc_criterium_file[,17])==cl[i],]
    temp <- temp[as.numeric(temp[,10]) > 0,1] # choose only drugs names that counteract DEGs
    if(length(temp)>0){
      
    } else {
      temp <- NA
    }
  } else {
    temp <- NA
  }
  return(cbind(temp, cl[i]))
}

colnames(z_d_fc_criteria) <- c("Drug", "Cluster")
z_d_fc_criteria <- table(z_d_fc_criteria[!is.na(z_d_fc_criteria[,1]),1])
temp <- names(z_d_fc_criteria)
temp <- temp %in% CD_drugs[,2]
names(z_d_fc_criteria) <- temp

bar_plot_data <- foreach(i = c(1:max(z_d_fc_criteria)), .combine = "rbind") %do% {
  temp <- z_d_fc_criteria[z_d_fc_criteria == i]
  return(c(sum(names(temp)==T), sum(names(temp)==F), i))
}
colnames(bar_plot_data) <- c("n_CD", "n_Other","n_Reoccuring")

#######################################
# MAKE BAR PLOT showing that reoccuring candidates are better than candidates occuring in only one cluster
#######################################
rownames(bar_plot_data) <- bar_plot_data[,3]
df <- melt(bar_plot_data[,c(1:2)])

p <- ggplot(df, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
  xlab("times predicted") +
  ylab("n drugs") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(colour = "black", size = 20),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = paste(dir, "PRECISION_for_known_CD_drugs_among_reoccuring_candidates.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")


####################################################################

#######################################
# MAKE PLOTS evaluating mean dc as aggregation parameter 
#######################################

# get drug candidates name and DrugBank ID 
fc_criterium_file <- read.table(file = paste(dir, "COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt",sep=""),sep="\t",stringsAsFactors = F,header = T, quote = "", fill = T)
z_d_fc_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  if(any(as.numeric(fc_criterium_file[,17])==cl[i])){
    temp <- fc_criterium_file[as.numeric(fc_criterium_file[,17])==cl[i],]
    temp <- temp[as.numeric(temp[,10]) > 0,] # choose only drugs names that counteract DEGs
    if(nrow(temp) == 0){
      temp <- rep(NA, times = ncol(fc_criterium_file))
    }
  } else {
    temp <- rep(NA, times = ncol(fc_criterium_file))
  }
  return(temp)
}
z_d_fc_criteria <- z_d_fc_criteria[!is.na(z_d_fc_criteria[,1]),] # all candidates
unique_drugs <- unique(z_d_fc_criteria[,1])
# DB data
drugbank <- as.matrix(read.table(file ="Input/DrugBank/all_drug_targets_drug_bank.txt",sep="\t", header = T, stringsAsFactors = F))
unique_drug_target_combs <- as.matrix(read.table(file ="Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt",sep="\t", header = T, stringsAsFactors = F))
unique_drugs <- unique(drugbank[drugbank[,2]%in%unique_drugs, 1:2])
unique_drugs <- cbind(unique_drugs, unique_drug_target_combs[match(unique_drugs[,1],unique_drug_target_combs[,2]),1])

# get dc from raw files 
file_pattern <- "drug-disease_closest_distances_vs_random_bin_adjusted__CD_DCA_MAST_DEGs_literature_PPI_Cluster_"
file_path <- "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/"
lf <- list.files(path = file_path, pattern = file_pattern)
lf1 <- lf[!grepl(lf, pattern = "_top_")]
lf2 <- lf[grepl(lf, pattern = "_top_1800")]
lf <- unique(c(lf1,lf2))
for(i in 1:length(cl)){
  if(any(grepl(lf,pattern = paste("Cluster_",cl[i], "_top_1800",sep="")))){
    lf <- lf[!grepl(lf, pattern = paste("Cluster_",cl[i],".txt",sep=""))]
  }
}

mean_dist <- matrix(NA, nrow = nrow(unique_drugs), ncol = length(cl))
rownames(mean_dist) <- unique_drugs[,2]
colnames(mean_dist) <- cl
for(i in 1:length(cl)){ # once for every cluster/column
  temp <- c(lf[grepl(lf, pattern = paste("Cluster_",cl[i], ".txt", sep=""))], lf[grepl(lf, pattern = paste("Cluster_",cl[i], "_top_", sep=""))])
  pred_mat <- as.matrix(read.table(file = paste(file_path, temp,sep=""), sep="\t", header = T, stringsAsFactors = F))  
  pred_mat <- pred_mat[match(unique_drugs[,3], pred_mat[,1]),]
  mean_dist[,i] <- as.numeric(pred_mat[,3])
}
colnames(mean_dist) <- paste("mean_dc_", colnames(mean_dist),sep="")
mean_dist <- cbind(mean_dist, rowMeans(mean_dist))
colnames(mean_dist)[ncol(mean_dist)] <- "mean_of_mean_dc"
mean_dist <- cbind(mean_dist, rownames(mean_dist)%in%CD_drugs[,2])
colnames(mean_dist)[ncol(mean_dist)] <- "CD_drug"


# get sum of eigenvector centrality
final <- as.matrix(read.table(file = "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_CD.txt", sep="\t", header = T, stringsAsFactors = F)[,c(1,63)])
mean_dist <- cbind(mean_dist, final[match(rownames(mean_dist), final[,1]),2])

# Histogram showing distribution of average dc (mean of mean dc)
hist(as.numeric(mean_dist[,29]))

df <- data.frame(mean_dc = as.numeric(mean_dist[,29]), 
                 CD_drug = as.logical(as.numeric(mean_dist[,30])),
                 eigenvector_sum= as.numeric(mean_dist[,31]))

# boxplot comparing dc between RA and non RA --> RA has higher dc
p <- ggplot(df, aes(y=mean_dc, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_DC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

# dotplot comparing dc depending on sum of eigenvector centrality
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_dc, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_DC_IN_RELATION_TO_EIGENVECTOR_SUM_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")


####################################################################

#######################################
# MAKE PLOTS evaluating n or % of counteracting drug targets as aggregation parameter 
#######################################

fc_criterium_file <- read.table(file = paste(dir, "COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt",sep=""),sep="\t",stringsAsFactors = F,header = T, quote = "", fill = T)
z_d_fc_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  if(any(as.numeric(fc_criterium_file[,17])==cl[i])){
    temp <- fc_criterium_file[as.numeric(fc_criterium_file[,17])==cl[i],]
    temp <- temp[as.numeric(temp[,10]) > 0,] # choose only drugs names that counteract DEGs
    if(nrow(temp) == 0){
      temp <- rep(NA, times = ncol(fc_criterium_file))
    }
  } else {
    temp <- rep(NA, times = ncol(fc_criterium_file))
  }
  return(temp)
}
z_d_fc_criteria <- z_d_fc_criteria[!is.na(z_d_fc_criteria[,1]),] # all candidates

unique_drugs <- cbind(unique(z_d_fc_criteria[,1]), NA, NA, NA,NA,NA,NA)
for(i in 1:nrow(unique_drugs)){
  if(sum(z_d_fc_criteria[,1]%in% unique_drugs[i,1])>1){
    temp <- as.matrix(z_d_fc_criteria[z_d_fc_criteria[,1]%in%unique_drugs[i,1],c(10,11,9,12:14)])
    mode(temp) <- "numeric"
    temp <- colMeans(temp)
  } else {
    temp <- z_d_fc_criteria[z_d_fc_criteria[,1] == unique_drugs[i,1],c(10,11,9,12:14)] #12:14]
  }
  unique_drugs[i,2:7] <- as.numeric(temp)
}
colnames(unique_drugs)[2:7] <- c("mean_n_counteracting", "mean_n_mimmicking", "mean_n_outside",
                                 "mean_%_counteracting", "mean_%_mimmicking", "mean_%_outside")
mean_dist <- cbind(mean_dist, unique_drugs[match(rownames(mean_dist),unique_drugs[,1]),2:7])

df <- data.frame(mean_per_counteracting = as.numeric(mean_dist[,35]), 
                 mean_per_mimmick = as.numeric(mean_dist[,36]), 
                 mean_per_outside = as.numeric(mean_dist[,37]),
                 mean_n_counteracting = as.numeric(mean_dist[,32]), 
                 mean_n_mimmick = as.numeric(mean_dist[,33]), 
                 mean_n_outside = as.numeric(mean_dist[,34]),
                 CD_drug = as.logical(as.numeric(mean_dist[,30])),
                 eigenvector_sum= as.numeric(gsub(",", ".", mean_dist[,31])))

# mean percent counteracting
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_per_counteracting, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_PERCENT_COUNTERACTING_FC_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")
p <- ggplot(df, aes(y=mean_per_counteracting, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_PERCENT_COUNTERACTING_FC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

# mean percent mimmicking
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_per_mimmick, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_PERCENT_MIMMICKING_FC_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")
p <- ggplot(df, aes(y=mean_per_mimmick, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_PERCENT_MIMMICKING_FC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

# mean percent outside
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_per_outside, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_PERCENT_OUTSIDE_FC_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")
p <- ggplot(df, aes(y=mean_per_outside, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_PERCENT_OUTSIDE_FC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

# mean number counteracting
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_n_counteracting, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_N_COUNTERACTING_FC_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")
p <- ggplot(df, aes(y=mean_n_counteracting, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_N_COUNTERACTING_FC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

# mean number mimmicking
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_n_mimmick, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_N_MIMMICKING_FC_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")
p <- ggplot(df, aes(y=mean_n_mimmick, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_N_MIMMICKING_FC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

# mean number outside
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_n_outside, color=CD_drug)) + geom_point()
ggsave(filename = paste(dir, "MEAN_N_OUTSIDE_FC_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")
p <- ggplot(df, aes(y=mean_n_outside, fill=CD_drug)) + geom_boxplot()
ggsave(filename = paste(dir, "MEAN_N_OUTSIDE_FC_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

setwd(fp)
rm(list = ls())
dev.off()

####################################################################

fp <- getwd()
setwd(paste(fp, "/..",sep=""))

#######################################
# MAKE PLOTS evaluating whether the centrality of targeted DEGS can be used as a aggregation parameter
#######################################

# LOAD DATA

# PPI
ppi <- read.table(file = "Input/literature_PPI/ppi.txt",sep="\t",header = T, stringsAsFactors = F)
unique_proteins <- unique(c(unique(ppi[,1]), unique(ppi[,2])))
#DEGs
degs <- read.table(file = "Input/CD GSE134809/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt",sep="\t", header = T, stringsAsFactors = F)
for(i in 1:ncol(degs)){
  temp <- degs[,i]
  degs[,i] <- NA
  temp <- temp[temp%in% as.numeric(unique_proteins)]
  if(length(temp)>0){
    degs[1:length(temp),i] <- temp
  }
}
degs <- degs[,colSums(!is.na(degs))>0]
degs <- degs[1:1800,] # restrict number of DEGs based on top number in network model
# Clusters
cl <- unlist(strsplit(colnames(degs), split = "Cluster_"))
cl <- cl[cl != ""]
# Drug candidate info
dir <- "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/"
fc_criterium_file <- read.table(file = paste(dir, "COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt",sep=""),sep="\t",stringsAsFactors = F,header = T, quote = "",fill=T)
fc_criterium_file <- as.matrix(fc_criterium_file[as.numeric(fc_criterium_file[,10]) > 0,]) # choose only drugs names that counteract DEGs
# Drug info
same_drugs <- read.table(file ="Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt",sep = "\t", header = T)
drug_targets <- read.table(file ="Input/Drugs/drug_targets_unique_literature_ppi.txt",sep = "\t", header = T)
drugbank <- as.matrix(read.table(file ="Input/DrugBank/all_drug_targets_drug_bank.txt",sep="\t", header = T, stringsAsFactors = F))
# CD drugs
CD_drugs <- as.matrix(read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", stringsAsFactors = F, header = T, quote = ""))
# get sum of eigenvector centrality
final <- as.matrix(read.table(file = paste(dir, "FINAL_DRUG_RANKING_CD.txt",sep=""), sep="\t", header = T, stringsAsFactors = F)[,c(1,63,65,33)])


# GET LCCs for disease neighboorhood of each cell type (neighboorhood = DEGs + proteins within d = 1) 

# create disease neighborhood LCCs for each cluster
source("drug_prediction_R/in_silico_drug_efficacy_screening_algorithm.R")
ppi <- data.frame(node1 = as.character(ppi[,1]), node2 = as.character(ppi[,2]))
ppi_graph <- graph_from_data_frame(ppi, directed = F) # igraph package

LCCs <- foreach(i = c(1:ncol(degs)), .combine = "cbind") %do% {
  V(ppi_graph)$disease <- 0
  neighborhood <- degs[!is.na(degs[,i]),i]
  V(ppi_graph)[V(ppi_graph)$name %in% neighborhood]$disease <- 1
  all_genes <- single_connected_subgraph_extraction(ppi_graph, plot = F) # only PPI proteins that are connected to other PPI proteins
  if(ncol(all_genes)>1){
    if(nrow(all_genes)>1){
      all_genes <- all_genes[rowSums(!is.na(all_genes))>0,] # removes empty rows
      if(sum(colSums(is.na(all_genes))==0)>1){
        all_genes <- as.vector(all_genes[, colSums(is.na(all_genes))==0][,1]) # 'all_genes' will be used later on as it includes all relevant unique proteins in the LCC
      } else {
        all_genes <- as.vector(all_genes[, colSums(is.na(all_genes))==0]) # 'all_genes' will be used later on as it includes all relevant unique proteins in the LCC
      }
    } else {
      all_genes <- NA # non-connected -> not possible to calculate centrality on only one gene
    }
  } else {
    if(nrow(all_genes)>1){
      all_genes <- as.vector(all_genes) # all genes connected
    } else {
      all_genes <- NA # non-connected -> not possible to calculate centrality on only one gene
    }
  }
  
  print(paste("I = ",i, " LCC includes % of DEGS:",sep=""))
  print(100*sum(degs[!is.na(degs[,i]),i] %in% as.numeric(all_genes))/sum(degs[,i]%in%c(ppi[,1],ppi[,2])))
  
  return(c(all_genes, rep(NA, times = length(unique(c(ppi[,1],ppi[,2])))-length(all_genes))))
}
colnames(LCCs) <- cl
print("cluster number and size of LCC")
print(colSums(!is.na(LCCs)))
LCCs <- LCCs[,colSums(!is.na(LCCs))>0]

library(igraph)
library(CINNA)

#cl <- unique(as.numeric(fc_criterium_file[,17]))

cl <- as.numeric(as.character(colnames(LCCs)))
cl <- cl[cl%in%fc_criterium_file[,17]]

gm_mean = function(a){prod(a)^(1/length(a))}

mean_target_centrality <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% { # for every cluster
  #for(i in 1:length(cl)){ 
  print(i)
  # select candidates targeting this cluster
  if(sum(as.numeric(fc_criterium_file[,17]) %in% cl[i])>1){
    cand <- fc_criterium_file[as.numeric(fc_criterium_file[,17]) %in% cl[i],1:18]
  } else {
    cand <- matrix(fc_criterium_file[as.numeric(fc_criterium_file[,17]) %in% cl[i],1:18], nrow = 1)
  }
  
  # Calculate centrality for proteins in disease neighborhood of this cluster
  ###########################################################################
  # make graph of disease neighboorhood LCC
  g <- ppi
  g <- g[g[,1]%in%LCCs[,colnames(LCCs) %in% cl[i]],]
  g <- g[g[,2]%in%LCCs[,colnames(LCCs) %in% cl[i]],]
  g <- graph_from_data_frame(g, directed = F)
  
  # Calculate centralities
  centrality_matrix <- calculate_centralities(g, include = c("eigenvector centralities",
                                                             #"K-core Decomposition", 
                                                             "Kleinberg's hub centrality scores"))#, 
  #"Laplacian Centrality", 
  #"Leverage Centrality", 
  #"Group Centrality", 
  #"Local Bridging Centrality"))
  centrality_matrix <- matrix(unlist(centrality_matrix), ncol = length(V(g)$name), byrow = T)
  rownames(centrality_matrix) <- c("eigenvector centralities",
                                   #"K-core Decomposition", 
                                   "Kleinberg's hub centrality scores")#,
  #"Laplacian Centrality", 
  #"Leverage Centrality", 
  #"Group Centrality", 
  #"Local Bridging Centrality")
  colnames(centrality_matrix) <- V(g)$name
  
  # Select only centralities of DEGs of this cluster
  #centrality_matrix <- centrality_matrix[,colnames(centrality_matrix)%in% degs[,colnames(degs) == paste("Cluster_",cl[i],sep="")]]
  
  cent <- vector()
  for(d in 1:nrow(cand)){ #for every drug candidate
    temp <- unique(drugbank[drugbank[,2] %in% cand[d,1],1]) # DrugBankID
    temp <- same_drugs[same_drugs[,2]%in%temp,1] # get unique drug target combo
    temp <- temp[temp%in%colnames(drug_targets)]
    if(length(temp)!=1){
      print(paste("PROBLEM: i = ",i," d = ", d, " temp = ", temp, sep=""))
    }
    temp <- as.vector(drug_targets[,colnames(drug_targets)%in% temp]) # drug targets in Entrez ID
    temp <- temp[!is.na(temp)]
    if(any(as.character(temp) %in% colnames(centrality_matrix))) {
      if(sum(as.character(temp) %in% colnames(centrality_matrix))>1){
        temp <- centrality_matrix[,colnames(centrality_matrix)%in% temp]
      } else {
        temp <- matrix(centrality_matrix[,colnames(centrality_matrix)%in% temp], ncol = 1)
      }
    } else {
      temp <- matrix(0, nrow = nrow(centrality_matrix),ncol = 1)
    }
    # Geometric mean of Eigenvector centrality of drug targets
    cent[d] <- gm_mean(temp[1,])#sum(temp[3,])
  }
  cand <- cbind(cand, cent)
  colnames(cand) <- NULL
  return(cand)
}
colnames(mean_target_centrality) <- c(colnames(fc_criterium_file)[1:18],"Centrality_sum")

# Now aggregate into mean centrality sum for all candidates
unique_drugs <- cbind(unique(mean_target_centrality[,1]), NA)
for(i in 1:nrow(unique_drugs)){
  if(sum(mean_target_centrality[,1]%in% unique_drugs[i,1])>1){
    temp <- mean(as.numeric(mean_target_centrality[mean_target_centrality[,1]%in%unique_drugs[i,1],19]))
  } else {
    temp <- as.numeric(mean_target_centrality[mean_target_centrality[,1] == unique_drugs[i,1],19])
  }
  unique_drugs[i,2] <- as.numeric(temp)
}

unique_drugs <- cbind(unique_drugs, unique_drugs[,1]%in%CD_drugs[,2])
unique_drugs <- cbind(unique_drugs, final[match(unique_drugs[,1], final[,1]),2:4])
colnames(unique_drugs)[2:6] <- c("mean_target_centrality", "CD_drug", "eigenvector_sum", "lit_search", "mean_dc")
write.table(unique_drugs, file = paste(dir, "ranking_based_on_average_drug_target_centrality_within_cell_type_specific_LCCs_of_DEGs.txt",sep=""),sep="\t", col.names = T, row.names = F)


df <- data.frame(mean_drug_target_centrality = as.numeric(unique_drugs[,2]),
                 CD_drug = as.logical(unique_drugs[,3]),
                 eigenvector_sum= as.numeric(gsub(",", ".", unique_drugs[,4])),
                 name = as.character(unique_drugs[,1]),
                 lit_search = as.character(unique_drugs[,5]),
                 mean_dc = as.numeric(gsub(",", ".", unique_drugs[,6])))

print(paste("Mean_CD: ",mean(df[df[,2]==T,1]), median(df[df[,2]==T,1]), sep=" "))
print(paste("Mean_non_CD: ",mean(df[df[,2]==F,1]),median(df[df[,2]==F,1]),sep=" "))

# Ranking better with cell type averaged geometric mean drug target centrality compared to only ranking by cell type eigenvector centrality sum? 
set.seed(1)
rand_lit_distribution <- rand_distribution <- vector()
semi_rand_lit_distribution <- semi_rand_distribution <- vector()
for(i in 1:1000000){ # one-million bootstrap repetitions for reproducibility
  df <- df[sample(c(1:nrow(df))),]
  rand_distribution[i] <- mean(which(df[,2]==T))
  rand_lit_distribution[i] <- mean(which(df[,5]=="Yes"))
  semi_rand_distribution[i] <- mean(which(df[order(df[,3], decreasing = T),2]==T))
  semi_rand_lit_distribution[i] <- mean(which(df[order(df[,3], decreasing = T),5]=="Yes"))
}
#hist(rand_distribution)
#hist(rand_lit_distribution)
#hist(semi_rand_distribution)
#hist(semi_rand_lit_distribution)
print("One-sided P-value showing that known RA drugs rank higher when using cell type eigenvector centrality rankingcompared to random:")
print(t.test(x = semi_rand_distribution, y = rand_distribution,alternative = "less"))
print("One-sided P-value showing that known RA drugs rank higher when using cell type eigenvector centrality rankingcompared to random:")
print(t.test(x = semi_rand_lit_distribution, y = rand_lit_distribution,alternative = "less"))

print("One-sided P-value showing that known RA drugs rank higher with drug target centrality compared to random:")
print(pnorm((mean(which(df[order(df[,1], decreasing = T),2]==T)) - mean(rand_distribution))/sd(rand_distribution)))
print("One-sided P-value showing that known RA drugs rank higher with drug target centrality compared to random:")
print(pnorm((mean(which(df[order(df[,1], decreasing = T),5]=="Yes")) - mean(rand_lit_distribution))/sd(rand_lit_distribution)))

print("One-sided P-value showing that known RA drugs rank higher with mean closest distance")
print(pnorm((mean(which(df[order(-df[,6], decreasing = T),2]==T)) - mean(rand_distribution))/sd(rand_distribution)))
print("One-sided P-value showing that known RA drugs rank higher with mean closest distance")
print(pnorm((mean(which(df[order(-df[,6], decreasing = T),5]=="Yes")) - mean(rand_lit_distribution))/sd(rand_lit_distribution)))

print("One-sided P-value showing that known RA drugs rank higher when combining cell type eigenvector centrality ranking with drug target centrality ranking instead of only cell type eigenvector centrality sum ranking:")
print(pnorm((mean(which(df[order(df[,3], df[,1], decreasing = T),2]==T)) - mean(semi_rand_distribution))/sd(semi_rand_distribution)))
print("One-sided P-value showing that known RA drugs rank higher with eigenvector centrality & mean target centrality")
print(pnorm((mean(which(df[order(df[,3], df[,1], decreasing = T),5]=="Yes")) - mean(semi_rand_lit_distribution))/sd(semi_rand_lit_distribution)))

print("One-sided P-value showing that known RA drugs rank higher with eigenvector centrality & mean closest distance")
print(pnorm((mean(which(df[order(df[,3], -df[,6], decreasing = T),2]==T)) - mean(semi_rand_distribution))/sd(semi_rand_distribution)))
print("One-sided P-value showing that known RA drugs rank higher with eigenvector centrality & mean closest distance")
print(pnorm((mean(which(df[order(df[,3], -df[,6], decreasing = T),5]=="Yes")) - mean(semi_rand_lit_distribution))/sd(semi_rand_lit_distribution)))

print("One-sided P-value showing that known RA drugs rank higher when combining cell type eigenvector centrality ranking with drug target centrality ranking instead of only cell type eigenvector centrality sum ranking:")
print(pnorm((mean(which(df[order(df[,3], df[,1], decreasing = T),2]==T)) - mean(rand_distribution))/sd(rand_distribution)))
print("One-sided P-value showing that known RA drugs rank higher with eigenvector centrality & mean target centrality")
print(pnorm((mean(which(df[order(df[,3], df[,1], decreasing = T),5]=="Yes")) - mean(rand_lit_distribution))/sd(rand_lit_distribution)))

print("One-sided P-value showing that known RA drugs rank higher with eigenvector centrality & mean closest distance")
print(pnorm((mean(which(df[order(df[,3], -df[,6], decreasing = T),2]==T)) - mean(rand_distribution))/sd(rand_distribution)))
print("One-sided P-value showing that known RA drugs rank higher with eigenvector centrality & mean closest distance")
print(pnorm((mean(which(df[order(df[,3], -df[,6], decreasing = T),5]=="Yes")) - mean(rand_lit_distribution))/sd(rand_lit_distribution)))


df <- df[order(df[,3],df[,1],decreasing = T),]
write.table(df, file = paste(dir,"ranking_by_cell_type_and_gm_mean_of_drug_target_centrality.txt",sep=""),sep="\t", col.names = T)

p <- ggplot(df, aes(x = 1:nrow(df), y = eigenvector_sum, color = CD_drug)) +
  geom_point(size = 1) + 
  geom_line(data = data.frame(y = c(0,6),x = c(50, 50)), aes(y=y, x=x),color = "black") +
  geom_line(data = data.frame(y = c(0,6),x = c(100, 100)), aes(y=y, x=x),color = "black")
ggsave(filename = paste(dir, "MEAN_DRUG_TARGET_CENTRALITY_and_EIGENCETOR_CENTRALITY_SCORE_SUM_for_ranking.pdf", sep=""), plot = p, device = "pdf", width = 200, height = 80, units = "mm")


# plots
p <- ggplot(df, aes(x = eigenvector_sum, y=mean_drug_target_centrality, color=CD_drug)) + 
  geom_point(size = 1+1*df$CD_drug) 
ggsave(filename = paste(dir, "MEAN_DRUG_TARGET_CENTRALITY_and_EIGENCETOR_CENTRALITY_SCORE_SUM_as_aggregation_parameter.pdf", sep=""), plot = p, device = "pdf", width = 100, height = 110, units = "mm")

p <- ggplot(df, aes(y=mean_drug_target_centrality, x=CD_drug, fill = CD_drug)) + 
  #geom_boxplot(outlier.shape = NA, color = "black", lwd = 0.25) + 
  ylim(0,1) + scale_y_continuous(trans = "log10") +
  geom_dotplot(color = NA, binaxis = "y", stackdir = "center", dotsize = 1, binwidth = 0.07, shape = 1) +
  #geom_jitter(position = position_jitter(0.2), size = 0.5) + #geom_point(size = 0.5) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(colour = "black", size = 20),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = paste(dir, "MEAN_DRUG_TARGET_CENTRALITY_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 160, height = 110, units = "mm")

# make categories for cell type eigenvector centrality sum 
cats <- seq(0,6,by = 0.5)
temp <- vector()
for(i in 2:length(cats)){
  if(any(df[,3]<cats[i])){
    temp <- rbind(temp, cbind(df[df[,3]<cats[i],1:4], paste(cats[i-1]," - ", cats[i],sep="")))
    df <- df[df[,3]>=cats[i],]
  }
}
df <- temp
colnames(df)[5] <- "cats"

p <- ggplot(df, aes(y=mean_drug_target_centrality, x=cats, fill = CD_drug)) + 
  #geom_boxplot(outlier.shape = NA, color = "black", lwd = 0.25) + 
  ylim(0,1) + scale_y_continuous(trans = "log10") +
  geom_dotplot(color = NA, binaxis = "y", stackdir = "center", dotsize = 1, binwidth = 0.03) +
  #geom_jitter(position = position_jitter(0.2), size = 0.5) + #geom_point(size = 0.5) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 20),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(colour = "black", size = 20),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = paste(dir, "MEAN_DRUG_TARGET_CENTRALITY_&_CELL_TYPE_CENTRALITY_AS_RANK_AGGREGATION_CRITERIA_for_known_CD_drugs.pdf", sep=""),
       plot = p, device = "pdf", width = 150, height = 150, units = "mm")

#rm(list = ls())
