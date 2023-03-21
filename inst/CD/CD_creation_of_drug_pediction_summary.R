#
# Get drug prediction data for all cell types and create table for quickly checking whether
# additional filtering criteria is fulfilled
#
# SAMUEL SCH?FER
#
##################################################################################################
library(doParallel)
registerDoParallel(cores = 1)

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

# Functions
##################################################################################################
additional_filtering <- function(summary, disease_model, drug_bank_target_data){
  print("Applying additional filtering criteria")
  
  # save colnames
  c_names <- colnames(summary)[c(1,3:ncol(summary),2)]
  
  if(sum(as.numeric(summary[,7])<0.05) > 0){
    if(sum(as.numeric(summary[,7])<0.05) > 1){
      summary <- summary[as.numeric(summary[,7])<0.05,] # significantly closer than random drug to random disease
    } else {
      summary <- summary[as.numeric(summary[,7])<0.05,] # significantly closer than random drug to random disease
      summary <- matrix(summary, nrow = 1)
    }
  } else {
    print("NO DRUGS matching additional filtering criteria: P < 0.05")
    return(NA)
  }
  if(sum(as.numeric(summary[,3])<1) > 0){
    if(sum(as.numeric(summary[,3])<1) > 1){
      summary <- summary[as.numeric(summary[,3])<1,] # significantly closer than random drug to random disease
    } else {
      summary <- summary[as.numeric(summary[,3])<1,] # significantly closer than random drug to random disease
      summary <- matrix(summary, nrow = 1)
    }
  } else {
    print("NO DRUGS matching additional filtering criteria: d(c) < 1")
    return(NA)
  }
  print(paste("matching n drugs:", nrow(summary), sep=" "))
  # Effect of drug on disease model?
  out <- foreach(i = 1:nrow(summary), .combine = rbind) %dopar% {
    temp <- summary[i,1]
    if(sum(drug_bank_target_data[,2]%in%temp)>1){
      temp <- drug_bank_target_data[grepl(drug_bank_target_data[,2], pattern = paste("\\Q",temp, "\\E", sep="")),]
    } else {
      temp <- as.matrix(t(drug_bank_target_data[drug_bank_target_data[,2] == temp,]))
    }
    temp <- temp[!is.na(temp[,6]),]
    if(is.vector(temp)){
      temp <- matrix(temp, nrow=1)
    }
    if(nrow(temp)==0){
      return(c(NA, NA, NA))
    }
    t <- paste(as.character(temp[,6]), " (", as.character(temp[,7]),")", sep="")
    t <- paste(t, sep = ", ")
    t2 <- cbind(as.character(temp[,6]), NA)
    t2[,2] <- disease_model[match(t2[,1],disease_model[,1]),3]
    t2[as.numeric(t2[,2])>0,2] <- "up"
    t2[as.numeric(t2[,2])<0,2] <- "down"
    t3 <- sum(is.na(t2[,2]))
    t2[is.na(t2[,2]),2] <- "not in disease model"
    t2 <- paste(as.character(t2[,1]), " (", t2[,2], ")", sep="")
    
    if(length(t)>1){ # more than one drug target
      for(j in 2:length(t)){
        t[1] <- paste(t[1], ", ", t[j], sep="")
        t2[1] <- paste(t2[1], ", ", t2[j], sep="")
      }
      t <- t[1]
      t2 <- t2[1]
    }
    return(c(t3, t2, t))
  }
  
  c_names <- c(c_names, "n_drug_targets_outside_disease_model", "drug_targets_in_disease", "drug_target_modification")
  if(is.vector(out)){
    summary <- as.vector(summary[1,])
    summary <- c(summary[c(1,3:length(summary),2)], out)
    summary <- matrix(summary, nrow = 1)
  } else {
    summary <- cbind(summary[,c(1,3:ncol(summary),2)], out)
  }
  colnames(summary) <- c_names
  return(summary)
}

# load input file
load_drug_prediction <- function(dir_network_pred_file){
  # load network prediction file
  in_file <- as.matrix(read.table(file = dir_network_pred_file, header = T, sep = "\t"))
  
  # order by Z-score
  in_file <- in_file[order(as.numeric(in_file[,4]), decreasing = F),]
  '
  # exclude Z-score > 0
  in_file <- in_file[as.numeric(in_file[,4]) <= 0,]
  
  # translate drug target combinations to drug names
  for(i in c(1:nrow(in_file))){
    temp <- in_file[i,1]
    temp <- c(temp, same_drugs[same_drugs[,1]%in%temp,2], same_drugs[same_drugs[,2]%in%temp,1])
    temp <- unique(as.character(all_drugs[match(temp, all_drugs[,1]),2]))
    if(length(temp)>1){
      for(j in 2:length(temp)){
        temp[1] <- paste(as.character(temp[1]), as.character(temp[j]), sep=", ")
      }
    }
    in_file[i,1] <- temp[1]
    rm(temp)
  }
  colnames(in_file) <- c("Drug_names", "n_drug_targets", "mean_network_distance_d(c)", "z_score_d(c)","mean_random_d(c)", "SD_d(c)_random", "P_based_on_z")
  '
  return(in_file)
}

# translate to individual drugs
translate_to_individual_drugs <- function(file, individual_drugs, CD_drugs){
  temp <- cbind(file[,1], 1:nrow(file))
  colnames(temp) <- c(colnames(file)[1], "Position")
  temp <- merge(x = temp, y = individual_drugs, by.x = colnames(temp)[1], by.y = colnames(individual_drugs)[1])
  out <- file[as.numeric(temp[,2]),]
  out[,1] <- temp[,3] # set individual drug names
  # known CD drugs?
  temp <- matrix(out[,1] %in% CD_drugs[,1], ncol = 1)
  colnames(temp) <- "known_CD_drug"
  out <- cbind(out, temp)
  return(out)
}
##################################################################################################


# LOAD files
##################################################################################################
# all drugs DrugBank
db <- all_drugs <- as.matrix(read.table(file = "Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", header = T))
all_drugs <- all_drugs[!is.na(all_drugs[,6]),] # gene symbol available
all_drugs <- all_drugs[grepl(all_drugs[,8], pattern = "Human"),] # for human
all_drugs <- all_drugs[grepl(all_drugs[,3], pattern = "approved"),] # only approved drugs
# Load unique drug target combinations
unique_drug_target_comb <- as.matrix(read.table(file = "Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T, stringsAsFactors = F))
# Prepare same_drugs
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = F)) 
unique_drugs <- unique(c(same_drugs[,1], same_drugs[,2]))
individual_drugs <- foreach(i = c(1:length(unique_drugs)), .combine = "rbind") %do% {
  temp <- c(same_drugs[same_drugs[,1]%in% unique_drugs[i],2],same_drugs[same_drugs[,2]%in% unique_drugs[i],1])
  temp <- unique(temp)
  return(cbind(unique_drugs[i], temp))
}
colnames(individual_drugs) <- c("Drug","Drug2")
individual_drugs <- individual_drugs[individual_drugs[,1]%in%colnames(unique_drug_target_comb),]
# known CD drugs
CD_drugs <- as.matrix(read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep="\t", header = T, quote = ""))
CD_drugs <- CD_drugs[CD_drugs[,1]!="",]
# Define biosimilars
biosimilar <- unique(c(CD_drugs[,1],same_drugs[same_drugs[,1]%in% CD_drugs[,1],2], same_drugs[same_drugs[,2]%in% CD_drugs[,1],1]))
# Define drug categories
unique_drugs <- cbind(all_drugs[,1:2], NA)
unique_drugs <- unique_drugs[duplicated(unique_drugs[,1]) != T, ]
unique_drugs[unique_drugs[,1]%in%biosimilar,3] <- "Biosimilar"
unique_drugs[unique_drugs[,1]%in%CD_drugs[,1],3] <- "CD_drug"

# Analysis
##################################################################################################

# CD DCA MAST DEGs
###################################

lf <- list.files(path = "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI")
lf <- lf[grepl(lf, pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__CD_DCA_MAST_DEGs_literature_PPI_Cluster_")]

dir.create(path = "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY")
dir.out <- "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/"
# create summary for all clusters
out <- vector()
for(i in 1:length(lf)){
  
  print(lf[i])
  
  dir <- paste("Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/", lf[i], sep="")
  # load prediction
  prediction <- load_drug_prediction(dir)
  # split prediction into individual drugs
  prediction <- translate_to_individual_drugs(prediction, individual_drugs, CD_drugs = CD_drugs)
  # translate DrugBank-ID to drug name
  drug_name <- cbind(db[,1:2], paste(db[,1], "_",db[,2],sep=""))
  drug_name <- drug_name[!duplicated(drug_name[,3]),1:2]
  prediction[,1] <- drug_name[match(prediction[,1],drug_name[,1]), 2]
  #save
  temp2 <- strsplit(lf[i], split = "\\__")[[1]][2]
  write.table(prediction, file = paste(dir.out, "INDIVIDUAL_DRUGS_prediction_based_on_", temp2, sep=""), sep="\t", col.names = T, row.names = F)
  
  # which cluster?
  temp2 <- strsplit(temp2, split="CD_DCA_MAST_DEGs_literature_PPI_Cluster_")[[1]][2]
  temp2 <- strsplit(temp2, split = "\\.txt")[[1]][1] # can output for example only cluster number "0", "1", or cluster description "0_top_10", "0_top25"
  
  # load translation matrix
  transl <- read.table(file = "Input/HGNC translation matrix 201108/transl.txt", sep="\t", header = T)
  # Get disease model data underlying drug prediction
  if(length(unlist(strsplit(temp2, split = "\\_")))>1){ # check whether all cluster genes were used or if only top X genes were used
    top_DEGs <- unlist(strsplit(temp2, split = "\\_"))[3]
    cluster <- unlist(strsplit(temp2, split = "\\_"))[1]
    disease_model <- read.table(file = paste("Input/CD/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/Cluster_", cluster, "_res=0.8_dims=32_k=15.txt", sep=""), sep="\t", header = T)
    #translate disease model to human gene symbols
    disease_model[,1] <- as.character(transl[match(disease_model[,1], transl[,7]),6]) # translate to NCBI gene symbols
    disease_model <- disease_model[!is.na(disease_model[,1]),]
    # select top DEGs
    disease_model <- disease_model[1:top_DEGs,]
    # translate only top DEGs to NCBI symbols
    disease_model[,1] <- as.character(transl[match(disease_model[,1], transl[,6]),2]) # translate to NCBI gene symbols
    rm(top_DEGs, cluster)
  } else { # all DEGs were used
    disease_model <- read.table(file = paste("Input/CD GSE134809/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/Cluster_", temp2, "_res=0.8_dims=32_k=15.txt", sep=""), sep="\t", header = T)
    #translate disease model to human gene symbols
   disease_model[,1] <- as.character(transl[match(disease_model[,1], transl[,7]),2]) # translate to NCBI gene symbols
  }
  disease_model <- disease_model[!is.na(disease_model[,1]),] # remove invalid mappings
  
  prediction <- additional_filtering(summary = prediction, disease_model = disease_model, drug_bank_target_data = db)
  
  if(!is.null(nrow(prediction))){
    if(nrow(prediction)>1){
      prediction <- prediction[order(as.numeric(prediction[,6]), decreasing = F),]
    }
    # only known CD drugs
    if(sum(prediction[,7] == TRUE)>0){
      if(sum(prediction[,7] == TRUE)==1){
        c_names <- c(colnames(prediction), "Cluster")
        temp3 <- prediction[prediction[,7] == TRUE,]
        temp3 <- matrix(temp3, nrow = 1)
        out <- rbind(out, cbind(temp3, rep(temp2, times = nrow(temp3))))
        colnames(out) <- c_names
        rm(c_names)
      } else {
        temp3 <- prediction[prediction[,7] == TRUE,]
        out <- rbind(out, cbind(temp3, rep(temp2, times = nrow(temp3))))
      }
    } else {
      out <- rbind(out, c(rep(NA, times = 11),temp2)) # times = 11 is resembling ncol(temp3)
    }
    # Specify what kind of known RA drug:
    prediction[prediction[,1]%in% unique_drugs[unique_drugs[,3]=="Biosimilar?",2],7] <- "Biosimilar?"
    prediction[prediction[,1]%in% unique_drugs[unique_drugs[,3] == "CD_drug",2],7] <- "CD_drug"
  } else {
    out <- rbind(out, c(rep(NA, times = 11),temp2)) # times = 11 is resembling ncol(temp3)
  }
  rm(temp3)
  
  # save drugs filtered with additional criteria
  write.table(prediction, file = paste(dir.out, "EVALUATION_predicted_drugs_MAST_DEGs_log(1.5)_threshold_cluster_", temp2,
                                 "_res=0.8_dims=32.txt", sep=""), sep="\t", col.names = T, row.names = F)
  
}
write.table(out, file = paste(dir.out, "EVALUATION_known_CD_drugs_in_all_clusters_log(1.5).txt", sep=""), sep="\t",
            col.names = c(colnames(prediction), "Cluster"), row.names = F)

setwd(fp)
rm(list = ls())

# Create precision-recall curves
###################################
fp <- getwd()
setwd(paste(fp, "/../Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY", sep=""))
library(viridis)

lf <- list.files()
lf <- lf[grepl(lf, pattern = "INDIVIDUAL_DRUGS_prediction_based_on_CD_DCA_MAST_DEGs_literature_PPI_Cluster_")]

# cut-offs used in attempt on threshold determination
#top_DEGs <- c(10,25,50,75,100,150,200,300,400,500,1000,1500,2000,2500,3000,3500)

# load DEGs
degs <- as.matrix(read.table(file = paste(fp, "/../Input/CD/DCA_MAST_DEGs_logFC_1.5_CD_ileal_cells/TRANSLATED_ENTREZ_ID_sig_MAST_DEGs_log(1.5)_all_clusters.txt",sep=""), sep="\t", header = T, stringsAsFactors = F))

out <- foreach(i = 1:length(lf), .combine = "rbind") %do% {
  print(lf[i])
  # LOAD infile
  in_file <- as.matrix(read.table(file = lf[i], sep="\t", header = T, stringsAsFactors = F))
  # number of CD drugs
  n_drugs <- sum(in_file[,8]==T)
  # select only drugs at z(c) < -1.64 // P < 0.05
  in_file <- in_file[as.numeric(in_file[,7]) < 0.05,]
  # calculate recall & precision for individual drugs
  recall <- sum(in_file[,8]==T) / n_drugs
  precision <- sum(in_file[,8]==T) / nrow(in_file)
  
  # which cluster? how many DEGs in calculation?
  cluster <- strsplit(lf[i], split = "INDIVIDUAL_DRUGS_prediction_based_on_CD_DCA_MAST_DEGs_literature_PPI_Cluster_")[[1]][2]
  cluster <- unlist(strsplit(cluster, split = "\\.txt"))
  if(grepl(cluster, pattern = "_top_")){
    n_degs <- strsplit(cluster, split = "_")[[1]][3]
    cluster <- strsplit(cluster, split = "_")[[1]][1]
  } else { # all DEGs used for drug prediction
    n_degs <- sum(!is.na(degs[,colnames(degs) == paste("Cluster_",cluster,sep="")]))
  }
  
  # save data
  return(c(cluster, n_degs, precision, recall))
}
colnames(out) <- c("Cluster", "n_degs", "precision", "recall")
out <- out[order(as.numeric(out[,2]), decreasing = F),]

# prepare for plotting
cluster <- unique(as.character(out[,1]))
colors <- viridis(n = length(cluster))

# Recall plot
pdf(file = "RECALL_at_thresholds_CD.pdf", paper = "a4r")
plot(x = c(min(as.numeric(out[,2])), max(as.numeric(out[,2]))),
     y = c(0,(max(as.numeric(out[,4]))+0.05)),type = "n", 
     xlab = "n top DEGs", ylab = "Recall", main = "RECALL for CD data")
for(i in 1:length(cluster)){
  x_temp <- as.numeric(out[out[,1]==cluster[i],2]) # n degs
  y_temp <- as.numeric(out[out[,1]==cluster[i],4]) # recall
  lines(x = x_temp, y = y_temp, type = "b", lwd = 3, col = colors[i], lty = 1)
}
legend("topright", cluster, cex=0.8, col=colors,  lty=1, lwd = 4, title="Clusters")
dev.off()

# Recall plot - log axis
pdf(file = "RECALL_at_thresholds_CD_log.pdf", paper = "a4r")
plot(x = c(min(as.numeric(out[,2])), max(as.numeric(out[,2]))),
     y = c(0,(max(as.numeric(out[,4]))+0.05)),type = "n", 
     xlab = "n top DEGs", ylab = "Recall", main = "RECALL for CD data", log = "x")
for(i in 1:length(cluster)){
  x_temp <- as.numeric(out[out[,1]==cluster[i],2]) # n degs
  y_temp <- as.numeric(out[out[,1]==cluster[i],4]) # recall
  lines(x = x_temp, y = y_temp, type = "b", lwd = 3, col = colors[i], lty = 1)
}
legend("topright", cluster, cex=0.8, col=colors,  lty=1, lwd = 4, title="Clusters")
dev.off()

# Precision plot
pdf(file = "PRECISION_at_thresholds_CD.pdf", paper = "a4r")
plot(x = c(min(as.numeric(out[,2])), max(as.numeric(out[,2]))),
     y = c(0,(max(as.numeric(out[,3]))+0.05)),type = "n", 
     xlab = "n top DEGs", ylab = "Precision", main = "PRECISION for CD data")
for(i in 1:length(cluster)){
  x_temp <- as.numeric(out[out[,1]==cluster[i],2]) # n degs
  y_temp <- as.numeric(out[out[,1]==cluster[i],3]) # precision
  lines(x = x_temp, y = y_temp, type = "b", lwd = 3, col = colors[i], lty = 1)
}
legend("topright", cluster, cex=0.8, col=colors,  lty=1, lwd = 4, title="Clusters")
dev.off()

# Precision plot
pdf(file = "PRECISION_at_thresholds_CD_log.pdf", paper = "a4r")
plot(x = c(min(as.numeric(out[,2])), max(as.numeric(out[,2]))),
     y = c(0,(max(as.numeric(out[,3]))+0.05)),type = "n", 
     xlab = "n top DEGs", ylab = "Precision", main = "PRECISION for CD data", log = "x")
for(i in 1:length(cluster)){
  x_temp <- as.numeric(out[out[,1]==cluster[i],2]) # n degs
  y_temp <- as.numeric(out[out[,1]==cluster[i],3]) # precision
  lines(x = x_temp, y = y_temp, type = "b", lwd = 3, col = colors[i], lty = 1)
}
legend("topright", cluster, cex=0.8, col=colors,  lty=1, lwd = 4, title="Clusters")
dev.off()

setwd(fp)
rm(list = ls())

# SUMMARIZE FOR LITERATURE REVIEW & FC CRITERIA CHECKING - using top 1800 DEGs
###################################
fp <- getwd()
setwd(paste(fp, "/..", sep=""))

# Select which files to include in summary
dir <- "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/"
lf <- list.files(path = dir)
lf <- lf[grepl(lf, pattern = "EVALUATION_predicted_drugs_MAST_DEGs_log\\(1.5)_threshold_")]
lf1 <- lf[grepl(pattern = "_top_1800_", x = lf)] # since top DEG cut-off was set to 1800 DEG
lf2 <- lf[!grepl(pattern = "_top_", x = lf)]
for(i in 1:length(lf1)){
  temp <- unlist(strsplit(lf1[i], split = "_")[[1]][9])
  lf2 <- lf2[!grepl(lf2, pattern = paste("EVALUATION_predicted_drugs_MAST_DEGs_log\\(1.5)_threshold_cluster_", temp, sep=""))]
}
lf <- c(lf1, lf2)
rm(lf1, lf2, i, temp)

# FUNCTION FOR LOADING INPUT FILES
load_drug_prediction <- function(dir_network_pred_file){
  # load network prediction file
  in_file <- as.matrix(read.table(file = dir_network_pred_file, header = T, sep = "\t", stringsAsFactors = F))
  if(ncol(in_file)==1){
    print(paste("NO DRUGS FOR: ",dir_network_pred_file, sep=""))
    return(NULL)
  } else {
    return(in_file)
  }
}

out <- foreach(i = 1:length(lf), .combine = "rbind") %do% {
  in_file <- load_drug_prediction(paste(dir, lf[i], sep=""))
  if(!is.null(in_file)){
    cluster <- unlist(strsplit(lf[i], split="_"))[9]
    return(cbind(in_file, cluster))
  } else {
    return(in_file)
  }
}
temp <- unique(out[,1])
temp2 <- vector()
for(i in 1:length(temp)){
  temp2[i] <- sum(out[,1]%in% temp[i])
}
temp <- cbind(temp, temp2)
rm(i, temp2)
out <- cbind(out, temp[match(out[,1], temp[,1]),2])
colnames(out)[ncol(out)] <- "in_n_clusters" 
out <- out[order(- as.numeric(out[,13]), out[,1], decreasing = F),]

write.table(out, file = paste(dir, "COMBINED_EVALUATION_FILES_for_checking_FC_criterion.txt", sep=""), sep="\t", col.names = T, row.names = F)





# LOAD centrality for CD clusters
centrality <- as.matrix(read.table(file = "Input/CD GSE134809/NicheNet/Cell-cell_centrality_summary_CD_ileum.txt", sep="\t", header = T, stringsAsFactors = F))
rownames(centrality) <- centrality[,1]
centrality <- centrality[,-1]


            
            
