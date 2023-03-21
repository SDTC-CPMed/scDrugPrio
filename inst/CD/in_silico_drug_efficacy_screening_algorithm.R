# Drug analysis script
# Memory heavy algorithm, for CPU heavy (slow) algorithm use the original python script
# R version of Barabasi Nature Communication 2016 algorithm

print("LOADING: in_silico_drug_efficacy_screening.R")

library(doParallel)
library(foreach)
library(igraph)
library(matrixStats)

# Drug prediction algorithm
###############################
# Variable input:
###############################
# ppi = matrix with 2 columns, interaction between proteins in the same row. Self-loops will be removed. LCC will be used as PPI.
# drug matrix = matrix specifying all drug targets for a given drug in one column, column name is equal with drug name
# disease_module = vector; specyfing genes that are disease associated
# min_bin_size = integer; minimum size for each bin. Recommended to be at least 100 nodes in a bin (based on Nat Com 2016 Barabasi et al publication)
# disease_module_name = character string indicating which groups of genes is tested, used for naming plots and files
# n_random_iterations = integer indicating n permutations
# cores = integer specifying n of cores that can be used for analysis (doParallel)
# unadj = boolean specifying if results not adjusted for node degree should be saved.
# bin_adj = boolean specyfing if results adjusted for node degree should be saved. If both unadj and bin_adj equal FALSE no results will be calculated.
# recycle = boolean; reuse:  ppi, ppi_distances, bins, all_genes   saved in "/temp" folder 
in_silico_drug_efficacy_screening <- function(ppi, drug_matrix, disease_module, min_bin_size = 100, disease_module_name = "no_name",
                                              n_random_iterations = 1000, cores = NA,  unadj =F, bin_adj = T,
                                              disease_module_lcc = F, recycle = F, out.dir){
  setwd(out.dir)
  if(!recycle){
    unlink("temp", recursive = T)
  }
  
  set.seed(35)
  
  if(all(c(unadj==F, bin_adj==F))){
    print("Both 'unadj' and 'bin_adj' ==  FALSE, no results will be calculated")
    return()
  }
  
  # set up parallel running
  if(is.na(cores)){
    cores <- detectCores()
    if(cores>1){
      cores <- cores -1
    }
  }
  registerDoParallel(cores = cores)
  
  # Load / format variables needed
  if(recycle){
    print(list.files(path= "temp"))
    ppi_dist <- read.table(file = "temp/ppi_distances.txt", sep="\t", header = F)
    colnames(ppi_dist) <- as.character(ppi_dist[1,])
    ppi_dist <- ppi_dist[-1,]
    rownames(ppi_dist) <- as.character(ppi_dist[,1])
    ppi_dist <- as.matrix(ppi_dist[,-1])
    print("ppi_dist LOADED")
    bins <- read.table(file = "temp/bins.txt", sep="\t", header = F)
    print("bins LOADED")
    all_genes <- as.vector(read.table(file = "temp/all_genes.txt", sep="\t", header = F)[,1])
    print("all_genes (included in LCC of PPI) LOADED")
    drug_matrix <- read.table(file = "temp/drug_matrix.txt", sep="\t", header = T)
    print("drug_matrix LOADED")
    print("")
    print("RECYCLING DONE")
    
  } else {
    
    dir.create(path = paste(getwd(), "/temp", sep=""))
    
    # Filtering PPI
    ppi <- ppi[!(ppi[,1]==ppi[,2]),] # no self loops
    ppi <- ppi[!is.na(ppi[,1]),] # no NA
    ppi <- ppi[!is.na(ppi[,2]),] # no NA
    
    # Setting up PPI
    ppi <- data.frame(node1 = ppi[,1], node2 = ppi[,2])
    ppi_graph <- graph_from_data_frame(ppi, directed = F) # igraph package
    V(ppi_graph)$disease <- 1
    all_genes <- single_connected_subgraph_extraction(ppi_graph, plot = F) # only PPI proteins that are connected to other PPI proteins
    all_genes <- all_genes[rowSums(is.na(all_genes))<ncol(all_genes),] # removes empty rows
    all_genes <- as.vector(all_genes[, colSums(is.na(all_genes))==0]) # 'all_genes' will be used later on as it includes all relevant unique proteins in the LCC
    ################ SAVE ALL_GENES ################
    write.table(matrix(all_genes, ncol=1), file = "temp/all_genes.txt", sep="\t", col.names = F, row.names = F)
    rm(ppi_graph)
    
    # Selecting only PPI proteins found in LCC
    ppi <- as.matrix(ppi)
    ppi <- ppi[ppi[,1]%in% all_genes,]
    ppi <- ppi[ppi[,2]%in% all_genes,]
    ################ SAVE PPI ################
    write.table(ppi, file = "temp/ppi_lcc.txt", sep="\t", col.names = T, row.names = F)
    ppi <- data.frame(node1 = ppi[,1], node2 = ppi[,2])
    ppi_graph <- graph_from_data_frame(ppi, directed = F) # igraph package
    print("PPI = READY")
    
    # Creating / Loading bins
    if(bin_adj){
      bins <- bin_creation_by_min_bin_size(ppi = ppi, min_bin_size = min_bin_size)
      ################ SAVE BINS ################
      write.table(bins, file = "temp/bins.txt", sep="\t", col.names = F, row.names = F)
      print("Degree binning = READY")
    }
    rm(ppi)
    
    # Calculate distances between all proteins in LCC
    ppi_dist <- distances(ppi_graph, v = all_genes, to = all_genes)
    ################ SAVE PPI_DISTANCES ################
    colnames(ppi_dist) <- as.character(colnames(ppi_dist))
    rownames(ppi_dist) <- as.character(rownames(ppi_dist))
    write.table(ppi_dist, file = "temp/ppi_distances.txt", sep="\t", col.names = NA, row.names = T)
    print("Distances between proteins = CALCULATED")
    rm(ppi_graph)
    
    # formatting drug matrix to only include drugs with at least one target found in PPI
    drug_matrix <- as.matrix(drug_matrix)
    for(i in c(1:ncol(drug_matrix))){
      temp <- as.character(drug_matrix[,i])
      drug_matrix[,i] <- NA
      if(any(!is.na(temp))){
        temp <- temp[!is.na(temp)]
      }
      temp <- temp[temp%in%all_genes]
      temp <- temp[!duplicated(temp)]
      if(length(temp)>0){
        drug_matrix[1:length(temp),i] <- temp
      }
    }
    rm(temp, i)
    drug_matrix <- drug_matrix[rowSums(!is.na(drug_matrix))>0,colSums(!is.na(drug_matrix))>0]
    
    # sorting drugs by n drug targets
    n_targets <- nrow(drug_matrix) - colSums(is.na(drug_matrix))
    drug_matrix <- drug_matrix[,order(n_targets, decreasing = T)]
    rm(n_targets)
    print("Drug matrix = FORMATTED")
    ################ SAVE DRUG_MATRIX ################
    write.table(drug_matrix, file = "temp/drug_matrix.txt", sep="\t", col.names = T, row.names = F)
  }
  
  # formatting disease_module
  disease_module <- as.vector(disease_module)
  disease_module <- disease_module[disease_module %in% all_genes] # must be in PPI
  
  if(disease_module_lcc){ # if TRUE calculate LCC formed by disease_module genes and use that instead
    ppi <- read.table(file = "temp/ppi_lcc.txt", sep="\t", header = T)
    ppi <- data.frame(node1 = ppi[,1], node2 = ppi[,2])
    ppi_graph <- graph_from_data_frame(ppi, directed = F)
    
    print("Identifying LCC for disease module genes")
    V(ppi_graph)$disease <- 0
    V(ppi_graph)[V(ppi_graph)$name %in% disease_module]$disease <- 1
    temp <- single_connected_subgraph_extraction(ppi_graph, plot = F) # only PPI proteins that are connected to other PPI proteins
    if(ncol(temp)>1){
      temp <- temp[rowSums(!is.na(temp))>0,]
      if(is.matrix(temp)){
        disease_module <- as.vector(temp[,colSums(is.na(temp))==0])
      } else {
        disease_module <- temp[!is.na(temp)]
      }
    } else {
      disease_module <- as.vector(temp[!is.na(temp[,1]),1])
    }
    
    rm(temp, ppi, ppi_graph)
  }
  
  if(length(disease_module)==0){
    print("NO DISEASE GENES")
    return()
  }
  
  #################################################
  # Unadjusted analysis
  # random drugs (NOT (!) the standard algorithm from Barabasi group)
  #################################################
  if(unadj){
    print("RUNNING unadjusted drug prediction based on closest distances")
    n_targets <- nrow(drug_matrix) - colSums(is.na(drug_matrix))
    n_targets <- unique(n_targets)
    random_drugs <- foreach(i = c(1:length(n_targets)), .combine = cbind) %dopar% {
      # returns vector with n_random_iterations elements, cotaining average shortest path distances between the given disease_module
      # and 'n_targets' random drug targets - not adjusted for node degree.
      temp <- random_drug_target_distances(all_genes, n_targets[i], disease_module, ppi_dist, n_random_iterations)
      temp <- c(n_targets[i], temp)
      return(temp)
    }
    colnames(random_drugs) <- random_drugs[1,]
    random_drugs <- random_drugs[-1,]
    print("Closest distance between random drugs and random disease modules = CALCULATED")
    # Statistics
    n_targets <- nrow(drug_matrix) - colSums(is.na(drug_matrix))
    
    # Analysis of closest distances
    out <- foreach(i = c(1:ncol(drug_matrix)), .combine = rbind) %dopar% {
      temp <- rep(0, times = 6)
      # n targets
      temp[1] <- n_targets[i]
      # Actual drug-disease distance
      drug_genes <- drug_matrix[,i]
      drug_genes <- unique(drug_genes[!is.na(drug_genes)])
      temp[2] <- average_closest_distance(ppi_dist, from = disease_module, to = drug_genes)
      # Mean random
      x <- as.numeric(random_drugs[,colnames(random_drugs)==n_targets[i]])
      temp[4] <- mean(x)
      # SD random
      temp[5] <- sd(x)
      # Z
      if(temp[2] !="NaN"){
        temp[3] <- (as.numeric(temp[2])-as.numeric(temp[4]))/as.numeric(temp[5])
      } else {
        temp[3] <- NA
        temp[4] <- NA
      }
      # P
      temp[6] <- pnorm(as.numeric(temp[3]))
      # drug name
      temp <- c(colnames(drug_matrix)[i],temp)
      return(temp)
    }
    
    colnames(out) <- c("Drug", "n_targets", "mean d(c)", "Z", "mean(random)", "SD(random)", "P")
    write.table(out, file = paste("drug-disease_closest_distances_vs_random_unadjusted__", disease_module_name,".txt", sep = ""), sep = "\t", col.names = T, row.names = F)
    
    '
    # Histogram
    n_targets <- unique(n_targets)
    for(i in c(1:length(n_targets))){
      x <- as.numeric(random_drugs[,colnames(random_drugs)==n_targets[i]])
      png(filename = paste("Histograms/Distances_random_drug_with_n_", n_targets[i], "_targets_to_",disease_module_name,".png", sep = ""), units = "px", pointsize = 12)
      h <- hist(x = x, main = paste("Histogram of distances random drugs w n = ", n_targets, " targets", sep =""), xlab = "d", ylab = "freq")
      plot(h)
      dev.off()
    }
    rm(random_drugs, out, i, n_targets)
    '
    print("Node degree unadjusted distance Z scores = CALCULATED")
  }
  
  #################################################
  # Adjusted analysis
  # random drugs (standard algorithm from Barabasi group)
  #################################################
  if(bin_adj){
    print("RUNNING bin-adjusted drug prediction based on closest distances")
    
    #
    # Maybe memory decreasing to remove bins that will not be needed for further analysis - should they exist.
    #
    
    # random_bin_drugs (1 column for every drug, 1 row for every n_random_iteration)
    random_bin_drugs <- foreach(i = c(1:n_random_iterations), .combine = rbind) %dopar% {
      # returns vector with n_random_iterations elements cotaining closest distances between disease and
      # (bin-adjusted) random drug tragets for a given drug
      random_drug_target_bin_adjusted_distances(bins, drug_matrix, disease_module, ppi_dist, i)
    }
    
    print("Closest distance between random drugs and random disease modules = CALCULATED")
    
    write.table(random_bin_drugs, file = paste("random_bin_drugs_",disease_module_name,".txt", sep=""), sep = "\t", col.names = colnames(drug_matrix), row.names = F)
    
    # Analysis of closest distances
    n_targets <- nrow(drug_matrix) - colSums(is.na(drug_matrix))
    out <- foreach(i = 1:ncol(drug_matrix), .combine = rbind) %dopar% {
      temp <- rep(NA, times = 6)
      # n targets
      temp[1] <- n_targets[i]
      # Actual drug-disease distance
      drug_genes <- drug_matrix[,i]
      drug_genes <- drug_genes[!is.na(drug_genes)]
      if(length(drug_genes)>0){
        drug_genes <- as.character(unique(drug_genes))
        temp[2] <- average_closest_distance(ppi_dist, from = disease_module, to = drug_genes)
      }
      # Mean random
      x <- as.numeric(random_bin_drugs[,i])
      temp[4] <- mean(x)
      # SD random
      temp[5] <- sd(x)
      # Z
      if(temp[2]!="NaN"){
        temp[3] <- (as.numeric(temp[2]) - as.numeric(temp[4]))/as.numeric(temp[5])
      } else {
        temp[3] <- NA
        temp[2] <- NA
      }
      # P
      temp[6] <- pnorm(as.numeric(temp[3]))
      # drug name
      temp <- c(colnames(drug_matrix)[i],temp)
      return(temp)
    }
    colnames(out) <- c("Drug", "n_targets", "mean_d(c)", "Z(c)", "mean(random)", "SD(random)", "P")
    write.table(out, file = paste("drug-disease_closest_distances_vs_random_bin_adjusted__",disease_module_name,".txt", sep=""), sep = "\t", col.names = T, row.names = F)
    rm(n_targets)
    print("Node degree adjusted closest distance Z scores = CALCULATED")
    '
    # Histograms
    for(i in c(1:ncol(drug_matrix))){
      x <- as.numeric(random_bin_drugs[,i])
      png(filename = paste("Histograms/", colnames(drug_matrix)[i], "_bin_adjusted_random_drug_distances_", disease_module_name,".png", sep = ""), units = "px", pointsize = 12)
        h <- hist(x = x, main = "Histogram of closest distances n random drug targets", xlab = "d", ylab = "freq")
        plot(h)
      dev.off()
    }
    '
    rm(random_bin_drugs)
  }
  print("Drug prediction = DONE")
  return(out)
}







# General LCC recognition algorithm
###############################literature_PPI
# Variable input:
###############################
# network = igraph object; protein-protein interaction network. V(network)$disease==1 if protein should be considered for LCC calculation.
# plot = boolean; returns plot of network
single_connected_subgraph_extraction <- function(network, plot=F){
  set.seed(35)
  
  extraction_para <- V(network)[V(network)$disease==1]
  subGr <- induced_subgraph(network, vids=extraction_para)
  
  if(plot){
    library(qgraph)
    e <- get.edgelist(subGr)
    e[,1] <- match(e[,1], V(subGr)$name)
    e[,2] <- match(e[,2], V(subGr)$name)
    mode(e) <- "numeric"
    l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subGr), area=8*(vcount(subGr)^2),repulse.rad=(vcount(subGr)^3.1))
    plot(subGr,layout=l,vertex.size=2,vertex.label=NA, vertex.color=c("grey","red")[1+(V(subGr)$disease==1)])
  }
  
  subGr <- decompose(subGr, mode = "weak")
  
  if(!length(subGr)==0){
    out <- matrix(data = NA, nrow = length(V(network)), ncol = length(subGr))
    for(i in c(1:ncol(out))){
      temp <- V(subGr[[i]])$name
      out[1:length(temp),i] <- temp
    }
    out <- out[rowSums(is.na(out))<ncol(out),]
    if(length(subGr)==1){
      out <- matrix(out, ncol = 1)
    }
    if(is.vector(out)){
      if(length(subGr)>1){
        out <- matrix(data = out, nrow = 1)
      }
    }
  } else {
    print("No disease associated genes found in graph")
    out <- 0
  }
  rm(subGr)
  return(out)
}

# Creates vector filled with average closest distances between x random drug targets and y random disease module genes 
###############################
# Variable input:
###############################
# all_genes = vector; all unique genes in PPI
# n_targets = integer; n drug targets
# disease_module = vector; specifying disease associated genes in PPI
# ppi_dist = matrix; specifying distances between all proteins in PPI
# n_random_iterations = integer; n permutations
random_drug_target_distances <- function(all_genes, n_targets, disease_module, ppi_dist, n_random_iterations){
  set.seed(n_random_iteration)
  out <- vector()# bins for drugs
  for(j in c(1:n_random_iterations)){
    temp <- sample(all_genes, size = n_targets, replace = F)
    random_disease_module <- sample(all_genes, size = length(disease_module), replace = F)
    out[j] <- average_closest_distance(ppi_dist, from = random_disease_module, to = temp)
  }
  return(out)
}

# Creates vector filled with average closest distances between x random drug targets (bin-adjusted for node degree) and y random disease module genes (bin-adjusted for node degree)
###############################
# Variable input:
###############################
# bins = matrix; one column per bin, each bin includes proteins with similar node-degree
# drugs = matrix; all drug targets for specific drug in a column
# disease_module = vector; specifying disease associated genes in PPI
# ppi_dist = matrix; specifying distances between all proteins in PPI
random_drug_target_bin_adjusted_distances <- function(bins, drugs, disease_module, ppi_dist, seed){
  set.seed(seed)
  # which bins needed for randomising drugs / randomising disease model?
  
  # which_bins_drugs -> columns = drugs, bins = rows
  which_bins_drugs <- foreach(k = c(1:ncol(drugs)), .combine = cbind) %dopar% {
    # bins for drugs
    temp <- drugs[,k]
    temp <- temp[!is.na(temp)]
    d.bins <- foreach(j = c(1:ncol(bins)), .combine = c) %do% {
      return(sum(temp%in%bins[,j]))
    }
    return(d.bins)
  }
  which_bins_disease <- vector()
  # bins for disease
  for(j in c(1:ncol(bins))){
    which_bins_disease <- c(which_bins_disease, sum(disease_module %in% bins[!is.na(bins[,j]),j]))
  }
  
  # randomize disease model
  random_disease <- vector()
  for(k in 1:ncol(bins)){
    if(which_bins_disease[k]>0){
      temp_bins <- bins[,k]
      temp_bins <- temp_bins[!is.na(temp_bins)]
      random_disease <-c(random_disease, sample(temp_bins, size = which_bins_disease[k], replace = F))
    }
  }
  if(length(random_disease)<length(disease_module[!is.na(disease_module)])){
    message("WARNING: Too few random disease genes!")
  }
  
  # randomize drug matrix & calculate average_closest_distance() between all 'random_drugs' and 'random_disease'
  if(length(disease_module)>1){
    ppi_dist <- ppi_dist[rownames(ppi_dist)%in%random_disease,]
  }
  out <- foreach(k = c(1:ncol(which_bins_drugs)), .combine = c) %do% { # for every drug
    # create random drug
    random_drug <- vector()
    for(j in 1:ncol(bins)){
      if(which_bins_drugs[j,k]>0){
        temp_bins <- bins[,j]
        temp_bins <- temp_bins[!is.na(temp_bins)]
        random_drug <- c(random_drug, sample(temp_bins, size = which_bins_drugs[j,k], replace = F))
      }
    }
    if(length(random_drug)<sum(which_bins_drugs[,k])){
      message("WARNING: Too few random drug targets!")
    }
    return(average_closest_distance(ppi_dist, from = random_disease, to = random_drug))
  }
  return(out)
}

# Calculate average closest path between drug targets and the nearest disease protein = d(c) in Nat Com 2016 Barabasi et al
###############################
# Variable input:
###############################
# ppi_dist = matrix; specifying distances between all proteins in PPI
# from & to = vectors; specifying genes between which distances should be calculated. E.g.: from = disease_module, to = drug
# to = drug targets
average_closest_distance <- function(ppi_dist, from, to){
  if(any(length(from)==0, length(to)==0)){
    print("Empty vector supplied for from / to of function average_closest_distance(ppi_dist, from, to)")
    return(NA)
  }
  if(length(from)>1){
    ppi_dist <- ppi_dist[rownames(ppi_dist)%in%from,]
    
    if(length(to)>1){ # > 1 diseas protein & > 1 target protein
      ppi_dist <- ppi_dist[,colnames(ppi_dist)%in%to]
      return(mean(colMins(ppi_dist))) # average closest path
      
    } else { # only 1 target protein
      ppi_dist <- ppi_dist[,colnames(ppi_dist)%in%to]
      return(min(ppi_dist)) # closest path - no averaging needed
    }
  } else { # only 1 disease protein!
    if(length(to)>1){
      ppi_dist <- ppi_dist[,colnames(ppi_dist) %in% to]
      ppi_dist <- ppi_dist[rownames(ppi_dist)%in%from,]
      return(mean(ppi_dist)) # average closest path
    } else { # only one drug target & 1 disease protein
      return(ppi_dist[rownames(ppi_dist)%in% from, colnames(ppi_dist) %in% to])
    }
  }
}

# Create node-degree bins for later adjustment
###############################
# Variable input:
###############################
# ppi = matrix with 2 columns, interaction between proteins in the same row.
# min_bin_size = integer; minimum size for each bin
bin_creation_by_min_bin_size <- function(ppi, min_bin_size){
  ppi <- ppi[,1:2]
  ppi <- c(as.character(ppi[,1]), as.character(ppi[,2]))
  
  # calculate node degree
  n_temp <- vector()
  temp <- unique(ppi)
  for(i in c(1:length(temp))){
    n_temp <- c(n_temp, sum(ppi == temp[i]))
  }
  names(n_temp) <- temp
  n_temp <- sort(n_temp)
  rm(temp)
  
  # create bins
  bins <- matrix(NA, nrow = length(n_temp), ncol = length(n_temp)/min_bin_size)
  bins_col <- vector()
  degree <- unique(n_temp)
  degree <- degree[order(degree, decreasing = F)]
  temp <- vector()
  b <- 1
  for(i in c(1:length(degree))){
    temp <- c(temp, names(n_temp[n_temp == degree[i]]))
    if(length(temp)>=min_bin_size){
      bins[1:length(temp),b] <- temp
      bins_col <- c(bins_col, i)
      b <- b+1
      temp <- vector()
    }
  }
  if(length(temp)>0 & length(temp)< min_bin_size){
    b <- b-1
    temp2 <- bins[,b]
    temp2 <- unique(temp2)
    temp2 <- temp2[!is.na(temp2)]
    temp2 <- c(temp2, temp)
    bins[1:length(temp2),b] <- temp2 
  }
  rm(temp, temp2, b, degree, bins_col)
  bins <- bins[,colSums(!is.na(bins))>0]
  bins <- bins[rowSums(!is.na(bins))>0,]
  return(bins)
}

# Precision and recall at different Z-score thresholds
###############################
# Variable input:
###############################
# predefined_list = vector; lists drug names that are known to have an effect on disease.
# results = matrix; output from in_silico_drug_efficacy_screening()
precision_and_recall_from_z_score <- function(predefined_list, results){
  predefined_list <- as.vector(predefined_list)
  min.z <- as.numeric(results[,4])
  min.z <- min.z[is.finite(min.z)]
  if(length(min.z)==0){
    print("No finite z-scores")
    return(NA)
  }
  min.z <- min(min.z)
  if(min.z > 0){
    print("No z-scores below 0")
    return(NA)
  }
  thresh <- seq(from=0, to = min.z, by = -0.1) # z-score thresholds for when a drug is considered "close" to a disease module
  precision <- vector()
  recall <- vector()
  proximit_predefined_drugs <- vector()
  proximit_drugs_n <- vector()
  for(i in c(1:length(thresh))){
    close <- results[as.numeric(results[,4]) <= thresh[i],1] # z-score less or equal to thresh is considered close 
    # recall
    proximit_predefined_drugs[i] <- sum(close%in%predefined_list)
    recall[i] <- as.numeric(proximit_predefined_drugs[i])/length(predefined_list) # true positived / (true positives + false negatives)
    # precision
    proximit_drugs_n[i] <- length(close)
    precision[i] <- as.numeric(proximit_predefined_drugs[i])/as.numeric(proximit_drugs_n[i])
  }
  out <- rbind(recall, precision, proximit_predefined_drugs, proximit_drugs_n)
  out <- cbind(c("Recall", "Precision", "n_proximit_predefined_drugs", "total_n_proximit_drugs"), out)
  colnames(out) <- c("Z-scores:", paste("z=", thresh, sep=""))
  return(out)
}



