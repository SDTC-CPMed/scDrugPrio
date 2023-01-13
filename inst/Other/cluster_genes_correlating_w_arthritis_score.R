# Correlation between gene expression and arthritis score in every single cluster - Oleg algorithm

#Assumptions: 
#   1 data set with column "Genes" and then other columns showing gene expression per mouse
#   2 data set with arthritis scores per mouse.

# geneExpr / geneX = data.frame, first column = factor of gene names, remaining columns = expression values for individual cells
# each expression value corresponding to one single cell that matches a certain cell type in a certain mouse
# E.g.: data frame for granulocytes
# gene I mouse1 I mouse2 I mouse3 I mouse4
#  A   I   0    I   1    I   NA   I   1
#  A   I   1    I   1    I   NA   I   0
#  A   I   2    I   NA   I   NA   I   0
#  B   I   4    I   7    I   NA   I   1
#  B   I   4    I   6    I   NA   I   1
#  B   I   5    I   NA   I   NA   I   1

# arthritis score = matrix. colnames correspond to mouse, row 1 contains arthritis scores for all mice.

library(doParallel)

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
# LOAD denoised DCA expression matrix 
geneExpr <- read.table(file = "Input/DCA_adjusted_matrix/mean.tsv", sep="\t", header = T)
rownames(geneExpr) <- geneExpr[,1]
geneExpr <- geneExpr[,-1]
clusters <- read.table(file = "Output/Seurat_clusters/Cell_identity.txt", sep="\t", header = T)
# rename cells to include cluster information
colnames(geneExpr) <- paste(colnames(geneExpr), "_cluster_", clusters[match(colnames(geneExpr),clusters[,1]),2], sep="")
# These are the actual values for our data set:
ArthScore <- matrix(c(0,0,0,0,1.5,1.5,2,2,2.5), nrow = 1)
colnames(ArthScore) <- c("Healthy_mouse_1", "Healthy_mouse_2", "Healthy_mouse_3", "Healthy_mouse_4", "Sick_mouse_1", "Sick_mouse_3", "Sick_mouse_4", "Sick_mouse_5", "Sick_mouse_6")
ArthScore=as.data.frame(ArthScore)

format_expression_matrix_for_bootstrap_correlation <- function(x, cell_type, sort_by = NULL){
  if(is.null(sort_by)){
    sort_by <- c(1,2,3,4,5,6,7,8,9)
  }
  x <- as.matrix(x[,grep(pattern = cell_type, x = colnames(x))]) # only 1 given cell type
  
  n_cells_mouse <- vector()
  for(i in 1:length(sort_by)){
    n_cells_mouse <- c(n_cells_mouse, length(grep(pattern = sort_by[i], x = colnames(x))))
  }
  print("Cells per mouse:")
  print(cbind(sort_by, n_cells_mouse))
  
  gene = sort(rep(rownames(x), times = max(n_cells_mouse))) # creating first col of data frame specifying gene
  x <- x[order(rownames(x)),]
  out <- matrix(NA, nrow = length(gene), ncol = 1+length(sort_by))
  out[,1] <- gene
  
  for(i in 1:length(sort_by)){
    if(n_cells_mouse[i]>3){ # removes values if less than 3 cells available for any given mouse
      temp <- x[,grep(pattern = sort_by[i], x = colnames(x))]
      for(j in 1:nrow(x)){
        temp2 <- as.numeric(temp[j,])
        row <- max(n_cells_mouse)*(j-1)+1
        out[row:(row-1+length(temp2)),(i+1)] <- temp2
      }
    }
  }
  print("DONE creation matrix")
  print("Now formatting into data frame")
  colnames(out) <- c("Genes", sort_by)
  out <- as.data.frame(out)
  out$Genes <- as.factor(out$Genes)
  return(out)
}


#first two parameters are arthritis data and gene data, 
# B is amount of bootstrap replications
# debug if percent of runtime execution is needed to be printet
testCorrelation<- function(artS, genX, dev=0.2, B=1000, debug=T, ncores = 1){
  
  if(ncores>1){
    registerDoParallel(cores = ncores)
  }
  
  fisher<- function(x) log((1+x)/(1-x))
  
  p=ncol(artS)
  n=nrow(genX)
  meanG=aggregate(.~Genes, data=genX, FUN=mean, na.rm = T, na.action = NULL)
  meanGM=as.matrix(t(meanG[,2:(p+1)]))
  
  
  scores=as.numeric(artS)
  cor0=as.numeric(cor(meanGM, scores,use = "pairwise.complete.obs"))
  
  #Boots=matrix(0, nrow=B, ncol=nlevels(genX$Genes)) #old
  ng=length(levels(genX$Genes))
  
  Boots <- foreach(b = 1:B, .combine = rbind) %dopar% { #new
  #for (b in 1:B){ #old
    meanM=matrix(NA, ncol=ng, nrow=p)
    for (d in 1:p){
      GM=matrix(genX[,d+1], ncol=ng)
      mode(GM) <- "numeric"
      ind2=which(!is.na(GM[,1]))
      ind3=sample(ind2, length(ind2), replace=T)
      genM=colMeans(GM[ind3,])
      if(!is.na(genM[1])) {
        meanM[d,]=genM
      }
    }
    colnames(meanM)<-levels(genX$Genes)
    
    artB=runif(p, min=pmax(scores-dev, numeric(p)), max=scores+dev )
    artBM=matrix(artB, ncol=1)
    
    cors=cor(meanM, artBM, use = "pairwise.complete.obs")
    
    
    #Boots[b,]=abs(as.numeric(cors)-cor0) #old
    if (b %%(B/10) == 0 & debug ) cat( b /(B/100), " percent done", "\n")
    return(abs(as.numeric(cors)-cor0)) #new
  }
  colnames(Boots)<-levels(genX$Genes)
  Boots1=fisher(Boots)
  Boots1SD=sqrt(colMeans(Boots1^2, na.rm = T))
  
  ASL=2-2*pnorm(fisher(abs(cor0)),0, Boots1SD)
  
  return(ASL)
}

#debug(testCorrelation)

setwd("Output/Correlation_genes") # output directory
cell_types <- paste("cluster_", unique(clusters[,2]),sep="")
for(i in 1:length(cell_types)){
  temp_cells <- format_expression_matrix_for_bootstrap_correlation(x = geneExpr, cell_type = cell_types[i], sort_by = colnames(ArthScore))
  print("Starting correlation analysis")
  p.vals <- testCorrelation(ArthScore, temp_cells, ncores = 20)
  write.table(data.frame(genes = rownames(geneExpr), P = p.vals), file = paste("Corr_P_", cell_types[i], ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  print(paste("Done with ", cell_types[i], sep=""))
}

# Create summary correlation genes
lf <- list.files()
out <- matrix(NA, nrow = nrow(geneExpr), ncol = length(lf))
for(i in 1:length(lf)){
  temp <- as.matrix(read.table(file = paste("Corr_P_cluster_",i,".txt",sep=""), sep="\t", header = T))
  x <- nrow(temp)
  temp <- temp[!is.na(temp[,2]),]
  print(paste("REMOVED ", x - nrow(temp), " genes with NAs instead of P values from 'Corr_P_cluster_",i,".txt'", sep=""))
  rm(x)
  
  if(any(temp[,2]<0.05)){
    if(sum(temp[,2]<0.05)==1){
      temp <- temp[temp[,2]<0.05,1]
    } else {
      temp <- temp[temp[,2]<0.05,]
      temp <- temp[order(temp[,2], decreasing = F),1] # most significant first
    }
    out[1:length(temp),i] <- temp
  }
  rm(temp)
}
colnames(out) <- paste("Cluster_",1:length(lf), sep="")
write.table(out, file = "Summary_Correlation_DCA_genes.txt", sep = "\t", col.names = T, row.names = F)
#out <- read.table(file = "Summary_Correlation_DCA_genes.txt", sep="\t", header = T, stringsAsFactors = F)

# Load mouse to human homologs
m2h <- as.matrix(read.table(file = "../../Input/Human-mouse_homologs/transl.txt",sep="\t", header = T, stringsAsFactors = F))

for(i in 1:ncol(out)){
  temp <- out[!is.na(out[,i]),i]
  #print(cbind(temp[1:10], m2h[match(temp[1:10],m2h[,3]),3:2]))
  temp <- as.numeric(m2h[match(temp,m2h[,3]),2])
  temp <- temp[!is.na(temp)]
  out[,i] <- NA
  if(length(temp)>0){
    out[1:length(temp),i] <- temp
  }
}
out <- out[,colSums(!is.na(out))>0]
write.table(out, file = "TRANSLATED_Summary_Correlation_DCA_genes_ordered.txt", sep = "\t", col.names = T, row.names = F)

# Remove everything from environment when done
rm(list= ls())

