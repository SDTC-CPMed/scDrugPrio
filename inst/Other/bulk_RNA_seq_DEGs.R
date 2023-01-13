# Create bulk-RNAseq DEGs

# Required functions:
##############################################################################

# Probe translation of micro array data
#########################################
# Input:
# GSE = GSE-ID for download from GEO
# x = array used
# gene = column name used for gene symbols in feature data
# split_gene_by = string seperation in "gene" column, commonly " /// ",  " // ",  ";",  "  " and so on.
# ref = some sort of reference vector containing all possible/interesting gene symbols
#
# Used for calling MicroArray_analysis
# g1 = group 1 GSE-IDs
# g2 = group 2 GSE-IDs
# stat_method = c("paired", "independent"); if paired than pairs are indicated by position in g1 and g2
#
# Output: selects useful probes for data set
MicroArray_probe_translation <- function(GSE, x = 1, gene, split_gene_by = "none", ref = NULL, g1, g2, stat_method = "independent", directory = NULL){
  library(BiocManager)
  library(GEOquery)
  
  # get data set
  gseX <- getGEO(as.character(GSE))
  
  # extract probes and genes
  transl <- cbind(gseX[[x]]@featureData@data[["ID"]], gseX[[x]]@featureData@data[[as.character(gene)]])
  print("Probes and genes loaded")
  
  # string split gene names
  if(!any(split_gene_by %in% "none")){
    #if probe ID corresponds to several entrez ID
    gene <- transl[,2]
    id <- transl[,1]
    
    r <- 0
    i <- 1
    while(i <= length(id)){
      gene2 <- gene[i]
      for(j in c(1:length(split_gene_by))){
        if(any(grepl(split_gene_by[j], gene2))){
          gene2 <- unlist(strsplit(gene2, split_gene_by[j]))
        }
      }
      if(!is.null(ref)){
        gene2 <- gene2[gene2%in%ref]
      }
      id2 <- rep(id[i], times=length(gene2))
      
      if(i == 1){
        id <- c(id2, id[2:length(id)])
        gene <- c(gene2, gene[2:length(gene)])
      } else if (i == length(id)){
        id <- c(id[1:(i-1)], id2)
        gene <- c(gene[1:(i-1)], gene2)
        break
      } else {
        id <- c(id[1:(i-1)], id2, id[(i+1):length(id)])
        gene <- c(gene[1:(i-1)], gene2, gene[(i+1):length(gene)])
      }
      
      # counting progress
      if(i > length(id)*r/10){
        print(paste("Probes ", r*10, "% processed", sep=""))
        r <- r+1
      }
      
      i <- i + length(gene2)
    }
    transl <- cbind(as.character(id), as.character(gene))
    transl <- transl[!duplicated(transl),]
    print("Probes translation done")
    rm(i, r, j, gene2, id2, gene, id)
  }
  
  # exclude probes mapping to multiple genes
  id <- unique(transl[,1])
  excl <- vector()
  for(i in c(1:length(id))){
    temp <- which(transl[,1]==id[i])
    if(length(temp)>1){
      excl <- c(excl, temp)
    }
  }
  transl <- transl[-excl,]
  print("Probes mapping to multiple genes excluded")
  rm(excl, temp, id)
  
  # find probe with highest FC if multiple probes for the same gene
  temp <- MicroArray_analysis(GSE, x, g1, g2, transl = cbind(transl[,1], transl[,1]), stat_method)
  transl <- cbind(rownames(temp), transl[match(transl[,1], rownames(temp)),2])
  transl <- transl[!duplicated(transl[,2]),]
  rm(temp)
  print("If multiple probes for the same gene, probes differing most between groups were selected!")
  
  if(is.null(directory)){
    write.table(transl, file = paste("selected_probe-gene_matches_", GSE, ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  } else {
    write.table(transl, file = paste(directory, "/selected_probe-gene_matches_", GSE, ".txt", sep=""), sep="\t", col.names = T, row.names = F)
  }
  
  return(transl)
}

# function for DEG calculation:
#########################################
# Input:
# GSE = GSE-ID for download from GEO
# x = array used
# g1 = group 1 GSE-IDs
# g2 = group 2 GSE-IDs
# transl = probe to Gene, probes = column1, genes = column2
# stat_method = c("paired", "independent"); if paired than pairs are indicated by position in g1 and g2
#
# Outputs FDR adjusted DEGs
MicroArray_analysis <- function(GSE, x = 1, g1, g2, transl, stat_method = "independent"){
  library(limma)
  library(BiocManager)
  library(GEOquery)
  
  # get data set
  gseX <- getGEO(as.character(GSE))
  expValues <- gseX[[x]]@assayData[["exprs"]]
  print("Expression values loaded")
  
  # remove unnecessary samples
  expValues <- expValues[,colnames(expValues)%in%c(g1,g2)]
  expValues <- expValues[rownames(expValues)%in%transl[,1],]
  rownames(expValues) <- transl[match(transl[,1],rownames(expValues)),2]
  class(expValues) <- "numeric"
  
  # build design table
  if(stat_method=="independent"){
    treat <- c(rep(1, times=length(g1)), rep(0, times = length(g2)))
    design <- model.matrix(~treat)
    rownames(design) <- c(g1, g2)
    rm(treat)
  }
  if(stat_method=="paired"){
    matrix_design <- matrix(data = NA, nrow=length(g1)*2, ncol = 3)
    colnames(matrix_design) <- c("sample", "pair", "treat")
    
    matrix_design[,1] <- c(g1, g2)
    matrix_design[1:length(g1),2] <- c(1:length(g1), 1:length(g2))
    matrix_design[1:length(g1),3] <- c(rep(1, times=length(g1)), rep(0, times = length(g2)))
    
    #order design table according to expression file columns
    matrix_design <- matrix_design[match(colnames(expValues), matrix_design[,1]),]
    
    #define design model.matrix
    group <- factor(matrix_design[,2])
    treat <- factor(matrix_design[,3], levels = c(0, 1))
    design <- model.matrix(~pair+treat)
    rownames(design) <- matrix_design[,1]
    rm(matrix_design, group, treat)
  }
  print("Design matrix built")
  
  #analysis
  fit <- lmFit(expValues, design)
  fit <- eBayes(fit)
  print("Limma package used")
  
  if(stat_method=="independent"){
    return(topTable(fit, number=nrow(fit), adjust.method = "BH", sort.by = "P"))
  } else {
    return(topTable(fit[,"treat1"], number=nrow(fit), adjust.method = "BH", sort.by = "P"))
  }
}

# bulk expression data for RA (microarray)
#########################################
# Synovial biopsies of rheumatoid arthritis and healthy controls
# PMID = 26711533
# Includes: 16 end-stage RA patients, 7 healthy controls. Unclear whether patients were treated or not.
# Platform: Affymetrix Human Genome U133 Plus 2.0 Array
'
GSE77298 <- function(){
 
  # Libraries
  library(BiocManager)
  library(limma)
  library(GEOquery)
  
  # Data
  gseX <- getGEO("GSE77298")
  
  # healthy and sick?
  temp <- gseX[["GSE77298_series_matrix.txt.gz"]]@phenoData@data
  healthy <- temp[1:7,2]
  sick <- temp[8:23,2]
  
  #translate (probe id to gene symbol)
  transl <- MicroArray_probe_translation(GSE = "GSE77298", x = 1, gene = "ENTREZ_GENE_ID", split_gene_by = " /// ", ref = NULL, g1 = sick, g2 = healthy,
                                         stat_method = "independent")
  
  # calculate DEGs
  degs <- MicroArray_analysis(GSE = "GSE77298", x= 1, g1 = sick, g2 = healthy, transl = transl, stat_method = "independent")
  return(degs)
}
'

# Synovial biopsies of rheumatoid arthritis and healthy controls
# PMID = 24690414
# Includes: 10 RA patients (undergoing joint replacement), 10 controls (early post mortem samples). Patients were treated.
# Platform: Affymetrix Human Genome U133A Array
GSE55235 <- function(){
  
  # Libraries
  library(BiocManager)
  library(limma)
  library(GEOquery)
  
  # Data
  gseX <- getGEO("GSE55235")
  
  # healthy and sick?
  # sick drug naive can be distinguished from sick-treated using Supplementary Data 5 in manuscript
  # here I only select sick drug naive and healthy drug naive from the same batch (batch: 1)
  temp <- gseX[["GSE55235_series_matrix.txt.gz"]]@phenoData@data
  
  healthy <- temp[1:10,2]
  sick <- temp[21:30,2]
  
  #translate (probe id to gene symbol)
  transl <- MicroArray_probe_translation(GSE = "GSE55235", x = 1, gene = "ENTREZ_GENE_ID", split_gene_by = " /// ", ref = NULL, g1 = sick, g2 = healthy,
                                         stat_method = "independent", directory = "Input/Bulk_data")
  
  # calculate DEGs
  degs <- MicroArray_analysis(GSE = "GSE55235", x= 1, g1 = sick, g2 = healthy, transl = transl, stat_method = "independent")
  return(degs)
}

# Whole blood gene expression of rheumatoid arthritis
# PMID = 30013029
# Includes: 45 drug naive RA patients, 21 treated RA patients and 35 healthy controls. Patients were untreated.
# Platform: Affymetrix Human Genome U133 Plus 2.0 Array

# We collect 30 healthy and 30 sick patients
GSE93272 <- function(){
  
  # Libraries
  library(BiocManager)
  library(limma)
  library(GEOquery)
  
  # Data
  gseX <- getGEO("GSE93272")
  
  # healthy and sick?
  # sick drug naive can be distinguished from sick-treated using Supplementary Data 5 in manuscript
  # here I only select sick drug naive and healthy drug naive from the same batch (batch: 1)
  temp <- gseX[["GSE93272_series_matrix.txt.gz"]]@phenoData@data
  
  healthy <- temp[grepl("Whole blood from healthy control",temp[,8]),]
  sick <- temp[grepl("Whole blood from rheumatoid arthritis",temp[,8]),]
  healthy <- healthy[healthy[,13]=="batch: 1",2]
  sick <- sick[sick[,13]=="batch: 1",2]
  
  #translate (probe id to gene symbol)
  transl <- MicroArray_probe_translation(GSE = "GSE93272", x = 1, gene = "ENTREZ_GENE_ID", split_gene_by = " /// ", ref = NULL, g1 = sick, g2 = healthy,
                                         stat_method = "independent", directory = "Input/Bulk_data")
  
  # calculate DEGs
  degs <- MicroArray_analysis(GSE = "GSE93272", x= 1, g1 = sick, g2 = healthy, transl = transl, stat_method = "independent")
  return(degs)
  
}

#########################################

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

degs <- GSE93272()
write.table(degs, file = "Input/Bulk_data/whole_blood_DEGs_GSE93272.txt", sep="\t", col.names = NA, row.names = T)
degs <- GSE55235()
write.table(degs, file = "Input/Bulk_data/synovium_DEGs_GSE55235.txt", sep="\t", col.names = NA, row.names = T)


rm(list = ls())
