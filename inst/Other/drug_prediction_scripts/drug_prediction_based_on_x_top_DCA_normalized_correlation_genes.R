# Drug prediction based on x top DCA normalized correlation genes

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

set.seed(35)

# x DEGs
x <- c(10,25,50,75,100,150,200,300,400,500,1000)
# Load DCA + MAST DEGs
degs <- read.table("Output/Correlation_genes/TRANSLATED_Summary_Correlation_DCA_genes.txt", sep="\t", header = T)
degs <- as.matrix(degs)
mode(degs) <- "character"
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/predefined_RA_drugs.txt", sep="\t", header = F)[,1])

# Set up output directories
dir.create("Output/Correlation_genes_predictions/literature_PPI")
dir.create("Output/Correlation_genes_predictions/HuRI_PPI")

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/Correlation_genes_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Analysis
save_name <- "DCA_normalized_correlating_genes_in_literature_PPI_1000"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){ # for every cluster
  for(j in 1:length(x)){ # using x top genes
    if(sum(!is.na(degs[,i]))>=x[j]){
      print(paste("CALCULATING: top ", x[j], " correlating genes from ", colnames(degs)[i], sep=""))
      out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[1:x[j],i]), disease_module_lcc = F, out.dir = out.dir,
                                               disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 32, recycle = rec, n_random_iterations = 1000)
      #rec <- T
      if(!(length(out)==1) & !all(is.na(out))){
        out <- precision_and_recall_from_z_score(predefined_list, out)
        print(out)
        write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], "_top_",x[j],"_correlating_genes.txt", sep=""), sep="\t", col.names = T, row.names = F)
      } else {
        print("No disease genes mapped to PPI")
      }
      rm(out)
    }
  }
}



# start drug prediction HuRI PPI
##################################################################################################################
print("START HuRI PPI network analysis")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/Correlation_genes_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists
# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/HuRI_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Analysis
save_name <- "DCA_normalized_correlating_genes_in_HuRI_PPI_1000"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){ # for every cluster
  for(j in 1:length(x)){ # using x top genes
    if(sum(!is.na(degs[,i]))>=x[j]){
      print(paste("CALCULATING: top ", x[j], " correlating genes from ", colnames(degs)[i], sep=""))
      out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(degs[1:x[j],i]), disease_module_lcc = F, out.dir = out.dir,
                                               disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 32, recycle = rec, n_random_iterations = 1000)
      #rec <- T
      if(!(length(out)==1) & !all(is.na(out))){
        out <- precision_and_recall_from_z_score(predefined_list, out)
        print(out)
        write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], "_top_",x[j],"_correlating_genes.txt", sep=""), sep="\t", col.names = T, row.names = F)
      } else {
        print("No disease genes mapped to PPI")
      }
      rm(out)
    }
  }
}


# remove everything when done
rm(list = ls())
