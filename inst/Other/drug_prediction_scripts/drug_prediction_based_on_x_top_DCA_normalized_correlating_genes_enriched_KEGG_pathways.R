# Drug prediction based on x top KEGG pathway enrichment with DCA normalized correlating genes

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

set.seed(35)

# x
x <- c(1,5,10,500) # n top pathways used in analysis; 500 deliberately bigger than n available pathways

# Load background genes (all genes measured with scRNAseq) - depending on that 
bg <- read.table("Input/Human-mouse_homologs/mouse_genes_from_scRNAseq.txt", sep="\t", header = F)
bg <- as.vector(bg[,1])

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

# Load KEGG pathways
kegg_pathways <- as.matrix(read.table(file = "Input/KEGG_pathways/TRANSLATED_kegg_genes_by_pathway.txt", sep="\t", header = T))

# Set up output directories
dir.create("Output/Correlation_genes_predictions/literature_PPI")
dir.create("Output/Correlation_genes_predictions/HuRI_PPI")
source(paste(fp, "/pathway_enrichment.R", sep=""))

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis")
out.dir <- paste(getwd(), "/Output/Correlation_genes_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = paste(fp,"/../Input/literature_PPI/ppi.txt",sep=""), sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])

# Analysis
save_name <- "DCA_normalized_correlating_gene_enriched_KEGG_pathways_in_literature_PPI_1000"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  # pathway enrichment
  sig_pathw <- pathway_enrichment_function(pathways = kegg_pathways, degs = degs[,i], bg = bg, cores = 16)
  
  if(!is.null(sig_pathw)){
    for(j in 1:length(x)){
      print(paste("CALCULATING: Top ", x[j], " pathways enriched with ", colnames(degs)[i], " correlating genes", sep=""))
      # drug prediction based on pathways
      if(x[j]< ncol(sig_pathw)){
        temp <- as.vector(sig_pathw[,1:x[j]])
        temp <- unique(temp[!is.na(temp)])
        save_name2 <- paste("_top_", x[j],"_significant_pathways",sep="")
      } else {
        temp <- as.vector(sig_pathw)
        temp <- unique(temp[!is.na(temp)])
        save_name2 <- "_all_significant_pathways"
      }
      out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(temp), disease_module_lcc = F, out.dir = out.dir,
                                               disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
      #rec <- T
      if(!(length(out)==1) & !all(is.na(out))){
        out <- precision_and_recall_from_z_score(predefined_list, out)
        print(out)
        write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], save_name2, ".txt", sep=""), sep="\t", col.names = T, row.names = F)
      } else {
        print("No disease genes")
      }
      rm(out, save_name2, temp)
    }
    rm(sig_pathw)  
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
save_name <- "DCA_normalized_correlating_gene_enriched_KEGG_pathways_in_literature_PPI_1000"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  # pathway enrichment
  sig_pathw <- pathway_enrichment_function(pathways = kegg_pathways, degs = degs[,i], bg = bg, cores = 16)
  
  if(!is.null(sig_pathw)){
    for(j in 1:length(x)){
      print(paste("CALCULATING: Top ", x[j], " pathways enriched with ", colnames(degs)[i], " correlating genes", sep=""))
      # drug prediction based on pathways
      if(x[j]< ncol(sig_pathw)){
        temp <- as.vector(sig_pathw[,1:x[j]])
        temp <- unique(temp[!is.na(temp)])
        save_name2 <- paste("_top_", x[j],"_significant_pathways",sep="")
      } else {
        temp <- as.vector(sig_pathw)
        temp <- unique(temp[!is.na(temp)])
        save_name2 <- "_all_significant_pathways"
      }
      out <- in_silico_drug_efficacy_screening(ppi = ppi, drug_matrix = drugs, disease_module = as.character(temp), disease_module_lcc = F, out.dir = out.dir,
                                               disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 16, recycle = rec, n_random_iterations = 1000)
      #rec <- T
      if(!(length(out)==1) & !all(is.na(out))){
        out <- precision_and_recall_from_z_score(predefined_list, out)
        print(out)
        write.table(out, file = paste(out.dir,"/precision_recall_matrix_", save_name, "_", colnames(degs)[i], save_name2, ".txt", sep=""), sep="\t", col.names = T, row.names = F)
      } else {
        print("No disease genes")
      }
      rm(out, save_name2, temp)
    }
    rm(sig_pathw)  
  }
}


# remove everything when done
rm(list = ls())
