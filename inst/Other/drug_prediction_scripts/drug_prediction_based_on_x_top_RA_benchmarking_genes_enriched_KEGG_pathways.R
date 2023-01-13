# Drug prediction based on x top KEGG pathway enrichment with benchmarking genes

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
set.seed(35)

# x
x <- c(1,5,10,500) # n top pathways used in analysis; 500 deliberately bigger than n available pathways

# Load background genes (all genes measured with scRNAseq) - depending on that 
bg <- read.table("Input/Human-mouse_homologs/mouse_genes_from_scRNAseq.txt", sep="\t", header = F)
bg <- as.vector(bg[,1])

# Load GWAS
gwas <- read.table(file = "Input/GWAScat/TRANSLATED_GWAS_RA_P<1e-8.txt", sep="\t", header = T)
gwas <- as.vector(gwas[,])
# Load OMIM
omim <- as.matrix(read.table(file = "Input/OMIM/OMIM.txt", sep="\t", header = T))
omim <- as.vector(omim)
# Load bulk RNA-seq data
bulk_synovium <- read.table("Input/Bulk_data/synovium_DEGs_GSE55235.txt", sep="\t", header = T)
bulk_synovium <- bulk_synovium[order(bulk_synovium[,6], decreasing = F),]
bulk_synovium <- bulk_synovium[bulk_synovium[,6]<0.05,1]

bulk_blood <- read.table("Input/Bulk_data/whole_blood_DEGs_GSE93272.txt", sep="\t", header = T)
bulk_blood <- bulk_blood[order(bulk_blood[,6], decreasing = F),]
bulk_blood <- bulk_blood[bulk_blood[,6]<0.05,1]

# Combine gene lists
degs <- cbind(c(gwas, rep(NA, times = 10000-length(gwas))), c(omim, rep(NA, times = 10000-length(omim))), c(gwas, omim, rep(NA, times = 10000-length(gwas)-length(omim))), 
              c(bulk_synovium, rep(NA, times = 10000-length(bulk_synovium))), c(bulk_blood, rep(NA, times = 10000-length(bulk_blood))))
colnames(degs) <- c("GWAS", "OMIM", "GWAS_&_OMIM", "Synovial_DEGs_bulk_RNA-seq", "Whole_blood_DEGs_bulk_RNA-seq")
degs <- degs[rowSums(!is.na(degs))>0,colSums(!is.na(degs))>0]
degs <- as.matrix(degs)
mode(degs) <- "character"
rm(gwas,omim, bulk_synovium, bulk_blood)

# Load KEGG pathways
kegg_pathways <- as.matrix(read.table(file = "Input/KEGG_pathways/TRANSLATED_kegg_genes_by_pathway.txt", sep="\t", header = T))

# Set up output directories
dir.create("Output/RA_benchmark_predictions")
dir.create("Output/RA_benchmark_predictions/literature_PPI")
dir.create("Output/RA_benchmark_predictions/HuRI_PPI")
source(paste(fp, "/pathway_enrichment.R", sep=""))

# start drug prediction literature PPI
##################################################################################################################
print("START literature PPI network analysis for benchmark enriched KEGG pathways")
out.dir <- paste(getwd(), "/Output/Other_predictions/literature_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = "Input/literature_PPI/ppi.txt", sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 16 unique drug combinations found in literature PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "benchmarking_genes_enriched_KEGG_pathways_in_literature_PPI"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  # pathway enrichment
  sig_pathw <- pathway_enrichment_function(pathways = kegg_pathways, degs = degs[,i], bg = bg, cores = 32)
  
  if(!is.null(sig_pathw)){
    for(j in 1:length(x)){
      print(paste("CALCULATING: Top ", x[j], " pathways enriched with ", colnames(degs)[i], sep=""))
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
print("START HuRI PPI network analysis for benchmark enriched KEGG pathways")
setwd(paste(fp, "/..", sep=""))
out.dir <- paste(getwd(), "/Output/Other_predictions/HuRI_PPI",sep="") # output folder for in_silico_drug_efficacy_screening()
source(paste(fp, "/in_silico_drug_efficacy_screening_algorithm.R",sep="")) # to remove "/temp" folder if exists

# LOAD PPI
ppi <- as.matrix(read.table(file = "Input/HuRI_PPI/ppi.txt", sep="\t", header = T))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])
# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_HuRI_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
# Load predefined drug list & format to vector
predefined_list <- as.vector(read.table("Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", header = T)[,1])
same_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_HuRI.txt", sep="\t", header = T, stringsAsFactors = T))
same_drugs <- same_drugs[same_drugs[,1]%in%predefined_list | same_drugs[,2]%in%predefined_list,]
predefined_list <- unique(c(same_drugs[,1], same_drugs[,2], predefined_list)) #29 drugs, 10 unique drug combinations found in HuRI PPI
predefined_list <- predefined_list[predefined_list %in% colnames(drugs)]
rm(same_drugs)

# Analysis
save_name <- "benchmarking_genes_enriched_KEGG_pathways_in_HuRI_PPI"
print(save_name)
rec <- F
for(i in 1:ncol(degs)){
  print(paste("CALCULATING: ", colnames(degs)[i], sep=""))
  # pathway enrichment
  sig_pathw <- pathway_enrichment_function(pathways = kegg_pathways, degs = degs[,i], bg = bg, cores = 16)
  
  if(!is.null(sig_pathw)){
    for(j in 1:length(x)){
      print(paste("CALCULATING: Top ", x[j], " pathways enriched with ", colnames(degs)[i], sep=""))
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
                                               disease_module_name = paste(save_name, "_", colnames(degs)[i],sep=""), cores = 20, recycle = rec, n_random_iterations = 1000)
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
