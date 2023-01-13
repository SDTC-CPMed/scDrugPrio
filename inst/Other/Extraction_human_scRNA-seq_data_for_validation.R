

########################################################################################################################################################################
# Extract DRESS data from h5 file
# *SKIN SAMPLES*
# PMID: 31959990 (NatMed Kim et al)
# GEO accession: GSE132802
########################################################################################################################################################################

'source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")'
library(rhdf5)
library(Seurat)
library(Matrix)
library(yarrr)
library(doParallel)
# Load my own h5 extrcation fucntion
source("Extract_h5_count_data_function.R")

# set wd
fp <- getwd()
setwd(paste(fp, "/..", sep=""))

# list objects in h5 file
h5ls(file = "Input/DRESS GSE132802/GSM3892569_Skin_Dress1_filtered_gene_bc_matrices_h5.h5")
# load h5 file
data <- h5read(file = "Input/DRESS GSE132802/GSM3892569_Skin_Dress1_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
# Apply function - returns sparse matrix
counts <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892569_Skin_DRESS")

gene_names <- cbind(data$gene_names, data$genes)
colnames(gene_names) <- c("Gene_symbol", "ENSEMBL_gene_ID")
write.table(gene_names, file = "Input/DRESS GSE132802/transl.txt", sep="\t", col.names = T, row.names = F)

########################################################################################
# Do this for all other samples as well and merge into one matrix
########################################################################################
# HV1 duplicate 1
data <- h5read(file = "Input/DRESS GSE132802/GSM3892577_SKIN_HV1_F1_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892577_SKIN_HV1_F1")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
# HV1 duplicate 2
data <- h5read(file = "Input/DRESS GSE132802/GSM3892578_SKIN_HV1_F2_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892578_SKIN_HV1_F2")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
# HV2
data <- h5read(file = "Input/DRESS GSE132802/GSM3892579_SKIN_HV2_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892579_SKIN_HV2")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
# HV3
data <- h5read(file = "Input/DRESS GSE132802/GSM3892580_SKIN_HV3_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892580_SKIN_HV3")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
# HV4
data <- h5read(file = "Input/DRESS GSE132802/GSM3892581_SKIN_HV4_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892581_SKIN_HV4")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
# HV5
data <- h5read(file = "Input/DRESS GSE132802/GSM3892582_SKIN_HV5_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892582_SKIN_HV5")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
rm(counts_1, data)

########################################################################################
# Save matrix as .txt
########################################################################################

#write.table(as.matrix(counts), file = "Input/DRESS GSE132802/Counts_for_skin_extracted.txt", sep="\t", col.names = NA, row.names = T)

########################################################################################
# Get baseline information
########################################################################################

temp <- as.matrix(counts)
n_exp_genes_per_cell <- colSums(temp>0)
n_reads_per_cell <- colSums(temp)
n_of_cells_in_which_gene_is_expressed <- rowSums(temp>0)
n_reads_per_gene <- rowSums(temp)
rm(temp)

samples <- c("DRESS", "HV1", "HV2", "HV3", "HV4", "HV5")
for(i in 1:length(samples)){
  # n cells per sample
  print(paste("n cells in:", samples[i], sum(grepl(pattern = samples[i], colnames(counts))), sep=" "))
  
  temp <- counts[,grepl(pattern = samples[i],colnames(counts))]
  
  # min n genes per sample
  out <- min(colSums(temp>0))
  print(paste("min genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- max(colSums(temp>0))
  print(paste("max genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- median(colSums(temp>0))
  print(paste("median genes in:", samples[i], out, sep=" "))
  
  # min n reads per sample
  out <- min(colSums(temp))
  print(paste("min reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- max(colSums(temp))
  print(paste("max reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- median(colSums(temp))
  print(paste("median reads in:", samples[i], out, sep=" "))
}


########################################################################################
# Quality control
########################################################################################

# Exclude genes that are expressed in less than 0.1% of cells
print(paste("n genes before filtering: ", nrow(counts), sep=""))
counts <- counts[n_of_cells_in_which_gene_is_expressed>(ncol(counts)*0.001),]
print(paste("n genes after filtering: ", nrow(counts), sep=""))

# recalculate
temp <- as.matrix(counts)
n_exp_genes_per_cell <- colSums(temp>0)
n_reads_per_cell <- colSums(temp)
n_of_cells_in_which_gene_is_expressed <- rowSums(temp>0)
n_reads_per_gene <- rowSums(temp)
rm(temp)

# Check for mitochondrial genes - labeled with "^MT-"
temp <- gene_names[grepl("^MT-" ,gene_names[,1]),2]
temp <- counts[rownames(counts)%in%temp,]
fraction_mitochondiral_genes_per_cell <- colSums(as.matrix(temp))/n_reads_per_cell
rm(temp)

temp <- cbind(fraction_mitochondiral_genes_per_cell, n_reads_per_cell, n_exp_genes_per_cell, colnames(counts))
temp1 <- temp[as.numeric(temp[,1])>0.20,] # cells with > 20% mitochondiral genes, excludes n = 303
temp <- temp[!as.numeric(temp[,1])>0.20,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,2])< 4000,]) # cells with less than 4000 reads, excludes 1452
temp <- temp[!as.numeric(temp[,2])<4000,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])<500,]) # cells with less than 500 genes, excludes n = 0
temp <- temp[!as.numeric(temp[,3])<500,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])>6000,]) # cells with more than 6000 genes, excludes n = 126
temp <- temp[!as.numeric(temp[,3])>6000,]

max <- (trunc(max(n_reads_per_cell)/1000)+1)*1000
t <- seq(1,10, by = 1)
labels <- vector()
for(i in 3:trunc(log10(max))){
  labels <- c(labels, t*10^i)
}
rm(t)

plot(x = 1, log = "x", type = "n", xlim = c(2000, max), ylim = c(0,1), pch = 19, 
     xlab = "n reads per cell", ylab = "fraction mitochondrial genes")
axis(side = 1, at = labels, labels = F)
points(pch=19, col = "black", x = as.numeric(temp[,2]), y = as.numeric(temp[,1]), cex = 0.3) # valid cells
points(pch = 19, col = "grey", x = as.numeric(temp1[,2]), y = as.numeric(temp1[,1]), cex = 0.3) # excluded cells

########################################################################################
# Apply Quality Criteria to data set & save
########################################################################################

counts <- counts[,colnames(counts)%in%temp[,4]]
write.table(as.matrix(counts), file = "Input/DRESS GSE132802/FILTERED_counts_for_skin_DRESS.txt", col.names = NA, row.names = T)

########################################################################################
# Get filtered information
########################################################################################

for(i in 1:length(samples)){
  # n cells per sample
  print(paste("n cells in:", samples[i], sum(grepl(pattern = samples[i], colnames(counts))), sep=" "))
  
  temp <- counts[,grepl(pattern = samples[i],colnames(counts))]
  
  # min n genes per sample
  out <- min(colSums(temp>0))
  print(paste("min genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- max(colSums(temp>0))
  print(paste("max genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- median(colSums(temp>0))
  print(paste("median genes in:", samples[i], out, sep=" "))
  
  # min n reads per sample
  out <- min(colSums(temp))
  print(paste("min reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- max(colSums(temp))
  print(paste("max reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- median(colSums(temp))
  print(paste("median reads in:", samples[i], out, sep=" "))
}

########################################################################################################################################################################
# Extract DRESS data from h5 file
# *PBMC SAMPLES*
# PMID: 31959990 (NatMed Kim et al)
# GEO accession: GSE132802
########################################################################################################################################################################

# list objects in h5 file
h5ls(file = "Input/DRESS GSE132802/GSM3892570_PBMC_DRESS1_filtered_gene_bc_matrices_h5.h5")
# DRESS PBMCs
data <- h5read(file = "Input/DRESS GSE132802/GSM3892570_PBMC_DRESS1_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892570_PBMC_DRESS")
# HV
data <- h5read(file = "Input/DRESS GSE132802/GSM3892571_PBMC_HV_filtered_gene_bc_matrices_h5.h5", name = "/GRCh38")
counts_1 <- extract_h5_count_data(my_h5_file = data, sample_name = "GSM3892571_PBMC_HV")
counts <- cbind(counts, counts_1[match(x = rownames(counts_1), table = rownames(counts)),])
rm(counts_1)

########################################################################################
# Get baseline information
########################################################################################

temp <- as.matrix(counts)
n_exp_genes_per_cell <- colSums(temp>0)
n_reads_per_cell <- colSums(temp)
n_of_cells_in_which_gene_is_expressed <- rowSums(temp>0)
n_reads_per_gene <- rowSums(temp)
rm(temp)

samples <- c("DRESS", "HV")
for(i in 1:length(samples)){
  # n cells per sample
  print(paste("n cells in:", samples[i], sum(grepl(pattern = samples[i], colnames(counts))), sep=" "))
  
  temp <- counts[,grepl(pattern = samples[i],colnames(counts))]
  
  # min n genes per sample
  out <- min(colSums(temp>0))
  print(paste("min genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- max(colSums(temp>0))
  print(paste("max genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- median(colSums(temp>0))
  print(paste("median genes in:", samples[i], out, sep=" "))
  
  # min n reads per sample
  out <- min(colSums(temp))
  print(paste("min reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- max(colSums(temp))
  print(paste("max reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- median(colSums(temp))
  print(paste("median reads in:", samples[i], out, sep=" "))
}


########################################################################################
# Quality control
########################################################################################

# Exclude genes that are expressed in less than 0.1% of cells
print(paste("n genes before filtering: ", nrow(counts), sep=""))
counts <- counts[n_of_cells_in_which_gene_is_expressed>(ncol(counts)*0.001),]
print(paste("n genes after filtering: ", nrow(counts), sep=""))

# recalculate
temp <- as.matrix(counts)
n_exp_genes_per_cell <- colSums(temp>0)
n_reads_per_cell <- colSums(temp)
n_of_cells_in_which_gene_is_expressed <- rowSums(temp>0)
n_reads_per_gene <- rowSums(temp)
rm(temp)

# Check for mitochondrial genes - labeled with "^MT-"
temp <- gene_names[grepl("^MT-" ,gene_names[,1]),2]
temp <- counts[rownames(counts)%in%temp,]
fraction_mitochondiral_genes_per_cell <- colSums(as.matrix(temp))/n_reads_per_cell
rm(temp)

temp <- cbind(fraction_mitochondiral_genes_per_cell, n_reads_per_cell, n_exp_genes_per_cell, colnames(counts))
temp1 <- temp[as.numeric(temp[,1])>0.20,] # cells with > 20% mitochondiral genes, excludes n = 303
temp <- temp[!as.numeric(temp[,1])>0.20,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,2])< 4000,]) # cells with less than 4000 reads, excludes 1452
temp <- temp[!as.numeric(temp[,2])<4000,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])<500,]) # cells with less than 500 genes, excludes n = 0
temp <- temp[!as.numeric(temp[,3])<500,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])>6000,]) # cells with more than 6000 genes, excludes n = 126
temp <- temp[!as.numeric(temp[,3])>6000,]

max <- (trunc(max(n_reads_per_cell)/1000)+1)*1000
t <- seq(1,10, by = 1)
labels <- vector()
for(i in 3:trunc(log10(max))){
  labels <- c(labels, t*10^i)
}
rm(t)

plot(x = 1, log = "x", type = "n", xlim = c(2000, max), ylim = c(0,1), pch = 19, 
     xlab = "n reads per cell", ylab = "fraction mitochondrial genes")
axis(side = 1, at = labels, labels = F)
points(pch=19, col = "black", x = as.numeric(temp[,2]), y = as.numeric(temp[,1]), cex = 0.3) # valid cells
points(pch = 19, col = "grey", x = as.numeric(temp1[,2]), y = as.numeric(temp1[,1]), cex = 0.3) # excluded cells

########################################################################################
# Apply Quality Criteria to data set & save
########################################################################################

counts <- counts[,colnames(counts)%in%temp[,4]]
write.table(as.matrix(counts), file = "Input/DRESS GSE132802/FILTERED_counts_for_PBMC_DRESS.txt", col.names = NA, row.names = T)

########################################################################################
# Get filtered information
########################################################################################

for(i in 1:length(samples)){
  # n cells per sample
  print(paste("n cells in:", samples[i], sum(grepl(pattern = samples[i], colnames(counts))), sep=" "))
  
  temp <- counts[,grepl(pattern = samples[i],colnames(counts))]
  
  # min n genes per sample
  out <- min(colSums(temp>0))
  print(paste("min genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- max(colSums(temp>0))
  print(paste("max genes in:", samples[i], out, sep=" "))
  # max n genes per sample
  out <- median(colSums(temp>0))
  print(paste("median genes in:", samples[i], out, sep=" "))
  
  # min n reads per sample
  out <- min(colSums(temp))
  print(paste("min reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- max(colSums(temp))
  print(paste("max reads in:", samples[i], out, sep=" "))
  # max n reads per sample
  out <- median(colSums(temp))
  print(paste("median reads in:", samples[i], out, sep=" "))
}



########################################################################################################################################################################
# Extract Crohn's Disease data
# * ILEUM SAMPLES*
# PMID: 31474370 (Cell, Martin et al)
# GEO accession: GSE134809
########################################################################################################################################################################

# detect all samples
ls <- list.files(path = "Input/CD GSE134809/GSE134809_RAW/")
ls <- ls[!grepl(pattern = ".gz", ls)]
ls_matrix <- ls[grepl(pattern = "_matrix.mtx", ls)]
ls_barcodes <- ls[grepl(pattern = "_barcodes.tsv", ls)]
ls_genes <- ls[grepl(pattern = "_genes.tsv", ls)]

# only relevant samples
r_samples <- c("GSM3972009", "GSM3972010", "GSM3972011", "GSM3972012", "GSM3972013", "GSM3972014", "GSM3972015", "GSM3972016", "GSM3972017", "GSM3972018", "GSM3972019",
               "GSM3972020", "GSM3972021", "GSM3972022", "GSM3972023", "GSM3972024", "GSM3972025", "GSM3972026", "GSM3972027", "GSM3972028", "GSM3972029", "GSM3972030") 
r_samples_type <- c("Involved", "Uninvolved","Involved", "Uninvolved","Involved", "Uninvolved","Uninvolved", "Involved", "Involved", "Uninvolved", "Uninvolved", "Involved",
                    "Uninvolved", "Involved", "Uninvolved", "Involved", "Uninvolved", "Involved", "Uninvolved", "Involved", "Uninvolved", "Involved")
samples <- cbind(r_samples, r_samples_type)

# Extract data from all included samples
counts <- foreach(i = 1:nrow(samples), .combine = "cbind") %do% {
  # Read from SparseMatrix format (MatrixMarket)
  data <- readMM(file = paste("Input/CD GSE134809/GSE134809_RAW/", ls_matrix[grepl(pattern = samples[i,1], ls_matrix)], sep=""))
  barcode <- as.matrix(read.table(file = paste("Input/CD GSE134809/GSE134809_RAW/", ls_barcodes[grepl(pattern = samples[i,1], ls_barcodes)], sep="")))
  barcode <- paste(as.character(samples[i,1]), "_", as.character(samples[i,2]), "_", as.character(barcode[,1]), sep="")
  genes <- as.matrix(read.table(file = paste("Input/CD GSE134809/GSE134809_RAW/", ls_genes[grepl(pattern = samples[i,1], ls_genes)], sep="")))
  rownames(data) <- genes[,1]
  colnames(data) <- barcode
  # Exclude obviously empty wells
  data <- data[,colSums(data)>0]
  return(data)
}

rm(ls, ls_barcodes, ls_matrix, ls_genes, i, barcode, data, r_samples, r_samples_type)

########################################################################################
# Get baseline information
########################################################################################

n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

out <- matrix(NA, nrow = nrow(samples), ncol = 8)
for(i in 1:nrow(samples)){
  # n cells per sample
  out[i,1] <- samples[i,1] # sample ID
  out[i,2] <- sum(grepl(pattern = samples[i,1], colnames(counts))) # n cells
  
  temp <- counts[,grepl(pattern = samples[i,1],colnames(counts))]
  
  # min n genes per sample
  out[i,3] <- min(colSums(temp>0))
  # max n genes per sample
  out[i,4] <- max(colSums(temp>0))
  # max n genes per sample
  out[i,5] <- median(colSums(temp>0))
  
  # min n reads per sample
  out[i,6] <- min(colSums(temp))
  # max n reads per sample
  out[i,7] <- max(colSums(temp))
  # max n reads per sample
  out[i,8] <- median(colSums(temp))
}
colnames(out) <- c("sampleID", "n_cells", "nGene_min", "nGene_max", "nGene_median", "nReads_min", "nReads_max", "nReads_median")
print(out)
write.table(out, file = "Input/CD GSE134809/Unfiltered_quality_statistics.txt", sep="\t", col.names = T, row.names = F)

########################################################################################
# Quality control
########################################################################################

# Exclude all cells with less than 1000 reads
temp <- counts[, colSums(counts)>=1000]

# Exclude genes that are expressed in less than 0.1% of cells that have >= 1000 reads
print(paste("n genes before filtering: ", nrow(counts), sep=""))
counts <- counts[rowSums(temp>0)>(ncol(temp)*0.001),]
print(paste("n genes after filtering: ", nrow(counts), sep=""))
rm(temp)

# recalculate
n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

# Check for mitochondrial genes - labeled with "^MT-"
temp <- genes[grepl("^MT-" ,genes[,2]),1]
temp <- counts[rownames(counts)%in%temp,]
fraction_mitochondiral_genes_per_cell <- colSums(as.matrix(temp))/n_reads_per_cell
rm(temp)

temp <- cbind(fraction_mitochondiral_genes_per_cell, n_reads_per_cell, n_exp_genes_per_cell, colnames(counts))
temp1 <- temp[as.numeric(temp[,2])< 1000,] # cells with less than 1000 reads, excludes n = 4323634
temp <- temp[!as.numeric(temp[,2])<1000,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,1])>0.20,]) # cells with > 20% mitochondiral genes, excludes n = 13396
temp <- temp[!as.numeric(temp[,1])>0.20,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])<250,]) # cells with less than 250 genes, excludes n = 493
temp <- temp[!as.numeric(temp[,3])<250,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])>6000,]) # cells with more than 6000 genes, excludes n = 21
temp <- temp[!as.numeric(temp[,3])>6000,]

max <- (trunc(max(n_reads_per_cell)/1000)+1)*1000
t <- seq(1,10, by = 1)
labels <- vector()
for(i in 3:trunc(log10(max))){
  labels <- c(labels, t*10^i)
}
rm(t)

plot(x = 1, log = "x", type = "n", xlim = c(1, max), ylim = c(0,1), pch = 19, 
     xlab = "n reads per cell", ylab = "fraction mitochondrial genes")
axis(side = 1, at = labels, labels = F)
points(pch=19, col = "black", x = as.numeric(temp[,2]), y = as.numeric(temp[,1]), cex = 0.2) # valid cells
points(pch = 19, col = "grey", x = as.numeric(temp1[,2]), y = as.numeric(temp1[,1]), cex = 0.2) # excluded cells

########################################################################################
# Apply Quality Criteria to data set & save
########################################################################################

counts <- counts[,colnames(counts)%in%temp[,4]]
write.table(as.matrix(counts), file = "Input/CD GSE134809/FILTERED_counts_for_CD_ileum_biopsies.txt", col.names = NA, row.names = T)

########################################################################################
# Get filtered information
########################################################################################

n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

out <- matrix(NA, nrow = nrow(samples), ncol = 8)
for(i in 1:nrow(samples)){
  # n cells per sample
  out[i,1] <- samples[i,1] # sample ID
  out[i,2] <- sum(grepl(pattern = samples[i,1], colnames(counts))) # n cells
  
  temp <- counts[,grepl(pattern = samples[i,1],colnames(counts))]
  
  # min n genes per sample
  out[i,3] <- min(colSums(temp>0))
  # max n genes per sample
  out[i,4] <- max(colSums(temp>0))
  # max n genes per sample
  out[i,5] <- median(colSums(temp>0))
  
  # min n reads per sample
  out[i,6] <- min(colSums(temp))
  # max n reads per sample
  out[i,7] <- max(colSums(temp))
  # max n reads per sample
  out[i,8] <- median(colSums(temp))
}
colnames(out) <- c("sampleID", "n_cells", "nGene_min", "nGene_max", "nGene_median", "nReads_min", "nReads_max", "nReads_median")
print(out)
write.table(out, file = "Input/CD GSE134809/Filtered_quality_statistics.txt", sep="\t", col.names = T, row.names = F)



########################################################################################################################################################################
# Extract Ulcerative Colitis data
# *COLON SAMPLES*
# PMID: 30814735 (Nature, Parikh et al)
# GEO accession: GSE116222
########################################################################################################################################################################

# Read in data
counts <- read.table(file = "Input/UC GSE116222/GSE116222_Expression_matrix.txt", sep="\t")
# data is formatted in tab-delimited UMI counts per 10,000 unique molecules detected according to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3214201

# Multiplication with 10 to avoid data loss when rounding the data
counts <- counts*10
counts <- counts+0.01
counts <- round(counts)

# Format sample names
temp <- colnames(counts)
sample_name <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
sample_group <- c("Healthy", "Non_inflamed", "Inflammed", "Healthy", "Non_inflamed", "Inflammed", "Healthy", "Non_inflamed", "Inflammed")
samples <- cbind(sample_name, sample_group)
rm(sample_group, sample_name)
for(i in 1:nrow(samples)){
  temp[grepl(pattern = paste(".", samples[i,1],sep=""), temp)] <- paste(temp[grepl(pattern = paste(".", samples[i,1],sep=""), temp)], samples[i,2], sep="_")
}
colnames(counts) <- temp

########################################################################################
# Get baseline information
########################################################################################

n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

out <- matrix(NA, nrow = nrow(samples), ncol = 8)
for(i in 1:nrow(samples)){
  # n cells per sample
  out[i,1] <- samples[i,1] # sample ID
  out[i,2] <- sum(grepl(pattern = samples[i,1], colnames(counts))) # n cells
  
  temp <- counts[,grepl(pattern = samples[i,1],colnames(counts))]
  
  # min n genes per sample
  out[i,3] <- min(colSums(temp>0))
  # max n genes per sample
  out[i,4] <- max(colSums(temp>0))
  # max n genes per sample
  out[i,5] <- median(colSums(temp>0))
  
  # min n reads per sample
  out[i,6] <- min(colSums(temp))
  # max n reads per sample
  out[i,7] <- max(colSums(temp))
  # max n reads per sample
  out[i,8] <- median(colSums(temp))
}
colnames(out) <- c("sampleID", "n_cells", "nGene_min", "nGene_max", "nGene_median", "nReads_min", "nReads_max", "nReads_median")
print(out)
write.table(out, file = "Input/UC GSE116222/Unfiltered_quality_statistics.txt", sep="\t", col.names = T, row.names = F)

########################################################################################
# Quality control
########################################################################################

# Exclude genes that are expressed in less than 0.1% of cells
print(paste("n genes before filtering: ", nrow(counts), sep=""))
counts <- counts[rowSums(counts>0)>(ncol(counts)*0.001),]
print(paste("n genes after filtering: ", nrow(counts), sep=""))

# recalculate
n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

# Check for mitochondrial genes - labeled with "^MT-"
temp <- rownames(counts)[grepl("^MT-" ,rownames(counts))]
temp <- counts[rownames(counts)%in%temp,]
fraction_mitochondiral_genes_per_cell <- colSums(as.matrix(temp))/n_reads_per_cell
rm(temp)
# min = 0.0002559181; max = 0.1385281

temp <- cbind(fraction_mitochondiral_genes_per_cell, n_reads_per_cell, n_exp_genes_per_cell, colnames(counts))
temp1 <- temp[as.numeric(temp[,3])<500,] # cells with less than 500 genes, excludes n = 50
temp <- temp[!as.numeric(temp[,3])<500,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])>6000,]) # cells with more than 6000 genes, excludes n = 0
temp <- temp[!as.numeric(temp[,3])>6000,]

max <- (trunc(max(n_reads_per_cell)/1000)+1)*1000
t <- seq(1,10, by = 1)
labels <- vector()
for(i in 3:trunc(log10(max))){
  labels <- c(labels, t*10^i)
}
rm(t)

plot(x = 1, log = "x", type = "n", xlim = c(1, max), ylim = c(0,1), pch = 19, 
     xlab = "n reads per cell", ylab = "fraction mitochondrial genes")
axis(side = 1, at = labels, labels = F)
points(pch=19, col = "black", x = as.numeric(temp[,2]), y = as.numeric(temp[,1]), cex = 0.2) # valid cells
points(pch = 19, col = "grey", x = as.numeric(temp1[,2]), y = as.numeric(temp1[,1]), cex = 0.2) # excluded cells

########################################################################################
# Apply Quality Criteria to data set & save
########################################################################################

counts <- counts[,colnames(counts)%in%temp[,4]]
write.table(counts, file = "Input/UC GSE116222/FILTERED_counts_for_UC_colon_biopsies.txt", sep="\t", col.names = T, row.names = T)

########################################################################################
# Get filtered information
########################################################################################

n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

out <- matrix(NA, nrow = nrow(samples), ncol = 8)
for(i in 1:nrow(samples)){
  # n cells per sample
  out[i,1] <- samples[i,1] # sample ID
  out[i,2] <- sum(grepl(pattern = samples[i,1], colnames(counts))) # n cells
  
  temp <- counts[,grepl(pattern = samples[i,1],colnames(counts))]
  
  # min n genes per sample
  out[i,3] <- min(colSums(temp>0))
  # max n genes per sample
  out[i,4] <- max(colSums(temp>0))
  # max n genes per sample
  out[i,5] <- median(colSums(temp>0))
  
  # min n reads per sample
  out[i,6] <- min(colSums(temp))
  # max n reads per sample
  out[i,7] <- max(colSums(temp))
  # max n reads per sample
  out[i,8] <- median(colSums(temp))
}
colnames(out) <- c("sampleID", "n_cells", "nGene_min", "nGene_max", "nGene_median", "nReads_min", "nReads_max", "nReads_median")
print(out)
write.table(out, file = "Input/UC GSE116222/Filtered_quality_statistics.txt", sep="\t", col.names = T, row.names = F)

########################################################################################
# Translate genes to Ensemble for Olegs method (needs to be excluded from code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
########################################################################################
'
# First try using NCBI

# human tax-id: 9606
info <- read.table(file = "Input/NCBI_annotation/gene_info", header = F, sep="\t", fill = T)
colnames(info) <- c("tax_id",	"GeneID",	"Symbol",	"LocusTag",	"Synonyms",	"dbXrefs",
                    "chromosome",	"map_location",	"description",	"type_of_gene",
                    "Symbol_from_nomenclature_authority",	"Full_name_from_nomenclature_authority",
                    "Nomenclature_status",	"Other_designations",	"Modification_date",
                    "Feature_type")
info <- info[info[,1]=="9606",]
info <- info[,1:3]
write.table(as.matrix(info), file = "Input/NCBI_annotation/Extracted_human_genes_from_geneinfo.txt", sep="\t", col.names = T, row.names = F)

transl <- merge(ensembl, info, by= "GeneID")
transl <- transl[,c(1,3,5)]
temp <- paste(transl[,1], transl[,2], transl[,3], sep="_")
transl <- transl[!duplicated(temp),]
rm(ensembl, info, temp)

write.table(as.matrix(transl), file = "Input/NCBI_annotation/transl.txt", sep="\t", col.names = T, row.names = F)
'
'
# Second try using NCBI

# human tax-id: 9606
# remove "#" in front of "taxid" in the file
ensembl <- read.table(file = "Input/NCBI_annotation/gene2ensembl", header = T, sep = "\t")
ensembl <- ensembl[ensembl[,1]=="9606",]
ensembl <- ensembl[,1:3]
temp <- paste(ensembl[,1], ensembl[,2], ensembl[,3], sep="_")
ensembl <- ensembl[!duplicated(temp),]
rm(temp)
write.table(as.matrix(ensembl), file = "Input/NCBI_annotation/Extracted_human_genes_from_gene2ensembl.txt", sep="\t", col.names = T, row.names = F)

genes <- rownames(counts)
print(length(genes)) # n = 18028
transl <- read.table(file = "Input/Human-mouse_homologs/transl.txt", sep="\t", header = T)
print(sum(genes %in% transl[,1])) # n = 10936

# merge transl with ensemble ID
print(length(unique((ensembl[,2])))) # n = 35026
print(sum(ensembl[,2]%in%transl[,2])) # n = 12711

transl <- merge(x = transl, y = ensembl, by = 2)
transl <- transl[,-5]
transl <- transl[!duplicated(transl[,2]),]
transl <- transl[!duplicated(transl[,5]),]
write.table(as.matrix(transl), file = "Input/NCBI_annotation/transl.txt", sep="\t", col.names = T, row.names = F)

'
# Using HGNC data base translation matrix
transl <- read.table(file = "Input/HGNC translation matrix 201108/downloaded_201108.txt", header = T, sep = "\t", fill = T, quote = "")
transl[transl[,7]=="", 7] <- NA
transl <- transl[- which(is.na(transl[,6]) & is.na(transl[,7])),] # exclude entries which do not match with a entrez ID nor a ensembl ID 
transl <- as.matrix(transl)
print(paste("Duplicated Entrez IDs? ", sum(duplicated(transl[,6])), sep=""))
print(paste("Duplicated Ensemble IDs? ", sum(duplicated(transl[,7])), sep=""))
print(transl[transl[,7] %in% transl[duplicated(transl[,7]),7],])
transl <- transl[transl[,2]!="SUGT1P4-STRA6LP",]
write.table(transl, file = "Input/HGNC translation matrix 201108/transl.txt", sep="\t", col.names = T, row.names = F)


print(sum(rownames(counts)%in%transl[,2]))
count_translation <- cbind(rownames(counts), NA)
colnames(count_translation) <- c("Annotation_data_set", "Official_HGNC_symbol_translation_matrix")
count_translation[,2] <- transl[match(count_translation[,1], transl[,2]),2]




#############################################
# checking if some symbols match to previous symbols or aliases
#############################################
print(colnames(transl))
other_symbols <- rbind(cbind(transl[,4], "Previous_symbol", transl[,c(2,6:7)]), cbind(transl[,5], "Alias_symbol", transl[,c(2,6:7)]))
other_symbols <- other_symbols[other_symbols[,1]!="",] # exclude entries for which no aliases/previous symbols exist
# split comma seperated aliases/previous symbols and add them as seperate rows
out <- vector()
pos <- vector()
for(i in 1:nrow(other_symbols)){
  temp <- other_symbols[i,1]
  if(grepl(pattern = ", ", x = temp)){
    temp <- unlist(strsplit(temp, split = ", "))
  }
  out <- c(out, temp)
  pos <- c(pos, rep(i, times = length(temp)))
}
other_symbols <- cbind(out, other_symbols[pos,2:5])
rm(out, pos)
colnames(other_symbols)[1:2] <- c("Other_symbol", "Type_of_symbol")
write.table(other_symbols, file = "Input/HGNC translation matrix 201108/other_gene_symbols_w_duplicates.txt", sep="\t", col.names = T, row.names = F)
# check if any additional symbols can be translated
temp <- count_translation[is.na(count_translation[,2]),1]
print(paste(length(temp)," gene symbols do not match to official HGNC symbol", sep=""))
print(paste(sum(temp%in%other_symbols[,1]), " gene symbols match to aliases/previous gene symbols", sep=""))

# handle duplicated mappings
#############################
other_symbols <- other_symbols[other_symbols[,1]%in% temp,]
# Exclude all repeats 
other_symbols <- other_symbols[!duplicated(paste(other_symbols[,1], "_", other_symbols[,3])),]
# one official gene symbol matching to several aliases/previous symbols?
# all gene mappings that mapped one alias towards several different official gene symbols were exlcuded (n = 30)
temp <- other_symbols[other_symbols[,1] %in% other_symbols[duplicated(other_symbols[,1]),1],]
print(temp)
temp <- cbind(temp, paste(temp[,1],"_", temp[,3], sep=""))
temp <- temp[! temp[,6]%in%temp[duplicated(temp[,6]),6],-6]
other_symbols <- other_symbols[!other_symbols[,1]%in%temp[,1],]
# matching to official gene symbol that is already included? (n = 10)
# if a gene symbol in data set directly matched with a official gene symbol no unmatched gene symbol might be translated to the same official gene symbol
other_symbols <- other_symbols[!other_symbols[,3]%in%count_translation[,2],]
# Quality control: Only one alias per official gene symbol? 
print(paste("Only one alias per official gene symbol? ", !any(any(duplicated(other_symbols[,3])), any(duplicated(other_symbols[,1]))), sep=""))

other_symbols <- other_symbols[match(count_translation[,1], other_symbols[,1]),]
count_translation[!is.na(other_symbols[,1]),2] <- other_symbols[!is.na(other_symbols[,1]),3]
count_translation <- cbind(count_translation, transl[match(count_translation[,2], transl[,2]),6:7])
write.table(count_translation, file = "Input/HGNC translation matrix 201108/transl_matrix_including_aliases_for_translation_of_UC_GSE116222_specifically.txt", sep="\t", col.names = T, row.names = F)
write.table(count_translation, file = "Input/UC GSE116222/transl_matrix_including_aliases_for_translation_of_UC_GSE116222_specifically.txt", sep="\t", col.names = T, row.names = F)

count_translation <- count_translation[!is.na(count_translation[,4]),]

# translate genes from scRNA-seq data set to ENSG
print(nrow(counts)) # n = 18028
counts <- counts[rownames(counts)%in%count_translation[,1],]
print(nrow(counts)) # n = 15297
rownames(counts) <- count_translation[match(rownames(counts), count_translation[,1]),4]

# Save
write.table(counts, file = "Input/UC GSE116222/FILTERED_&_TRANSLATED_counts_for_UC_colon_biopsies.txt", sep="\t", col.names = T, row.names = T)


########################################################################################################################################################################
# Extract Ulcerative Colitis data --------------------------------------------------------- NOT USED ---------------------------------------------------------
# *COLON SAMPLES*
# PMID: 32826341 (Science Immunology, Boland et al)
# GEO accession: GSE125527
########################################################################################################################################################################
'
# Load relevant files
data <- read.csv(file = "Input/UC GSE125527/GSE125527_UMI_cell_table_sparse.csv", header= T)
genes <- as.vector(read.csv(file = "Input/UC GSE125527/GSE125527_gene_id_rownames.csv", header = F)[,1])
cells <- as.vector(read.csv(file = "Input/UC GSE125527/GSE125527_cell_id_colnames.csv", header = F)[,1])
# Extract matrix from spare matrix format
counts <- sparseMatrix(i = data[,1], j = data[,2], x = data[,3], dims = c(length(genes), length(cells)), dimnames = list(genes, cells))


transl <- read.table(file = "Input/NCBI_annotation/transl.txt", sep="\t", header = T)
print(sum(genes %in% transl[,3])) #n = 4853
transl <- read.table(file = "Input/NCBI_annotation/Extracted_human_genes_from_geneinfo.txt", sep="\t", header = T)
print(sum(genes %in% transl[,3])) #n = 4878
transl <- read.table(file = "Input/Human-mouse_homologs/transl.txt", sep="\t", header = T)
print(sum(genes %in% transl[,1])) #n = 8600
'

########################################################################################################################################################################
# Extract Multiple sclerosis data
# * CSF SAMPLES *
# PMID: 31937773 (NatCom, Schafflick et al)
# GEO accession: GSE138266
########################################################################################################################################################################

# detect all samples
ls <- list.files(path = "Input/MS GSE138266/GSE138266_RAW/")
ls <- ls[!grepl(pattern = ".gz", ls)]
ls_matrix <- ls[grepl(pattern = "_matrix.mtx", ls)]
ls_barcodes <- ls[grepl(pattern = "_barcodes.tsv", ls)]
ls_genes <- ls[grepl(pattern = "_genes.tsv", ls)]

# only relevant samples
r_samples <- c("GSM4104122","GSM4104123","GSM4104124","GSM4104125","GSM4104126","GSM4104127","GSM4104128","GSM4104129","GSM4104130","GSM4104131","GSM4104132","GSM4104133") 
r_samples_type <- c(rep("MS", times = 6), rep("IIH", times = 6))
samples <- cbind(r_samples, r_samples_type)

# Extract data from all included samples
counts <- foreach(i = 1:nrow(samples), .combine = "cbind") %do% {
  # Read from SparseMatrix format (MatrixMarket)
  data <- readMM(file = paste("Input/MS GSE138266/GSE138266_RAW/", ls_matrix[grepl(pattern = samples[i,1], ls_matrix)], sep=""))
  barcode <- as.matrix(read.table(file = paste("Input/MS GSE138266/GSE138266_RAW/", ls_barcodes[grepl(pattern = samples[i,1], ls_barcodes)], sep="")))
  barcode <- paste(as.character(samples[i,1]), "_", as.character(samples[i,2]), "_", as.character(barcode[,1]), sep="")
  genes <- as.matrix(read.table(file = paste("Input/MS GSE138266/GSE138266_RAW/", ls_genes[grepl(pattern = samples[i,1], ls_genes)], sep="")))
  rownames(data) <- genes[,1]
  colnames(data) <- barcode
  # Exclude obviously empty wells
  data <- data[,colSums(data)>0]
  return(data)
}

rm(ls, ls_barcodes, ls_matrix, ls_genes, i, barcode, data, r_samples, r_samples_type)

########################################################################################
# Get baseline information
########################################################################################

n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

out <- matrix(NA, nrow = nrow(samples), ncol = 8)
for(i in 1:nrow(samples)){
  # n cells per sample
  out[i,1] <- samples[i,1] # sample ID
  out[i,2] <- sum(grepl(pattern = samples[i,1], colnames(counts))) # n cells
  
  temp <- counts[,grepl(pattern = samples[i,1],colnames(counts))]
  
  # min n genes per sample
  out[i,3] <- min(colSums(temp>0))
  # max n genes per sample
  out[i,4] <- max(colSums(temp>0))
  # max n genes per sample
  out[i,5] <- median(colSums(temp>0))
  
  # min n reads per sample
  out[i,6] <- min(colSums(temp))
  # max n reads per sample
  out[i,7] <- max(colSums(temp))
  # max n reads per sample
  out[i,8] <- median(colSums(temp))
}
colnames(out) <- c("sampleID", "n_cells", "nGene_min", "nGene_max", "nGene_median", "nReads_min", "nReads_max", "nReads_median")
print(out)
write.table(out, file = "Input/MS GSE138266/Unfiltered_quality_statistics.txt", sep="\t", col.names = T, row.names = F)

########################################################################################
# Quality control
########################################################################################

# Exclude all cells with less than 1000 reads
temp <- counts[, colSums(counts)>=1000]

# Exclude genes that are expressed in less than 0.1% of cells that have >= 1000 reads
print(paste("n genes before filtering: ", nrow(counts), sep=""))
counts <- counts[rowSums(temp>0)>(ncol(temp)*0.001),]
print(paste("n genes after filtering: ", nrow(counts), sep=""))
rm(temp)

# recalculate
n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

# Check for mitochondrial genes - labeled with "^MT-"
temp <- genes[grepl("^MT-" ,genes[,2]),1]
temp <- counts[rownames(counts)%in%temp,]
fraction_mitochondiral_genes_per_cell <- colSums(as.matrix(temp))/n_reads_per_cell
rm(temp)

temp <- cbind(fraction_mitochondiral_genes_per_cell, n_reads_per_cell, n_exp_genes_per_cell, colnames(counts))
temp1 <- temp[as.numeric(temp[,2])< 1000,] # cells with less than 1000 reads, excludes n = 4323634
temp <- temp[!as.numeric(temp[,2])<1000,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,1])>0.20,]) # cells with > 20% mitochondiral genes, excludes n = 13396
temp <- temp[!as.numeric(temp[,1])>0.20,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])<500,]) # cells with less than 500 genes, excludes n = 493
temp <- temp[!as.numeric(temp[,3])<500,]
temp1 <- rbind(temp1, temp[as.numeric(temp[,3])>6000,]) # cells with more than 6000 genes, excludes n = 21
temp <- temp[!as.numeric(temp[,3])>6000,]

max <- (trunc(max(n_reads_per_cell)/1000)+1)*1000
t <- seq(1,10, by = 1)
labels <- vector()
for(i in 3:trunc(log10(max))){
  labels <- c(labels, t*10^i)
}
rm(t)

plot(x = 1, log = "x", type = "n", xlim = c(1, max), ylim = c(0,1), pch = 19, 
     xlab = "n reads per cell", ylab = "fraction mitochondrial genes")
axis(side = 1, at = labels, labels = F)
points(pch=19, col = "black", x = as.numeric(temp[,2]), y = as.numeric(temp[,1]), cex = 0.2) # valid cells
points(pch = 19, col = "grey", x = as.numeric(temp1[,2]), y = as.numeric(temp1[,1]), cex = 0.2) # excluded cells

########################################################################################
# Apply Quality Criteria to data set & save
########################################################################################

counts <- counts[,colnames(counts)%in%temp[,4]]
write.table(as.matrix(counts), file = "Input/MS GSE138266/FILTERED_counts_for_MS_CSF_samples.txt", col.names = NA, row.names = T)

########################################################################################
# Get filtered information
########################################################################################

n_exp_genes_per_cell <- colSums(counts>0)
n_reads_per_cell <- colSums(counts)
n_of_cells_in_which_gene_is_expressed <- rowSums(counts>0)
n_reads_per_gene <- rowSums(counts)

out <- matrix(NA, nrow = nrow(samples), ncol = 8)
for(i in 1:nrow(samples)){
  # n cells per sample
  out[i,1] <- samples[i,1] # sample ID
  out[i,2] <- sum(grepl(pattern = samples[i,1], colnames(counts))) # n cells
  
  temp <- counts[,grepl(pattern = samples[i,1],colnames(counts))]
  
  # min n genes per sample
  out[i,3] <- min(colSums(temp>0))
  # max n genes per sample
  out[i,4] <- max(colSums(temp>0))
  # max n genes per sample
  out[i,5] <- median(colSums(temp>0))
  
  # min n reads per sample
  out[i,6] <- min(colSums(temp))
  # max n reads per sample
  out[i,7] <- max(colSums(temp))
  # max n reads per sample
  out[i,8] <- median(colSums(temp))
}
colnames(out) <- c("sampleID", "n_cells", "nGene_min", "nGene_max", "nGene_median", "nReads_min", "nReads_max", "nReads_median")
print(out)
write.table(out, file = "Input/MS GSE138266/Filtered_quality_statistics.txt", sep="\t", col.names = T, row.names = F)

















