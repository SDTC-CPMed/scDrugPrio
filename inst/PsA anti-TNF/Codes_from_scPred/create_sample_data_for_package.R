#
#
# Create example data
#
# BY SAMUEL SCHAEFER
# 2022-01-24
#
################################
setwd("~/Library/CloudStorage/OneDrive-Linköpingsuniversitet/Old RA mouse project - Barabasi/scPred_R_package/new_220809/scPred/data")
set.seed(1)

# Drug Bank data
drug_bank_example_data <- read.table(file = "drug_targets_drug_bank.txt",sep="\t", header = T, stringsAsFactors = F)
drug_bank_example_data <- as.matrix(drug_bank_example_data)

temp1 <- drug_bank_example_data[!grepl(pattern = "Humans", x = drug_bank_example_data[,8]),] # human
drug_bank_example_data <- drug_bank_example_data[grepl(pattern = "Humans", x = drug_bank_example_data[,8]),] # human

temp2 <- drug_bank_example_data[!grepl(pattern = "approved", x = drug_bank_example_data[,3]),] # approved
drug_bank_example_data <- drug_bank_example_data[grepl(pattern = "approved", x = drug_bank_example_data[,3]),] # approved

temp3 <- drug_bank_example_data[is.na(drug_bank_example_data[,6]),] # valid drug target
drug_bank_example_data <- drug_bank_example_data[!is.na(drug_bank_example_data[,6]),] # valid drug target

# select 100 valid drugs
drug_bank_example_data <- drug_bank_example_data[drug_bank_example_data[,1] %in% unique(drug_bank_example_data[,1])[1:100],]
# add some non-valid drugs
drug_bank_example_data <- rbind(drug_bank_example_data,
                                temp1[temp1[,1] %in% unique(temp1[,1])[1:25],],
                                temp2[temp2[,1] %in% unique(temp2[,1])[1:17],],
                                temp3[temp3[,1] %in% unique(temp3[,1])[1:8],])

save(x = drug_bank_example_data, file = "drug_bank_example_data.rda")
rm(drug_bank_example_data, temp1, temp2, temp3)

# translation between mouse, human and symbols vs id's.
translation_mouse_human <- read.table(file = "transl.txt", sep="\t", header = T, stringsAsFactors = F)
translation_mouse_human <- as.matrix(translation_mouse_human)
translation_mouse_human[,2] <- as.character(as.numeric(translation_mouse_human[,2]))
translation_mouse_human[,4] <- as.character(as.numeric(translation_mouse_human[,4]))
save(x = translation_mouse_human, file = "translation_mouse_human.rda")
rm(translation_mouse_human)

# literature curated PPIN
lit_ppi <- read.table(file = "lit_ppi.txt", head = T, stringsAsFactors = F, sep="\t")
lit_ppi <- as.matrix(lit_ppi)
save(x = lit_ppi, file = "lit_ppi.rda")
rm(lit_ppi)


# DCA data
set.seed(2)
denoised_DCA <- read.table("mean.tsv", sep = "\t", header = F, stringsAsFactors = F)
rownames(denoised_DCA) <- denoised_DCA[,1]
denoised_DCA <- denoised_DCA[,-1]
colnames(denoised_DCA) <- denoised_DCA[1,]
denoised_DCA <- denoised_DCA[-1,]
denoised_DCA <- as.matrix(denoised_DCA)

latent_DCA <- read.table("latent.tsv", sep = "\t", header = F, stringsAsFactors = F)
rownames(latent_DCA) <- latent_DCA[,1]
latent_DCA <- latent_DCA[,-1]
colnames(latent_DCA) <- paste("PC_",1:ncol(latent_DCA),sep="")
latent_DCA <- t(latent_DCA)
latent_DCA <- as.matrix(latent_DCA)

#random_selection <- sample(x = c(1:ncol(latent_DCA)), size = 1500)
#denoised_DCA <- denoised_DCA[,random_selection]

cell_id <- read.table(file = "Cell_identity.txt", sep ="\t", stringsAsFactors = F, head = T)
denoised_DCA <- denoised_DCA[, colnames(denoised_DCA) %in% cell_id[cell_id[,2] %in% c(14,4,9,0,1,12,5,6,7),1]]

temp_h <- denoised_DCA[,grepl("Healthy",colnames(denoised_DCA))]
temp_s <- denoised_DCA[,!(colnames(denoised_DCA) %in% colnames(temp_h))]

cell_type <- rbind(c(4,9,14), c(0,12,1), c(5,6,7))
for(i in 1:nrow(cell_type)){
  temp <- temp_h[, colnames(temp_h) %in% cell_id[cell_id[,2] %in% cell_type[i,],1]]
  if(ncol(temp)>250){
    temp <- temp[,sample(x = c(1:ncol(temp)), size = 250)]
    temp_h <- temp_h[,!(colnames(temp_h) %in% cell_id[cell_id[,2] %in% cell_type[i,],1])]
    temp_h <- cbind(temp, temp_h)
  }
}
for(i in 1:nrow(cell_type)){
  temp <- temp_s[, colnames(temp_s) %in% cell_id[cell_id[,2] %in% cell_type[i,],1]]
  if(ncol(temp)>250){
    temp <- temp[,sample(x = c(1:ncol(temp)), size = 250)]
    temp_s <- temp_s[,!(colnames(temp_s) %in% cell_id[cell_id[,2] %in% cell_type[i,],1])]
    temp_s <- cbind(temp, temp_s)
  }
}
denoised_DCA <- cbind(temp_s, temp_h)
latent_DCA <- latent_DCA[,match(colnames(denoised_DCA), colnames(latent_DCA))]

save(x = denoised_DCA, file = "denoised_DCA.rda")
save(x = latent_DCA, file = "latent_DCA.rda")

# somewhat equally distributed selection
#table(cell_id[cell_id[,1] %in% colnames(denoised_DCA),2])


# Marker genes for cell typing
marker_genes <- read.table(file = "Markers_for_cell_typing_heatmap_RA.txt", sep = "\t",stringsAsFactors = F, header = T)
marker_genes[marker_genes == ""] <- NA
marker_genes <- marker_genes[!is.na(marker_genes[,3]),]
marker_genes <- marker_genes[order(marker_genes[,3], decreasing = F),1] # n = 62
marker_genes <- marker_genes[c(19,23, 29,31, 41,46:50,52)]
save(x = marker_genes, file = "marker_genes.rda")

