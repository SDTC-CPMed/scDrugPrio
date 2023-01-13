# Extract GWAS

####################################################################################
# GWAS catalog downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads -> "All associations v1.0" 31/1 2019
####################################################################################
library(readxl)

fp <- getwd()
setwd(paste(fp, "/..", sep=""))

GWAScat <- as.matrix(read_xlsx("Input/GWAScat/gwas_catalog_v1.0.2-associations_e93_r2019-01-31.tsv.xlsx", sheet = 1, col_names = TRUE))
print("GWAS loaded")

disease <- unique(GWAScat[grepl("Rheumatoid arthritis", GWAScat[,8]),8])
print(disease)

GWAScat <- GWAScat[GWAScat[,8]%in% disease,] # Choosing only "Rheumatoid arthritis *"
GWAScat <- GWAScat[GWAScat[,28] < 1e-8,] # GWAS has to be significant with P < 1e-8

# all genes belonging to genic SNPs + nearest gene to intergenic SNPs
temp <- GWAScat[GWAScat[,26]!="0",] #intergenic
GWAScat <- GWAScat[GWAScat[,26]=="0",] #genic

# Map intergenic genes to closest upstream or downstream gene
temp[,19] <- as.numeric(temp[,19])
temp[,20] <- as.numeric(temp[,20])
print("Start mapping intergenic genes to closest upstream or downstream gene")
for(i in c(1:nrow(temp))){
  # col 15 = Mapped_Gene
  # col 19 = Upstream gene distance
  # col 20 = Downstream gene distance
  if(!any(is.na(temp[i,19:20]))){
    if(temp[i,19]> temp[i,20]){
      temp[i,15] <- strsplit(temp[i,15], split = " - ")[[1]][2] #downstream gene is closer
    } else {
      temp[i,15] <- strsplit(temp[i,15], split = " - ")[[1]][1] #upstream gene is closer
    }
  } 
}

temp <- temp[!is.na(temp[,24]),] #exclude SNPs without ID, as this is an analysis of SNP1 x SNP2
temp <- temp[!is.na(temp[,15]),] #exclude SNPs without mapped genes

GWAScat <- rbind(GWAScat, temp)
GWAScat <- GWAScat[,c(8,15)]
rm(temp)

#now split all gene names
print("Splitting gene names")
gwas <- unlist(strsplit(GWAScat[,2], split = ", ", fixed = T))
gwas <- unlist(strsplit(temp, split="; ", fixed = T))
gwas <- matrix(gwas, ncol = 1)
colnames(gwas) <- "Rheumatoid_athritis_GWAS"
print("GWAS genes extracted")

write.table(gwas, file = "Input/GWAScat/GWAS_RA_P<1e-8.txt", sep="\t", col.names = T, row.names = F)

rm(list = ls())
