# Translate mouse results to human genes


# Functions
#############################################################################################################
human_mouse_homologs_function <- function(){
  human_mouse_homolog <- read.table(file = "Input/Human-mouse_homologs/HOM_MouseHumanSequence.txt", sep = "\t", header = T)
  print("Unique human mouse homologs: (n)")
  print(length(unique(human_mouse_homolog[,1])))
  print("Between:")
  print(paste(unique(human_mouse_homolog[,2]), sep=""))
  
  temp1 <- human_mouse_homolog[human_mouse_homolog[,2]=="human",c(1,4:5)] # human
  temp2 <- human_mouse_homolog[human_mouse_homolog[,2]=="mouse, laboratory",c(1,4:5)]
  transl <- merge(temp1, temp2, by=colnames(human_mouse_homolog)[1])
  colnames(transl) <- c("HomologGene.ID", "human_gene_symbol", "human_entrez_ID", "mouse_gene_symbol", "mouse_entrez_ID")
  transl <- as.matrix(transl)
  
  print(paste("n of duplicate mappings to same HomologGene.ID: ", sum(duplicated(transl[,1])), sep=""))
  print("Deleting genes not found in scRNAseq data")
  temp.transl <- transl[transl[,4] %in% mouse_genes,]
  # Check for incomplete matchings
  temp <- mouse_genes[!(mouse_genes%in%temp.transl[,4])]
  potential_match <- c("HomologGene.ID", "human_gene_symbol", "human_entrez_ID", "mouse_gene_symbol", "mouse_entrez_ID", "potential_match")
  for(i in 1:length(temp)){
    if(any(grepl(temp[i], transl[,4]))){
      pos <- which(grepl(temp[i], transl[,4]))
      if(length(pos)>1){
        t <- matrix(NA, ncol = 6, nrow = length(pos))
        for(j in 1:length(pos)){
          t[j,1:5] <- transl[pos[j],]
        }
        t[,6] <- temp[i]
        potential_match <- rbind(potential_match, t)
        rm(t,j)
      } else {
        potential_match <- rbind(potential_match, c(transl[pos,], temp[i]))
      }
      rm(pos)
    }
  }
  potential_match <- potential_match[c(24,25,36,39,50,54,58,59,61,70,72,73,74,75,77,79,84,86,88,103,194,241),] # manual determined potential matches
  potential_match <- potential_match[!(potential_match[,4]%in%temp.transl[,4]),]
  rm(temp.transl)
  # incorporate incomplete matchings
  for(i in 1:nrow(potential_match)){
    transl[transl[,4] == potential_match[i,4],4] <- potential_match[i,6]
  }
  write.table(potential_match, file = "Input/Human-mouse_homologs/incorporated_incomplete_matchings.txt", sep="\t", col.names = T, row.names = F)
  rm(potential_match)
  transl <- transl[transl[,4] %in% mouse_genes,]
  print(paste("n of duplicate mappings several mouse genes to same human gene: ", sum(duplicated(transl[,4])), sep=""))
  
  # Deal with duplicates
  duplicates <- transl[duplicated(transl[,4]),4]
  duplicates <- transl[transl[,4]%in% duplicates,]
  t <- 1:nrow(duplicates)
  t <- t[!(t%in%c(1,3,6,7,10,11,13,16,17,19,21,24,25,27,30,32,34,36,39,40,42,45,46,48,50,52,54,57,59,61,67,70,71,74,77,79,80,83,85,88,91))] # unwanted duplicates 
  duplicates <- duplicates[t,]
  transl <- transl[!(paste(transl[,2], "_", transl[,4],sep="")%in%paste(duplicates[,2], "_", duplicates[,4], sep="")),]
  rm(t, duplicates)
  
  return(transl)
}

# in.obj = matrix
# translation = matrix; column 1 = gene names as in in.obj, column 2 = gene names to be used instead
translate_by_column <- function(in.obj, translation){
  for(i in 1:ncol(in.obj)){
    temp <- in.obj[,i]
    in.obj[,i] <- NA
    temp <- translation[translation[,1]%in%temp,2] # mouse gene symbol to human entrez ID
    temp <- unique(temp)
    if(length(temp)>0){
      in.obj[1:length(temp),i] <- temp
    }
  }
  in.obj <- in.obj[rowSums(!is.na(in.obj))>0, colSums(!is.na(in.obj))>0]
  return(in.obj)
}

# Analysis
#############################################################################################################

fp <- paste(getwd(),"/..",sep="")
setwd(fp)

# LOAD: mouse genes from scRNA-seq
'
scRNAseq_data <- read.table(file = "Input/DCA_adjusted_matrix/mean.tsv", sep="\t", header = T)
mouse_genes <- scRNAseq_data[,1]
mouse_genes <- as.character(mouse_genes[!is.na(mouse_genes)])
write.table(matrix(mouse_genes, ncol = 1), file = "Input/Human-mouse_homologs/mouse_genes_from_scRNAseq.txt", sep="\t", col.names = F, row.names = F)
rm(scRNAseq_data)
'
mouse_genes <- read.table(file = "Input/Human-mouse_homologs/mouse_genes_from_scRNAseq.txt", sep="\t", header = F)
mouse_genes <- as.vector(mouse_genes[,1])

# LOAD: human-mouse homologs
human_mouse_homologs <- human_mouse_homologs_function()
human_mouse_homologs <- human_mouse_homologs[,-1]
write.table(human_mouse_homologs, file = "Input/Human-mouse_homologs/transl.txt", col.names = T, row.names = F, sep="\t")

print(paste("Total n mouse genes (from scRNA-seq): ", length(mouse_genes), sep=""))
mouse_genes <- mouse_genes[mouse_genes %in% human_mouse_homologs[,3]]
print(paste("Translatable mouse genes:", length(mouse_genes), sep=" "))

# LOAD: DCA MAST DEGs and translate
degs <- read.table(file = "Output/DCA_MAST_DEGs/Summary_sig_adj_MAST_DEGs_log(1.5)_all_clusters_res=0.6_dims=32.txt", sep="\t", header = T)
degs <- translate_by_column(degs, translation = human_mouse_homologs[,3:2]) # mouse gene symbol to human entrez ID
write.table(degs, file = "Output/DCA_MAST_DEGs/TRANSLATED_Summary_sig_adj_MAST_DEGs_all_clusters_res=0.6_dims=32.txt", sep="\t", col.names =T, row.names = F)
# LOAD: DCA MAST DEGs (high threshold) and translate
degs <- read.table(file = "Output/DCA_MAST_DEGs/Summary_sig_adj_MAST_DEGs_high_threshold_all_clusters_res=0.6_dims=32.txt", sep="\t", header = T)
degs <- translate_by_column(degs, translation = human_mouse_homologs[,3:2]) # mouse gene symbol to human entrez ID
write.table(degs, file = "Output/DCA_MAST_DEGs/TRANSLATED_Summary_sig_adj_MAST_DEGs_high_threshold_all_clusters_res=0.6_dims=32.txt", sep="\t", col.names =T, row.names = F)
# LOAD: Correlation DCA genes and translate
degs <- read.table(file = "Output/Correlation_genes/Summary_Correlation_DCA_genes.txt", sep="\t", header = T)
degs <- translate_by_column(degs, translation = human_mouse_homologs[,3:2]) # mouse gene symbol to human entrez ID
write.table(degs, file = "Output/Correlation_genes/TRANSLATED_Summary_Correlation_DCA_genes.txt", sep="\t", col.names =T, row.names = F)
# LOAD: KEGG pathway genes and translate
degs <- as.matrix(read.table(file = "Input/KEGG_pathways/kegg_genes_by_pathway.txt", sep="\t", header = T))
colnames(degs) <- as.character(degs[1,])
degs <- degs[-1,]
degs <- translate_by_column(degs, translation = human_mouse_homologs[,1:2]) # human gene symbol to human entrez ID
write.table(degs, file = "Input/KEGG_pathways/TRANSLATED_kegg_genes_by_pathway.txt", sep="\t", col.names =T, row.names = F)
# LOAD: GWAS and translate
degs <- read.table(file = "Input/GWAScat/GWAS_RA_P<1e-8.txt", sep="\t", header = T)
degs <- translate_by_column(degs, translation = human_mouse_homologs[,1:2]) # human gene symbol to human entrez ID
write.table(degs, file = "Input/GWAScat/TRANSLATED_GWAS_RA_P<1e-8.txt", sep="\t", col.names =T, row.names = F)


rm(list = ls())