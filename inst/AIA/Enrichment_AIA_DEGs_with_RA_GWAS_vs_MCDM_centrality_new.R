#
#
# Validation of centrality by GWAS enrichment
#
# BY: SAMUEL SCHAEFER
###################################################

library(scales)
library(ggplot2)

# Load data
gwas <- as.matrix(read.table(file = "../Input/GWAScat/GWAS_RA_P_smaller_1e-8.txt",sep="\t", header = T))
gwas <- as.vector(gwas)
degs <- as.matrix(read.table(file = "../Output/DCA_MAST_DEGs/TRANLATED_SUMMARY_sig_adj_MAST_DEGs_log(1.5)_all_clusters_entrez_symbol.txt",sep="\t",header = T))
bg_genes <- as.matrix(read.table(file = "../Output/NicheNet/background_genes.txt",sep="\t",header = T))

centrality <- as.matrix(read.table(file = "../Output/NicheNet/Cell_type_centrality_summary.txt",sep="\t",header = T, stringsAsFactors = F))
rownames(centrality) <- centrality[,1]
centrality <- centrality[,-1]
mode(centrality) <- "numeric"
centrality <- t(centrality)

temp <- rownames(centrality)
temp <- unlist(strsplit(temp, "X\\."))
temp <- unlist(strsplit(temp, "X"))
temp <- temp[temp != ""]
rownames(centrality) <- temp

colors <- as.matrix(read.table(file = "../Output/Clustering/Cluster_colors.txt",sep="\t", header = T))

# GWAS enrichment
##################################
fisher_enrichment_function <- function(gene_sets, genes_of_interest, bg, return_p_all_gene_sets = T, return_genes_sig_gene_sets = F, FDR_adj = T, cores = 1){
  
  library(doParallel)
  registerDoParallel(cores=cores)
  
  print("Calculating Fisher enrichment")
  
  ################################ gene_set = yes  # gene_set = no
  #################################################################
  # genes_of_interest = yes      #   f[1,1]      #   f[1,2]
  #################################################################
  # genes_of_interest = no       #   f[2,1]      #   f[2,2]
  
  out <- foreach(i = 1:ncol(gene_sets), .combine = c) %dopar% {
    
    f <- array(dim = c(2,2))
    f[1,1] <- length(which(!is.na(intersect(gene_sets[,i], genes_of_interest)))) #intersect genes_of_interest and gene_sets, NA is removed!
    f[1,2] <- length(which(!is.na(genes_of_interest)))-f[1,1] #rest genes_of_interest
    f[2,1] <- sum(!is.na(gene_sets[,i])) - f[1,1] #rest gene_set genes 
    f[2,2] <- sum(!is.na(bg)) - f[1,1] - f[1,2] - f[2,1]  #rest background genes
    
    return(fisher.test(f, alternative = "greater")$p.value)  # alternative="greater" --> one-sided test for "significant enrichment"
    
  }
  
  # Bonferroni adjustment
  if(FDR_adj){
    print("Applying FDR adjustment")
    out <- p.adjust(out, method = "fdr")
  }
  
  if(return_p_all_gene_sets){
    temp <- cbind(colnames(gene_sets), out)
    if(FDR_adj){
      colnames(temp) <- c("gene_set", "FDR_adj_P")
    } else {
      colnames(temp) <- c("gene_set", "unadj_P")
    }
    return(temp)
  }
  
  if(return_genes_sig_gene_sets){
    temp <- cbind(colnames(gene_sets), out)
    if(sum(temp[,2]<0.05)>0){
      if(sum(temp[,2]<0.05)>1){
        temp <- temp[temp[,2]<0.05,]
        temp[order(temp[,2], decreasing = F),]
        return(gene_sets[,temp[,1]])
      } else { # only one sig enriched gene_set
        temp <- temp[temp[,2]<0.05,]
        return(as.matrix(gene_sets[,colnames(gene_sets)==temp[1]]))
      }
    } else {
      print("No significant enrichment")
      return()
    }
  }
  
}

result <- fisher_enrichment_function(gene_sets = degs, genes_of_interest = gwas, bg = bg_genes)

# Correlation enrichment P vs centrality
##################################

col <- colors[,c(4,1)]
col[,2] <- paste("Cluster_", as.numeric(col[,2]),sep="")

result <- cbind(result, centrality[match(result[,1],rownames(centrality)),5], col[match(result[,1],col[,2]),1])
colnames(result) <- c("cluster", "P", "centrality", "color")
result <- as.data.frame(result)
result$P <- -log10(as.numeric(as.character(result$P)))
result$centrality <- as.numeric(as.character(result$centrality))

pdf(file = "../Output/NicheNet/Correlation_between_GWAS_enrichment_and_centrality_old_RA_own_palette_221026.pdf", width = 3.5, height = 3)

ggplot(result, aes(x = P,y = centrality)) + 
  geom_point(size = 4, color = result$color, alpha = 1) + 
  ylim(c(0,1.5)) +
  xlab("-log10(adj P value)") + ylab("eigenvector centrality") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 1.5),
        axis.ticks = element_line(color = "black", size = 1.5,linetype = 1),
        axis.text = element_text(color = "black", size = 14, family = "Helvetica"),
        axis.title = element_text(color = "black", size = 14, family = "Helvetica"),
        axis.title.y = element_text(margin = margin(0,8,0,0)), 
        panel.grid.major = element_line(color = "lightgrey", size = 0.75, linetype = 3),
        axis.ticks.length = unit(0.06,"inch")) +
  scale_x_continuous(expand = c(0.05,0)) + scale_y_continuous(expand = c(0,0)) +
  stat_smooth(method = "lm", se = T, alpha = 0.25, fill = "lightgrey")


dev.off()

write.table(result, file = "../Output/NicheNet/Correlation_between_GWAS_enrichment_and_centrality_old_RA_own_palette_221026.txt", sep="\t", col.names = NA, row.names = T)
cor.test(x = result$P, y = result$centrality, method = "pearson")

# Fisher exact tests for central cell types precision compared to random

t <- matrix(0, nrow = 2, ncol = 2)
colnames(t) <- c("RA", "Other")
rownames(t) <- c("predicted candidate", "not-predicted")

# Cluster 10
t[1,] <- c(29,226-29)
t[2,1] <- 57-t[1,1]
t[2,2] <- 1840-sum(t[1,])-t[2,1]

fisher.test(t, alternative = "greater", conf.int = T)

# Cluster 12
t[1,] <- c(22,193-22)
t[2,1] <- 57-t[1,1]
t[2,2] <- 1840-sum(t[1,])-t[2,1]

fisher.test(t, alternative = "greater", conf.int = T)

# Cluster 1
t[1,] <- c(1,16-1)
t[2,1] <- 57-t[1,1]
t[2,2] <- 1840-sum(t[1,])-t[2,1]

fisher.test(t, alternative = "greater", conf.int = T)

