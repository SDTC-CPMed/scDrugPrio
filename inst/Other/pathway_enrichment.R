# Pathway enrichment (Fisher exact - one sided)

###########################
# Variable input
###########################
# pathways = matrix; one pathway per column
# degs = vector; genes of interest
# bg = vector; all measured genes (used as background)
# return_p_all_pathways = Booleam; if T returns a matrix listing all pathways in column 1 and enrichment P values in column 2
# return_genes_sig_pathways = Boolean; if T returns a matrix corresponding to "pathways" ordered by enrichment significant, excluding all unsignificant pathways
# FDR_adj = Boolean; if T p values will be adjusted using Benjamini Hochberg correction
###########################
# NOTE: uses one-sided fisher exact test
pathway_enrichment_function <- function(pathways, degs, bg, return_p_all_pathways = F, return_genes_sig_pathways = T, FDR_adj = T, cores = 1){
  
  # can not return p values for all pathways and pathway genes for significant pathways
  stopifnot(any(c(return_p_all_pathways, return_genes_sig_pathways)==F))
  library(doParallel)
  registerDoParallel(cores=cores)

  print("Calculating Fisher enrichment")
  
  ################### pathway = yes  # pathway = no
  ################################################
  # DEGs = yes      #   f[1,1]      #   f[1,2]
  ################################################
  # DEGs = no       #   f[2,1]      #   f[2,2]
  
  out <- foreach(i = 1:ncol(pathways), .combine = c) %dopar% {
    
    f <- array(dim = c(2,2))
    f[1,1] <- length(which(!is.na(intersect(pathways[,i], degs)))) #intersect DEGs and pathways, NA is removed!
    f[1,2] <- length(which(!is.na(degs)))-f[1,1] #rest DEGs
    f[2,1] <- sum(!is.na(pathways[,i])) - f[1,1] #rest pathway genes 
    f[2,2] <- sum(!is.na(bg)) - f[1,1] - f[1,2] - f[2,1]  #rest background genes
    
    return(fisher.test(f, alternative = "greater")$p.value)  # alternative="greater" --> one-sided test for "significant enrichment"
    
  }
  
  # Bonferroni adjustment
  if(FDR_adj){
    print("Applying Bonferroni adjustment")
    out <- p.adjust(out, method = "fdr")
  }
  
  if(return_p_all_pathways){
    temp <- cbind(colnames(pathways), out)
    if(FDR_adj){
      colnames(temp) <- c("Pathway", "FDR_adj_P")
    } else {
      colnames(temp) <- c("Pathway", "unadj_P")
    }
    return(temp)
  }
  
  if(return_genes_sig_pathways){
    temp <- cbind(colnames(pathways), out)
    if(sum(temp[,2]<0.05)>0){
      if(sum(temp[,2]<0.05)>1){
        temp <- temp[temp[,2]<0.05,]
        temp[order(temp[,2], decreasing = F),]
        return(pathways[,temp[,1]])
      } else { # only one sig enriched pathway
        temp <- temp[temp[,2]<0.05,]
        return(as.matrix(pathways[,colnames(pathways)==temp[1]]))
      }
    } else {
      print("No significant enrichment")
      return()
    }
  }

}


