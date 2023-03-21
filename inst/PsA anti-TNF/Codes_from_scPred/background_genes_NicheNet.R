#' Calculation of background genes for NicheNet ligand activity analysis
#'
#' @param data Matrix. DCA denoised, log10 transformed expression matrix. Gene names as row names. Cell IDs as column names.
#' @param cell_IDs Vector or list. Cell IDs for the cell population (e.g. the cluster / cell type) of interest. If list, every list object contains all cell labels for a given cell type.
#' @param exp_cutoff Integer or float. Expression cut-off to be considered a background gene. The standard setting corresponds to a mean log2 transformed gene expression of >= 0.1 among the cell population of interest. This cut-off was chosen to generate ca 10 000 background genes per cell population as recommended by the authors of NicheNet. The value might need to be adapted to your individual data set. The same cut-off should be applied to all cell populations within a data set.
#'
#' @return A vector of the length nrow(data) that specifies all genes passing the expression cut-off. If cell_IDs is a matrix, specifying all
#' @export
#'
background_genes_NicheNet <- function(data, cell_IDs, exp_cutoff = 0.2){

  if(is.list(cell_IDs)){
    out <- foreach(i = c(1:length(cell_IDs)), .combine = "cbind") %do% {
      temp <- cell_IDs[[i]]
      # reverse log10 tranformation
      expressed_genes <- 10**data[,colnames(data) %in% temp]

      # apply NicheNet recommended expression threshold
      temp <- log2(rowMeans(expressed_genes)) >= exp_cutoff
      expressed_genes <- rownames(expressed_genes)[temp]

      print(paste("n genes expressed", length(expressed_genes), sep=" "))

      exp_genes <- c(expressed_genes, rep(NA, times = nrow(data)-length(expressed_genes)))
      return(exp_genes)
    }
    exp_genes <- out
    colnames(exp_genes) <- names(cell_IDs)
  } else {
    # reverse log10 tranformation
    expressed_genes <- 10**data[,colnames(data) %in% cell_IDs]

    # apply NicheNet recommended expression threshold
    temp <- log2(rowMeans(expressed_genes)) >= exp_cutoff
    expressed_genes <- rownames(expressed_genes)[temp]

    print(paste("n genes expressed", length(expressed_genes), sep=" "))

    exp_genes <- c(expressed_genes, rep(NA, times = nrow(data)-length(expressed_genes)))
  }
  return(exp_genes)
}
