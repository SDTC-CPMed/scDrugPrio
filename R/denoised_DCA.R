#' Deep count autoencoder denoised scRNA-seq expression matrix. 
#'  
#' Gene annotation is mouse Entrez gene symbols.
#' 
#' This expression matrix is only sample data that was created by randomly sampling the DCA adjusted expression matrix for 250 healthy and 250 sick cells   
#' from B cell clusters, T cell clusters an  myeloid clusters respectivly (as identified in the original publication).
#' 
#' 
#' @format A matrix with 16751 genes (rows) and 1500 cells (columns). 
#' 
"denoised_DCA"