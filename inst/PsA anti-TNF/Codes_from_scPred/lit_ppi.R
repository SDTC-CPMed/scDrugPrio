#' Literature curated, undirected, human protein-protein interaction (PPI) network
#' 
#' Gene annotation is Entrez gene IDs  
#' The human interactome was assembled from 16 databases as described by 
#' do Valle et al. Genes were mapped to their Entrez ID based on the 
#' National Centre for Biotechnology Information (NCBI, July, 2017) database 
#' as well as their official gene symbols. The resulting interactome includes 
#' 351,444 protein-protein interactions (PPIs) connecting 17,706 unique proteins. 
#' 
#' The largest connected component (LCC) includes 351,393 PPIs and 17,651 proteins. 
#' The LCC has to be extracted in a additional step before use of this network.
#'
#' @format A matrix with 351444 rows and 2 variables:
#' \describe{
#'   \item{Protein_A}{Connected to Protein_B}
#'   \item{Protein_B}{Connected to Protein_A}
#' }
#' @source do Valle, I. F. et al. Network medicine framework shows that proximity 
#' of polyphenol targets and disease proteins predicts therapeutic effects of 
#' polyphenols. Nature Food 2, 143-155, doi:10.1038/s43016-021-00243-7 (2021).
"lit_ppi"