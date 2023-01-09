#' Load sample data
#' Downsized data sets for learning the structure and uses of the R package
#'
#' @export
#' @return
#' drug_bank_example_data = DrugBank data for 146 drugs
#' transl = translation matrix for human and mouse gene annotations as derived from NCBI
#' lit_ppi = literature derived protein-protein interaction network
#' denoised_DCA = DCA denoised scRNA-seq matrix for a subsample of the original scRNA-seq data set
#' latent_DCA = DCA derived latent features for thecells in denoised_DCA
#' marker_genes = genes used for cell typing of the original data set
#' RA_drugs = DrugBank IDs for known drugs against rheumathoid arthritis
#' fc_evaluation_done = file in which counteraction of fold change of drug targets has been evaluated
#'
load_sample_data <- function(){

  # get DrugBank matrix sample
  #load(file = "drug_bank_example_data.rda")
  data("drug_bank_example_data", package = "scDrugPrio")
  assign("drug_bank_example_data",drug_bank_example_data, .GlobalEnv)

  #transl
  #load(file = "translation_mouse_human.rda")
  data("translation_mouse_human", package = "scDrugPrio")
  assign("transl",translation_mouse_human, .GlobalEnv)

  # lit PPIN
  #load(file = "lit_ppi.rda")
  data("lit_ppi", package = "scDrugPrio")
  assign("lit_ppi",lit_ppi, .GlobalEnv)

  # scRNA-seq data after DCA adjustment
  #load(file = "denoised_DCA.rda")
  data("denoised_DCA_1", package = "scDrugPrio")
  data("denoised_DCA_2", package = "scDrugPrio")
  data("denoised_DCA_3", package = "scDrugPrio")
  data("denoised_DCA_4", package = "scDrugPrio")
  data("denoised_DCA_5", package = "scDrugPrio")
  data("denoised_DCA_6", package = "scDrugPrio")
  
  denoised_DCA <- cbind(denoised_DCA_1, denoised_DCA_2, denoised_DCA_3, denoised_DCA_4, denoised_DCA_5, denoised_DCA_6)
  assign("denoised_DCA", denoised_DCA, .GlobalEnv)
  
  #load(file = "latent_DCA.rda")
  data("latent_DCA", package = "scDrugPrio")
  assign("latent_DCA", latent_DCA, .GlobalEnv)

  #load(file = "marker_genes.rda")
  data("marker_genes", package = "scDrugPrio")
  assign("marker_genes", marker_genes, .GlobalEnv)

  # load RA drugs
  #load(file = "RA_drugs.rda")
  data("RA_drugs", package = "scDrugPrio")
  assign("RA_drugs", RA_drugs, .GlobalEnv)

  # Load FC_evaluation outcome file
  #load(file = "fc_evaluation_done.rda")
  data("fc_evaluation_done", package = "scDrugPrio")
  assign("fc_evaluation_done", fc_evaluation_done, .GlobalEnv)
  
  rm(list = c(grep(ls(pos = .GlobalEnv), pattern = "denoised_DCA_", value = T)), pos = .GlobalEnv)
}
