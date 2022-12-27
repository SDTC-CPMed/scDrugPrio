#' Load sample data
#' Downsized data sets for learning the structure and uses of the R package
#'
load_sample_data <- function(){

  # get DrugBank matrix sample
  #load(file = "drug_bank_example_data.rda")
  data("drug_bank_example_data", package = "scDrugPrio")
  assign("drug_bank_example_data",drug_bank_example_data, .GlobalEnv)

  #transl
  #load(file = "translation_mouse_human.rda")
  data("translation_mouse_human", package = "scDrugPrio")
  assign("transl",transl, .GlobalEnv)

  # lit PPIN
  #load(file = "lit_ppi.rda")
  data("lit_ppi", package = "scDrugPrio")
  assign("lit_ppi",lit_ppi, .GlobalEnv)

  # scRNA-seq data after DCA adjustment
  #load(file = "denoised_DCA.rda")
  data("denoised_DCA", package = "scDrugPrio")
  assign("denoised_DCA", denoised_DCA, .GlobalEnv)
  
  #load(file = "latent_DCA.rda")
  data("latent_DCA", package = "scDrugPrio")
  assign("latent_DCA", latent_DCA, .GlobalEnv)
  
  #load(file = "marker_genes.rda")
  data("marker_genes", package = "scDrugPrio")
  marker_genes <- as.vector(marker_genes)
  marker_genes <- marker_genes[!is.na(marker_genes)]
  assign("marker_genes", marker_genes, .GlobalEnv)
  
  # load RA drugs
  #load(file = "RA_drugs.rda")
  data("RA_drugs", package = "scDrugPrio")
  marker_genes <- as.vector(marker_genes)
  marker_genes <- marker_genes[!is.na(marker_genes)]
  assign("RA_drugs", marker_genes, .GlobalEnv)

  # Load FC_evaluation outcome file
  #load(file = "fc_evaluation_done.rda")
  data("fc_evaluation_done", package = "scDrugPrio")
  marker_genes <- as.vector(marker_genes)
  marker_genes <- marker_genes[!is.na(marker_genes)]
  assign("fc_evaluation_done", marker_genes, .GlobalEnv)
  
}
