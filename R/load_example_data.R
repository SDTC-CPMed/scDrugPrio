#' Load scPred example data used for demonstrations of functions
#'
#' @return Several variables to the global environment
#'
#' @export
#' @examples load_example_data()
#'
load_example_data <- function(){

  lit_ppi <<- as.matrix(read.table(file = "data-raw/lit_ppi.txt", sep="\t", header = T, stringsAsFactors = T))

  drug_target_matrix <<- as.matrix(read.table(file = "data/DrugBank_drug_target_matrix.txt",sep="\t", header = T, stringsAsFactors = F))

}
