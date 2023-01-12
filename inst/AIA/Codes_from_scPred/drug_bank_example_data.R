#' A sample dataset containing the information for 150 drugs in DrugBank. 
#'
#' @format A matrix with 758 rows and 8 variables:
#' \describe{
#'   \item{drugID}{DrugBank ID}
#'   \item{status}{Status of the drug, 'approved' representing that the drug has in at least one jurisdiction anywhere, at some point in time has been approved.}
#'   \item{gene_symbol}{Entrez gene symbols of drug targets}
#'   \item{drug_action}{The effect if of the drug in the drug target}
#'   ...
#' }
#' @source Free DrugBank data, downloaded July 2019. \url{https://go.drugbank.com}
"drug_bank_example_data"