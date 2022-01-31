#' Combination of evaluation files for checking manualy if drugs counteract targeted DEGs fold change
#'
#' @param files character vector; specifying full file names including paths for all evaluation files that should be combined into one output file.
#' @param label_for_files character vector; same length as 'files'. Specifies the label for each file that should be used in output file.
#' @param output_file_name character string; specifies name of output file (including output directory path and ".txt" ending)
#' @param cores number of cores used for computation.
#'
#' @return Saves a combined file in the chosen directory.
#' @export
#'
combine_evaluation_files <- function(files, label_for_files, output_file_name, cores = 1){

  # CHECK IF BASIC REQUIREMENTS ARE FULFILLED
  ######################################################################

  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop(
      "Package \'doParallel\' must be installed to use this function",
      call. = FALSE
    )
  }
  library(doParallel)

  # SET UP ENVIRONMENT
  ######################################################################

  # set up parallel computing
  if(cores > 1){
    max_cores <- detectCores()
    if(cores <= max_cores){
      registerDoParallel(cores = cores)
    } else {
      stop(
        "Variable \'cores\' higher than available max available",
        call. = FALSE
      )
    }
  }

  # ANALYSIS
  ########################################################

  load_files <- function(file = files[file]){
    # load network prediction file
    in_file <- as.matrix(read.table(file = file, header = T, sep = "\t", stringsAsFactors = F))
    if(ncol(in_file)==1){
      print(paste("NO DRUGS FOR: ",file, sep=""))
      return(NULL)
    } else {
      return(in_file)
    }
  }

  out <- foreach(file = c(1:length(files)), .combine = "rbind") %dopar% {
    in_file <- load_drug_prediction(files[file])
    if(!is.null(in_file)){
      return(cbind(in_file, output_file_name[file]))
    } else {
      return(in_file)
    }
  }

  # SAVE FOR LITERATURE REVIEW & FC CRITERIA CHECKING
  ###################################
  write.table(out, file = output_file_name, sep="\t", col.names = T, row.names = F)
}


