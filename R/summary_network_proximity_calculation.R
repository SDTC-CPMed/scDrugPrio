#' Evaluation of network screening
#'
#' Creation of summary files used for evaluation of drug effect on the targeted genes.
#' Counteraction of differential expressed genes by drug targets is later used for drug candidate selection.
#'
#' @param files character vector; full path and file names for files from ```seperate_unique_drug_target_combinations_into_individual_drugs()```.
#' @param deg_files character vector same length as `files`; specifying full path and file name for the DEG file that should be used for a given element in `files`.
#' @param ppin matrix with 2 columns representing all interactions in the protein-protein interaction network (PPIN) used for network calculations.
#' @param top_degs either integer vector of same length as 'files' specifying the number of top DEGs used for calculation of network distances or a character string  specifying "all", in case all DEGs were used. If 'top_DEGs' is supplied as a numeric vector, a value in top_DEGs that is higher than the available number of genes in 'deg_files' will result in the use of all DEGs.
#' @param drug_effect_data matrix, first column = DrugBank ID, second = drug name, third = drug target, fourth = pharmacological effect on target.
#' @param transl matrix; 3 columns. Column 1 matches gene annotation in `ppin`, column 2 matches gene annotation in `deg_files`, column 3 matches gene annotation in `DrugBank_info`. Not required if genes are annotated similarly.
#' @param file_names vector same length as `files`; specifying output file names without file path. If `file_names` is equal to 'files' and `out_dir` is getwd() this function will override network distance calculation files.
#' @param out_dir character string; specifying output directory.
#' @param cores number of cores used for computing.
#'
#' @export Creates evaluation files for every element in `files` that can be used for manual evaluation of whether drugs counteract the fold change of targeted DEGs.
#'
summary_network_proximity_calculation <- function(files,
                                                  deg_files,
                                                  ppin,
                                                  top_degs = "all",
                                                  drug_effect_data,
                                                  file_names,
                                                  transl = NULL,
                                                  out_dir = getwd(),
                                                  cores = 1){

  # CHECK IF BASIC REQUIREMENTS ARE FULFILLED
  ######################################################################

  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop(
      "Package \'doParallel\' must be installed to use this function",
      call. = FALSE
    )
  }
  library(doParallel)

  if (length(top_degs) > 1) {
    if(length(top_degs) != length(files)){
      stop(
        "'top_degs' has wrong dimensions, please ensure that 'top_degs' has is a numeric vector with same length as files or that 'top_degs' is a character string of length = 1.",
        call. = FALSE
      )
    } else if (!is.numeric(top_degs)) {
      stop(
        "'top_degs' has wrong format, please ensure that 'top_degs' is a numeric vector.",
        call. = FALSE
      )
    }
    if(any(is.na(top_degs))){
      stop(
        "'top_degs' includes NA values, please make sure that one integer is supplied for every element in 'files'.",
        call. = FALSE
      )
    }
  }

  # SET UP ENVIRONMENT
  ######################################################################

  dir.create(out_dir, showWarnings = F)

  # set up parallel running
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

  ppin <- unique(c(unique(ppin[,1]), unique(ppin)))
  # make sure to exclude all genes that were not in PPIN
  #transl <- transl[transl[,1]%in% ppin, ]
  drug_effect_data <- drug_effect_data[drug_effect_data[,3] %in% transl[,3],] # keep only human targets with matchable gene names
  drug_effect_data <- drug_effect_data[!is.na(drug_effect_data[,3]),]

  # ANALYSIS
  ########################################################

  foreach(file = c(1:length(files))) %dopar% { # for every file

    # Load files
    in_file <- as.matrix(read.table(file = files[file], sep="\t", header = T, stringsAsFactors = F))
    degs <- as.matrix(read.table(file = deg_files[file], sep="\t", header = T, stringsAsFactors = F))

    if(!is.null(transl)){
      degs <- degs[degs[,1]%in% transl[,2],]
      degs[,1] <- transl[match(degs[,1], transl[,2]),2] # traslate to gene annotation used in network distance calculation to select the genes used there.
      degs <- degs[!is.na(degs[,1]),] # remove invalid mappings
    }

    degs <- degs[as.numeric(degs[,6]) < 0.05,] # adj P value < 0.05
    if(!(length(top_degs) == 1 && top_degs %in% c("all", "All", "ALL"))){
      if(top_degs[file] <= nrow(degs)){
        degs <- degs[1:top_degs[file],]
      }
    }
    degs[,1] <- as.character(transl[match(degs[,1], transl[,2]),3]) # translate to gene annotation of DrugBank for output file (NCBI symbols)
    degs <- degs[!is.na(degs[,1]),] # remove invalid mappings


    # select drugs with dc < 1 and zc < -1.64 (that is P < 0.05 for network proximity)
    in_file <- drug_candidate_selection(input= in_file)

    if(!is.null(in_file)){
      # examine pharmacological effect on drug targets as well as fold change direction of targeted DEGs
      in_file <- drug_effect_on_targeted_degs(input = in_file, degs = degs, drug_effect = drug_effect_data[,-2])
    }

    if(!is.null(in_file)){
      # save drugs filtered with dc < 1 and zc < -1.64
      write.table(in_file, file = paste(out_dir, file_names[file], sep=""), sep="\t", col.names = T, row.names = F)
    }
  }
  print("All files processed!")
}

# USED AS EXAMPLE
# files = paste("data/AIA/DCA_MAST_DEGs_predictions/SUMMARY/", list.files("data/AIA/DCA_MAST_DEGs_predictions/SUMMARY/", pattern = "Cluster_14"),sep="") # vector; full path and file name for files from seperate_unique_drug_target_combinations_into_individual_drugs().
# files <- files[!grepl(files, pattern = "top")]
# deg_files = rep("data/AIA/DCA_MAST_DEGs/Cluster_14_res=0.6_dims=32_k=20.txt", times = length(files)) # vector same length as 'files'; specifying full path and file name for the DEG file that should be used for a given element in 'files'.
# drug_effect_data = as.matrix(read.table(file = "data-raw/all_drug_targets_drug_bank.txt",sep="\t", header = T, stringsAsFactors = F))[,c(1,2,6,7)]  # first column = DrugBank ID, second = drug name, third = drug target, fourth = pharmacological effect on target
# ppin = as.matrix(read.table(file = "data-raw/lit_ppi.txt",sep="\t", header = T, stringsAsFactors = F)) # PPIN used for network calculations
# top_degs = "all" # either integer vector of same length as 'files' specifying the number of top DEGs used for calculation of network distances or a character string of length 1 specyfing "all" if all DEGs were used. If 'top_DEGs' is supplied as a numeric vector, a value in top_DEGs that is higher than the available number of genes in 'deg_files' will result in the use of all DEGs.
# file_names = paste("EVALUATION_",files,sep="") # vector same length as 'files'; specifying output file names without file path.
# out_dir = "data/AIA/DCA_MAST_DEGs_predictions/SUMMARY/"
# cores = 5
# transl <- as.matrix(read.table(file = "data-raw/Human-mouse_homologs/transl.txt",sep="\t", header = T, stringsAsFactors = F))[,c(2,3,1)]
# transl2 <- as.matrix(read.table(file = "data-raw/HGNC translation matrix 201108/transl.txt",sep="\t", header = T, stringsAsFactors = F))
# transl[,1] <- as.numeric(transl[,1])
# file_names = list.files("data/AIA/DCA_MAST_DEGs_predictions/SUMMARY/", pattern = "Cluster_14")
#
# summary_network_proximity_calculation(files = files,
#                                       deg_files = deg_files,
#                                       ppin = ppin,
#                                       top_degs = top_degs,
#                                       drug_effect_data = drug_effect_data,
#                                       transl = transl,
#                                       file_names = list.files("data/AIA/DCA_MAST_DEGs_predictions/SUMMARY/", pattern = "Cluster_14"),
#                                       out_dir = out_dir,
#                                       cores = 5)









