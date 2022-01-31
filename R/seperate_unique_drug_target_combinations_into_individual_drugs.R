#' Seperation of unique drug targets to individual drugs.
#'
#' @param files character vector; specifying file names for ```average_closest_distance_network_screening()``` output files that this function shall be applied to.
#' @param same_drugs matrix with 2 columns; specifying mapping of individual drugs to unique drug target combinations created by ```prepare_drug_target_matrix_for_network_distance_calculation()```. First column indicating the DrugBankID used for the unique drug target combination in average_closest_distance_network_screening() and second column indicating the mapping to individual drugs.
#' @param disease_specific_drugs character vector; specifying DrugBankID of disease specific drugs (drug names for known drug-disease pairs).
#' @param in_dir character string; specifying path to directory of 'files'. Too long directory paths can cause problems for ```read.table()``` in this function (at least on Windows).
#' @param out_dir character string; specifying path to directory for the output files
#' @param output_file_name_add_on character string; file name add on for output files. Set to "" if you do not wish to change filenames. If out_dir = in_dir and output_file_name_add_on = "" this function will overwrite 'files'.
#' @param remove_part_of_file_name_for_output character string; part of file name in 'files' that is to be removed from file name of output file.
#' @param cores integer; number of cores used for parallel computing.
#'
#' @return Creates a number of files equal to length(files) in out_dir that includes network distances for individual drugs instead of unique drug target combinations.
#' @export
#'
#' @examples
#' same_drugs = as.matrix(read.table(file = "data/SAME_DRUG_TARGETS_in_lit_ppin.txt", sep="\t", header = T, stringsAsFactors = F))
#' disease_specific_drugs = as.matrix(read.table(file = "data-raw/RA_drugs_from_DrugBank_200214.txt",sep="\t", header = T, stringsAsFactors = F))
#' seperate_unique_drug_target_combinations_into_individual_drugs(files = "data/test.txt", same_drugs = same_drugs, disease_specific_drugs = disease_specific_drugs)
#'
seperate_unique_drug_target_combinations_into_individual_drugs <- function(files,
                                                                           same_drugs,
                                                                           disease_specific_drugs,
                                                                           in_dir = getwd(),
                                                                           out_dir = getwd(),
                                                                           output_file_name_add_on = "INDIVIDUAL_DRUGS_",
                                                                           remove_part_of_file_name_for_output = "drug-disease_closest_distances_vs_random_bin_adjusted\\__",
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

  # Seperation of files
  ########################################################
  foreach(file = c(1:length(files)), .export = NULL) %dopar% { # for every file

    # Load file
    in_file <- as.matrix(read.table(file = paste(in_dir, files[file],sep=""), header = T, sep = "\t", stringsAsFactors = F))
    # Check for special case in which randomization led to SD = 0 meaning that Z & P value could not be calculated
    if(any(is.na(in_file[,7]))){
      temp <- in_file[is.na(in_file[,7]),]
      in_file <- in_file[!is.na(in_file[,7]),]
      for(i in 1:nrow(temp)){
        if(as.numeric(temp[i,3])==as.numeric(temp[i,5]) && as.numeric(temp[i,6])==0){ # actual dc is equal to mean dc of random computation & SD = 0 -> P = 1 and Z = 0
          temp[i,c(4,7)] <- c(0,1)
        }
      }
      in_file <- rbind(in_file, temp)
    }
    if(exists("temp")){
      rm(temp)
    }

    # translate to individual drugs
    if(any(duplicated(in_file[,1]))){
      stop("Output files of average_closest_distance_network_screening() can not be uniquely linked to drugs in column 1 of 'same_drugs'",
           call. = FALSE)
    }

    same_drugs <- as.matrix(same_drugs)
    mode(same_drugs) <- "character"
    same_drugs <- same_drugs[same_drugs[,1] %in% in_file[,1],]
    out <- merge(x = same_drugs, y = in_file, by.x = colnames(same_drugs)[1], by.y = colnames(in_file)[1])
    out <- out[,-1]

    # disease specific drugs?
    temp <- matrix(out[,1] %in% disease_specific_drugs)
    out <- cbind(out[,1], temp, out[,2:ncol(out)])
    colnames(out)[1:2] <- c("DrugBankID", "known_disease_specific_drug")
    rm(temp)

    #save
    temp <- strsplit(files[file], split = remove_part_of_file_name_for_output)[[1]][2]
    write.table(out, file = paste(out_dir, output_file_name_add_on, temp, ".txt", sep=""), sep="\t", col.names = T, row.names = F)

    return(paste(out_dir, output_file_name_add_on, temp, ".txt", sep=""))
  }
  print("Drug target combinations successfully seperated into individual drugs")
}


# USED FOR TESTING
# seperate_unique_drug_target_combinations_into_individual_drugs(
#   disease_specific_drugs = as.matrix(read.table(file = "data-raw/RA_drugs_from_DrugBank_200214.txt",sep="\t", stringsAsFactors = F, header = T, quote = ""))[,1],
#   same_drugs = as.matrix(read.table(file = "data/SAME_DRUG_TARGETS_in_lit_ppin.txt",sep="\t", stringsAsFactors = F, header = T)),
#   in_dir = "data/AIA/DCA_MAST_DEGs_predictions/",
#   out_dir = "data/AIA/DCA_MAST_DEGs_predictions/SUMMARY/",
#   files = list.files("data/AIA/DCA_MAST_DEGs_predictions/", pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__"),
#   cores = 5)
