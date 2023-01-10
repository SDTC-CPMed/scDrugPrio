#' NicheNet ligand activity analysis for all cell types
#' Will test NicheNet ligand activity for all possible cell type interactions (inclduign self-interactions). Genes are formatted as Entrez Symbols in order to be compatible with the standard ligand target interactions and ligand receptor network.
#'
#' @param degs Matrix. Every column represents one cell types DEGs. Column names indicate the cell type.
#' @param background_genes Matrix. Every column represents one cell types background genes. Same column names as in `degs`. Background genes for cell types can be derived from `background_genes_NicheNet()`
#' @param out_dir String. Indicates the full path to the folder in which NicheNet ligand activity meta-data and results should be stored.
#' @param ligand_target Matrix. Targets in rows, ligands in columns. If NULL, the standard information from NicheNet will be downloaded: https://zenodo.org/record/3260758/files/ligand_target_matrix.rds
#' @param lr_network Matrix. If NULL, the standard information from NicheNet will be downloaded: https://zenodo.org/record/3260758/files/lr_network.rds
#' @param cores Integer. Number of cores to be used in parallel computing. Defaults to number of cores available - 1.
#' @param only_pos Boolean. By default this function returns only ligand interactions with positive Pearsson coefficients, given that these symbolize a interaction between the ligand and the DEGs in the downstream cell type. Negative Pearsson coefficients instead symbolize a interaction between the ligand and the background gene of the downstream cell type.
#' @return Returns a matrix for all ligand interactions between all possible cell types / clusters.
#' @export
#'
#'
NicheNet_ligand_activity_analysis <- function(degs,
                                              background_genes,
                                              out_dir = "",
                                              ligand_target = NULL,
                                              lr_network = NULL,
                                              cores = 1,
                                              only_pos = T){

  set.seed(35)

  # CHECK IF BASIC REQUIREMENTS ARE FULFILLED
  ######################################################################

  if(out_dir != ""){
    if(substr(x = out_dir, start = nchar(out_dir), stop = nchar(out_dir)) != "/"){
      out_dir <- paste(out_dir, "/",sep="")
    }
  }

  if(is.null(degs)){
    stop(
      "Variable \'ppin\' not specified.",
      call. = FALSE
    )
  }

  if(is.null(background_genes)){
    stop(
      "Variable \'background_genes\' not specified.",
      call. = FALSE
    )
  }

  if(is.null(ligand_target)){
    options(timeout = 5*60)
    ligand_target <- as.matrix(readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))) # targets = rows, ligands = columns
  }

  if(is.null(lr_network)){
    options(timeout = 5*60)
    lr_network = as.matrix(readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds")))
  }

  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    stop(
      "Package \'nichenetr\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(nichenetr)

  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop(
      "Package \'doParallel\' must be installed to use this function.",
      call. = FALSE
    )
  }
  library(doParallel)
  registerDoParallel(cores = cores)

  # Only include cell types that express DEGs
  ######################################################################
  degs <- degs[colSums(!is.na(degs))>0,]
  background_genes <- background_genes[, colnames(background_genes) %in% colnames(degs)]

  # Loop through all possible cell type - cell type combinations
  ######################################################################
  
  all_ligand_activity <- foreach(i = c(1:ncol(degs)), .combine = "rbind", .packages = c("nichenetr", "doParallel")) %dopar% {

    # sender DEGs
    sender_degs <- degs[!is.na(degs[,i]),i]
    sender_degs <- sender_degs[sender_degs %in% colnames(ligand_target)]
    sender_degs <- sender_degs[!is.na(sender_degs)]

    if(length(sender_degs)>0){
      
      out <- foreach(j = c(1:ncol(background_genes)), .combine = "rbind", .packages = "nichenetr") %do% {

        # background genes
        target_background <- background_genes[!is.na(background_genes[,j]),j]

        # target DEGs
        target_DEGs <- degs[,colnames(degs) == colnames(background_genes)[j]]
        target_DEGs <- target_DEGs[!is.na(target_DEGs)]

        # NicheNet Ligand Analysis
        ligand_activity <- predict_ligand_activities(geneset = target_DEGs,
                                                     background_expressed_genes = target_background,
                                                     ligand_target_matrix = ligand_target,
                                                     potential_ligands = sender_degs,
                                                     single = T)

        #rank ligands based on ligand activity (pearson correlation coefficient)
        if(exists("ligand_activity")){
          ligand_activity <- as.matrix(ligand_activity)
          if(is.vector(ligand_activity)){
            temp_names <- names(ligand_activity)
            ligand_activity <- matrix(ligand_activity, nrow = 1)
            colnames(ligand_activity) <- temp_names
            rm(temp_names)
          } else if (nrow(ligand_activity)>1){
            ligand_activity <- ligand_activity[order(ligand_activity[,4], decreasing = T),]
          }
          write.table(ligand_activity, file = paste(out_dir, "ligand_activity_",colnames(degs)[i], "_to_",colnames(background_genes)[j],".txt", sep=""),
                      sep="\t", col.names = T, row.names = F)
        } else {
          ligand_activity <- matrix(NA, nrow = 1, ncol = 4)
        }
        temp <- cbind(ligand_activity, colnames(degs)[i], colnames(background_genes)[j])
        rm(ligand_activity)
        return(temp)
      }
      
    } else {
      out <- matrix(NA, nrow = 1, ncol = 6)
      colnames(out) <- c("test_ligand", "auroc", "aupr", "pearson", "", "")
    }
    return(out)
  }
  
  all_ligand_activity <- all_ligand_activity[!is.na(all_ligand_activity[,1]),]
  colnames(all_ligand_activity)[(ncol(all_ligand_activity)-1):ncol(all_ligand_activity)] <- c("Sender", "Target")
  if(only_pos){
    all_ligand_activity <- all_ligand_activity[as.numeric(all_ligand_activity[,4])>0,]
  }
  write.table(all_ligand_activity, file = paste(out_dir, "all_ligand_activity.txt", sep=""), sep="\t", col.names = T, row.names = F)

  # Summarize all ligand interactions
  ######################################################################
  for(i in 1:length(unique(all_ligand_activity[,5]))){
    if(i == 1 && exists("overview_ligand_activity")){
      rm(overview_ligand_activity)
    }
    for(j in 1:length(unique(all_ligand_activity[,6]))){
      t_i <- unique(all_ligand_activity[,5])[i]
      t_j <- unique(all_ligand_activity[,6])[j]

      temp <- all_ligand_activity

      if(sum(temp[,5]== t_i)>1){
        temp <- temp[temp[,5]== t_i,]
        if(sum(temp[,6]== t_j)>1){
          temp <- temp[temp[,6]== t_j,]
          temp <- nrow(temp)
        } else {
          if(sum(temp[,6]== t_j)==1){
            temp <- 1
          } else {
            temp <- 0
          }
        }
      } else {
        if(sum(temp[,5]== t_i)==1){
          if(which(temp[,5]==t_i) %in% which(temp[,6] == t_j)){
            temp <- 1
          } else {
            temp <- 0
          }
        } else {
          temp <- 0
        }
      }

      if(!exists("overview_ligand_activity")){
        overview_ligand_activity <- matrix(c(t_i, t_j, temp), nrow = 1)
        colnames(overview_ligand_activity) <- c("Sender", "Target", "n_interactions")
      } else {
        overview_ligand_activity <- rbind(overview_ligand_activity, c(t_i, t_j, temp))
      }
    }
  }
  rm(t_i, t_j, i, j, temp)
  write.table(overview_ligand_activity, file = paste(out_dir, "overview_interactions.txt", sep=""), sep="\t", col.names = T, row.names = F)

  return(all_ligand_activity)
}

