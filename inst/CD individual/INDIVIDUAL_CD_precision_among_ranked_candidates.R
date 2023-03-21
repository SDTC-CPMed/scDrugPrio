
#######################################################################################
# Final drug ranking evaluation INDIVIDUAL CD 
#
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
library(doParallel)
library(readxl)
library(data.table)
library(ggplot2)


# Add literature search outcomes to evaluation files
#######################################################################################
'
lit_search <- as.matrix(read_xlsx(path = "CD_batch_corrected/Output/Final_ranking/FINAL_drug_ranking_evaluated.xlsx",sheet = 1, col_names = T))
lit_search <- lit_search[,c(1:2,10:12)]
lit_search <- lit_search[!is.na(lit_search[,3]),]
lit_search2 <- as.matrix(read_xlsx(path = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/COMBINED_FINAL_TOP_100_INDIVIDUAL_PATIENTS_EXCEPT_1_&_10_evaluated.xlsx",sheet = 1, col_names = T))
lit_search2 <- lit_search2[,c(1:2,11:13)]
lit_search2 <- lit_search2[!is.na(lit_search2[,3]),]
lit_search3 <- as.matrix(read_xlsx(path = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_INDIVIDUAL_CD_PATIENT_1.xlsx",sheet = 1, col_names = T))
lit_search3 <- lit_search3[,c(1:2,30:32)]
lit_search3 <- lit_search3[!is.na(lit_search3[,3]),]
lit_search4 <- as.matrix(read_xlsx(path = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_INDIVIDUAL_CD_PATIENT_10.xlsx",sheet = 1, col_names = T))
lit_search4 <- lit_search4[,c(1:2,22:24)]
lit_search4 <- lit_search4[!is.na(lit_search4[,3]),]
lit_search5 <- as.matrix(read_xlsx(path = "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_CD.xlsx",sheet = 1, col_names = T))
lit_search5 <- lit_search5[,c(1:2,67:69)]
lit_search5 <- lit_search5[!is.na(lit_search5[,3]),]
lit_search <-rbind(lit_search, lit_search2, lit_search3, lit_search4, lit_search5)
rm(lit_search2, lit_search3, lit_search4, lit_search5)
lit_search <- lit_search[!duplicated(lit_search[,1]),]

in_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"

patient <- paste("Patient_",c(1:11),sep="")

for(pat in patient){
  temp <- read.table(file = paste(in_path, pat, "/Drug_ranking/FINAL_drug_ranking.txt", sep=""), sep = "\t", header = T, stringsAsFactors = F)
  colnames(temp)[1:2] <- c("Drug", "DrugBank")
    temp <- cbind(temp, lit_search[match(temp[,1], lit_search[,1]),3:5])
  
  write.table(temp, file = paste(in_path, pat, "/Drug_ranking/FINAL_drug_ranking_evaluated.txt", sep=""), sep = "\t", col.names = T, row.names = F)
}

# Additional literature search needed?
#######################################################################################

out <- foreach(i = c(1:11), .combine = "rbind") %do%{
  pat <- patient[i]
  temp <- read.table(file = paste(in_path, pat, "/Drug_ranking/FINAL_drug_ranking_evaluated.txt", sep=""), sep = "\t", header = T, stringsAsFactors = F)
  temp <- temp[temp[,9] <= 100,]
  temp <- temp[is.na(temp[,10]),]
  
  return(temp)
}
out <- out[!duplicated(out[,1]),]
write.table(out, file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Additional_lit_search_needed.txt", sep="\t", col.names = T, row.names = F)
rm(out)

#
#
#
# MANUAL WORK
#
#
#

# Load outcome for additional search
lit_search2 <- as.matrix(read_xlsx(path = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Additional_lit_search_needed_CD_evaluated.xlsx",sheet = 1, col_names = T))
lit_search2 <- lit_search2[,c(1:2,10:12)]
lit_search2 <- lit_search2[!is.na(lit_search2[,3]),]

lit_search <-rbind(lit_search, lit_search2)
rm(lit_search2)
lit_search <- lit_search[!duplicated(lit_search[,1]),]


for(pat in patient){
  temp <- read.table(file = paste(in_path, pat, "/Drug_ranking/FINAL_drug_ranking.txt", sep=""), sep = "\t", header = T, stringsAsFactors = F)
  colnames(temp)[1:2] <- c("Drug", "DrugBank")
  temp <- cbind(temp, lit_search[match(temp[,1], lit_search[,1]),3:5])
  
  write.table(temp, file = paste(in_path, pat, "/Drug_ranking/FINAL_drug_ranking_evaluated.txt", sep=""), sep = "\t", col.names = T, row.names = F)
} '

# Precision plots among ranked candidates
#######################################################################################

cd_drug <- read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt",sep="\t", header = T, stringsAsFactors = F, quote = "")
cd_drug <- cd_drug[cd_drug[,1] != "",]
n_cd_drug <- nrow(cd_drug)

in_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"

patient <- paste("Patient_",1:11,sep="")

for(pat in patient){
  
  #INPUT
  ###########################################################
  ranking <- as.matrix(read.delim(file = paste(in_path, pat, "/Drug_ranking/FINAL_drug_ranking_evaluated.txt", sep=""), sep="\t", header = T))
  
  #DATA PREPARATION
  ###########################################################
  
  # APPROVED DRUGS
  known_drugs <- ranking[,3] == "CD_drug"
  if(nrow(ranking)>10){
    intervall <- seq(from = 10, to = nrow(ranking), by = 10)
    if(intervall[length(intervall)]<nrow(ranking)){
      intervall <- c(intervall, nrow(ranking))
    }
  } else {
    intervall <- nrow(ranking)
  }
  
  temp <- temp2 <- temp3 <- vector()
    for(i in 1:length(intervall)){
    temp[i] <- intervall[i]
    temp2[i] <- sum(known_drugs[as.numeric(ranking[,9]) <= intervall[i]])
    temp3[i] <- sum(as.numeric(ranking[,9]) <= intervall[i])
  }
  known_drugs <- cbind("n_rank" = temp, "n_drugs" = temp3, "n_known" = temp2, "precision" = (temp2*100)/temp3)
  
  # LITERATURE DRUGS
  lit_drugs <- grepl(ranking[as.numeric(ranking[,9]) <= 100,11], pattern = "Yes")
  if(length(lit_drugs)>10){
    intervall <- seq(from = 10, to = length(lit_drugs), by = 10)
    if(intervall[length(intervall)]<length(lit_drugs)){
      intervall <- c(intervall, length(lit_drugs))
    }
  } else {
    intervall <- length(lit_drugs)
  }
  
  
  temp <- temp2 <- temp3 <- vector()
  for(i in 1:length(intervall)){
    temp[i] <- intervall[i]
    temp2[i] <- sum(lit_drugs[as.numeric(ranking[,9]) <= intervall[i]])
    temp3[i] <- sum(as.numeric(ranking[,9]) <= intervall[i])
  }
  lit_drugs <- cbind("n_rank" = temp, "n_drugs" = temp3, "n_lit" = temp2, "precision" = (temp2*100)/temp3)
  
  plot <- rbind(cbind(known_drugs, type = "known"), cbind(lit_drugs, type = "lit"))
  plot <- as.data.frame(plot)
  plot$n_drugs <- as.numeric(as.character(plot$n_drugs))
  plot$precision <- as.numeric(as.character(plot$precision))
  
  write.table(plot, file = paste(in_path, pat, "/Precision_among_ranked_candidates.txt",sep=""),sep="\t", col.names = T, row.names = F)
  
  
  plot <- rbind(cbind(known_drugs, type = "known"), cbind(lit_drugs, type = "lit"))
  plot <- as.data.frame(plot)
  plot$n_drugs <- as.numeric(as.character(plot$n_drugs))
  plot$precision <- as.numeric(as.character(plot$precision))
  
  p <- ggplot(data = plot, aes(x = n_drugs, y = precision, group = type)) + 
    geom_line(aes(linetype = type)) + geom_point(aes(shape = type)) +
    xlab("n of top candidates") +
    ylab("Precision (%)") +
    scale_y_continuous(breaks = seq(0,max(plot$precision)+10, by = 10), limits = c(0,40)) +
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10),
          legend.text = element_text(size = 10),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(filename = paste(in_path, pat, "/Drug_ranking/PRECISION_for_CD_drugs_among_ranked_candidates.pdf", sep=""), plot = p, device = "pdf", width = 4, height = 2.5, units = "in")
}

rm(list = ls())







