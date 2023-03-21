
#######################################################################################
# Final drug ranking evaluation PsA anti-IL17 
#
# BY SAMUEL SCHAEFER
#######################################################################################

library(ggplot2)

psa_drug <- read.table(file = "../Input/Drugs/PsA_drugs_from_DrugBank.txt",sep="\t", header = T, stringsAsFactors = F, quote = "")
psa_drug <- psa_drug[psa_drug[,1] != "",]
n_psa_drug <- nrow(psa_drug)

in_path <- "../Output/Final_ranking/"

patient <- c("Responder", "Non_responder")

for(pat in patient){
  
  #INPUT
  ###########################################################
  ranking <- as.matrix(read.delim(file = paste(in_path, pat, "/FINAL_drug_ranking_",pat,"_evaluated.txt", sep=""), sep="\t", header = T))
  
  #DATA PREPARATION
  ###########################################################
  
  # APPROVED DRUGS
  known_drugs <- ranking[,3] == T
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
  
  write.table(plot, file = paste(in_path, pat, "/Precision_among_ranked_",pat,"candidates.txt",sep=""),sep="\t", col.names = T, row.names = F)
  
  
  plot <- rbind(cbind(known_drugs, type = "known"), cbind(lit_drugs, type = "lit"))
  plot <- as.data.frame(plot)
  plot$n_drugs <- as.numeric(as.character(plot$n_drugs))
  plot$precision <- as.numeric(as.character(plot$precision))
  
  p <- ggplot(data = plot, aes(x = n_drugs, y = precision, group = type)) + 
    geom_line(aes(linetype = type)) + geom_point(aes(shape = type)) +
    xlab("n of top candidates") +
    ylab("Precision (%)") +
    scale_y_continuous(breaks = seq(0,20, by = 10), limits = c(0,20)) +
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
  
  ggsave(filename = paste(in_path, pat, "/PRECISION_for_PSA_drugs_among_",pat,"_ranked_candidates.pdf", sep=""), plot = p, device = "pdf", width = 4, height = 2.5, units = "in")
}

rm(list = ls())
