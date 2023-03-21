
#######################################################################################
# Final drug ranking evaluation INDIVIDUAL CD 
#
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/..", sep=""))
library(doParallel)
#library(readxl)
#library(data.table)
library(ggplot2)
library(reshape2)


# INPUT
#######################################################################################

in_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"

patient <- paste("Patient_",c(1:11),sep="")

cd_drug <- read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt",sep="\t", header = T, stringsAsFactors = F, quote = "")
cd_drug <- cd_drug[cd_drug[,1] != "",]
n_cd_drug <- nrow(cd_drug)

out <- foreach(i = c(1:length(patient)), .combine = "rbind") %do% {
  temp <- as.matrix(read.delim(file = paste(in_path, patient[i], "/Drug_ranking/FINAL_drug_ranking_evaluated.txt", sep=""), sep="\t", header = T))
  temp <- temp[,c(1,3,8,11)]
  temp <- cbind(temp, rank = rank( - as.numeric(temp[,3])), patient = patient[i])
  return(temp)
}

#DATA PREPARATION
###########################################################

precision <- foreach(x = c(1:length(patient)), .combine = "rbind") %do% {
  ranking <- out[out[,6] %in% patient[x],]
  
  # APPROVED DRUGS
  known_drugs <- ranking[,2] == "CD_drug"
  temp <- temp2 <- temp3 <- vector()
  if(nrow(ranking)>10){
    intervall <- seq(from = 10, to = nrow(ranking), by = 10)
    if(intervall[length(intervall)]<nrow(ranking)){
      intervall <- c(intervall, nrow(ranking))
    }
  } else {
    intervall <- nrow(ranking)
  }

  for(i in 1:length(intervall)){
    temp[i] <- intervall[i]
    temp2[i] <- sum(known_drugs[as.numeric(ranking[,5]) <= intervall[i]])
    temp3[i] <- sum(as.numeric(ranking[,5]) <= intervall[i])
  }
  known_drugs <- cbind("n_rank" = temp, "n_drugs" = temp3, "n_known" = temp2, "precision" = (temp2*100)/temp3)

  # LITERATURE EVIDENCE
  lit_drugs <- grepl(ranking[as.numeric(ranking[,5]) <= 100,4], pattern = "Yes")
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
    temp2[i] <- sum(lit_drugs[as.numeric(ranking[,5]) <= intervall[i]])
    temp3[i] <- sum(as.numeric(ranking[,5]) <= intervall[i])
  }
  lit_drugs <- cbind("n_rank" = temp, "n_drugs" = temp3, "n_lit" = temp2, "precision" = (temp2*100)/temp3)

  plot <- rbind(cbind(known_drugs, type = "known"), cbind(lit_drugs, type = "lit"))
  plot <- as.data.frame(plot)
  plot$n_drugs <- as.numeric(as.character(plot$n_drugs))
  plot$precision <- as.numeric(as.character(plot$precision))
  plot$pat <- as.character(patient[x])
  
  return(plot)
}

write.table(precision, file = paste(in_path, "Precision_among_ranked_candidates.txt",sep=""),sep="\t", col.names = T, row.names = F)


precision$pat <- factor(x = precision$pat, levels = patient)




# col_palette2 <- c("#660000","#800000", "#FF0000", "#e7772b", "#f49728", "#FF7F50", # Reds
#                   "#e8c124", "#f6ed31", "#F0E68C", # Yellows
#                   
#                   "#000000", "#808080", "#A9A9A9", "#D3D3D3", # Black & greys 4
#                   
#                   "#000080", "#2005A5", "#0000FF", "#2b58a3", "#1F65CC", "#3686D3", "#4197D9","#4AA4DE",  # Blues 8
#                   "#4169E1", "#6495ED", "#1E90FF", "#00BFFF", "#87CEFA", "#ADD8E6", "#2f9da7", "#008080", # Blues 8
#                   "#553c93","#82318f","#882865","#4B0082", # Purples
#                   
#                   "#006400", "#228B22", "#556B2F","#4F7942","#6B8E23","#9ACD32","#BFE890","#AFE1AF", # Greens 8
#                   "#32CD32", "#2E8B57", "#3CB371", "#66CDAA", "#7CFC00", "#00FF00", "#8FBC8F", "#ADFF2F") # Greens 8
#
#show_col(col_palette2)

my_col <- c("#6B8E23", "#f49728", "#1E90FF", "#006400","#D3D3D3", "#228B22","#1F65CC","#82318f","#87CEFA","#660000","#FF0000")
#show_col(my_col)
names(my_col) <- patient



# PLOT
###########################################################
p <- ggplot(data = precision, aes(x = n_drugs, y = precision, group = type)) + 
  
  geom_rect(data = subset(precision,pat == "Patient_1"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.2) +
  geom_rect(data = subset(precision,pat == "Patient_2"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.2) +
  geom_rect(data = subset(precision,pat == "Patient_3"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.1) +
  geom_rect(data = subset(precision,pat == "Patient_4"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.1) +
  geom_rect(data = subset(precision,pat == "Patient_5"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.01) +
  geom_rect(data = subset(precision,pat == "Patient_6"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.1) +
  geom_rect(data = subset(precision,pat == "Patient_7"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.1) +
  geom_rect(data = subset(precision,pat == "Patient_8"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.05) +
  geom_rect(data = subset(precision,pat == "Patient_9"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.5) +
  geom_rect(data = subset(precision,pat == "Patient_10"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.03) +
  geom_rect(data = subset(precision,pat == "Patient_11"), aes(fill = pat),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = my_col) +
  
  geom_line(aes(linetype = type)) + 
  geom_point(aes(shape = type), size = 2) +
  facet_grid(pat ~ ., scales = "free_y") +
  xlab("n of top candidates") +
  ylab("Precision (%)") +
  scale_y_continuous(breaks = seq(0,60,by=20), minor_breaks = NULL) +
  ylim(0,60) +
  geom_hline(yintercept = 0, color = "black") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), #element_rect(fill = "lightgrey"), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(), #element_line(colour = alpha("black",1), linewidth = 0.5, ), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
        )

ggsave(filename = paste(in_path, "ALL_INDIVIDUAL_CD_PRECISION_for_CD_drugs_among_ranked_candidates.pdf", sep=""), plot = p, device = "pdf", width = 10, height = 12, units = "in")

ggsave(filename = paste(in_path, "ALL_INDIVIDUAL_CD_PRECISION_for_CD_drugs_among_ranked_candidates.png", sep=""), plot = p, device = "png", width = 10, height = 10, units = "in")




p <- ggplot(data = precision, aes(x = n_drugs, y = precision, group = type)) + 
  geom_line(aes(linetype = type)) + 
  geom_point(aes(shape = type), size = 2) +
  facet_grid(pat ~ ., scales = "free_y") +
  xlab("n of top candidates") +
  ylab("Precision (%)") +
  scale_y_continuous(breaks = seq(0,60,by=20), minor_breaks = NULL) +
  ylim(0,60) +
  geom_hline(yintercept = 0, color = "black") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), #element_rect(fill = "lightgrey"), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(), #element_line(colour = alpha("black",1), linewidth = 0.5, ), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave(filename = paste(in_path, "ALL_INDIVIDUAL_CD_PRECISION_for_CD_drugs_among_ranked_candidates_no_background_colors.pdf", sep=""), plot = p, device = "pdf", width = 10, height = 12, units = "in")

ggsave(filename = paste(in_path, "ALL_INDIVIDUAL_CD_PRECISION_for_CD_drugs_among_ranked_candidates_no_background_colors.png", sep=""), plot = p, device = "png", width = 10, height = 10, units = "in")


