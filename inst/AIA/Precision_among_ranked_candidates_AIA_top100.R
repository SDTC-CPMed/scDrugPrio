#
#
# Make precision plots for ranked predictions
#
# BY SAMUEL SCHAEFER
#
######################################################################

library(readxl)
library(data.table)
library(ggplot2)


#INPUT
###########################################################
ranking <- as.matrix(read_xlsx(path = "../Output/Final_ranking/FINAL_drug_ranking_evaluated.xlsx", sheet = 1))
ranking <- ranking[as.numeric(ranking[,9]) <= 100,]


ra_drug <- read.table(file = "../Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt",sep="\t", header = T, stringsAsFactors = F, quote = "")
n_ra_drug <- nrow(ra_drug)


#DATA PREPARATION
###########################################################

dir <- "../Output/Final_ranking/"

# APPROVED DRUGS
known_drugs <- ranking[,3] == T
temp <- temp2 <- temp3 <- vector()
intervall <- seq(from = 10, to = nrow(ranking), by = 10)
if(intervall[length(intervall)]<nrow(ranking)){
  intervall <- c(intervall, nrow(ranking))
}

for(i in 1:length(intervall)){
  temp[i] <- intervall[i]
  temp2[i] <- sum(known_drugs[as.numeric(ranking[,9]) <= intervall[i]])
  temp3[i] <- sum(as.numeric(ranking[,9]) <= intervall[i])
}
known_drugs <- cbind("n_rank" = temp, "n_drugs" = temp3, "n_known" = temp2, "precision" = (temp2*100)/temp3)

# LITERATURE EVIDENCE
lit_drugs <- grepl(ranking[as.numeric(ranking[,9]) <= 100,11], pattern = "Yes")
intervall <- seq(from = 10, to = length(lit_drugs), by = 10)
if(intervall[length(intervall)]<length(lit_drugs)){
  intervall <- c(intervall, length(lit_drugs))
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


p <- ggplot(data = plot, aes(x = n_drugs, y = precision, group = type)) + 
  geom_line(aes(linetype = type)) + geom_point(aes(shape = type)) +
  xlab("n of top candidates") +
  ylab("Precision (%)") +
  scale_y_continuous(breaks = seq(0,max(plot$precision)+10, by = 10), limits = c(0,80)) +
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

ggsave(filename = paste(dir, "PRECISION_for_AIA_drugs_among_ranked_candidates_top100.pdf", sep=""), plot = p, device = "pdf", width = 4, height = 2.5, units = "in")

rm(list = ls())



