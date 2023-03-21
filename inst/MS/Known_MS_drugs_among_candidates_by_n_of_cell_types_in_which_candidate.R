
#######################################################################################
# How many known MS drugs found in more than one cluster among final candidates?
#######################################################################################

#install.packages("ggplot2")
library(ggplot2)
library(scales)
library(doParallel)

out <- read.table(file = "../Output/Final_ranking/FINAL_drug_ranking.txt", sep="\t", header = T, stringsAsFactors = F)

out <- out[,c(1,3,6)]

n_cl <- unique(out[,3])

tabl_in <- foreach(i =c(1:max(n_cl)), .combine = "rbind") %do% {
  if(sum(out[,3]==i)>0){
    if(sum(out[,3]==i)>1){
      temp <- out[out[,3] == i,] 
      return(rbind(c(sum(temp[,2]!= "MS_drug"), "nonMS", i), c(sum(temp[,2]=="MS_drug"), "MS", i)))
    } else {
      temp <- as.vector(as.matrix(out[out[,3]==i,]))
      return(rbind(c(sum(temp[2]!= "MS_drug"), "nonMS", i), c(sum(temp[2]=="MS_drug"), "MS", i)))
    }
    return()
  } else {
    return(rbind(c(0,"nonMS", i), c(0, "MS", i)))
  }
}
df <- data.frame(n_drugs = as.numeric(tabl_in[,1]), MS = factor(tabl_in[,2], levels = c("MS","nonMS"), ordered = T), n_cl = as.numeric(tabl_in[,3]))



p <- ggplot(df, aes(x = n_cl, y=n_drugs, fill = MS)) + 
  geom_col(width = 0.85) + 
  scale_fill_manual(values = hue_pal()(2)) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(family = "sans", colour = "black", size = 6),
        axis.title = element_text(family = "sans",colour = "black", size = 6),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "in n clusters", breaks = seq(from = 1, to = max(df$n_cl), by = 1))

ggsave(filename = "../Output/Final_ranking/known_MS_drugs_among_final_candidates_vs_n_clusters_in_which_candidates.pdf", plot = p, device = "pdf", units = "mm", width = 70, height = 60)
#dev.off()

rm(list = ls())

