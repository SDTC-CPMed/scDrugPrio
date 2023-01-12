#
#
# Make plots establishing the value of rank aggregation parameters.
########################################################################################################

library(doParallel)
library(ggplot2)
library(reshape2)
fp <- getwd()
setwd(paste(fp, "/..", sep=""))

##############################################################################################################################################
# old RA mouse
##############################################################################################################################################

RA_drugs <- as.matrix(read.table(file = "Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt", sep="\t", stringsAsFactors = F, header = T, quote = ""))
dir <- "Output/Network_distances/"
lf <- list.files(path = dir)
lf <- lf[grepl(pattern = "INDIVIDUAL_DRUGS_Cluster_", x = lf)]

# Read in silico drug efficacy output
out <- foreach(i = (1:length(lf)), .combine = "rbind") %do% {
  infile <- as.matrix(read.table(file = paste(dir, lf[i], sep=""), sep="\t",header = T, stringsAsFactors = F))
  temp <- unlist(strsplit(strsplit(lf[i], split = "_")[[1]][4], split = ".txt")) # get cluster number
  temp <- temp[temp != ""]
  infile <- cbind(infile, temp)
  return(infile)
}
colnames(out)[9] <- "Cluster"
cl <- unique(out[,9]) # clusters
n_RA <- length(unique(out[grepl(out[,2],pattern = "TRUE"),1]))

fc_criterium_file <- read.delim(file = "Output/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt",sep="\t",stringsAsFactors = F,header = T)
z_d_fc_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  if(any(fc_criterium_file[,19] %in% paste("FC_INDIVIDUAL_DRUGS_Cluster_",cl[i],"_res=0.6_dims=32_k=15",sep=""))){
    temp <- fc_criterium_file[fc_criterium_file[,19] %in% paste("FC_INDIVIDUAL_DRUGS_Cluster_",cl[i],"_res=0.6_dims=32_k=15",sep=""),]
    temp <- temp[as.numeric(temp[,12]) > 0,1] # choose only drugs names that counteract DEGs
    if(length(temp)>0){
      
    } else {
      temp <- NA
    }
  } else {
    temp <- NA
  }
  return(cbind(temp, cl[i]))
}

colnames(z_d_fc_criteria) <- c("Drug", "Cluster")
z_d_fc_criteria <- table(z_d_fc_criteria[!is.na(z_d_fc_criteria[,1]),1])
temp <- names(z_d_fc_criteria)
temp <- temp %in% RA_drugs[,2]
names(z_d_fc_criteria) <- temp

bar_plot_data <- foreach(i = c(1:max(z_d_fc_criteria)), .combine = "rbind") %do% {
  temp <- z_d_fc_criteria[z_d_fc_criteria == i]
  return(c(sum(names(temp)==T), sum(names(temp)==F), i))
}
colnames(bar_plot_data) <- c("n_RA", "n_Other","n_Reoccuring")

#######################################
# MAKE BAR PLOT showing that reoccuring candidates are better than candidates occuring in only one cluster
#######################################
rownames(bar_plot_data) <- bar_plot_data[,3]
df <- melt(bar_plot_data[,c(1:2)])

p <- ggplot(df, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
  xlab("times predicted") +
  ylab("n drugs") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 20),
        axis.title = element_text(colour = "black", size = 20),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "Output/Final_ranking/PRECISION_for_known_RA_drugs_among_reoccuring_candidates.pdf", plot = p, device = "pdf", width = 100, height = 110, units = "mm")
