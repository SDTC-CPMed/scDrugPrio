#
#
# Investigate which criteria could be used to select drugs from network drug prediction
#
################################################################################################

library(doParallel)
library(vctrs)
library(ggplot2)
library(reshape2)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)



# INPUT
################################################################################################
dir <- "../Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/"

# Read network distances
lf <- list.files(path = dir)
lf <- lf[grepl(pattern = "INDIVIDUAL_DRUGS_prediction_based_on_CD_DCA_MAST_DEGs_literature_PPI_Cluster_", x = lf)]
out <- foreach(i = (1:length(lf)), .combine = "rbind") %do% {
  infile <- as.matrix(read.table(file = paste(dir, lf[i], sep=""), sep="\t",header = T, stringsAsFactors = F))
  temp <- strsplit(lf[i], split = "_")[[1]][13] # get cluster number
  infile <- cbind(infile, temp)
  return(infile)
}
colnames(out)[9] <- "Cluster"


clust_colors <- read.table(file = "../Output/CD/Cluster_colors.txt",sep="\t", stringsAsFactors = F, header = T)
clust_colors <- clust_colors[,c(1,2,5,4,3)]
fc_criterium_file <- read.table(file = "../Output/FC_criteria_checking/SUMMARY_only_drugs_passing_network_crit_with_lit_search_added_evaluated.txt",sep="\t",stringsAsFactors = F,header = T, quote = "")


# PLOT FOR z(d(c))
################################################################################################
z_plot <- out[,c(8,4,9)]
z_plot[z_plot[,1] == " TRUE",1] <- "TRUE" 
z_plot <- as.data.frame(z_plot)
z_plot$Cluster <- factor(x = clust_colors[match(z_plot$Cluster, clust_colors[,1]),3], levels = clust_colors[,3], ordered = T)
z_plot$Z.c. <- as.numeric(as.character(z_plot$Z.c.))
z_plot$known_CD_drug <- factor(x = z_plot$known_CD_drug, levels = c("FALSE","TRUE"), ordered = T)

print(max(z_plot$Z.c.))
print(min(z_plot$Z.c.))

# Basic dot plot
p <- ggplot(z_plot, aes(x=Cluster, y=Z.c., fill=known_CD_drug)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.035, stackratio = 0.8, position = position_dodge((0.7)), colour = NA) +
  scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey")) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size = 20),
        axis.text.x = element_text(angle = 90, color = "black", size = 10, hjust = 0.95, vjust = 0.2),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(-12,12) + scale_y_continuous(breaks = seq(-12,12,2))
ggsave(filename = paste(dir, "z_scores_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 300, height = 200, units = "mm")

p <- ggplot(z_plot, aes(x=Cluster, y=Z.c., fill=known_CD_drug)) +
  geom_violin(aes(fill = known_CD_drug), scale = 3) +
  #scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey"))+
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size = 20),
        axis.text.x = element_text(angle = 90, color = "black", size = 10, hjust = 0.95, vjust = 0.2),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(-12,12) + scale_y_continuous(breaks = seq(-12,12,2))
ggsave(filename = paste(dir, "z_scores_for_known_CD_drugs_violin.pdf", sep=""), plot = p, device = "pdf", width = 300, height = 200, units = "mm")

# PLOT FOR d(c)
################################################################################################
d_plot <- out[,c(8,3,9)]
d_plot[d_plot[,1] == " TRUE",1] <- "TRUE" 
d_plot <- as.data.frame(d_plot)
d_plot$Cluster <- factor(x = clust_colors[match(d_plot$Cluster, clust_colors[,1]),3], levels = clust_colors[,3], ordered = T)
d_plot$mean_d.c. <- as.numeric(as.character(d_plot$mean_d.c.))
d_plot$known_CD_drug <- factor(d_plot$known_CD_drug, levels = c("FALSE","TRUE"), ordered = T)

# Basic dot plot
p <- ggplot(d_plot, aes(x=Cluster, y=mean_d.c., fill=known_CD_drug)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, stackratio = 0.8, position = position_dodge((0.7)), colour = NA) +
  #scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey")) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size = 20),
        axis.text.x = element_text(angle = 90, color = "black", size = 10, hjust = 0.95, vjust = 0.2),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(0,5) + scale_y_continuous(breaks = seq(0,5,1))
ggsave(filename = paste(dir, "mean_closest_distance_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 300, height = 200, units = "mm")

# Violin plot
p <- ggplot(d_plot, aes(x=Cluster, y=mean_d.c., fill=known_CD_drug)) +
  geom_violin(aes(fill = known_CD_drug), scale = 3) +
  #scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey"))+
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size = 20),
        axis.text.x = element_text(angle = 90, color = "black", size = 10, hjust = 0.95, vjust = 0.2),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(0,5) + scale_y_continuous(breaks = seq(0,5,1))
ggsave(filename = paste(dir, "mean_closest_distance_for_known_CD_drugs_violin.pdf", sep=""), plot = p, device = "pdf", width = 300, height = 200, units = "mm")


# PLOT FOR PRECISION AT z(d(c)) cut-offs
################################################################################################
cl <- unique(out[,9]) # clusters
n_MS <- length(unique(out[grepl(out[,8],pattern = "TRUE"),1]))
prec_rec_matrix <- foreach(c = c(1:length(cl)), .combine = "rbind") %do% {
  temp <- out[out[,9]==cl[c],]
  cut_offs <- seq(ceiling(min(as.numeric(temp[,4]))*10)/10,0,0.1)
  prec_rec <- foreach(cut = c(1:length(cut_offs)), .combine = "rbind") %do% {
    temp2 <- temp[as.numeric(temp[,4]) <= cut_offs[cut],8]
    prec <- sum(grepl(temp2, pattern = "TRUE"))/length(temp2)
    recall <- sum(grepl(temp2, pattern = "TRUE"))/n_MS
    return(c(cut_offs[cut],prec,recall))
  }
  prec_rec <- cbind(prec_rec, cl[c])
  return(prec_rec)
}
colnames(prec_rec_matrix) <- c("cut_off", "precision", "recall", "Cluster")
prec_rec_matrix <- as.data.frame(prec_rec_matrix)
prec_rec_matrix$cut_off <- as.numeric(as.character(prec_rec_matrix$cut_off))
prec_rec_matrix$precision <- as.numeric(as.character(prec_rec_matrix$precision))
prec_rec_matrix$recall <- as.numeric(as.character(prec_rec_matrix$recall))
prec_rec_matrix$Cluster <- factor(x = clust_colors[match(prec_rec_matrix$Cluster, clust_colors[,1]),3], levels = clust_colors[,3], ordered = T)

# Precision
p <- ggplot(prec_rec_matrix, aes(x = cut_off, y = precision, color = Cluster)) + 
  geom_line(size = 1) +
  scale_size(range= c(0,1.5),guide = F) +
  scale_color_manual(values = clust_colors[clust_colors[,3] %in% as.character(prec_rec_matrix$Cluster),4], aesthetics = "color") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black",size = 12, face = "bold"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1.5),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        legend.key.size = unit(0.7,"lines"), 
        #legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"))
ggsave(filename = paste(dir, "PRECISION_CURVES_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 200, height = 80, units = "mm")

# RECALL
p <- ggplot(prec_rec_matrix, aes(x = cut_off, y = recall, color = Cluster)) + 
  geom_line(size = 1) + ylim(0,1) +
  scale_size(range= c(0,1.5),guide = F) +
  scale_color_manual(values = clust_colors[clust_colors[,3] %in% as.character(prec_rec_matrix$Cluster),4], aesthetics = "color") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black",size = 12, face = "bold"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1.5),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        legend.key.size = unit(0.7,"lines"),
        legend.title = element_text(face = "bold", size = 10, color = "black"))
ggsave(filename = paste(dir, "RECALL_CURVES_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 200, height = 80, units = "mm")

# RECALL vs PRECISION
p <- ggplot(prec_rec_matrix, aes(x = recall, y = precision, color = Cluster)) + 
  geom_line(size = 1) +
  scale_size(range= c(0,1.5),guide = F) +
  scale_color_manual(values = clust_colors[clust_colors[,3] %in% as.character(prec_rec_matrix$Cluster),4], aesthetics = "color") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black",size = 12, face = "bold"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1.5),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        legend.key.size = unit(0.7,"lines"),
        legend.title = element_text(face = "bold", size = 10, color = "black"))
ggsave(filename = paste(dir, "RECALL_PRECISION_CURVES_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 200, height = 80, units = "mm")


# Selection criteria precision and recall plots
################################################################################################

z_0.15_criterium <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
  temp <- temp[as.numeric(temp[,4]) < -0.15,8]
  prec <- sum(grepl(temp, pattern = "TRUE")) / length(temp)
  rec <- sum(grepl(temp, pattern = "TRUE")) / n_MS
  return(c(prec, rec, cl[i], "z_0.15_criterium"))
}

z_criterium <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
  temp <- temp[intersect(which(as.numeric(temp[,7]) < 0.05),which(as.numeric(temp[,4]) < 0)),8]
  prec <- sum(grepl(temp, pattern = "TRUE")) / length(temp)
  rec <- sum(grepl(temp, pattern = "TRUE")) / n_MS
  return(c(prec, rec, cl[i], "z_criterium"))
}

foreach(i = c(1:length(cl))) %do% {
  temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
  temp <- temp[intersect(which(as.numeric(temp[,7]) < 0.05),which(as.numeric(temp[,4]) < 0)),8]
  print(paste("n of candidates at z(c) < -1.64: ", length(temp), " in cluster ", cl[i],sep=""))
}

z_and_d_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
  temp <- temp[intersect(which(as.numeric(temp[,7]) < 0.05),which(as.numeric(temp[,4]) < 0)),]
  temp <- temp[as.numeric(temp[,3])<1,8]
  if(length(temp)>0){
    prec <- sum(grepl(temp, pattern = "TRUE")) / length(temp)
    rec <- sum(grepl(temp, pattern = "TRUE")) / n_MS
  } else {
    prec <- rec <- 0
  }
  return(c(prec, rec, cl[i], "z_and_d_criteria"))
}

z_d_fc_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
  if(any(fc_criterium_file[,20]==paste("FC_INDIVIDUAL_DRUGS_Cluster_",cl[i],"_res=0.8_dims=32_k=15",sep=""))){
    temp <- fc_criterium_file[fc_criterium_file[,20]==paste("FC_INDIVIDUAL_DRUGS_Cluster_",cl[i],"_res=0.8_dims=32_k=15",sep=""),]
    temp <- temp[as.numeric(temp[,12]) > 0,3] # choose only drugs that counteract DEGs
    if(length(temp)>0){
      prec <- sum(temp == "TRUE") / length(temp)
      rec <- sum(temp == "TRUE") / n_MS
    } else {
      prec <- rec <- 0
    }
  } else {
    prec <- rec <- 0
  }
  return(c(prec, rec, cl[i], "z_d_fc_criteria"))
}

foreach(i = c(1:length(cl))) %do% {
  temp <- fc_criterium_file[fc_criterium_file[,20]==paste("FC_INDIVIDUAL_DRUGS_Cluster_",cl[i],"_res=0.8_dims=32_k=15",sep=""),]
  temp <- temp[as.numeric(temp[,12]) > 0,3] # choose only drugs that counteract DEGs
  print(paste("n of candidates when applying all criteria: ", length(temp), " in cluster ", cl[i],sep=""))
}

combined_rec_prec <- rbind(z_0.15_criterium,z_criterium, z_and_d_criteria, z_d_fc_criteria)
colnames(combined_rec_prec) <- c("Precision", "Recall", "Cluster", "Criteria")
combined_rec_prec <- as.data.frame(combined_rec_prec)
combined_rec_prec$Precision <- as.numeric(as.character(combined_rec_prec$Precision))*100
combined_rec_prec$Recall <- as.numeric(as.character(combined_rec_prec$Recall))*100
combined_rec_prec$Cluster <- factor(x = clust_colors[match(combined_rec_prec$Cluster, clust_colors[,1]),3], levels = clust_colors[,3], ordered = T)
combined_rec_prec$Criteria <- factor(x = combined_rec_prec$Criteria, levels = c("z_0.15_criterium","z_criterium","z_and_d_criteria", "z_d_fc_criteria"), ordered = T)


# NEW PRECISION PLOT
p <- ggplot(combined_rec_prec, aes(x = Cluster, y= Precision, fill = Cluster, color = Cluster)) +
  geom_bar_pattern(stat = "Identity",
                   aes(pattern = Criteria),
                   position = position_dodge2(preserve = "single", padding = 0.2),
                   #color = "white", 
                   pattern_fill = "white",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.005,
                   #width = 0.01,
                   #pattern_key_scale_factor = 0.01
  ) + 
  scale_pattern_manual(values = c(z_0.15_criterium = "none", z_criterium = "circle", z_and_d_criteria = "stripe", z_d_fc_criteria = "crosshatch"),
                       guide = guide_legend(override.aes = list(fill = "lightgray"))) +
  scale_fill_manual(values = clust_colors[clust_colors[,3]%in%combined_rec_prec$Cluster,4], aesthetics = "color", guide = F) +
  scale_fill_manual(values = clust_colors[clust_colors[,3]%in%combined_rec_prec$Cluster,4], aesthetics = "fill",
                    #guide = guide_legend(override.aes = list(pattern = "none"))
                    guide = F) +
  ylab("Precision (%)") + xlab("") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black",size = 15, face = "bold"),
        axis.text.x = element_text(angle = 90, color = "black", size = 8, hjust = 0.95, vjust = 0.2, face = "plain"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.line = element_line(color = "black", size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.background = element_blank(),
        #legend.key = element_blank(), 
        legend.key.size = unit(1.5,"lines"),
        legend.title = element_text(face = "bold", size = 15, color = "black"))
ggsave(filename = paste(dir, "PRECISION_for_known_CD_drugs_SELECTION_CRITERIA.pdf", sep=""), plot = p, device = "pdf", width = 300, height = 175, units = "mm")

# NEW RECALL PLOT
p <- ggplot(combined_rec_prec, aes(x = Cluster, y= Recall, fill = Cluster, color = Cluster)) +
  geom_bar_pattern(stat = "Identity",
                   aes(pattern = Criteria),
                   position = position_dodge2(preserve = "single", padding = 0.2),
                   #color = "white", 
                   pattern_fill = "white",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.005,
                   #width = 0.01,
                   #pattern_key_scale_factor = 0.01
  ) + 
  scale_pattern_manual(values = c(z_0.15_criterium = "none", z_criterium = "circle", z_and_d_criteria = "stripe", z_d_fc_criteria = "crosshatch"),
                       guide = guide_legend(override.aes = list(fill = "lightgray"))) +
  scale_fill_manual(values = clust_colors[clust_colors[,3]%in%combined_rec_prec$Cluster,4], aesthetics = "color", guide = F) +
  scale_fill_manual(values = clust_colors[clust_colors[,3]%in%combined_rec_prec$Cluster,4], aesthetics = "fill",
                    #guide = guide_legend(override.aes = list(pattern = "none"))
                    guide = F) +
  ylab("Recall (%)") + xlab("") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black",size = 15, face = "bold"),
        axis.text.x = element_text(angle = 90, color = "black", size = 8, hjust = 0.95, vjust = 0.2, face = "plain"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.line = element_line(color = "black", size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.background = element_blank(),
        #legend.key = element_blank(), 
        legend.key.size = unit(1.5,"lines"),
        legend.title = element_text(face = "bold", size = 15, color = "black"))
ggsave(filename = paste(dir, "RECALL_for_known_CD_drugs_SELECTION_CRITERIA.pdf", sep=""), plot = p, device = "pdf", width = 300, height = 175, units = "mm")


# Correlation between precision and centrality
################################################################################################
centrality <- as.matrix(read.table(file = "../Output/NicheNet/Cell_type_centrality_summary.txt",sep="\t", stringsAsFactors = F, header = T))
centrality <- centrality[5,-1]
centrality <- cbind(centrality, unlist(strsplit(unlist(strsplit(names(centrality),split = "Cluster_"))[seq(2,length(centrality)*2,2)], split = "_MS_P1")))
colnames(centrality) <- c("Centrality", "Cluster")
centrality <- as.data.frame(centrality)
centrality$Centrality <- as.numeric(as.character(centrality$Centrality))
centrality$Cluster <- clust_colors[match(centrality$Cluster, clust_colors[,1]),3]
rownames(combined_rec_prec) <- NULL
rownames(centrality) <- NULL
combined_rec_prec <- cbind(combined_rec_prec, centrality[match(as.character(combined_rec_prec$Cluster), as.character(centrality$Cluster)),1])
combined_rec_prec$Centrality <- combined_rec_prec$`centrality[match(as.character(combined_rec_prec$Cluster), as.character(centrality$Cluster)), `

p <- ggplot(combined_rec_prec, aes(x = Centrality, y = Precision, fill = Criteria)) +
  geom_point(size = 2, aes(colour = Criteria))
ggsave(filename = paste(dir, "CORRELATION_PRECISION_CENTRALITY_for_CD.pdf", sep=""), plot = p, device = "pdf", width = 200, height = 100, units = "mm")

temp <- combined_rec_prec[combined_rec_prec$Criteria == "z_d_fc_criteria",]
print("Pearson correlation for z_d_fc_criteria")
print(cor.test(temp$Precision, temp$Centrality, method = "pearson"))

