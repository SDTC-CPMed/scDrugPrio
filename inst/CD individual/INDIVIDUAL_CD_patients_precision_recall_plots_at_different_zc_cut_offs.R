#
#
# Investigate which criteria could be used to select drugs from network drug prediction
#
################################################################################################

library(doParallel)
library(ggplot2)
library(reshape2)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

# INPUT
##############################################################################################################################################

dir2 <- "../Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
pat <- paste("Patient_",1:11,sep="")

foreach(x = c(1:11)) %do% {
  
  print(pat[x])
  
  # INPUT
  ##############################################################################################################################################
  
  dir <- paste(dir2, pat[x],"/literature_PPI/SUMMARY/",sep="")
  
  lf <- list.files(path = dir)
  lf <- lf[grepl(pattern = "INDIVIDUAL_DRUGS_prediction_based_on_CD_DCA_MAST_DEGs_literature_PPI_Cluster_", x = lf)]
  lf1 <- lf[!grepl(pattern = "top",x = lf)] # selects predictions based on all DEGs
  lf2 <- lf[grepl(pattern = "top_3500",x = lf)] # selects predictions based on only 3500 top DEGs
  temp <- unlist(strsplit(lf2, split = "_top_3500.txt"))
  temp <- temp[temp != ""]
  temp <- paste(temp, ".txt",sep="")
  lf1 <- lf1[!(lf1 %in% temp)]
  lf <- c(lf1, lf2)
  rm(lf1, lf2, temp)
  
  if(pat[x] %in% c("Patient_1", "Patient_10")){
    cl_col <- read.table(file = paste("../Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/",pat[x],"/Cluster_colors.txt",sep=""),sep="\t",header = T, stringsAsFactors = F)[, - c(2:3)]
  } else{
    cl_col <- read.table(file = paste("../Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/",pat[x],"/Automatic_cluster_colors.txt",sep=""),sep="\t",header = T, stringsAsFactors = F)
    temp <- unlist(strsplit(x = cl_col[,1], split = paste(pat[x],"__",sep="")))
    temp <- temp[temp != ""]
    cl_col[,1] <- temp
  }
  
  # Read in silico drug efficacy output
  out <- foreach(i = (1:length(lf)), .combine = "rbind") %do% {
    infile <- as.matrix(read.table(file = paste(dir, lf[i], sep=""), sep="\t",header = T, stringsAsFactors = F))
    temp <- unlist(strsplit(strsplit(lf[i], split = "_")[[1]][13], split = ".txt")) # get cluster number
    #infile <- cbind(infile, paste("Cluster_",temp, sep=""))
    infile <- cbind(infile, temp)
    return(infile)
  }
  colnames(out)[9] <- "Cluster"
  out <- out[!is.na(out[,4]),] 
  
  # PLOT FOR z(d(c))
  ################################################################################################
  z_plot <- out[,c(8,4,9)]
  z_plot[z_plot[,1] == " TRUE",1] <- "TRUE" 
  z_plot <- as.data.frame(z_plot)
  z_plot$Cluster <- factor(x = z_plot$Cluster)
  z_plot$Z.c. <- as.numeric(as.character(z_plot$Z.c.))
  z_plot$known_CD_drug <- factor(x = z_plot$known_CD_drug, levels = c("FALSE","TRUE"), ordered = T)
  
  # Basic dot plot
  p <- ggplot(z_plot, aes(x=Cluster, y=Z.c., fill=known_CD_drug)) +
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.03, stackratio = 0.8, position = position_dodge((0.7)), colour = NA) +
    scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey")) +
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black",size = 20),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(-8,8) + scale_y_continuous(breaks = seq(-8,8,2))
  ggsave(filename = paste(dir, "z_scores_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 220, height = 110, units = "mm")
  
  p <- ggplot(z_plot, aes(x=Cluster, y=Z.c., fill=known_CD_drug)) +
    geom_violin(aes(fill = known_CD_drug), scale = 3) +
    scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey"))+
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 20),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(-8,8) + scale_y_continuous(breaks = seq(-8,8,2))
  ggsave(filename = paste(dir, "z_scores_for_known_CD_drugs_violin.pdf", sep=""), plot = p, device = "pdf", width = 220, height = 110, units = "mm")
  
  
  
  # PLOT FOR d(c)
  ################################################################################################
  d_plot <- out[,c(8,3,9)]
  d_plot[d_plot[,1] == " TRUE",1] <- "TRUE" 
  d_plot <- as.data.frame(d_plot)
  d_plot$Cluster <- factor(x = d_plot$Cluster)
  d_plot$mean_d.c. <- as.numeric(as.character(d_plot$mean_d.c.))
  d_plot$known_CD_drug <- factor(d_plot$known_CD_drug, levels = c("FALSE","TRUE"), ordered = T)
  
  # Basic dot plot
  p <- ggplot(d_plot, aes(x=Cluster, y=mean_d.c., fill=known_CD_drug)) +
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.01, stackratio = 0.8, position = position_dodge((0.7)), colour = NA) +
    scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey")) +
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 20),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(0,5) + scale_y_continuous(breaks = seq(0,5,1))
  ggsave(filename = paste(dir, "mean_closest_distance_for_known_CD_drugs.pdf", sep=""), plot = p, device = "pdf", width = 220, height = 110, units = "mm")
  
  # Violin plot
  p <- ggplot(d_plot, aes(x=Cluster, y=mean_d.c., fill=known_CD_drug)) +
    geom_violin(aes(fill = known_CD_drug), scale = 3) +
    scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c("black", "darkgrey"))+
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 20),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(0,5) + scale_y_continuous(breaks = seq(0,5,1))
  ggsave(filename = paste(dir, "mean_closest_distance_for_known_CD_drugs_violin.pdf", sep=""), plot = p, device = "pdf", width = 220, height = 110, units = "mm")
  
  
  # PLOT FOR PRECISION AT z(d(c)) cut-offs
  ################################################################################################
  cl <- unique(out[,9]) # clusters
  n_CD <- length(unique(out[grepl(out[,8],pattern = "TRUE"),1]))
  prec_rec_matrix <- foreach(c = c(1:length(cl)), .combine = "rbind") %do% {
    temp <- out[out[,9]==cl[c],]
    cut_offs <- seq(ceiling(min(as.numeric(temp[,4]))*10)/10,0,0.1)
    prec_rec <- foreach(cut = c(1:length(cut_offs)), .combine = "rbind") %do% {
      temp2 <- temp[as.numeric(temp[,4]) <= cut_offs[cut],8]
      prec <- sum(grepl(temp2, pattern = "TRUE"))/length(temp2)
      recall <- sum(grepl(temp2, pattern = "TRUE"))/n_CD
      return(c(cut_offs[cut],prec,recall))
    }
    prec_rec <- cbind(prec_rec, cl[c])
    return(prec_rec)
  }
  colnames(prec_rec_matrix) <- c("cut_off", "precision", "recall", "Cluster")
  prec_rec_matrix <- as.data.frame(prec_rec_matrix)
  prec_rec_matrix$cut_off <- as.numeric(as.character(prec_rec_matrix$cut_off))
  prec_rec_matrix$precision <- 100*as.numeric(as.character(prec_rec_matrix$precision))
  prec_rec_matrix$recall <- 100*as.numeric(as.character(prec_rec_matrix$recall))
  cl_order <- as.character(c(1, 0, 2, 7, 6, 5, 4, 3))
  prec_rec_matrix$Cluster <- factor(x = as.numeric(as.character(prec_rec_matrix$Cluster)), levels = cl_order[cl_order %in% z_plot$Cluster], ordered = T)
  
  # Precision
  p <- ggplot(prec_rec_matrix, aes(x = cut_off, y = precision, group = Cluster)) + 
    geom_line(aes(color = Cluster, size = 1)) +
    scale_color_manual(values = cl_col[match(as.numeric(as.character(levels(prec_rec_matrix$Cluster))),cl_col[,1]),2]) +
    scale_size(range= c(0,1.5),guide = F) +
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 20),
          legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(filename = paste(dir, "PRECISION_CURVES_for_known_CD_drugs_",pat[x],".pdf", sep=""), plot = p, device = "pdf", width = 150, height = 100, units = "mm")
  
  # RECALL
  p <- ggplot(prec_rec_matrix, aes(x = cut_off, y = recall, group = Cluster)) + 
    geom_line(aes(color = Cluster, size = 1)) +
    scale_color_manual(values = cl_col[match(as.numeric(as.character(levels(prec_rec_matrix$Cluster))),cl_col[,1]),2]) +
    scale_size(range= c(0,1.5),guide = F) +
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 20),
          legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(filename = paste(dir, "RECALL_CURVES_for_known_CD_drugs_",pat[x],".pdf", sep=""), plot = p, device = "pdf", width = 150, height = 100, units = "mm")
  
  # RECALL vs PRECISION
  p <- ggplot(prec_rec_matrix, aes(x = recall, y = precision, group = Cluster)) + 
    geom_line(aes(color = Cluster, size = 1)) +
    scale_color_manual(values = cl_col[match(as.numeric(as.character(levels(prec_rec_matrix$Cluster))),cl_col[,1]),2]) +
    scale_size(range= c(0,1.5),guide = F) +
    theme(plot.background = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 20),
          legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(filename = paste(dir, "RECALL_PRECISION_CURVES_for_known_CD_drugs_",pat[x],".pdf", sep=""), plot = p, device = "pdf", width = 150, height = 100, units = "mm")

  
  
  
  
  # Selection criteria precision and recall plots
  ################################################################################################
  
  z_0.15_criterium <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
    temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
    temp <- temp[as.numeric(temp[,4]) < -0.15,8]
    prec <- sum(grepl(temp, pattern = "TRUE")) / length(temp)
    rec <- sum(grepl(temp, pattern = "TRUE")) / n_CD
    return(c(prec, rec, cl[i], "z_0.15_criterium"))
  }
  
  z_criterium <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
    temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
    temp <- temp[intersect(which(as.numeric(temp[,7]) < 0.05),which(as.numeric(temp[,4]) < 0)),8]
    prec <- sum(grepl(temp, pattern = "TRUE")) / length(temp)
    rec <- sum(grepl(temp, pattern = "TRUE")) / n_CD
    return(c(prec, rec, cl[i], "z_criterium"))
  }
  
  z_and_d_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
    temp <- out[as.numeric(out[,9])==as.numeric(cl[i]),]
    temp <- temp[intersect(which(as.numeric(temp[,7]) < 0.05),which(as.numeric(temp[,4]) < 0)),]
    temp <- temp[as.numeric(temp[,3])<1,8]
    if(length(temp)>0){
      prec <- sum(grepl(temp, pattern = "TRUE")) / length(temp)
      rec <- sum(grepl(temp, pattern = "TRUE")) / n_CD
    } else {
      prec <- rec <- 0
    }
    return(c(prec, rec, cl[i], "z_and_d_criteria"))
  }
  
  fc_criterium_file <- read.table(file = paste(dir, "../../",pat[x],"_FC_criteria_evaluated.txt",sep=""),sep="\t",stringsAsFactors = F,header = T)
  
  z_d_fc_criteria <- foreach(i = c(1:length(cl)), .combine = "rbind") %do% {
    if(any(fc_criterium_file[,18]==cl[i])){
      temp <- fc_criterium_file[fc_criterium_file[,18]==cl[i],]
      temp <- temp[as.numeric(temp[,11]) > 0,8] # choose only drugs that counteract DEGs
      if(length(temp)>0){
        prec <- sum(grepl(temp, pattern = "CD_drug")) / length(temp)
        rec <- sum(grepl(temp, pattern = "CD_drug")) / n_CD
      } else {
        prec <- rec <- 0
      }
    } else {
      prec <- rec <- 0
    }
    return(c(prec, rec, cl[i], "z_d_fc_criteria"))
  }
  
  combined_rec_prec <- rbind(z_0.15_criterium,z_criterium, z_and_d_criteria, z_d_fc_criteria)
  colnames(combined_rec_prec) <- c("Precision", "Recall", "Cluster", "Criteria")
  rownames(combined_rec_prec) <- NULL
  combined_rec_prec <- as.data.frame(combined_rec_prec)
  combined_rec_prec$Precision <- as.numeric(as.character(combined_rec_prec$Precision))*100
  combined_rec_prec$Recall <- as.numeric(as.character(combined_rec_prec$Recall))*100
  combined_rec_prec$Cluster <- factor(x = cl_col[match(as.numeric(combined_rec_prec$Cluster), as.numeric(cl_col[,1])),3], levels = cl_col[,3], ordered = T)
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
                     pattern_spacing = 0.015,
                     #width = 0.01,
                     #pattern_key_scale_factor = 0.01
    ) + 
    scale_pattern_manual(values = c(z_0.15_criterium = "none", z_criterium = "circle", z_and_d_criteria = "stripe", z_d_fc_criteria = "crosshatch"),
                         guide = guide_legend(override.aes = list(fill = "lightgray"))) +
    scale_fill_manual(values = cl_col[cl_col[,3]%in%combined_rec_prec$Cluster,2], aesthetics = "color", guide = "none") +
    scale_fill_manual(values = cl_col[cl_col[,3]%in%combined_rec_prec$Cluster,2], aesthetics = "fill", guide = "none") +
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
  ggsave(filename = paste(dir, "PRECISION_for_known_CD_drugs_SELECTION_CRITERIA_",pat[x],".pdf", sep=""), plot = p, device = "pdf", width = 200, height = 100, units = "mm")
  
  # NEW RECALL PLOT
  p <- ggplot(combined_rec_prec, aes(x = Cluster, y= Recall, fill = Cluster, color = Cluster)) +
    geom_bar_pattern(stat = "Identity",
                     aes(pattern = Criteria),
                     position = position_dodge2(preserve = "single", padding = 0.2),
                     #color = "white", 
                     pattern_fill = "white",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.015,
                     #width = 0.01,
                     #pattern_key_scale_factor = 0.01
    ) + 
    scale_pattern_manual(values = c(z_0.15_criterium = "none", z_criterium = "circle", z_and_d_criteria = "stripe", z_d_fc_criteria = "crosshatch"),
                         guide = guide_legend(override.aes = list(fill = "lightgray"))) +
    scale_fill_manual(values = cl_col[cl_col[,3]%in%combined_rec_prec$Cluster,2], aesthetics = "color", guide = F) +
    scale_fill_manual(values = cl_col[cl_col[,3]%in%combined_rec_prec$Cluster,2], aesthetics = "fill",
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
  ggsave(filename = paste(dir, "RECALL_for_known_CD_drugs_SELECTION_CRITERIA_patient_",pat[x],".pdf", sep=""), plot = p, device = "pdf", width = 200, height = 100, units = "mm")
  
}






