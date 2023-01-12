#
#
# FC vs P value plots for OLD RA DATA cluster 4, 9, 14
#
# SAMUEL SCHAEFER
#################################################################
library(dplyr)
library(ggplot2)
library(ggrepel)
set.seed(54)

#######################################
# Load data
#######################################

# load known RA drug targets
transl <- read.table(file = "../Input/Human-mouse_homologs/transl.txt",sep="\t", header = T, stringsAsFactors = F)
targets <- as.matrix(read.table(file = "../Input/Drugs/drug_targets_unique_literature_ppi.txt",sep="\t",header = T, stringsAsFactors = F))
same_drugs <- read.table(file = "../Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt",sep="\t",header = T, stringsAsFactors = F)
RA_drugs <- as.matrix(read.table(file = "../Input/Drugs/drugs_w_indication_RA_DrugBank_20200214.txt",sep="\t", header = T, stringsAsFactors = F, quote = ""))[,1]
RA_drugs <- unique(c(RA_drugs, same_drugs[same_drugs[,1]%in%RA_drugs,2], same_drugs[same_drugs[,2]%in%RA_drugs,1]))
RA_targets <- unique(as.vector(targets[,colnames(targets)%in% RA_drugs]))
RA_targets <- RA_targets[!is.na(RA_targets)]
RA_targets <- as.character(transl[as.numeric(transl[,2])%in%as.numeric(RA_targets),3])
print(RA_targets)
RA_targets <- RA_targets[order(RA_targets, decreasing = T)]

outdir <- "../Output/DCA_MAST_DEGs/"
files <- list.files(path = outdir)
files <- files[grepl(files, pattern = "_res=0.6_dims=32_k=15")]
files <- files[grepl(files, pattern = ".txt")]


force <- c(30, 40, 50)
y_limit <- 50

'for(i in 1:length(files)){
  degs <- read.delim(file = paste(outdir, "/", files[i], sep=""), sep="\t", header = T)
  head(degs)
  degs <- degs[,c(1,3,6)]
  col <- rep("grey", times = nrow(degs))
  col[intersect(which(degs[,2] < - log(1.5)),which(degs[,3] < 0.05))] <- "blue"
  col[intersect(which(degs[,2] > log(1.5)),which(degs[,3] < 0.05))] <- "red"
  col <- factor(x = col, levels = c("blue", "grey", "red"),ordered = T)
  degs[,3] <- -log10(degs[,3])
  degs <- as.data.frame(degs)
  degs <- cbind(degs, col)
  degs$avg_log2FC <- as.numeric(degs$avg_log2FC)
  
  
  fc_p_plot <- ggplot(data = degs, aes(x = avg_log2FC, y = p_val_adj)) + 
    geom_point(aes(color = col), size = 0.1) + 
    scale_color_manual(values =  c("blue","grey","red")[c("blue","grey","red") %in% as.character(unique(degs$col))]) +
    labs(x = "logFC", y = "-log10(adj P value)") +
    xlim(-2,2) + ylim(0,y_limit) +
    theme(rect = element_blank(),plot.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggrepel::geom_text_repel(data = degs_labeled, aes(label = X), segment.size = 0.5, size = 2, force = force[i],
    min.segment.length = unit(0, "lines"))
  
  ggsave(filename = paste(outdir, "/FC_vs_P_plot_all_RA_targets_", strsplit(files[i], split = ".txt")[[1]][1], ".pdf", sep=""),plot = fc_p_plot, device = "pdf", width = 3, height = 2)
}'

RA_targets <- c("Jak1","Jak2","Jak3","Il6ra","Il1","Il1r1","Il1r2","Il6","Il17","Tnf","Tnfrsf1b","Bcl2","Anxa1","Rxra","Cd80", "Ms4a1", "Nfkbia", "Cd86", "Ccl2", "Icam1", "Vcam1")
RA_targets <- RA_targets[order(RA_targets, decreasing = T)]
force <- 20
y_limit <- 40

for(i in 1:length(files)){
  degs <- read.delim(file = paste(outdir, "/", files[i], sep=""), sep="\t", header = T)
  degs <- degs[,c(1,3,6)]
  col <- rep("grey", times = nrow(degs))
  col[intersect(which(degs[,2] < - log(1.5)),which(degs[,3] < 0.05))] <- "blue"
  col[intersect(which(degs[,2] > log(1.5)),which(degs[,3] < 0.05))] <- "red"
  col <- factor(x = col, levels = c("blue", "grey", "red"),ordered = T)
  degs[,3] <- -log10(degs[,3])
  degs <- as.data.frame(degs)
  degs <- cbind(degs, col)
  degs$avg_log2FC <- as.numeric(degs$avg_log2FC)
  
  degs_labeled <- degs[degs[,1]%in%RA_targets,]
  degs_labeled <- degs_labeled[order(degs_labeled[,1],decreasing = T),]
  degs_labeled <- degs_labeled[degs_labeled[,4] != "grey",]
  
  fc_p_plot <- ggplot(data = degs, aes(x = avg_log2FC, y = p_val_adj)) + 
    geom_point(aes(colour = col),size = 0.1) + 
    scale_color_manual(values =  c("blue","grey","red")[c("blue","grey","red") %in% as.character(unique(degs$col))]) +
    labs(x = "logFC", y = "-log10(adj P value)") +
    xlim(-2,2) + ylim(0,y_limit) +
    theme(rect = element_blank(),plot.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggrepel::geom_text_repel(data = degs_labeled, aes(label = X), segment.size = 0.5, size = 2, force = force,
                             min.segment.length = unit(0, "lines"))
  
  ggsave(filename = paste(outdir, "/FC_vs_P_plot_selected_RA_targets_", strsplit(files[i], split = ".txt")[[1]][1], ".pdf", sep=""),plot = fc_p_plot, device = "pdf", width = 3, height = 2)
}

rm(list = ls())



