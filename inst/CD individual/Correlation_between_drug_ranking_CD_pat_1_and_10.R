#
# Comparison of drug ranking for CD patient 1 vs CD patient 10
#
# BY SAMUEL SCHAEFER
#
#################################################

fp <- getwd()
setwd(paste(fp, "/..",sep=""))

# LOAD DATA
dir <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
cd1 <- as.matrix(read.table(file = paste(dir,"Patient_1/Drug_ranking/FINAL_drug_ranking_evaluated.txt",sep=""),sep="\t",header = T))
cd10 <- as.matrix(read.table(file = paste(dir,"Patient_10/Drug_ranking/FINAL_drug_ranking_evaluated.txt",sep=""),sep="\t",header = T))

# VENN DIAGRAM of overlap
library(VennDiagram)

cd1_list <- cd1[,1]
cd10_list <- cd10[,1]

venn.diagram(x = list(cd1_list, cd10_list),
             category.names = c("CD patient 1", "CD patient 10"),
             filename = paste(dir,"Overlap_of_drug_candidates_patient_1_and_10.png",sep=""),
             output = T)


library(ggvenn)

temp <- list(pat1 = cd1[,1],
             pat10 = cd10[,1])
p <- ggvenn(temp, fill_color = c("white", "lightgrey"), stroke_size = 1, set_name_size = 4, show_percentage = F)
ggsave(filename = paste(dir, "Overlap_drug_candidates_pat1_and_pat10.pdf",sep=""),plot = p,device = "pdf",width = 3, height = 3)

# Rank correlation scatter plot
library(ggpubr)

cd1 <- cd1[,c(1,3,9)]
cd10 <- cd10[,c(1,3,9)]

cd1 <- cd1[order(as.numeric(cd1[,3]),decreasing = F),]

cd <- merge(cd1, cd10, by = "Drug")
colnames(cd)[c(3,5,2)] <- c("PAT1", "PAT10", "known_CD_drug")
cd <- as.data.frame(cd)
cd$PAT1 <- as.numeric(x = cd$PAT1)
cd$PAT10 <- as.numeric(cd$PAT10)
cd$known_CD_drug <- factor(cd$known_CD_drug, levels = c("FALSE", "CD_drug"))
cd$size <- rep(1.5,times = length(cd$known_CD_drug))
cd$size[cd$known_CD_drug == "CD_drug"] <- 3


# PLOT
cor_plot <- ggscatter(cd, x = "PAT1", y = "PAT10", color = "known_CD_drug", shape = "known_CD_drug",
          label = "Drug", repel = T, label.select = cd$Drug[cd$known_CD_drug == "CD_drug"],
          size = "size", palette = c("black", "#ff0000"))

ggsave(filename = paste(dir,"Correlation_of_drug_rank_for_candidates_CD_pat_1_and_pat_10.pdf",sep=""),plot = cor_plot, device = "pdf", width = 5, height = 4)

