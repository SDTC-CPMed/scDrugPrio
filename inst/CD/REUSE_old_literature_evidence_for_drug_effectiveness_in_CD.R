
# RE-USE OLD literature evidence for final prediction
# Run in R 4.0.4
######################################################################

old_lit <- read.table(file = "../../Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_CD_evaluated.txt", sep="\t", header = T, quote = "")
old_lit <- old_lit[,c(1,2,67:69)]

final_drug <- read.table(file = "../Output/Final_ranking/FINAL_drug_ranking.txt",sep="\t", header = T)

final_drug <- cbind(final_drug, old_lit[match(final_drug[,2], old_lit[,2]),3:5])

write.table(final_drug, file = "../Output/Final_ranking/FINAL_drug_ranking_evaluated.txt",sep="\t", col.names = T, row.names = F)