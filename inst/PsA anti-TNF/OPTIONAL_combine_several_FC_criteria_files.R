#
#
# COMBINE SEVERAL PARTS FOR FC CHECKING
#
# BY SAMUEL SCHAEFER
##########################################################################


library(readxl)

# INPUT
part1 <- read.delim(file = "../Output/FC_criteria_checking/Non-responder/SUMMARY_only_drugs_passing_network_criteria_part_1_evaluated.txt",header = T, sep = "\t")[,1:20]
full <- read.table(file = "../Output/FC_criteria_checking/Non-responder/SUMMARY_only_drugs_passing_network_criteria.txt", sep="\t", header = T)

# Merge files
full <- cbind(full[,1:11], NA,NA,NA,NA,NA, full[,12:13], NA, full[,14])
colnames(full) <- colnames(part1)
full <- full[!(full[,20]%in% part1[,20]),]
full <- rbind(part1, full)

write.table(full, file = "../Output/FC_criteria_checking/Non-responder/SUMMARY_only_drugs_passing_network_criteria_evaluated.txt",sep="\t", col.names = T, row.names = F)
