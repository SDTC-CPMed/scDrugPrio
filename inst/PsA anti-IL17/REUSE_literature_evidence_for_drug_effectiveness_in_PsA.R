
# RE-USE OLD literature evidence for final prediction
# Run in R 4.0.4
######################################################################

library(readxl)

# anti-TNF responders
old_lit <- as.matrix(read_xlsx(path = "../../Doctis_v9/Output/Final_ranking/Responder/FINAL_drug_ranking_Responder_evaluated.xlsx", sheet = 1, col_names = T))
old_lit <- old_lit[,c(1,2,10:12)]


# APPLY TO anti-IL17 RESPONDERS
final_drug <- read.table(file = "../Output/Final_ranking/Responder/FINAL_drug_ranking_Responder.txt",sep="\t", header = T)
final_drug <- cbind(final_drug, old_lit[match(final_drug[,2], old_lit[,2]),3:5])
for(i in 10:12){
  final_drug[is.na(final_drug[,i]),i] <- ""
}
write.table(final_drug, file = "../Output/Final_ranking/Responder/FINAL_drug_ranking_Responder_evaluated.txt",sep="\t", col.names = T, row.names = F)

# APPLY TO anti-IL17 NON-RESPONDERS
old_lit2 <- as.matrix(read_xlsx(path = "../Output/Final_ranking/Responder/FINAL_drug_ranking_Responder_evaluated.xlsx", sheet = 1, col_names = T))
old_lit2 <- old_lit2[,c(1,2,10:12)]
old_lit <- rbind(old_lit2[!is.na(old_lit2[,3]),], old_lit[!is.na(old_lit2[,3]),])
old_lit <- old_lit[! duplicated(old_lit[,1]),]

final_drug <- read.table(file = "../Output/Final_ranking/Non_responder/FINAL_drug_ranking_Non_responder.txt",sep="\t", header = T)
final_drug <- cbind(final_drug, old_lit[match(final_drug[,2], old_lit[,2]),3:5])
for(i in 10:12){
  final_drug[is.na(final_drug[,i]),i] <- ""
}
write.table(final_drug, file = "../Output/Final_ranking/Non_responder/FINAL_drug_ranking_Non_responder_evaluated.txt",sep="\t", col.names = T, row.names = F)
