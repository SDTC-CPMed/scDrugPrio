
#######################################################################################
# # REUSE POOLED CD BATCH CORRECTION LITERATURE SEARCH FOR INDIVIDUAL PATIENTS
# 
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/../", sep=""))

install.packages("openxlsx")
library(openxlsx)

#######################################################################################
# Input data
#######################################################################################

cd_pooled <- as.matrix(read_xlsx(path = "CD_batch_corrected/Output/Final_ranking/FINAL_drug_ranking_evaluated.xlsx", sheet = 1))
cd_pooled <- cd_pooled[,c(1,2,10:12)]

# Conflicts?
unique_drugs <- unique(cd_pooled[,1])
for(i in 1:length(unique_drugs)){
  if(sum(cd_pooled[,1] %in% unique_drugs[i])>1){
    temp <- cd_pooled[cd_pooled[,1]%in% unique_drugs[i],]
    for(j in 2:nrow(temp)){
      if(temp[j,3] != temp[1,3]){
        print(paste(temp[1,1], " has a conflict!", sep =""))
      }
    }
  }
}

cd_pooled <- cd_pooled[!duplicated(cd_pooled[,1]),]
cd_pooled <- cd_pooled[!is.na(cd_pooled[,3]),]

in_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"
final_rank <- as.matrix(read_xlsx(path = paste(in_path, "COMBINED_FINAL_TOP_100_INDIVIDUAL_PATIENTS_EXCEPT_1_&_10_evaluated.xlsx", sep=""), sheet = 1))

#######################################################################################
# Make sure that above is consistent with final_rank
#######################################################################################

final_rank2 <- final_rank
final_rank2[final_rank2[,1] %in% cd_pooled[,1],11:13] <- cd_pooled[match(final_rank2[final_rank2[,1] %in% cd_pooled[,1],1], cd_pooled[,1]),3:5]
final_rank2 <- as.data.frame(final_rank2)

write.xlsx(final_rank2, file = paste(in_path, "COMBINED_FINAL_TOP_100_INDIVIDUAL_PATIENTS_EXCEPT_1_&_10_evaluated_checked.xlsx",sep=""), 
           colNames = T, rowNames = F, append = F, sheetName = "Sheet1")




