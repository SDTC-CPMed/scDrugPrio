
#######################################################################################
# Combine literature evidence for CD drugs
# FOR CD DATA SET - all individuals patients bu patient 1 and 10 that were processed earlier
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/../", sep=""))

library(readxl)
library(foreach)

#######################################################################################
# Input data
#######################################################################################

in_path <- "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/"

cd_1 <- as.matrix(read_xlsx(path = paste(in_path, "Patient_1/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_INDIVIDUAL_CD_PATIENT_1.xlsx", sep =""), sheet = 1))
cd_10 <- as.matrix(read_xlsx(path = paste(in_path, "Patient_10/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_INDIVIDUAL_CD_PATIENT_10.xlsx", sep =""), sheet = 1))
cd_general <- as.matrix(read_xlsx(path = "Output/CD/DCA_MAST_DEGs_predictions/literature_PPI/SUMMARY/FINAL_DRUG_RANKING_CD.xlsx", sheet = 1))

current_literature <- rbind(cd_1[,c(1:2,30:32)], 
                            cd_10[,c(1:2,22:24)],
                            cd_general[,c(1:2,67:69)])
current_literature <- cbind(current_literature, c(rep("Pat_1", times = nrow(cd_1)), rep("Pat_10", times = nrow(cd_10)), rep("General", times = nrow(cd_general))))

current_literature <- current_literature[!is.na(current_literature[,3]),]
current_literature <- current_literature[order(current_literature[,1]),]

# Any conflicts? 
unique_drugs <- unique(current_literature[,1])
for(i in 1:length(unique_drugs)){
  if(sum(current_literature[,1] %in% unique_drugs[i])>1){
    temp <- current_literature[current_literature[,1]%in% unique_drugs[i],]
    for(j in 2:nrow(temp)){
      if(temp[j,3] != temp[1,3]){
        print(paste(temp[1,1], " has a conflict!", sep =""))
      }
    }
  }
}

# If no output, move on!


#######################################################################################
# Combine other patients top 100 in one file
#######################################################################################

pat <- paste("Patient_",1:11,sep="")
pat <- pat[!(pat %in% paste("Patient_",c(1,10),sep=""))]

out <- foreach(i = c(1:length(pat)), .combine = "rbind") %do% {
  temp <- read.table(file = paste(in_path, pat[i], "/Drug_ranking/FINAL_drug_ranking.txt", sep = ""), sep ="\t", header = T, stringsAsFactors = F)
  temp <- temp[temp[,9] <= 100,]
  temp <- cbind(temp, Patient = pat[i])
  return(temp)
}

out <- cbind(out, current_literature[match(out[,1], current_literature[,1]),3:5])
out[is.na(out[,ncol(out)-2]), ncol(out)-2] <- ""
out[is.na(out[,ncol(out)-1]), ncol(out)-1] <- ""
out[is.na(out[,ncol(out)]), ncol(out)] <- ""

write.table(out, file = paste(in_path, "COMBINED_FINAL_TOP_100_INDIVIDUAL_PATIENTS_EXCEPT_1_&_10.txt",sep=""),sep="\t", col.names = T, row.names = F)
