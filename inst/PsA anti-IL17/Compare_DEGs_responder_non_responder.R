
# COMPARE DEGS
# Run in R 4.0.4
####################################################################################

setwd("/data/samsc76/Doctis_data_v2/R")

# LOAD FILES
degs_resp <- read.table(file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_R.txt", sep = "\t", header = T, stringsAsFactors = F)
degs_nonresp <- read.table(file = "../Output/DCA_MAST_DEGs/SUMMARY_MAST_DEGs_all_clusters_NR.txt", sep = "\t", header = T, stringsAsFactors = F)

for(i in 1:ncol(degs_resp)){
  print(paste(colnames(degs_resp)[i], "  n DEGs:", sum(!is.na(degs_resp[,i])), "Responder,", sum(!is.na(degs_nonresp[,i])), "Non-Resp. Shared (% of Responder):", 
              100* sum(degs_resp[!is.na(degs_resp[,i]),i] %in% degs_nonresp[!is.na(degs_nonresp[,i]),i]) / sum(!is.na(degs_resp[,i])), sep = " "))
}