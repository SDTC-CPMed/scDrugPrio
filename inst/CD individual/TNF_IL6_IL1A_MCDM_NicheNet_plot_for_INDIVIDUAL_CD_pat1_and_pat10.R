#
#
# VINOVA FIGURE FOR IBD PATIENTS
#
# By SAMUEL SCHAEFER
#
##############################################

fp <- getwd()
setwd(paste(fp, "/..",sep=""))

# INPUT FILES

# anti-TNF non responder 
pat1 <- as.matrix(read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/NicheNet/all_ligand_activity.txt",sep="\t", header = T))
pat1 <- pat1[as.numeric(pat1[,4])>0,]
# TNF responder
pat10 <- as.matrix(read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/NicheNet/all_ligand_activity.txt",sep="\t", header = T))
pat10 <- pat10[as.numeric(pat10[,4])>0,]

# Set Pearson correlation threshold:
# print(max(as.numeric(pat1[,4])))
# print(min(as.numeric(pat1[,4])))
# print(sum(as.numeric(pat1[,4])>0.01) / nrow(pat1))
# 
# print(max(as.numeric(pat10[,4])))
# print(min(as.numeric(pat10[,4])))
# print(sum(as.numeric(pat10[,4])>0.01) / nrow(pat10))
# 
# # Threshold set to +0.01
# pat1 <- pat1[as.numeric(pat1[,4])>0.01,]
# pat10 <- pat10[as.numeric(pat10[,4])>0.01,]
# 
# # number of interactions per patient?
# print(length(table(pat1[,1]))) # 412 different ligands
# print(max(table(pat1[,1]))) # max = 36
# print(min(table(pat1[,1]))) # min = 2
# print(table(pat1[pat1[,1]=="TNF",1])) # TNF = 12 interactions
# print(table(pat1[pat1[,1]=="IL1A",1])) # IL1A = 20 interactions
# print(table(pat1[pat1[,1]=="IL1B",1])) # IL1B = 24 interactions
# print(table(pat1[pat1[,1]=="IL6",1])) # IL6 = 20 interactions
# 
# 
# print(length(table(pat10[,1]))) # 319 different ligands
# print(max(table(pat10[,1]))) # max = 15
# print(min(table(pat10[,1]))) # min = 1
# print(table(pat10[pat10[,1]=="TNF",1])) # TNF = 10 interactions
# print(table(pat10[pat10[,1]=="IL1A",1])) # IL1A = 4 interactions
# print(table(pat10[pat10[,1]=="IL1B",1])) # IL1B = 15 interactions
# print(table(pat10[pat10[,1]=="IL6",1])) # IL6 = 8 interactions

# Set Ligands of interest
lig_of_interest <- c("TNF", "IL6", "IL1A")
pat1 <- pat1[pat1[,1]%in%lig_of_interest,]
pat10 <- pat10[pat10[,1]%in%lig_of_interest,]

# Threshold set to +0.02
# pat1 <- pat1[as.numeric(pat1[,4])>0.02,]
# pat10 <- pat10[as.numeric(pat10[,4])>0.02,]
# 
# # Prepare matrix
# pat1[!(pat1[,1] %in% lig_of_interest),1] <- "Other"
# pat10[!(pat10[,1] %in% lig_of_interest),1] <- "Other"

# Set colors
cl_col_pat1 <- read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/Cluster_colors.txt",sep="\t", header = T, stringsAsFactors = F, comment.char = "")
pat1[,4] <- as.numeric(pat1[,4])
temp <- matrix(NA, ncol = ncol(pat1), nrow = sum(!(paste("Cluster_",cl_col_pat1[,1],sep="") %in% pat1[,6])))
temp[,6] <- paste("Cluster_",cl_col_pat1[!(paste("Cluster_",cl_col_pat1[,1],sep="") %in% pat1[,6]),1],sep="")
pat1 <- rbind(pat1, temp)
pat1 <- cbind(pat1,color1 = NA, color2 = NA)
pat1[,8:9] <- cl_col_pat1[match(pat1[,6], paste("Cluster_",cl_col_pat1[,1],sep="")),4]

cl_col_pat10 <- read.table(file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/Cluster_colors.txt",sep="\t", header = T, stringsAsFactors = F, comment.char = "")
pat10[,4] <- as.numeric(pat10[,4])
temp <- matrix(NA, ncol = ncol(pat10), nrow = sum(!(paste("Cluster_",cl_col_pat10[,1],sep="") %in% pat10[,6])))
temp[,6] <- paste("Cluster_",cl_col_pat10[!(paste("Cluster_",cl_col_pat10[,1],sep="") %in% pat10[,6]),1],sep="")
pat10 <- rbind(pat10, temp)
pat10 <- cbind(pat10,color1 = NA, color2 = NA)
pat10[,8:9] <- cl_col_pat10[match(pat10[,6], paste("Cluster_",cl_col_pat10[,1],sep="")),4]

# Write files
write.table(pat1, file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_1/NicheNet/Processed_NicheNet_lig_for_cytoscape_pat1.txt",sep="\t", col.names = T, row.names = F)
write.table(pat10, file = "Output/CD/Individual_patients/DCA_MAST_DEGs_predictions/Patient_10/NicheNet/Processed_NicheNet_lig_for_cytoscape_pat10.txt",sep="\t", col.names = T, row.names = F)

# write.table(pat1[pat1[,1] != "Other",], file = "Output/CD/Vinova figure for IBD patients/Processed_NicheNet_lig_for_cytoscape_pat1_no_other.txt",sep="\t", col.names = T, row.names = F)
# write.table(pat10[pat10[,1] != "Other",], file = "Output/CD/Vinova figure for IBD patients/Processed_NicheNet_lig_for_cytoscape_pat10_no_other.txt",sep="\t", col.names = T, row.names = F)

# Sum up pearsons of "other"
# sum_up_other <- function(input){
#   input <- input[!is.na(input[,6]),]
#   
#   comb <- unique(paste(input[,5], "___", input[,6],sep=""))
#   out <- vector()
#   for(i in 1:length(comb)){
#     send <- strsplit(comb[i], split = "___")[[1]][1]
#     targ <- strsplit(comb[i], split = "___")[[1]][2]
#     col <- unique(input[intersect(which(input[,5] %in% send),which(input[,6] %in% targ)),7])
#     if(length(intersect(which(input[,5] %in% send),which(input[,6] %in% targ)))==1){
#       temp <- input[intersect(which(input[,5] %in% send),which(input[,6] %in% targ)),]
#       temp <- matrix(temp, nrow = 1)
#     } else {
#       temp <- input[intersect(which(input[,5] %in% send),which(input[,6] %in% targ)),]
#     }
#     if(any(temp[,1]%in% "Other")){
#       temp <- sum(as.numeric(temp[temp[,1] %in% "Other",4]))
#       out <- rbind(out, c("Other", NA, NA, temp, send, targ, col))
#     }
#   }
#   input <- input[input[,1] != "Other",]
#   input <- rbind(input, out)
#   return(input)
# }
# 
# pat1x <- sum_up_other(pat1)
# pat10x <- sum_up_other(pat10)
# 
# write.table(pat1x, file = "Output/CD/Vinova figure for IBD patients/Processed_NicheNet_lig_for_cytoscape_pat1_summed_up_other.txt",sep="\t", col.names = T, row.names = F)
# write.table(pat10x, file = "Output/CD/Vinova figure for IBD patients/Processed_NicheNet_lig_for_cytoscape_pat10_summed_up_other.txt",sep="\t", col.names = T, row.names = F)
# 

'
# Match clusters of individual patients with clusters from pooled patients for quick & dirty cell typing
pat1_cell <- as.matrix(read.table(file = "Input/CD GSE134809/Individual_patients/Cell_identities_patient_1.txt",sep="\t", header = T))
pat10_cell <- as.matrix(read.table(file = "Input/CD GSE134809/Individual_patients/Cell_identities_patient_10.txt",sep="\t", header = T))
pooled_cell <- as.matrix(read.table(file = "Input/CD GSE134809/CD_all/Cell_ID_to_cluster_ID_k=15_res=0.8.txt",sep="\t", header = T))
cell_typing <- as.matrix(read.table(file = "Input/CD GSE134809/Marker_genes/Cell_typing_by_marker_expression.txt",sep="\t", header = T))

# change cell ty?pe labeling
cell_type <- vector()
cell_type[1] <- cell_typing[1,3]
for(i in 2:nrow(cell_typing)){
  if(cell_typing[i,3] %in% cell_typing[1:(i-1),3]){
    n <- 1 + sum(cell_typing[1:(i-1),3] %in% cell_typing[i,3])
    cell_type[i] <- paste(cell_typing[i,3], "_",n,sep="")
  } else {
    cell_type[i] <- cell_typing[i,3]
  }
}
cell_typing[,3] <- cell_type


pat1_cell <- cbind(pat1_cell, pooled_cell[match(pat1_cell[,1],pooled_cell[,1]),2])
pat10_cell <- cbind(pat10_cell, pooled_cell[match(pat10_cell[,1],pooled_cell[,1]),2])

sort_in_table <- function(pat_cell, cell_typing){
  table <- matrix(0, nrow = length(unique(pat_cell[,2])), ncol = length(unique(pat_cell[,3])))
  colnames(table) <- unique(pat_cell[,3])
  rownames(table) <- unique(pat_cell[,2])
  
  for(r in 1:nrow(table)){
    temp <- pat_cell[pat_cell[,2] %in% unique(pat_cell[,2])[r],]
    temp <- table(temp[,3])
    #print(paste("r:",r,sep=" "))
    for(c in 1:ncol(table)){
      #print(paste("c:",c,sep=" "))
      if(any(names(temp)%in%colnames(table)[c])){
        table[r,c] <- temp[(names(temp) %in% unique(pat_cell[,3])[c])]
      }
    }
  }
  
  temp <- colnames(table)
  temp <- cell_typing[match(as.numeric(temp), as.numeric(cell_typing[,1])),3]
  colnames(table) <- temp
  temp <- cell_typing[cell_typing[,3]%in%colnames(table),]
  
  table <- table[,match(temp[,3],colnames(table))]
  
  return(table)
}

pat1_sorted <- sort_in_table(pat1_cell, cell_typing)
pat10_sorted <- sort_in_table(pat10_cell, cell_typing)
write.table(pat1_sorted, file = "Output/CD/Vinova figure for IBD patients/Pat_1_cell_type_inference.txt",sep="\t", col.names = NA, row.names = T)
write.table(pat10_sorted, file = "Output/CD/Vinova figure for IBD patients/Pat_10_cell_type_inference.txt",sep="\t", col.names = NA, row.names = T)
'


