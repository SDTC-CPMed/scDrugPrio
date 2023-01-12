#
#
# Bar graph for cell type centrality in MCDM
########################################################################################################

# INPUT
cl_col <- read.table(file ="../Output/Clustering/Cluster_colors.txt", sep="\t", header = T, stringsAsFactors = F)

centr <- read.table(file = "../Output/NicheNet/Cell_type_centrality_summary.txt", sep="\t", header = T, stringsAsFactors = F)
rownames(centr) <- centr[,1]
centr <- t(centr[,-1])

# PREPARE DATA
cl_col <- cbind(cl_col, centrality = centr[match(paste("Cluster_",cl_col[,1],sep=""), rownames(centr)) ,5])
cl_col <- cl_col[!is.na(cl_col[,5]),]

# BAR PLOT
pdf(file = "../Output/NicheNet/Cell_type_centrality_unsorted.pdf")

barplot(height = cl_col[,5], col = cl_col[,4], font.axis = 1, cex.axis = 1.5)

dev.off()

# BAR PLOT SORTED
cl_col <- cl_col[order(cl_col[,5], decreasing = T),]

# BAR PLOT
pdf(file = "../Output/NicheNet/Cell_type_centrality.pdf")

barplot(height = cl_col[,5], col = cl_col[,4], font.axis = 1, cex.axis = 1.5)

dev.off()