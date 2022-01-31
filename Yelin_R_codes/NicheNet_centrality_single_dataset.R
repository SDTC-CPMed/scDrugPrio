

##########################################################################################################
# Calculation of cell-cell centralities
#
# Assuming NicheNet_analysis.R has been run
##########################################################################################################
#install.packages("igraph")
#install.packages("CINNA")

starttime<-Sys.time()

library(igraph)
library(CINNA)
 

#########################################################
# Input
#########################################################

# LOAD NicheNet ligand interactions
ligands_interaction <- read.table(file = paste(outputdir,"/all_ligand_activity_",datasource,".txt",sep=""), header = T) # cell-cell interactions + ligand information
ligands_interaction <- as.matrix(ligands_interaction) # cell-cell interactions + ligand information

# Create graph
g <- graph_from_edgelist(el = ligands_interaction[,5:6], directed = T)
                 
out <- matrix(NA, nrow = 4, ncol = length(V(g)$name))

rownames(out) <- c("node_degree_all", "node_degree_in", "node_degree_out", "closeness")
colnames(out) <- paste("Cluster_", V(g)$name, sep="")

# Calculate centrality degree
out[1,] <- centr_degree(g, mode = "all")$res # all
out[2,] <- centr_degree(g, mode = "in")$res # in degree
out[3,] <- centr_degree(g, mode = "out")$res # out degree

# Calculate closeness
out[4,] <- closeness(g, mode = "all")

# Calculate centralities

#pr_cent<-proper_centralities(g)
#calculate_centralities(g, include = pr_cent[1:20])  %>% pca_centralities(scale.unit = TRUE)
centrality_matrix <- calculate_centralities(g, include = c("eigenvector centralities", "K-core Decomposition", "Kleinberg's hub centrality scores", 
                                                           "Laplacian Centrality", "Leverage Centrality", "Group Centrality", "Local Bridging Centrality"))

centrality_matrix <- matrix(unlist(centrality_matrix), ncol = length(V(g)$name), byrow = T)
colnames(centrality_matrix) <- paste("Cluster_", V(g)$name, sep="")
rownames(centrality_matrix) <- c("eigenvector centralities", "K-core Decomposition", "Kleinberg's hub centrality scores", 
                                 "Laplacian Centrality", "Leverage Centrality", "Group Centrality", "Local Bridging Centrality")

out <- rbind(out, centrality_matrix)

write.table(out, file=paste(outputdir,"/Cell-cell_centrality_summary_",datasource,".txt",sep=""), sep="\t", col.names = NA, row.names = T)


endtime <-Sys.time()
print(starttime) 
print(endtime)


