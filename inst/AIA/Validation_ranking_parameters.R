#
#
# Make plots establishing the value of rank aggregation parameters.
########################################################################################################

library(doParallel)
library(ggplot2)
library(reshape2)

# INPUT
################################################################

final_drug_rank <- read.table(file ="../Output/Final_ranking/FINAL_drug_ranking_evaluated.txt", sep="\t", header = T)

df <- data.frame(mean_drug_target_centrality = as.numeric(final_drug_rank[,7]),
                 RA_drug = as.logical(final_drug_rank[,3]),
                 eigenvector_sum= as.numeric(final_drug_rank[,6]),
                 name = as.character(final_drug_rank[,1]),
                 combined_cent = as.numeric(final_drug_rank[,8]),
                 lit_search = as.character(final_drug_rank[,11]))


# rank eigenvector centrality
df <- df[order(df$eigenvector_sum, decreasing = T),]
df$rank_eigenvector_sum <- rank( - df$eigenvector_sum)
# rank drug target centrality
df <- df[order(df$mean_drug_target_centrality, decreasing = T),]
df$rank_mean_drug_target_cent <- rank(- df$mean_drug_target_centrality)
# total rank based on intercellular centrality and drug traget centrality
df$combined_cent <- df$eigenvector_sum - df$mean_drug_target_centrality*0.1
df <- df[order(df$combined_cent, decreasing = T),]
df$rank_combined_cent <- rank( - df$combined_cent)


# Ranking better with cell type averaged geometric mean drug target centrality compared to only ranking by cell type eigenvector centrality sum? 
set.seed(1)
rand_distribution <- drug_rand_distribution <-  semi_rand_distribution_drug_centr <- semi_rand_distribution_drug_centr_lit <- vector()
for(i in 1:1000){ # bootstrap 
  # random sampling of eigenvector centrality rank
  temp <- df
  temp$rank_eigenvector_sum <- sample(df$rank_eigenvector_sum)
  rand_distribution[i] <- mean(temp[temp$RA_drug==T,]$rank_eigenvector_sum)
  
  # random sampling of mean drug centrality rank
  temp <- df
  temp$rank_mean_drug_target_cent <- sample(df$rank_mean_drug_target_cent)
  drug_rand_distribution[i] <- mean(temp[temp$RA_drug==T,]$rank_mean_drug_target_cent)
  
  # random sampling of mean drug target centrality, while using original eigenvector centrality 
  temp <- df[order(df$combined_cent),]
  temp$mean_drug_target_centrality <- sample(temp$mean_drug_target_centrality)  
  temp$combined_cent <- temp$eigenvector_sum - temp$mean_drug_target_centrality*0.1
  # calculate random combined rank from intercellular eigenvector sum and random mean drug target centrality
  temp$rank_combined_cent <- rank(- temp$combined_cent)
  semi_rand_distribution_drug_centr[i] <- mean(temp[temp$RA_drug==T,]$rank_combined_cent)
}
# hist(df$rank_eigenvector_sum[df$RA_drug == T])
# hist(rand_distribution)
# hist(drug_rand_distribution)
# hist(df$rank_combined_cent[df$RA_drug == T])
# hist(semi_rand_distribution_drug_centr)
print("One-sided P-value showing that known RA drugs rank higher when using intercellular centrality ranking compared to random:")
print(paste("P = ", pnorm((mean(df$rank_eigenvector_sum[df$RA_drug == T]) - mean(rand_distribution))/sd(rand_distribution)),
             " mean real = ", mean(df$rank_eigenvector_sum[df$RA_drug == T]), " mean rand = ", mean(rand_distribution), " SD rand = ", sd(rand_distribution), sep = " "))
 
print("One-sided P-value showing that known RA drugs rank higher when using intracellular centrality compared to random:")
print(paste("P = ", pnorm((mean(df$rank_mean_drug_target_cent[df$RA_drug == T]) - mean(drug_rand_distribution))/sd(drug_rand_distribution)),
             " mean real = ", mean(df$rank_mean_drug_target_cent[df$RA_drug == T]), " mean rand = ", mean(drug_rand_distribution), " SD rand = ", sd(drug_rand_distribution), sep = " "))
 
print("One-sided P-value showing that known RA drugs rank higher when using intercellular & intracellular centrality compared to only intercellular  eigenvector centrality sum:")
print(paste("P = ", pnorm((mean(df$rank_combined_cent[df$RA_drug == T]) - mean(semi_rand_distribution_drug_centr))/sd(semi_rand_distribution_drug_centr)),
            " mean real = ", mean(df$rank_combined_cent[df$RA_drug == T]), " mean rand = ", mean(semi_rand_distribution_drug_centr), " SD rand = ", sd(semi_rand_distribution_drug_centr), sep = " "))


print("Correlation between ranking using inter- and intracellular centrality respectively")
print(cor.test(x = df[,7], y = df[,8]))



