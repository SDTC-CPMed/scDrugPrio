

#######################################################################################
# Script for checking node degree of drugs for rebuttal
# 2022-10-10 
# BY SAMUEL SCHAEFER
#######################################################################################

fp <- getwd()
setwd(paste(fp, "/../", sep=""))

library(doParallel)
library(igraph)
library(ggplot2)
library(ggridges)

#######################################################################################
# LOAD files
#######################################################################################

# Load drug_matrix 
drugs <- read.table("Input/Drugs/drug_targets_unique_literature_ppi.txt", sep="\t", header = T)
drugs <- as.matrix(drugs)
mode(drugs) <- "character"
ident_drugs <- read.table(file = "Input/Drugs/drugs_with_identical_drug_target_combinations_literature_PPI.txt", sep ="\t", stringsAsFactors = F, header = T)
# Load predefined drug list & format to vector
ra_drug <- as.vector(read.table("Input/Drugs/predefined_RA_drugs.txt", sep="\t", header = F)[,1])
ms_drug <- as.vector(read.table(file = "Input/Drugs/MS_drugs_from_DrugBank.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
cd_drug <- as.vector(read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
psa_drug <- as.vector(read.table(file = "Input/Drugs/PsA_drugs_from_DrugBank.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
# Load DrugBank data
drugbank <- as.matrix(read.table(file = "Input/Drugs/all_drug_targets_drug_bank.txt", sep="\t", header = T, stringsAsFactors = F))
temp <- unique(drugbank[endsWith(drugbank[,2], "ab"),1:2])
write.table(temp, file ="Input/Drugs/antibody_drugs_for_network_degrees.txt",sep="\t", col.names =  T)
antibody_drug <- unique(drugbank[endsWith(drugbank[,2], "ab"),1])
print(length(antibody_drug)) # n = 316

# LOAD PPI
ppi <- as.matrix(read.table(file = "Input/literature_PPI/ppi.txt", sep="\t", header = T, stringsAsFactors = F))
ppi[,1] <- as.character(ppi[,1])
ppi[,2] <- as.character(ppi[,2])

#######################################################################################
# Get node degree
#######################################################################################
ppin <- graph_from_edgelist(el = ppi, directed = F)

degrees <- degree(ppin)
degrees <- cbind(gene = names(degrees), degree = degrees)

drug_degree <- foreach(i = c(1:ncol(drugs)), .combine = "cbind") %do% {
  out <- drugs[,i]
  out <- as.numeric(degrees[match(out, degrees[,1]),2])
  return(out)
}
colnames(drug_degree) <- colnames(drugs)

# for all individual drugs
ident_drugs <- ident_drugs[ident_drugs[,1]%in% colnames(drug_degree),]
drug_degree <- drug_degree[,match(ident_drugs[,1], colnames(drug_degree))]
colnames(drug_degree) <- ident_drugs[,2]

#######################################################################################
# sum of all drug target's degrees
#######################################################################################

all_drugs <- colSums(drug_degree, na.rm = T)
all_drugs <- log10(all_drugs)
hist(all_drugs)

background <- as.numeric(degrees[,2])
background <- log10(background)
hist(background)

ra_drug <- all_drugs[names(all_drugs) %in% ra_drug]
hist(ra_drug)
ms_drug <- all_drugs[names(all_drugs) %in% ms_drug]
hist(ms_drug)
cd_drug <- all_drugs[names(all_drugs) %in% cd_drug]
hist(cd_drug)
psa_drug <- all_drugs[names(all_drugs) %in% psa_drug]
hist(psa_drug)

plot_data <- rbind(cbind(degree = background, group = "Background PPIN"),
                   cbind(degree = all_drugs, group = "All DrugBank drugs"),
                   cbind(degree = ra_drug, group = "Approved RA drugs"),
                   cbind(degree = ms_drug, group = "Approved MS drugs"),
                   cbind(degree = cd_drug, group = "Approved CD drugs"),
                   cbind(degree = psa_drug, group = "Approved PsA drugs"))
plot_data <- as.data.frame(plot_data)
plot_data$degree <- as.numeric(plot_data$degree)
plot_data$group <- factor(x = as.character(plot_data$group), 
                          levels = c("Background PPIN", "All DrugBank drugs", "Approved RA drugs", "Approved MS drugs", "Approved CD drugs", "Approved PsA drugs"))

my_colors <- c("#553c93", "#2b58a3", "#2f9da7", "#4F7942", "#AFE1AF", "#e7772b")

# Plot
dir.create(path = "Output/Drug_target_centrality/",showWarnings = F)

p <- ggplot(plot_data, aes(degree, group)) + 
  geom_density_ridges(aes(fill = group), alpha = 0.8) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(limits = c(0,5)) +
  scale_fill_manual(values = my_colors) +
  labs(x="log10 network degree", y="Counts") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p

ggsave(filename = "Output/Drug_target_centrality/Summed_degree_of_drugs_targets.pdf", plot = p, device = "pdf", width = 7, height = 4.5)


#######################################################################################
# sum of all drug target's degrees WITH ANTIBODY DRUGS
#######################################################################################

all_drugs <- colSums(drug_degree, na.rm = T)
all_drugs <- log10(all_drugs)
hist(all_drugs)

antibody_drugs <- all_drugs[names(all_drugs) %in% antibody_drug]
hist(antibody_drugs)

plot_data <- rbind(cbind(degree = background, group = "Background PPIN"),
                   cbind(degree = all_drugs, group = "All DrugBank drugs"),
                   cbind(degree = ra_drug, group = "Approved RA drugs"),
                   cbind(degree = ms_drug, group = "Approved MS drugs"),
                   cbind(degree = cd_drug, group = "Approved CD drugs"),
                   cbind(degree = psa_drug, group = "Approved PsA drugs"),
                   cbind(degree = antibody_drugs, group = "Antibody drugs"))
plot_data <- as.data.frame(plot_data)
plot_data$degree <- as.numeric(plot_data$degree)
plot_data$group <- factor(x = as.character(plot_data$group), 
                          levels = c("Background PPIN", "All DrugBank drugs", "Approved RA drugs", "Approved MS drugs", "Approved CD drugs", "Approved PsA drugs", "Antibody drugs"))

my_colors <- c("#553c93", "#2b58a3", "#2f9da7", "#4F7942", "#AFE1AF", "#e7772b", "#A9A9A9")

# Plot
dir.create(path = "Output/Drug_target_centrality/",showWarnings = F)

p <- ggplot(plot_data, aes(degree, group)) + 
  geom_density_ridges(aes(fill = group), alpha = 0.8) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(limits = c(0,5)) +
  scale_fill_manual(values = my_colors) +
  labs(x="log10 network degree", y="Counts") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p

ggsave(filename = "Output/Drug_target_centrality/Summed_degree_of_drugs_targets_with_antibody_drugs.pdf", plot = p, device = "pdf", width = 7, height = 5.5)

#######################################################################################
# mean of all drug target's degrees
#######################################################################################

# Load predefined drug list & format to vector
ra_drug <- as.vector(read.table("Input/Drugs/predefined_RA_drugs.txt", sep="\t", header = F)[,1])
ms_drug <- as.vector(read.table(file = "Input/Drugs/MS_drugs_from_DrugBank.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
cd_drug <- as.vector(read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
psa_drug <- as.vector(read.table(file = "Input/Drugs/PsA_drugs_from_DrugBank.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])

all_drugs <- colMeans(drug_degree, na.rm = T)
all_drugs <- log10(all_drugs)
hist(all_drugs)

background <- as.numeric(degrees[,2])
background <- log10(background)
hist(background)

ra_drug <- all_drugs[names(all_drugs) %in% ra_drug]
hist(ra_drug)
ms_drug <- all_drugs[names(all_drugs) %in% ms_drug]
hist(ms_drug)
cd_drug <- all_drugs[names(all_drugs) %in% cd_drug]
hist(cd_drug)
psa_drug <- all_drugs[names(all_drugs) %in% psa_drug]
hist(psa_drug)

plot_data <- rbind(cbind(degree = background, group = "Background PPIN"),
                   cbind(degree = all_drugs, group = "All DrugBank drugs"),
                   cbind(degree = ra_drug, group = "Approved RA drugs"),
                   cbind(degree = ms_drug, group = "Approved MS drugs"),
                   cbind(degree = cd_drug, group = "Approved CD drugs"),
                   cbind(degree = psa_drug, group = "Approved PsA drugs"))
plot_data <- as.data.frame(plot_data)
plot_data$degree <- as.numeric(plot_data$degree)
plot_data$group <- factor(x = as.character(plot_data$group), 
                          levels = c("Background PPIN", "All DrugBank drugs", "Approved RA drugs", "Approved MS drugs", "Approved CD drugs", "Approved PsA drugs"))

my_colors <- c("#553c93", "#2b58a3", "#2f9da7", "#4F7942", "#AFE1AF", "#e7772b")

# Plot

p <- ggplot(plot_data, aes(degree, group)) + 
  geom_density_ridges(aes(fill = group), alpha = 0.8) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(limits = c(0,5)) +
  scale_fill_manual(values = my_colors) +
  labs(x="log10 network degree", y="Counts") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p

ggsave(filename = "Output/Drug_target_centrality/Mean_degree_of_drugs_targets.pdf", plot = p, device = "pdf", width = 7, height = 4.5)



#######################################################################################
# mean of all drug target's degrees WITH ANTIBODY DRUGS
#######################################################################################

# Load predefined drug list & format to vector
antibody_drugs <- all_drugs[names(all_drugs) %in% antibody_drug]
hist(antibody_drugs)

plot_data <- rbind(cbind(degree = background, group = "Background PPIN"),
                   cbind(degree = all_drugs, group = "All DrugBank drugs"),
                   cbind(degree = ra_drug, group = "Approved RA drugs"),
                   cbind(degree = ms_drug, group = "Approved MS drugs"),
                   cbind(degree = cd_drug, group = "Approved CD drugs"),
                   cbind(degree = psa_drug, group = "Approved PsA drugs"),
                   cbind(degree = antibody_drugs, group = "Antibody drugs"))
plot_data <- as.data.frame(plot_data)
plot_data$degree <- as.numeric(plot_data$degree)
plot_data$group <- factor(x = as.character(plot_data$group), 
                          levels = c("Background PPIN", "All DrugBank drugs", "Approved RA drugs", "Approved MS drugs", "Approved CD drugs", "Approved PsA drugs", "Antibody drugs"))

my_colors <- c("#553c93", "#2b58a3", "#2f9da7", "#4F7942", "#AFE1AF", "#e7772b", "#A9A9A9")

# Plot

p <- ggplot(plot_data, aes(degree, group)) + 
  geom_density_ridges(aes(fill = group), alpha = 0.8) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(limits = c(0,5)) +
  scale_fill_manual(values = my_colors) +
  labs(x="log10 network degree", y="Counts") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p

ggsave(filename = "Output/Drug_target_centrality/Mean_degree_of_drugs_targets_with_antibody_drugs.pdf", plot = p, device = "pdf", width = 7, height = 5.5)


#######################################################################################
# max of all drug target's degrees
#######################################################################################

# Load predefined drug list & format to vector
ra_drug <- as.vector(read.table("Input/Drugs/predefined_RA_drugs.txt", sep="\t", header = F)[,1])
ms_drug <- as.vector(read.table(file = "Input/Drugs/MS_drugs_from_DrugBank.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
cd_drug <- as.vector(read.table(file = "Input/Drugs/CD_drugs_from_DrugBank_201015.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])
psa_drug <- as.vector(read.table(file = "Input/Drugs/PsA_drugs_from_DrugBank.txt", sep = "\t", stringsAsFactors = F, quote = "")[,1])

all_drugs <- foreach(i = c(1:ncol(drug_degree))) %do% {
  temp <- drug_degree[,i]
  temp <- temp[!is.na(temp)]
  return(max(temp))
}
all_drugs <- unlist(all_drugs)
names(all_drugs) <- colnames(drug_degree)
all_drugs <- log10(all_drugs)
hist(all_drugs)

background <- as.numeric(degrees[,2])
background <- log10(background)
hist(background)

ra_drug <- all_drugs[names(all_drugs) %in% ra_drug]
hist(ra_drug)
ms_drug <- all_drugs[names(all_drugs) %in% ms_drug]
hist(ms_drug)
cd_drug <- all_drugs[names(all_drugs) %in% cd_drug]
hist(cd_drug)
psa_drug <- all_drugs[names(all_drugs) %in% psa_drug]
hist(psa_drug)

plot_data <- rbind(cbind(degree = background, group = "Background PPIN"),
                   cbind(degree = all_drugs, group = "All DrugBank drugs"),
                   cbind(degree = ra_drug, group = "Approved RA drugs"),
                   cbind(degree = ms_drug, group = "Approved MS drugs"),
                   cbind(degree = cd_drug, group = "Approved CD drugs"),
                   cbind(degree = psa_drug, group = "Approved PsA drugs"))
plot_data <- as.data.frame(plot_data)
plot_data$degree <- as.numeric(plot_data$degree)
plot_data$group <- factor(x = as.character(plot_data$group), 
                          levels = c("Background PPIN", "All DrugBank drugs", "Approved RA drugs", "Approved MS drugs", "Approved CD drugs", "Approved PsA drugs"))

# Plot

p <- ggplot(plot_data, aes(degree, group)) + 
  geom_density_ridges(aes(fill = group), alpha = 0.8) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(limits = c(0,5)) +
  scale_fill_manual(values = my_colors) +
  labs(x="log10 network degree", y="Counts") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p

ggsave(filename = "Output/Drug_target_centrality/Max_degree_of_drugs_targets.pdf", plot = p, device = "pdf", width = 7, height = 4.5)


#######################################################################################
# max of all drug target's degrees WITH ANTIBODY DRUGS
#######################################################################################

# Load predefined drug list & format to vector
antibody_drugs <- all_drugs[names(all_drugs) %in% antibody_drug]
hist(antibody_drugs)

plot_data <- rbind(cbind(degree = background, group = "Background PPIN"),
                   cbind(degree = all_drugs, group = "All DrugBank drugs"),
                   cbind(degree = ra_drug, group = "Approved RA drugs"),
                   cbind(degree = ms_drug, group = "Approved MS drugs"),
                   cbind(degree = cd_drug, group = "Approved CD drugs"),
                   cbind(degree = psa_drug, group = "Approved PsA drugs"),
                   cbind(degree = antibody_drugs, group = "Antibody drugs"))
plot_data <- as.data.frame(plot_data)
plot_data$degree <- as.numeric(plot_data$degree)
plot_data$group <- factor(x = as.character(plot_data$group), 
                          levels = c("Background PPIN", "All DrugBank drugs", "Approved RA drugs", "Approved MS drugs", "Approved CD drugs", "Approved PsA drugs", "Antibody drugs"))

my_colors <- c("#553c93", "#2b58a3", "#2f9da7", "#4F7942", "#AFE1AF", "#e7772b", "#A9A9A9")

# Plot

p <- ggplot(plot_data, aes(degree, group)) + 
  geom_density_ridges(aes(fill = group), alpha = 0.8) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(limits = c(0,5)) +
  scale_fill_manual(values = my_colors) +
  labs(x="log10 network degree", y="Counts") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

p

ggsave(filename = "Output/Drug_target_centrality/Max_degree_of_drugs_targets_with_antibody_drugs.pdf", plot = p, device = "pdf", width = 7, height = 5.5)




