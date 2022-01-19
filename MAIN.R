# R script for coordinating all scripts

fp <- getwd()
# Setup directories % R environment
source(paste(fp, "/SETUP.R", sep=""))


# Create / format data
##################################################

# Apply DCA (deep count autoencoder) to scRNA-seq matrix (python)

# Calculate clusters based on latent DCA features & cell type clusters
# Calculate DEGs between healthy & sick cells within each cluster using MAST algorithm
# Cell type clusters
source(paste(fp, "/", sep=""))
# Calculate randomized LASSO, dependent variable = arthritis scores (for each cluster)
source(paste(fp, "/random_LASSO_on_clustered_DCA_genes.R")) # not done yet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Identify genes correlating with arthritis scores (in each cluster)
source(paste(fp, "/cluster_genes_correlating_w_arthritis_score.R", sep="")) # done
# Identify genes multidimensional threshold genes
source(paste(fp, "/", sep=""))
# Create drug matrix
source(paste(fp, "/", sep=""))
# Create bulk-RNAseq DEGs
source(paste(fp, "/bulk_RNA_seq_DEGs.R", sep="")) # done
# Extract GWAS genes from GWAS cataloge
source(paste(fp, "/gwas_extraction.R", sep="")) # done
# Download & format KEGG pathways
source(paste(fp, "/kegg_pathways_formatting.R",sep="")) # done


# translate mouse results to human
source(paste(fp, "/translate_mouse_to_human.R",sep=""))

# Analyze data
##################################################
# Screening DCA MAST DEGs
source(paste(fp,"/drug_prediction_based_on_DCA_MAST_DEGs.R", sep="")) # done
# Screening x top DCA MAST DEGs
source(paste(fp, "/drug_prediction_based_on_x_top_DCA_MAST_DEGs.R", sep="")) # done
# Screening LCCs formed by DCA MAST DEGs
source(paste(fp, "/drug_prediction_based_on_LCCs_from_DCA_MAST_DEGs.R", sep="")) # done
# Screening overlap between cluster-specific DCA MAST DEGs
source(paste(fp, "/drug_prediction_based_on_DCA_MAST_DEGs_overlapping_between_x_clusters.R", sep="")) # done
# Screening overlap between x LCCs formed by DCA MAST DEGs
source(paste(fp, "/drug_prediction_based_on_LCCs_of_overlapping_DCA_MAST_DEGs.R")) # done
# Screening x top pathways enriched by DCA MAST DEGs
source(paste(fp, "/drug_prediction_based_on_x_top_DCA_MAST_DEGs_enriched_KEGG_pathways.R", sep="")) # done
# Screening overlap of pathways enriched by x clusters DCA MAST DEGs
source(paste(fp, "/drug_prediction_based_on_overlapping_pathways.R", sep="")) # done


# CORRELATING GENES (DCA genes correlating to arthritis score - Oleg algorithm)
######################
# Screening correlation genes
source(paste(fp, "/drug_prediction_based_on_correlation_genes_DCA_normalized.R", sep="")) # done
# Screening x top DCA normalized correlating genes
source(paste(fp, "/drug_prediction_based_on_x_top_DCA_normalized_correlation_genes.R", sep="")) # done
# Screening LCCs formed by DCA normalized correlating genes
source(paste(fp, "/drug_prediction_based_on_LCCs_from_DCA_normalized_correlating_genes.R", sep="")) # done
# Screening overlap between cluster-specific correlating genes
source(paste(fp, "/drug_prediction_based_on_DCA_normaized_correlating_genes_overlapping_between_x_clusters.R", sep="")) # done
# Screening LCCs formed by overlapping correlating genes
source(paste(fp, "/drug_prediction_based_on_LCCs_of_overlapping_DCA_normalized_correlating_genes.R", sep="")) # done
# Screening x top pathways enriched by correlating genes
source(paste(fp, "/drug_prediction_based_on_x_top_DCA_normalized_correlating_genes_enriched_KEGG_pathways.R", sep="")) # done
# Screening overlap of pathways enriched by x clusters correlating genes
source(paste(fp, "/drug_prediction_based_on_overlapping_pathways_enriched_with_DCA_normalized_correlating_genes.R", sep="")) # done


# DCA MULTIDIMENSIONAL THRESHOLDS
######################
source()

# RANDOMIZED LASSO GENES: dependent variable= arthritis score
######################
source()

# Screening GWAS, OMIM, GWAS & OMIM, bulk-RNAseq DEGs
######################
# Screening benchmark genes
source(paste(fp, "/drug_prediction_based_on_benchmark_genes.R", sep="")) # done
# Screening x top bulk RNA-seq DEGs
source(paste(fp, "/drug_prediction_based_on_x_top_bulk_RNA-seq_DEGs", sep="")) # done
# Screening LCCs formed by benchmark genes
source(paste(fp, "/drug_prediction_based_on_LCCs_formed_by_benchmarking_genes.R", sep="")) # done
# Screening overlap between bulk RNA-seq data sets
source(paste(fp, "/drug_prediction_based_on_bulk_RNA-seq_DEGs_overlapping_between_data_sets.R", sep="")) # done
# Screening LCCs formed by overlapping bulk genes
source(paste(fp, "/drug_prediction_based_on_LCCs_from_overlapping_bulk_RNA-seq_DEGs.R", sep="")) # done
# Screening enriched pathways overlapping between bulk data sets
source(paste(fp, "/drug_prediction_based_on_KEGG_pathways_overlapping_between_bulk_RNA-seq_data_sets.R", sep="")) # done
# Screening x top pathways enriched by benchmark genes
source(paste(fp, "/drug_prediction_based_on_x_top_KEGG_pathways_enriched_with_benchmarking_genes.R", sep="")) # done


# Summarize analysis
##################################################
# Recalculate precision / recall for RA drug target combinations & individual drugs using new RA drug list as well as EULAR drug list
# Also creates precision-recall plots using "summary_plots_drug_predictions.R"
source(paste(fp, "/OPTIONAL_recalculate_precision_recall_matrix_for_individual_RA_drugs_and_unique_drug_target_combinations", sep=""))
# Apply additional criteria to drug predictions & create summaries
source(paste(fp, "/OPTIONAL_creation_of_drug_prediction_summary.R", sep=""))
# Overlap between best predictions?
source(paste(fp, "/OPTIONAL_overlap_between_best_predictions.R", sep=""))


# Investigate overlap between summaries




