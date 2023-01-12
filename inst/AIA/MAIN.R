
# This script indicates the order in which scripts were run to attain results.
# If needed, this script will download the necessary meta-data from figshare.


# FUNDAMENTAL RESULTS
############################################################################
############################################################################

# Quality control and filtering of AIA scRNA-seq data was performed prior to the analysis. Cut-offs are specified in the Supplementary files. 

# Required input data need to be downloaded, including the DCA de-noised scRNA-seq matrix. 
# As this is not possible without authenticating to figshare, we encourage the user to manually download the required files. 
# Upon unzipping the downloaded input and output files, codes should be deposit in the "your_file_path/R" directory, which should be made the working directory.

# When this is done, codes can be run in the following order.

source("Clustering_and_DEGs.R")

source("NicheNet_analysis.R")

source("Intracellular_centrality.R")

source("Network_predictions.R")

source("FC_evaluation.R")

# Here the FC criteria evaluation was performed manually. This file will need to be downloaded before continuing.

source("Final_drug_prioritization.R")

# Again a literature search was carried out for the top 100 ranking drugs. This file will need to be downloaded before continuing.

source("Validation_ranking_parameters.R")

source("Precision_among_ranked_candidates_AIA.R")

source("Precision_among_ranked_candidates_AIA_top100.R.R")


# ADDITIONAL VALIDATIONS
############################################################################
############################################################################

# Establishing Importance of MCDM interactions

source("Establishing_impoartance_of_MCDM_interactions.R")

# Enrichment AIA DEGs with RA GWAS genes

source("Enrichment_AIA_DEGs_with_RA_GWAS_vs_MCDM_centrality_new.R")

# Enrichment AIA DEGs with RA GWAS genes

source("Enrichement_B_cell_DEGs_w_RA_drug_targets.R")

# Heterogenity in latent space

source("Comparison_of_heterogeneity_in_latent_space.R")

# Precision and recall at different zc cut-offs

source("AIA_drug_selection_criteria_and_precision_at_different_zc.R")


# Figures
source("Bar_graph_cell_type_centrality.R")
source("Bar_graph_RA_drugs_among_reoccuring_candidates.R")
source("AIA_B_cell_clusters_FC_vs_P_plots.R")

# OTHER
############################################################################
############################################################################

source("Summary_mouse_RA_DEGS_FOR_SUPPLEMENT.R")


