<!--  
#By SAMUEL SCHÄFER
#2022-01-21 
-->

# scPred: Network analyses of single cell-based digital twins for personalized treatment of inflammatory diseases

### scPred: General information

The goal of scPred is the creation of multicellular network-based disease models 
from single-cell RNA-sequencing (scRNA-seq) data for drug repositioning.
scPred was developed using mouse data and validated in in vitro experiments 
and by application to human scRNA-seq data consistently showing a high 
precision for already approved drugs for the investigated diseases.

Analysis was mainly conducted in R 3.4.4.
Input data can be downloaded as described in PMID: XXXXX 

### Introduction to scPred

<br>
<br> 
<img src="vignettes/Overview fig v4.png" width="800" />
<br>
<br>

scPred was developed based on scRNA-seq data from an antigen-induced arthritis mouse model and healthy controls (1), initial data analysis included denoising (2), clustering (3), cell typing and DEG calculation (4). The DEGs derived by comparison cells from healthy and sick mice were then used to create inter- and intracellular disease models. For intercellular disease models, interaction of differentially expressed upstream ligands and downstream genes were inferred using NicheNet ligand activity analysis (5). Based on the resulting cell-cell interaction network the centrality of each cell type, serving as a proxy for therapeutic importance, in the disease could be calculated (6). For intracellular disease models, the DEGs of a cell type alongside the drug targets of known drugs were mapped onto the human protein-protein interaction network and drug candidates for a cell type were selected based on network distance (7). Furthermore, as a proxy for the relevance of drug’s targets, centrality of the drug’s targets in the largest connected component (LCC) formed by a cell type’s DEGs was calculated (8). To rank all drugs, the combined drug target and intercellular centralities of drugs were calculated (9). Combined intercellular centrality score was the sum of the intercellular centrality scores of cell types in which the drug was a candidate. Combined drug target centrality was calculated as the geometric mean of the mean drug target centrality in all cell types in which this drug was selected as a candidate. Crossed out values in the table presenting the drug target and intercellular centrality indicate that a drug was not selected as a candidate in this cell type.

### Creation of 3D network visualization for interaction between plasma cell DEGs and drug candidates

In order to better understand the interactions between DEGs and drug candidates, selected based on zc < -1.64 and dc < 1, we created a 3D [visualisation](https://scpred.shinyapps.io/3D_network/) for the most central cell type (plasma cells) in the antigen induced arthritis mouse data. Interactions between DEGs (blue) are representing protein-protein interactions (PPI) described in the literature-curated PPI network by do Valle et al.[^1]. DEGs node size is based on fold change. Drug candidates are connected to their respective gene drug targets by edges. Potential drug candidates are shown in red. Established drugs for human rheumatoid arthritis are represented in yellow. The higher the absolute value of a drugs Y-axis value the higher the drug rank. Drug candidates that counteracted at least one DEGs fold change received positive Y-axis values while drug candidates that did not counteract the fold change of any targeted DEG received negative Y-axis values.  The visualization was created in R version 4.1.1. The following R packages were used: igraph (1.2.6), plotly (4.10.0) and shiny (1.6.0).


## References

[^1]: do Valle, I. F. et al. Network medicine framework shows that proximity of polyphenol targets and disease proteins predicts therapeutic effects of polyphenols. Nature Food 2, 143-155, doi:10.1038/s43016-021-00243-7 (2021).
