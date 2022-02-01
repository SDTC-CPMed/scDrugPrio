## Overview

scPred presents a strategy for drug repositioning based on scRNA-seq
based, multidimensional and multicellular disease models that
incorporate the key biological and pharmacological properties. scPred
was primarily developed using data from a mouse model (antigen induced
arthritis, AIA) of rheumatoid arthritis (RA) and validated in *in vitro*
experiments and by prediction precision for known drug-disease pairs in
AIA and two human scRNA-seq data sets (multiple sclerosis, MS, and
Crohn’s disease, CD).

In this vignette we will go through all steps required to format the
input data, create intercellular and intracellular disease models,
select drug candidates for each cell type and create a final ranking for
all the candidates.

Workflow of scPred: 1. Setup - Installation of scPred 2. Data formatting
- Make drug target matrix - Extract largest connected component (LCC) of
protein-protein interaction network (PPIN) - Select unique drug target
combinations found in LCC of PPIN 3. Deep count autoencoder denoising of
expression data 4. Formatting of single cell RNA-sequencing data -
Seurat for clustering - MAST framework for calculation of differentially
expressed genes - Cell typing using marker genes 5.

This vignette guides you in detail through all these steps. As example
data, we are using the scRNA-seq data for the antigen induced arthritis
(AIA) model.

## Setup and data formatting

### Installation of scPred

``` r
# From GitHub:
install.packages("devtools")
devtools::install_github("SDTC-CPMed/scPred")

# or local installation:
install.packages("scPred_1.0.0.tar.gz", repos = NULL, type = "source")
library(scPred)
```

Should difficulties occur when installing the dependencies for scPred
try `source("scPred/inst/SETUP.R")`.

## Data formatting

scPred comes equipped with some raw data files in the directory
‘data-raw’. Additionally, more raw data files as well as processed data
files are made available at [zenodo](zenodo.com).

``` r
setwd("scPred")
```

``` r
lf <- list.files(path = "data-raw/")
print(lf)
#>  [1] "AIA"                               "all_drug_targets_drug_bank.txt"   
#>  [3] "CD"                                "CD_drugs_from_DrugBank_201015.txt"
#>  [5] "DCA_adjusted_matrix"               "HGNC translation matrix 201108"   
#>  [7] "HGNC_transl.txt"                   "Human-mouse_homologs"             
#>  [9] "HuRI_ppi.txt"                      "INDIVIDUAL_CD"                    
#> [11] "ligand_target_NicheNet.txt"        "lit_ppi.txt"                      
#> [13] "MS"                                "MS_drugs_from_DrugBank_201214.txt"
#> [15] "RA_drugs_from_DrugBank_200214.txt"
```

### Make drug target matrix

Format files provided in data-raw for use in analysis.

``` r
DrugBank_info <- read.table(file = "data-raw/all_drug_targets_drug_bank.txt",sep="\t", header = T, stringsAsFactors = F)
head(DrugBank_info)
#>    drugID drug_name   status Swiss_target_ID
#> 1 DB00001 Lepirudin approved          P00734
#> 2 DB00002 Cetuximab approved          P00533
#> 3 DB00002 Cetuximab approved          O75015
#> 4 DB00002 Cetuximab approved          P00736
#> 5 DB00002 Cetuximab approved          P02745
#> 6 DB00002 Cetuximab approved          P02746
#>                                                  target_name gene_symbol
#> 1                                                Prothrombin          F2
#> 2                           Epidermal growth factor receptor        EGFR
#> 3 Low affinity immunoglobulin gamma Fc region receptor III-B      FCGR3B
#> 4                                Complement C1r subcomponent         C1R
#> 5                      Complement C1q subcomponent subunit A        C1QA
#> 6                      Complement C1q subcomponent subunit B        C1QB
#>   target_organism drug_action
#> 1       inhibitor      Humans
#> 2      antagonist      Humans
#> 3            <NA>      Humans
#> 4            <NA>      Humans
#> 5            <NA>      Humans
#> 6            <NA>      Humans
```

``` r
# Summarize drug targets for all drugs, translates them to Entrez Gene IDs and saves matrix in 'data/'
# Every column includes the drug targets of one drug.
data_analysis_and_formatting()
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
drug_target_matrix <- as.matrix(read.table(file = "data/DrugBank_drug_target_matrix.txt",sep="\t", header = T, stringsAsFactors = F))
drug_target_matrix[1:5,1:6]
#>      DB00001 DB00002 DB00004 DB00005 DB00006 DB00007
#> [1,]    2147    1956    3559    7124    2147    2798
#> [2,]      NA    2215    3560    7133      NA      NA
#> [3,]      NA     715    3561    2209      NA      NA
#> [4,]      NA     712      NA    2214      NA      NA
#> [5,]      NA     713      NA    2212      NA      NA
```

### Selection of the largest connected component (LCC) of the protein-protein interaction network (PPIN)

Exclusion of nodes that are not connected to the largest connected
component (LCC) of the protein-protein interaction network (PPIN) is
crucial for network distance calculations so that path lengths have
finite values. In order to use only the LCC of the chosen PPIN in
further calculations we apply the following function:

``` r
# Full PPIN (literature curated) - genes annotated as Entrez Gene IDs
lit_ppi <- as.matrix(read.table(file = "data-raw/lit_ppi.txt",sep="\t", header = T, stringsAsFactors = F))
head(lit_ppi)
#>      Protein_A Protein_B
#> [1,]         1       310
#> [2,]         1       368
#> [3,]         1      1026
#> [4,]         1      2232
#> [5,]         1      2886
#> [6,]         1      3172
dim(lit_ppi)
#> [1] 351444      2

# n unique proteins in PPIN
length(unique(as.vector(lit_ppi)))
#> [1] 17706

# Select LCC
ppin <- ppin_formatting(as.matrix(lit_ppi))
#> [1] "n of unique proteins/genes in PPIN: n =  17651"
dim(ppin)
#> [1] 346324      2

# n unique proteins in LCC of PPIN
length(unique(as.vector(ppin)))
#> [1] 17651
```

As can be seen above the LCC of the literature curated PPIN includes
17651 unique proteins, compared to the 17706 proteins in the full
literature PPIN. Based on the difference in dimension we can also tell
that exclusion of the 55 proteins that were not connected to the LCC
resulted reduced the number of protein-protein interaction by 351444 -
346324 = 5120.

### Select unique drug target combinations found in LCC of PPIN

``` r
# number of drug targets in drug_target_matrix
sum(!is.na(drug_target_matrix))
#> [1] 7924

# set drug targets that are not found in PPIN to NA
drug_target_matrix[!(drug_target_matrix %in% unique(as.vector(ppin)))] <- NA

# number of drug targets in the drug_target_matrix that were also found in PPIN
sum(!is.na(drug_target_matrix))
#> [1] 7801

# are there drugs with no drug targets in PPIN?
sum(colSums(!is.na(drug_target_matrix)) == 0)
#> [1] 4

# are there drugs with identical drug targets?
sum(duplicated(t(drug_target_matrix)))
#> [1] 569
```

As seen above, the raw version of `drug_target_matrix` is not in optimal
shape for network distance calculation yet. It still includes drug
targets that are not found in the PPIN as well as 569 duplicated drug
target entries. To speed up computations and decrease file size we will
apply the
function`prepare_drug_target_matrix_for_network_distance_calculation()`
that will exclude all drug targets not found in the PPIN as well as
excludes duplicated sets of drug targets. Only unique drug target
combinations will remain. As we after network distance calculations plan
to map results for unique drug target combinations back to individual
drugs, `prepare_drug_target_matrix_for_network_distance_calculation()`
prepares a file titled
`paste("SAME_DRUG_TARGETS_", file_name, ".txt", sep="")` that notes
which drugs had identical targets in a given PPIN. This is especially
important given that two drugs might have the same drug targets and
hence the same network distances in the PPIN, but yet different
pharmacological actions on their target.

``` r
#Preparation of drug target matrix for network distance calculation.
drug_target_matrix <- as.matrix(read.table(file = "data/DrugBank_drug_target_matrix.txt",sep="\t", header = T, stringsAsFactors = F))
drug_target_matrix <- prepare_drug_target_matrix_for_network_distance_calculation(ppin, drug_target_matrix, file_name = "in_lit_ppin", out_dir = "data")
#> [1] "n of unique proteins/genes in PPIN: n =  17651"
#> [1] "n unique drugs included in analysis: n = 1840"
#> [1] "n drugs with unique drug target combinations: n = 1204"
#> [1] "n unique drugs that can be extracted from UNIQUE_DRUG_TARGET_COMBINATION_in_lit_ppin.txt using SAME_DRUG_TARGETS_in_lit_ppin.txt: n = 1840"
#> [1] "Drug matrix = FORMATTED"

# File noting which drugs had identical drug tragets in the PPIN
# Left column represents the DrugBank ID of a drug that will represent a unique set of drug targets in the network distance calculation
# Right column represents the DrugBank IDs of all individual drugs that have this drug target combination 
same_drugs <- read.table(file = "data/SAME_DRUG_TARGETS_in_lit_ppin.txt", sep="\t", header = T, stringsAsFactors = F)
head(same_drugs)
#>     Drug1   Drug2
#> 1 DB00001 DB00001
#> 2 DB00001 DB00006
#> 3 DB00001 DB00278
#> 4 DB00001 DB01123
#> 5 DB00001 DB04898
#> 6 DB00001 DB06695
```

## Deep count autoencoder denoising of expression data

Having prepared drug data and protein interaction networks for analysis,
we need to format the scRNA-seq data set for calculation. Initially,
quality criteria is applied to the raw single cell data aiming to
exclude cells of poor quality and genes that are only expressed in very
few cells. After this, we applied Deep Count Autoencoder (DCA)
denoising, as described by [Gökcen et
al.](https://www.nature.com/articles/s41467-018-07931-2), following the
recommendations at [theislab/dca](https://github.com/theislab/dca).
Briefly that meant:

    # Installation
    $ conda install -c bioconda dca

    # Application to single cell matrix (represented by matrix.csv)
    $ dca matrix.csv results

DCA creates several output files documenting the denoising process,
namely: dispersion.tsv, dropout.tsv, latent.tsv and mean.tsv

mean.tsv corresponds to the DCA denoised gene expression matrix. DCA
also outputs a representation of the original single-cell data in the
latent space. This representation has fewer features than the original
data and represents intercellular expression differences generally
better than purely linear PCA models.The latent space representation is
also corrected for single-cell data artefacts such as dropouts and
varying library sizes.

## Formatting of single cell RNA-sequencing data

### Seurat for clustering

Using Seurat 3.1.0, we cluster the data based on DCA-derived latent
features for optimal cell typing. A detailed vignette on clustering
using Seurat can be found [here](https://github.com/satijalab/seurat).
In the code below `mean` and `latent` represent the DCA derived .tsv
files. In the single cell RNA-sequencing based gene expression matrix
each row represents a gene and each column an individual cell. Observe
that some functions might have been altered slightly in case you are
using Seurat 4.1.0.

``` r
#install.packages("Seurat")
library(Seurat)
set.seed(12)

# Downloads DCA adjusted data from zenodo.com if needed.
if(!file.exists("data-raw/AIA/mean.tsv")){
    if (!requireNamespace("RCurl", quietly = TRUE)) {
    stop(
      "Package \'RCurl\' must be installed to download DCA denoised data",
      call. = FALSE
    )
  }
  library(RCurl)
  #download.file("URL",destfile="data-raw/AIA/mean.tsv",method="libcurl") # INSERT REAL LINK!!!!!!!!!!!!!!!!!!!!!!!
  #download.file("URL",destfile="data-raw/AIA/latent.tsv",method="libcurl") # INSERT REAL LINK!!!!!!!!!!!!!!!!!!!!!!!
}

# Load data
denoised_DCA <- read.table("data-raw/DCA_adjusted_matrix/mean.tsv", sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(denoised_DCA) <- denoised_DCA[,1]
denoised_DCA <- denoised_DCA[,-1]
colnames(denoised_DCA) <- denoised_DCA[1,]
denoised_DCA <- denoised_DCA[-1,]
dim(denoised_DCA)
denoised_DCA[1:5,1:3]

lat <- read.table("data-raw/DCA_adjusted_matrix/latent.tsv", sep = "\t", header = F, stringsAsFactors = FALSE)
rownames(lat) <- lat[,1]
lat <- lat[,-1]
colnames(lat) <- paste("PC_",1:ncol(lat),sep="")
lat <- t(lat)
dim(lat)
lat[1:5,1:5]

# Create Seurat object using denoised expression values
denoised_DCA_clustered <- CreateSeuratObject(counts = denoised_DCA, project = "DCA_denoised", min.cells = 1, min.features = 1)


# Dimensional reduction using DCA derived latent features instead of linear principal components
pca <- new("DimReduc", cell.embeddings = t(lat), assay.used = "RNA", key = "PC_")
denoised_DCA_clustered@reductions$pca <- pca

# Clustering
denoised_DCA_clustered <- FindNeighbors(denoised_DCA_clustered, dims = 1:nrow(lat), k = 20) # shared nearest neighbor graph + Jaccard index
denoised_DCA_clustered <- FindClusters(denoised_DCA_clustered, resolution = 0.6) # Louvain algorithm
```

    #> Loading required package: Matrix
    #>               Joint_Healthy_mouse_1_CAGATGTTGGAT
    #> 0610007P14Rik                           1.498564
    #> 0610009B22Rik                           1.337379
    #> 0610009E02Rik                           0.532906
    #> 0610009L18Rik                           0.726244
    #> 0610009O20Rik                           1.186423
    #>               Joint_Healthy_mouse_1_TTTCCTCATCCN
    #> 0610007P14Rik                           0.960031
    #> 0610009B22Rik                           0.855918
    #> 0610009E02Rik                           0.649258
    #> 0610009L18Rik                           0.565916
    #> 0610009O20Rik                           0.755993
    #>               Joint_Healthy_mouse_1_GGATTTCTAGTT
    #> 0610007P14Rik                           1.058407
    #> 0610009B22Rik                           0.631603
    #> 0610009E02Rik                           0.268703
    #> 0610009L18Rik                           0.333690
    #> 0610009O20Rik                           0.537956
    #> [1]    32 16751
    #>      Joint_Healthy_mouse_1_CAGATGTTGGAT Joint_Healthy_mouse_1_TTTCCTCATCCN
    #> PC_1                          -0.951603                          -1.043974
    #> PC_2                           4.576371                          10.440075
    #> PC_3                          -1.135241                         -15.394433
    #> PC_4                          -6.364090                          -5.961215
    #> PC_5                          14.607316                          22.068670
    #>      Joint_Healthy_mouse_1_GGATTTCTAGTT Joint_Healthy_mouse_1_TAACAGCTATCG
    #> PC_1                         -35.096344                          -4.729716
    #> PC_2                         -11.036093                           8.775923
    #> PC_3                           0.572823                         -15.965678
    #> PC_4                         -22.459824                          -3.946763
    #> PC_5                          11.196144                          21.643744
    #>      Joint_Healthy_mouse_1_CACCGGAGTTCT
    #> PC_1                           1.838377
    #> PC_2                          12.284655
    #> PC_3                         -15.080997
    #> PC_4                          -8.692765
    #> PC_5                          30.427040

Above code results in clustering of the scRNA-seq data. In
`CreateSeuratObject()` we can set min.cells = 1 and min.features = 1, as
we previously have filtered the raw scRNA-seq expression matrix and
applied quality criteria. Application of the `FindNeighbors()`and
`FindClusters()`function with standard parameters will result in Seurat
constructing a shared nearest-neighbour graph followed by application of
the Louvain algorithm to the shared nearest neighbour graph to identify
clusters. In `FindNeighbors()`, specifying `k = 20`leads to Seurat using
the 20 nearest neighbours as defined by Jaccard index. In
`FindClusters()` the resolution parameter will affect the number of
identified clusters. The value of the resolution parameter will depend
on the size of the data set (number of cells). Arbitrary optimization of
`k` and `resolution` will be required for new data sets in order to
identify biologically valuable clusters, though there have been recent
suggestions by
[satijalab](https://satijalab.org/seurat/articles/get_started.html) on a
more systematic optimization process in which data is initially
overclustered and nodes og the cluster tree are merged until the ideal
threshold is reached.

Having now clustered the data we will need to save cluster outcomes.

    #> Attaching SeuratObject

``` r
# Save the clustering outcome
id <- Idents(denoised_DCA_clustered)
id <- cbind(as.character(names(id)), as.character(id))
colnames(id) <- c("Cell_ID", "Cluster_ID")
write.table(id, file = "Cell_identity.txt", sep="\t", col.names = T, row.names = F)

# n cells per cluster
table(id[,2])
#> 
#>    0    1   10   11   12   13   14   15   16   17    2    3    4    5    6    7 
#> 1874 1477  838  746  724  696  624  613  305  274 1365 1214 1140 1101  981  934 
#>    8    9 
#>  932  913
```

### MAST framework for calculation of differentially expressed genes

DEG calculation was performed using the MAST framework which was called
on through Seurat. The MAST framework deploys a scRNA-seq tailored
hurdle model for evaluation of expression changes between groups. Here
we chose to calculate DEGs between cells from sick and healthy samples
within the same cluster. Before we start DEG calculations we will need
to define groups of cells between we wish to calculate DEGs. While we
believe the DEG calculation between sick and healthy cells within the
same cluster to be crucial for scPred, there are many ways in which a
user can make Seurat understand which group of cells it should use. Our
approach below is merely a suggestion.

``` r
id[grepl("Healthy", id[,1]),2] <- as.numeric(id[grepl("Healthy", id[,1]),2]) + 0.1
table(id[,2])
#> 
#>    0  0.1    1  1.1   10 10.1   11 11.1   12 12.1   13 13.1   14 14.1   15 15.1 
#> 1308  566 1274  203  694  144  667   79  595  129  606   90  367  257  449  164 
#>   16 16.1   17 17.1    2  2.1    3  3.1    4  4.1    5  5.1    6  6.1    7  7.1 
#>  273   32  200   74  870  495 1051  163  208  932  957  144  874  107  822  112 
#>    8  8.1    9  9.1 
#>  715  217  270  643
# --> we have > 3 healthy and sick cells in each cluster enabling DEG calculation for all clusters

# Set new identities in Seurat object
temp_id <- as.factor(as.numeric(id[,2]))
denoised_DCA_clustered@active.ident <- temp_id

# unique clusters
temp_id_unique <- sort(as.numeric(unique(id[,2])))
temp_id_unique
#>  [1]  0.0  0.1  1.0  1.1  2.0  2.1  3.0  3.1  4.0  4.1  5.0  5.1  6.0  6.1  7.0
#> [16]  7.1  8.0  8.1  9.0  9.1 10.0 10.1 11.0 11.1 12.0 12.1 13.0 13.1 14.0 14.1
#> [31] 15.0 15.1 16.0 16.1 17.0 17.1
```

Now that we have defined groups for DEG calculation we calculate DEGs
between healthy and sick cells in each cluster individually.

``` r
# Running MAST DEG calculation on DCA denoised data
calculate_DEGs <- function(temp_id_unique, outdir, seurat_obj, pseudo.count = 0, logfc = log(1.5)){
  # Running MAST DEGs analysis
  for(i in 1:(length(temp_id_unique)/2)){
    pos <- i*2-1
    temp <- FindMarkers(seurat_obj,
                        slot = "counts", # DCA denoised, log10 transformed counts
                        test.use = "MAST",
                        ident.1 = as.character(temp_id_unique[pos]), # important that this is a factor vector with cell labels as names
                        ident.2 = as.character(temp_id_unique[pos+1]),
                        only.pos = F,
                        random.seed = 3,
                        pseudocount.use = pseudo.count, # not needed as DCA is already log10 transformed 
                        logfc.threshold = logfc,
                        min.pct = 0.1,
                        min.cells.group = 3)
      if(exists("temp")){
        write.table(temp, file = paste(outdir, "/Cluster_", temp_id_unique[pos], "_res=0.6_dims=32_k=20.txt",sep=""),sep="\t",col.names=NA,row.names=T)
        rm(temp)
      }
  }
}
calculate_DEGs(temp_id_unique, outdir = "data/AIA", seurat_obj = denoised_DCA_clustered)

# Example output (plasma B cell DEGs):
head(read.table(file = "data/AIA/Cluster_14_res=0.6_dims=32_k=20.txt", sep="\t", header = T, stringsAsFactors = F))
```

    #>          X        p_val  avg_logFC pct.1 pct.2    p_val_adj
    #> 1  Mir6236 4.786502e-51 -0.7921202     1     1 1.162928e-46
    #> 2 AF357399 2.686115e-34 -0.4910407     1     1 6.526185e-30
    #> 3     Rtp4 7.752321e-31  0.8546260     1     1 1.883504e-26
    #> 4     Irf7 3.479412e-29  0.8098554     1     1 8.453580e-25
    #> 5   Ifitm6 4.021597e-29  0.9002376     1     1 9.770873e-25
    #> 6     Lrg1 3.646098e-27  0.8390006     1     1 8.858559e-23

Next we will make a summary of all clusters DEGs for easy access by
future functions.

### Cell typing using marker genes

## Network distance calculation

### Calculation of average closest network distance compared

### Summary of network distance calculations

### Evaluating whether drug candidates counteracted fold change of targeted DEGs

### Final drug candidate ranking
