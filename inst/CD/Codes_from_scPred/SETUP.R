# SETUP on R 3.4.4

# R script for package setup!
# Used to build up package library from empty 3.4 folder
###########################################################################################


# Quite a lot of loaded DLLs - ca 100
# Depending on local environment needed to increase R_MAX_NUM_DLLS=150
#
# Open terminal (bash) and run 'nano ~/.Renviron'
# Type 'R_MAX_NUM_DLLS=150'
# Close nano using crtl + X
# Save the modified .Renviron file
# DONE


if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)


if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

# INSTALL SEURAT 3.1.0 ON R 3.4.4
if (!requireNamespace("Seurat", quietly = TRUE)) {
  print("INSTALL SEURAT 3.1.0")
  devtools::install_version(package = 'plotrix', version = package_version('3.6'), dependencies = T)
  library(plotrix)
  devtools::install_version(package = 'metap', version = package_version('1.1'), dependencies = T)
  library(metap)
  install.packages("https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-8.tar.gz", repos=NULL, dependencies = T)
  library(mvtnorm)
  devtools::install_version(package = 'Seurat', version = package_version('3.1.0'), dependencies = T) # non-zero exit status

  # Install unavailable SEURAT dependencies (and their dependencies):
  devtools::install_version(package = 'cowplot', version = package_version('0.9.3'), dependencies = T)
  library(cowplot)
  devtools::install_version(package = 'fitdistrplus', version = package_version('1.0-14'), dependencies = T)
  library(fitdistrplus)
  devtools::install_version(package = 'irlba', version = package_version('2.3.2'), dependencies = T)
  library(irlba)
  devtools::install_version(package = 'caTools', version = package_version('1.17.1.3'), dependencies = T)
  # install ROCR
  library(caTools)
  devtools::install_version(package = 'gplots', version = package_version('3.0.1.1'), dependencies = T)
  library(gplots)
  devtools::install_version(package = 'ROCR', version = package_version('1.0.7'), dependencies = F)
  library(ROCR)

  devtools::install_version(package = 'rsvd', version = package_version('1.0.2'), dependencies = T)
  library(rsvd)
  devtools::install_version(package = 'sctransform', version = package_version('0.2.1'), dependencies = T)
  library(sctransform)
  devtools::install_version(package = 'SDMTools', version = package_version('1.1-221.1'), dependencies = T) # update none
  library(SDMTools)
  # install uwot
  devtools::install_version(package = 'RcppParallel', version = package_version('4.4.4'), dependencies = T)
  library(RcppParallel)
  devtools::install_version(package = 'uwot', version = package_version('0.1.4'), dependencies = F)
  devtools::install_github('jlmelville/uwot')
  library(uwot)
  devtools::install_version(package = 'Matrix', version = package_version('1.2.16'), dependencies = F)
  library(Matrix)
  install.packages("rlang") # update = TRUE
  library(rlang)
  devtools::install_version(package = 'Seurat', version = package_version('3.1.0'), dependencies = F)
}
library(Seurat)
print("Seurat working!")


if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}
library(doParallel)


if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
library(igraph)


if (!requireNamespace("GEOquery", quietly = TRUE)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GEOquery")
}
library(GEOquery)


if (!requireNamespace("R.filesets", quietly = TRUE)) {
  devtools::install_version("R.filesets", version = "2.14.0")
}
library(R.filesets)


if (!requireNamespace("matrixStats", quietly = TRUE)) {
  devtools::install_version("matrixStats", version = "0.55.0")
}
library(matrixStats)


if (!requireNamespace("nichenetr", quietly = TRUE)) {
  devtools::install_version("latticeExtra", version = "0.6.28")
  devtools::install_version("Hmisc", version = "4.2.0")
  devtools::install_version("progressr", version = "0.4.0")
  devtools::install_version("lava", version = "1.6.5")
  devtools::install_version("ipred", version = "0.9.9") # no updates
  devtools::install_version("caret", version = "6.0.81") # no updates
  devtools::install_version("misc3d", version = "0.8.4")
  devtools::install_version("plot3D", version = "1.1") # no updates
  devtools::install_version("smoof", version = "1.5.1") # no updates
  devtools::install_version("mlrMBO", version = "1.1.2") # no updates
  devtools::install_version("ggpubr", version = "0.2")
  devtools::install_version("rjson", version = "0.2.20")
  devtools::install_version("ggpubr", version = "0.1.7")
  source("https://bioconductor.org/biocLite.R")
  biocLite("ComplexHeatmap") # no updates
  #devtools::install_version("car", version = "3.0.3")
  devtools::install_github("saeyslab/nichenetr") # no updates
}
library(nichenetr)

if (!requireNamespace("CINNA", quietly = TRUE)) {
  devtools::install_version("network", version = "1.14.377") # no updates
  devtools::install_version("nloptr", version = "1.2.1")
  devtools::install_github("lme4/lme4")
  devtools::install_version("pbkrtest", version = "0.4.7") # no updates
  devtools::install_version("carData", version = "3.0.3") # no updates
  devtools::install_version("quantreg", version = "5.54") # no updates
  devtools::install_version("car", version = "3.0.2") # no updates
  devtools::install_version("rstatix", version = "0.3.1") # no updates
  devtools::install_version("intergraph", version = "2.0.1") # no updates
  devtools::install_version("FactoMineR", version = "1.42") # no updates
  devtools::install_version("factoextra", version = "1.0.6") # no updates
  devtools::install_version("statnet.common", version = "4.1.4")
  devtools::install_version("sna", version = "2.4") # no updates
  devtools::install_version("CINNA", version = "1.1.51") # no updates
}
library(CINNA)

'
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("KEGGREST", version = "3.6")
}
library(KEGGREST)

if (!requireNamespace("KEGGgraph", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("KEGGgraph", version = "3.6")
}
library(KEGGgraph)'







print("INSTALL stats")
install.packages("stats")
library(stats)
print("dplyr working!")

print("INSTALL limma")
install.packages("limma")
library(limma)
print("limma working!")

print("INSTALL gespeR")
#install.packages("gespeR")
print("gespeR working!")
BiocManager::install("gespeR", version = "3.6")
library(gespeR)

print("INSTALL jaccard")
install.packages("jaccard")
library(jaccard)
print("jaccard working!")

print("INSTALL corrplot")
install.packages("corrplot")
library(corrplot)
print("corrplot working!")

print("INSTALL CINNA")
install.packages("CINNA")
library(CINNA)
print("CINNA working!")

print("INSTALL igraph")
install.packages("igraph")
library(igraph)
print("igraph working!")

print("INSTALL hrbrthemes")
install.packages("hrbrthemes")
library(hrbrthemes)
print("hrbrthemes working!")

print("INSTALL CINNA")
devtools::install_version("nloptr", version = "1.2.2.1", repos = "http://cran.us.r-project.org", dependencies = F)
#devtools::install_version("lme4", version = "1.1-21", repos = "http://cran.us.r-project.org", dependencies = F)
#devtools::install_version("carData", version = "3.0-2", repos = "http://cran.us.r-project.org", dependencies = F)
#devtools::install_version("quantreg", version = "3.0-2", repos = "http://cran.us.r-project.org", dependencies = F)
devtools::install_version("network", version = "1.16.0", repos = "http://cran.us.r-project.org", dependencies = T)
devtools::install_version("statnet.common", version = "4.1.4", repos = "http://cran.us.r-project.org", dependencies = T)
library(nloptr)
library(network)
library(statnet.common)





devtools::install_version("CINNA", version = "1.1.53", dependencies = T)
library(CINNA)
print("CINNA working!")

print("Unloading all installed packages")
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
rm(list=ls())
