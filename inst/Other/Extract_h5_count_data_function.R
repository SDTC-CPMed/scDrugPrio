#
#
# Read cellranger hdf5 .h5 files in R
#
###############################################################

library(Matrix)
library(rhdf5)

extract_h5_count_data <- function(my_h5_file, sample_name) {
  # Each cell has a barcode nucleotide sequence.
  my_barcodes <- my_h5_file$barcodes
  # Merge barcode with sample name
  my_barcodes <- paste(sample_name, "_", my_barcodes, sep="")
  # Read the data into a sparse matrix.
  counts <- sparseMatrix(
    dims = my_h5_file$shape,
    i = as.numeric(my_h5_file$indices),
    p = as.numeric(my_h5_file$indptr),
    x = as.numeric(my_h5_file$data),
    index1 = FALSE
  )
  #counts <- as.matrix(counts)
  colnames(counts) <- my_barcodes
  rownames(counts) <- my_h5_file$genes
  return(counts)
}


