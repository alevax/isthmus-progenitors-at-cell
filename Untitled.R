# Luca Zanella
# 11.28.2023

# Converting TE001 samples to h5ad for usage with Python (for scanpy and CellRank2 analysis)
# INPUTS (MANDATORY): 
# rdsFolder: path to where .rds Seurat data are stored (scratch)
# h5adFolder: path to where .h5ad Seurat data will be saved (scratch)
#
rm(list = ls()) # cleans workspace
gc()
cat("\014")

# Set paths
{
  cat("YOU MUST SET 'rdsFolder' and 'h5adFolder'.\n")
  rdsFolder <- "/Volumes/ac_lab_scratch/lz2841/isc-rebuttal/TE001/"
  h5adFolder <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001-h5ad" 
}



