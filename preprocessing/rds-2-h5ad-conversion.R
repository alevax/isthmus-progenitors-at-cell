# Luca Zanella
# 12.21.2023

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
  rdsFolder <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001/"
  h5adFolder <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001-h5ad" 
}

# Libraries
{
  suppressWarnings(suppressMessages(library(Seurat)))
  suppressWarnings(suppressMessages(library(SeuratData)))
  suppressWarnings(suppressMessages(library(SeuratDisk)))
}




# Reading in the .rds Seurat objects and converting them to .h5ad 
{
  
  cpm <- t(readRDS(file.path(rdsFolder, "TE001-cpm.rds"))) # cells x genes
  
  counts_Seurat <- readRDS(file.path(rdsFolder, "TE001-seurat-analysis-data.rds"))
  
  vp_Seurat <- readRDS(file.path(rdsFolder, "TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-with-paneth.rds"))

}

# Filter out Seurat object with more stringent thresholds
{
  counts_Seurat <- subset(counts_Seurat, subset = nFeature_RNA > 1500)
}



# Export files to .tsv and .h5ad for use with Python
{
  # cpm as .tsv
  write.table(x=cpm, file = file.path(h5adFolder, "TE001-cpm.tsv"), sep="\t")
  
  # counts to AnnData 
  SaveH5Seurat(counts_Seurat, filename=file.path(h5adFolder,"TE001-counts.h5Seurat"))
  Convert(source=file.path(h5adFolder,"TE001-counts.h5Seurat"), dest="h5ad", assay="RNA")
  
  
  # SCT assay to AnnData
  SaveH5Seurat(counts_Seurat, filename=file.path(h5adFolder,"TE001-sct.h5Seurat"))
  Convert(source=file.path(h5adFolder,"TE001-sct.h5Seurat"), dest="h5ad", assay="SCT")
  
  
  
  # viper to AnnData
  SaveH5Seurat(vp_Seurat, filename=file.path(h5adFolder,"TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-with-paneth.h5Seurat"))
  Convert(source=file.path(h5adFolder,"TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-with-paneth.h5Seurat"), dest="h5ad")

}