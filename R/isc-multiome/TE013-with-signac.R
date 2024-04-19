# https://www.youtube.com/watch?app=desktop&v=e2396GKFMRY&feature=youtu.be
# ------------------------------------------------------------------------
# system.time({source("sources/isc-multiome/TE013-with-signac.R")})

source("../vaxtools/R/utils.R")
my_sample_id <- "TE013"
create_workspace( isWorkspaceToClean = FALSE , 
									experiments_dir = "experiments" , 
									run_dir = paste0(my_sample_id,"-signac"))

library(logger)
# Using logger library
log_threshold(DEBUG)
log_threshold(DEBUG,index=2)

log_filename <- file.path(reports.dir,"log.txt")
log_appender(appender = appender_file(log_filename),index=2)
# log_formatter(formatter_glue)
log_layout(layout_glue_colors)
log_warn('--- Start of Analysis ----')
# demo(colors, package = 'logger', echo = FALSE)	

log_info("Working on Sample: " , my_sample_id)

library(tidyverse)
library(hdf5r)
library(Signac)
library(Seurat)

library(future)
options(future.globals.maxSize = 100 * 1024 ^ 3) # for 50 Gb RAM
plan()
library(parallel)
plan("multicore", workers = 8)

log_info(">>> Loading Peaks")
x_counts <- Read10X_h5("~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/peaks/filtered_feature_bc_matrix.h5")
peaks_counts <- x_counts$Peaks
dim(peaks_counts)
head(peaks_counts)

log_info(">>> Loading Metadata")
my_meta <- read.csv2("~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/per_barcode_metrics.csv",sep = ",")
rownames(my_meta) <- my_meta$barcode

log_info(">>> Loading Fragments")
fragments.file <- "~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/fragments/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay( counts = peaks_counts ,
																		 sep = c(":","-") ,
																		 genome = "mm10" ,
																		 fragments = fragments.file ,
																		 min.cells = 10 , 
																		 min.features = 200 )
log_info(">>> Creating Seurat Object")
seurat_obj <- CreateSeuratObject( counts = chrom_assay ,
																	assay = "peaks" , 
																	meta.data = my_meta )
# seurat_obj[[]] %>% head()


library(EnsDb.Mmusculus.v79) # BiocManager::install("EnsDb.Mmusculus.v79")
library(GenomeInfoDb) # BiocManager::install("GenomeInfoDb")
library(biovizBase) # BiocManager::install("biovizBase")

annotations <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"

Annotation(seurat_obj) <- annotations
log_info(">>> Computing Nucleosome Signal")
seurat_obj <- NucleosomeSignal(seurat_obj)
log_info(">>> Computing TSS Enrichment")
seurat_obj <- TSSEnrichment(seurat_obj)

VlnPlot( seurat_obj , features = c("atac_peak_region_fragments" , "nucleosome_signal" , "TSS.enrichment") , pt.size = 0.1 , ncol = 5 )

print(seurat_obj)
log_info(">>> >> # Cells before filtering: " , ncol(seurat_obj))
seurat_obj_subset <- seurat_obj %>% subset( subset = atac_peak_region_fragments > 5000 & TSS.enrichment > 4 )
print(seurat_obj_subset)
log_info(">>> >> # Cells after filtering: " , ncol(seurat_obj_subset))

log_info(">>> Running RunTFIDF")
seurat_obj_subset <- RunTFIDF(seurat_obj_subset)
log_info(">>> Running FindTopFeatures")
seurat_obj_subset <- FindTopFeatures(seurat_obj_subset,min.cutoff = 'q0')
log_info(">>> Running RunSVD")
seurat_obj_subset <- RunSVD(seurat_obj_subset)

DepthCor(seurat_obj_subset) # The first component correlates with depth, so we remove it

log_info(">>> Running RunUMAP")
seurat_obj_subset <- RunUMAP(seurat_obj_subset,reduction = 'lsi',dims=2:30)
log_info(">>> Running FindNeighbors")
seurat_obj_subset <- FindNeighbors(seurat_obj_subset,reduction = 'lsi',dims=2:30)
log_info(">>> Running FindClusters")
seurat_obj_subset <- FindClusters(seurat_obj_subset,algorithm = 3)	
DimPlot(seurat_obj_subset,label = TRUE) + NoLegend()

log_info(">>> Running GeneActivity")
gene_activities <- GeneActivity(seurat_obj_subset)

seurat_obj_subset[['gene_activities']] <- CreateAssayObject( gene_activities )

log_info(">>> Running NormalizeData")
seurat_obj_subset <- NormalizeData( seurat_obj_subset ,
																		assay = 'gene_activities' ,
																		normalization.method = 'LogNormalize' ,
																		scale.factor = median(seurat_obj_subset$nCount_gene_activities ) )

DefaultAssay(seurat_obj_subset) <- 'gene_activities'
FeaturePlot( seurat_obj_subset, features = c('Mki67' , 'Lgr5' , 'Dclk1' , 'Dll1' , 'Lyz1')  , max.cutoff = 'q95' )
DefaultAssay(seurat_obj_subset) <- 'peaks'	

DefaultAssay(seurat_obj_subset) <- 'gene_activities'
log_info(">>> Running FindMarkers on Gene Activities")
da_genes <- FindMarkers(seurat_obj_subset , ident.1 = '8' , min.pct = 0.05 , test.use = 'LR')

DefaultAssay(seurat_obj_subset) <- 'peaks'
log_info(">>> Running FindMarkers on Peaks")
da_peaks <- FindMarkers(seurat_obj_subset , ident.1 = '8' , min.pct = 0.05 , test.use = 'LR')

log_info(">>> Saving Signac Data Objects")
saveRDS(seurat_obj,"~/Clouds/Dropbox/Data/isc/TE013/TE013-signac.rds")
saveRDS(seurat_obj_subset,"~/Clouds/Dropbox/Data/isc/TE013/TE013-signac-subset.rds")

da_peaks$distance <- ClosestFeature(seurat_obj_subset,regions = rownames(da_peaks))$distance
da_peaks$closest_gene <- ClosestFeature(seurat_obj_subset,regions = rownames(da_peaks))$closest_gene

CoveragePlot( seurat_obj_subset , 
							region = rownames(da_peaks)[1] ,
							extend.upstream =  10000 , extend.downstream = 5000 , group.by = 'seurat_clusters' )


cyto_te013 <- read_delim("~/Clouds/Dropbox/Data/isc/TE013/TE013-cytotrace.csv",delim = ",")

peaks_counts.mat <- as.matrix(peaks_counts)
index <- rowSums(peaks_counts.mat) > 0
sum(index)
my_peaks <- colSums(peaks_counts.mat > 0)

peaks_tbl <- tibble( cell_id = names(my_peaks) ,
										 open_peaks = my_peaks )

my_tibble <- left_join( cyto_te013 , peaks_tbl )

ggplot(my_tibble, aes(x=cytotrace_score,y=open_peaks)) +
	geom_point(size=0.05) +
	xlim(c(1,0)) +
	theme_minimal()

cor(my_tibble$cytotrace_score,my_tibble$open_peaks,method = "pea")

# save(list = ls(),file="~/Desktop/TE013-R-Session-for-rebuttal")
# save.image()
