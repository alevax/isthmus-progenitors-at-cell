##
# Weighted Nearest Neighbor Analysis
# ----------------------------------
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# -----------------------------------------------------------------------------
# system.time({ source("sources/isc-multiome/TE013-wnn-analysis.R") })

source("../vaxtools/R/utils.R")
my_seed <- 42
set.seed(my_seed)

# sample_id <- "MJ019"
# counts_path <- paste0("~/Clouds/Dropbox/Data/mouse-organoids/new-data/",sample_id,"/filtered_feature_bc_matrix/")
# atac_path <- paste0("~/Clouds/Dropbox/Data/mouse-organoids/ATAC-Seq/",sample_id)
sample_id <- "TE013"
counts_path <- "~/Clouds/Dropbox/Data/isc/TE013/filtered_feature_bc_matrix/"
atac_path <- "~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/fragments/"

# gene_list <- c("Ar","Vim","Nr3c1","Ascl1","Ascl2","Nsd2",
# 							 "Foxm1","Atoh1","Zeb1","Zeb2","Snai1",
# 							 "Snai2","Twist1","Twist2","Sox2","Sox11")
gene_list <- c("Dclk1","Atoh1","Olfm4","Mki67","Lgr5",
"Ascl2","Sox2","Krt19")

create_workspace( paste0("wnn-analysis-",sample_id))

library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)

## Data Input ----

# load the RNA and ATAC data
my_10x_data <- Read10X(counts_path)
# extract RNA and ATAC data
rna_counts <- my_10x_data$`Gene Expression`
atac_counts <- my_10x_data$Peaks

# Create Seurat object
my_seurat.obj <- CreateSeuratObject(counts = rna_counts)
my_seurat.obj[["percent.mt"]] <- PercentageFeatureSet(my_seurat.obj, pattern = "^Mt-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- file.path(atac_path,"/atac_fragments.tsv.gz")
chrom_assay <- CreateChromatinAssay(
	counts = atac_counts,
	sep = c(":", "-"),
	genome = 'mm10',
	fragments = frag.file,
	min.cells = 10,
	annotation = annotations
)
my_seurat.obj[["ATAC"]] <- chrom_assay

VlnPlot(my_seurat.obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
				log = TRUE, pt.size = 0) + NoLegend()

my_seurat.obj <- subset(
	x = my_seurat.obj,
	subset = nCount_ATAC < 7e4 &
		nCount_ATAC > 5e3 &
		nCount_RNA < 25000 &
		nCount_RNA > 1000 &
		percent.mt < 20
)

## Filtering out CD45+ cells ----
index_CD45 <- grep('^Ptprc$', rownames(my_seurat.obj), value = FALSE)
index_Epcam <- grep('^Epcam$', rownames(my_seurat.obj), value = FALSE)
rownames(my_seurat.obj)[index_CD45]
index_CD45 <- my_seurat.obj@assays$RNA@data[index_CD45,] > 0
sum(index_CD45)
index_Epcam <- my_seurat.obj@assays$RNA@data[index_Epcam,] > 0
sum(index_Epcam)
sum( index_CD45 & index_Epcam )
ncol(my_seurat.obj)
my_seurat.obj <- my_seurat.obj %>% subset( cells = colnames(my_seurat.obj)[!index_CD45] )
ncol(my_seurat.obj)

# RNA analysis ----
DefaultAssay(my_seurat.obj) <- "RNA"
my_seurat.obj <- SCTransform(my_seurat.obj, verbose = FALSE) %>% 
	RunPCA() %>% 
	RunUMAP(dims = 1:50, 
					reduction.name = 'umap.rna', 
					reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(my_seurat.obj) <- "ATAC"
my_seurat.obj <- RunTFIDF(my_seurat.obj)
my_seurat.obj <- FindTopFeatures(my_seurat.obj, min.cutoff = 'q0')
my_seurat.obj <- RunSVD(my_seurat.obj)
my_seurat.obj <- RunUMAP(my_seurat.obj, 
												 reduction = 'lsi', dims = 2:50, 
												 reduction.name = "umap.atac", 
												 reduction.key = "atacUMAP_")


## calculate a WNN graph, ----
# representing a weighted combination of RNA and ATAC-seq modalitie
my_seurat.obj <- FindMultiModalNeighbors(my_seurat.obj, 
																				 reduction.list = list("pca", "lsi"), 
																				 dims.list = list(1:50, 2:50))
my_seurat.obj <- RunUMAP(my_seurat.obj, nn.name = "weighted.nn", 
												 reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
my_seurat.obj <- FindClusters(my_seurat.obj, graph.name = "wsnn", 
															resolution = 0.1 ,
															algorithm = 3, verbose = FALSE)
Idents(my_seurat.obj)

DefaultAssay(my_seurat.obj) <- "RNA"
cluster3.markers <- FindMarkers(my_seurat.obj, ident.1 = 3 )
head(cluster3.markers, n = 5)

# # perform sub-clustering on cluster 6 to find additional structure
# my_seurat.obj <- FindSubCluster(my_seurat.obj, cluster = 6, graph.name = "wsnn", algorithm = 3)
# Idents(my_seurat.obj) <- "sub.cluster"

print_msg_info(">>> Loading VIPER Clusters from ingest Analysis")
{
	my_seurat.obj@meta.data$cell_id <- colnames(my_seurat.obj)
	
	my_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/multiome-metadata-ingest.csv")
	vp_metadata <- read_csv(my_data_dir)

	table(vp_metadata$cluster_id)
	table(vp_metadata$iter_cluster_id)
	# vp_metadata$cell_id_new <- paste0("TE013#",vp_metadata$cell_id)
	cells_in_common <- intersect(vp_metadata$cell_id,my_seurat.obj$cell_id)

	# TE013_ATAC.project <- TE013_ATAC.project[ cells_in_common ]
	# 	vp_data$cell_id_new <- paste0("TE013#",vp_data$cell_id)

	my_seurat.obj <- my_seurat.obj %>% subset( cells = cells_in_common )
	
	index <- match( my_seurat.obj$cell_id , vp_metadata$cell_id )

	vp_metadata <- vp_metadata[index,]
	dim(vp_metadata)

	stopifnot(identical( as.character(my_seurat.obj$cell_id) , vp_metadata$cell_id ))

	# TE013_ATAC.project@cellColData$ClustersVIPER <- as.character(vp_metadata$cluster_id)
	my_seurat.obj@meta.data$ClustersVIPER <- as.character(vp_metadata$iter_cluster_id)

}

p1 <- DimPlot(my_seurat.obj, reduction = "umap.rna", 
							group.by = "seurat_clusters", label = TRUE, label.size = 4, 
							repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(my_seurat.obj, reduction = "umap.atac", 
							group.by = "seurat_clusters", label = TRUE, label.size = 4, 
							repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(my_seurat.obj, reduction = "wnn.umap", group.by = "seurat_clusters", 
							label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")

p4 <- DimPlot(my_seurat.obj, reduction = "umap.rna", group.by = "ClustersVIPER", 
							label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ClustersVIPER")
p4 <- p4 + scale_color_brewer(palette = "Set2")

p5 <- DimPlot(my_seurat.obj, reduction = "umap.atac", group.by = "ClustersVIPER", 
							label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ClustersVIPER")
p5 <- p5 + scale_color_brewer(palette = "Set2")

p6 <- DimPlot(my_seurat.obj, reduction = "wnn.umap", group.by = "ClustersVIPER", 
							label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ClustersVIPER")
p6 <- p6 + scale_color_brewer(palette = "Set2")

p <- p1 + p2 + p3 + p4 + p5 + p6 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
pdf( file.path(reports.dir,"umap.pdf") , width = 12 , height = 8 )
	print(p)
dev.off()

## to make the visualization easier, subset T cell clusters
# celltype.names <- levels(my_seurat.obj)
# tcell.names <- grep("CD4|CD8|Treg", celltype.names,value = TRUE)
# tcells <- subset(pbmc, idents = tcell.names)

## Color Palettes ----
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
# tab:blue : #1f77b4
# tab:orange : #ff7f0e
# tab:green : #2ca02c
# tab:red : #d62728
# tab:purple : #9467bd
# tab:brown : #8c564b
# tab:pink : #e377c2
# tab:gray : #7f7f7f
# tab:olive : #bcbd22
# tab:cyan : #17becf

Idents(my_seurat.obj) <- my_seurat.obj@meta.data$ClustersVIPER

for ( a_gene in gene_list )
{
	if ( !(a_gene %in% rownames(my_seurat.obj@assays$SCT@data)) )
	{
		next ;
	}
		
	p <- CoveragePlot(my_seurat.obj, 
										region = a_gene ,
										features = a_gene, 
										# ranges.group.by = "ClustersVIPER" ,
										extend.upstream = 50 , extend.downstream = 50 ,
										assay = 'ATAC', 
										expression.assay = 'SCT', peaks = TRUE)
	

	my_palette <- tab10_palette
	names(my_palette) <- names(table(my_seurat.obj@meta.data$ClustersVIPER))
	p <- p & scale_fill_manual(values = my_palette)
	
	filename <- file.path(reports.dir,paste0("coverage-plot-",a_gene,".pdf"))
	pdf(filename , width = 8 , height = 4.5)
		print(p)
	dev.off()	
	
}


## Run ChromVAR Analysis ----
# BiocManager::install(c("chromVAR","JASPAR2020","TFBSTools","motifmatchr"))
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)

# library(devtools)
# install_github('immunogenomics/presto')

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(my_seurat.obj) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, 
												# opts = list(species = 9606, all_versions = FALSE))
												opts = list(collection="CORE", 
																		species = 10090,
																		# species = 9606,
																		all_versions = FALSE))
# sapply(pwm_set,function(x) x@name )
motif.matrix <- CreateMotifMatrix(features = granges(my_seurat.obj), 
																	pwm = pwm_set, genome = 'mm10', 
																	use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
my_seurat.obj <- SetAssayData(my_seurat.obj, assay = 'ATAC', 
															slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
system.time({ my_seurat.obj <- RunChromVAR(
	object = my_seurat.obj,
	genome = BSgenome.Mmusculus.UCSC.mm10
)})

for ( a_gene in gene_list )
{
	motif.name <- ConvertMotifID(my_seurat.obj, name = a_gene)
	if(is.na(motif.name)) next ;
	if ((a_gene %in% rownames(my_seurat.obj@assays$SCT)) == FALSE) next ;
	
	gene_plot <- FeaturePlot(my_seurat.obj, features = paste0("sct_",a_gene), order = TRUE ,
													 reduction = 'wnn.umap')
	motif_plot <- FeaturePlot(my_seurat.obj, features = motif.name, 
														order = TRUE ,
														min.cutoff = 0, cols = c("lightgrey", "darkred"), 
														reduction = 'wnn.umap')
	p <- gene_plot | motif_plot
	filename <- file.path(reports.dir,paste0("motif-plot-",a_gene,".pdf"))
	pdf(filename,width = 8,height = 4)
		print(p)
	dev.off()	
}

print_msg_info(">>> Plotting RNA-Seq expression of some genes")
gene_list.subset <- gene_list[ gene_list %in% rownames(my_seurat.obj@assays$SCT) ]
p <- FeaturePlot(my_seurat.obj, 
								 features = paste0("sct_",gene_list.subset), 
								 col = c("white","red") ,
								 order = TRUE , reduction = 'wnn.umap' )
filename <- file.path(reports.dir,"features-plot.pdf")
pdf(filename,width = 16,height = 16)
	print(p)
dev.off()	

print_msg_info(">>> Plotting RNA-Seq expression of some genes")
gene_list.subset <- gene_list[ gene_list %in% rownames(my_seurat.obj@assays$SCT) ]
p <- FeaturePlot(my_seurat.obj, 
								 features = paste0("sct_",gene_list.subset), 
								 col = c("white","red") , keep.scale = "all" ,
								 order = TRUE , reduction = 'wnn.umap' )

filename <- file.path(reports.dir,"features-plot-cutoff.pdf")
pdf(filename,width = 16,height = 16)
print(p)
dev.off()	

require(scales)
p <- DotPlot(my_seurat.obj, 
						 assay = "RNA" ,
						 cluster.idents = FALSE , 
						 scale = FALSE ,
						 features = paste0("rna_",gene_list.subset) ) + RotatedAxis() +
	scale_colour_viridis_c(option = "D") +
	theme(axis.text.x = element_text(size = 8))
filename <- file.path(reports.dir,"dot-plot.pdf")
pdf(filename,width = 5,height = 4)
	print(p)
dev.off()	

p <- DotPlot(my_seurat.obj, 
						 assay = "SCT" ,
						 cluster.idents = FALSE , 
						 scale = TRUE ,
						 features = paste0("sct_",gene_list.subset) ) + RotatedAxis() +
scale_colour_gradient2(low = muted("blue") , mid = "white" , 
											 high = muted("red"), midpoint = 0) +
	theme(axis.text.x = element_text(size = 8))
filename <- file.path(reports.dir,"dot-plot-scaled.pdf")
pdf(filename,width = 5,height = 4)
	print(p)
dev.off()	

# 
# # library(presto)
# markers_rna <- presto:::wilcoxauc.Seurat(X = my_seurat.obj, 
# 																				 group_by = 'seurat_clusters', assay = 'data', seurat_assay = 'SCT')
# markers_motifs <- presto:::wilcoxauc.Seurat(X = my_seurat.obj, 
# 																						group_by = 'seurat_clusters', assay = 'data', seurat_assay = 'chromvar')
# motif.names <- markers_motifs$feature
# colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
# colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
# markers_rna$gene <- markers_rna$RNA.feature
# markers_motifs$gene <- ConvertMotifID(my_seurat.obj, id = motif.names)
# 
# # a simple function to implement the procedure above
# topTFs <- function(celltype, padj.cutoff = 1e-2) {
# 	ctmarkers_rna <- dplyr::filter(
# 		markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
# 		arrange(-RNA.auc)
# 	ctmarkers_motif <- dplyr::filter(
# 		markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
# 		arrange(-motif.auc)
# 	top_tfs <- inner_join(
# 		x = ctmarkers_rna[, c(2, 11, 6, 7)], 
# 		y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
# 	)
# 	top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
# 	top_tfs <- arrange(top_tfs, -avg_auc)
# 	return(top_tfs)
# }
# 
# head(topTFs("1"), 5)
# motif.name <- ConvertMotifID(my_seurat.obj, name = 'Rarg')
# gene_plot <- FeaturePlot(my_seurat.obj, features = "sct_Rarg", reduction = 'wnn.umap' , order = TRUE )
# motif_plot <- FeaturePlot(my_seurat.obj, features = motif.name, min.cutoff = 0, 
# 													order = TRUE , ,cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
# gene_plot | motif_plot
# 
# 
