##
# Analysis with ArchR
# -------------------
# system.time({source("sources/isc-multiome/TE013-archR-analysis.R")})

RHDF5_USE_FILE_LOCKING=FALSE

setwd("~/Workspace/isc-project/")
source("../vaxtools/R/utils.R")

print_msg_info(">>> Install or Load Libraries")
{
	# remotes::install_version("RSQLite", version = "2.2.5") To fix a Bug in 2.2.6 version and org.Mm.eg.db
	# remotes::install_github('daroczig/logger')
	# https://daroczig.github.io/logger/articles/Intro.html
	
	# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
	# devtools::install_github("GreenleafLab/ArchR",
													 # ref="dev_groupCoverage", repos = BiocManager::repositories())
	
	library(ArchR)
	# ArchR::installExtraPackages()
	library(tidyverse)
	library(ggplot2)
	library(ggrepel)
	library(parallel) # mclapply
}

## Setting Parameters ----
print_msg_info(">>> Setting Parameters")
{
	sample_id <- "TE013"
	
	my_seed <- 42
	N_threads <- 1
	set.seed(my_seed)
	addArchRThreads(threads = N_threads,force = TRUE)
	addArchRGenome("mm10") # Working with Mouse Data
	
	my_dir <- paste0(sample_id,"-with-ArchR")
	create_workspace( my_dir )
	# reports.dir <- "/Users/av2729/Workspace/jia-experiment/experiments/2021-03-26-TE013-with-ArchR"	
	
	print_msg_info(">>> Setting Logger")
	{
		library(log4r)
		
		my_logfile = file.path(reports.dir,paste0(sample_id,"-atac-log.txt") )
		my_console_appender = console_appender(layout = default_log_layout())
		my_file_appender = file_appender(my_logfile, append = TRUE, 
																		 layout = default_log_layout())
		
		my_logger <- log4r::logger(threshold = "INFO", 
															 appenders= list(my_console_appender,my_file_appender))
		
		log4r_info <- function(msg) {
			log4r::info(my_logger, paste0("INFO" , " - " , msg ) )
		}
		
		log4r_error <- function() {
			log4r::error(my_logger, "ERROR")
		}
		
		log4r_debug <- function() {
			log4r::debug(my_logger, "DEBUG")
		}		
		
		log4r_info("Starting Analysis") 
		
	}
}


# ## Creating ArrowFiles ----
# print_msg_info("Creating ArrowFiles ...")
# {
# 	addArchRGenome("mm10")
# 	base_dir <- "~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq"
# 	# base_dir <- "/Volumes/Data Transfer/data-dumps/TE013-multiome-data"
# 	output_dir <- file.path(base_dir,"ArchR-analysis")
# 	# output_dir <- reports.dir
# 
# 	dir.create(output_dir)
# 	print_msg_warn("*** Changing Working Directory to: " , output_dir )
# 	setwd(output_dir)
# 
# 	filtered_fragments.filename <- file.path(base_dir,"fragments/atac_fragments.tsv.gz")
# 	names(filtered_fragments.filename) <- sample_id
# 	ArrowFiles <- createArrowFiles(
# 		inputFiles = filtered_fragments.filename,
# 		QCDir = file.path(output_dir,"quality-control-analysis"),
# 		logFile = file.path(output_dir,"logfile.txt") ,
# 		outputNames = names(filtered_fragments.filename) ,
# 		sampleNames = names(filtered_fragments.filename) ,
# 		# validBarcodes = rna_meta$sample_id ,
# 		minTSS = 4, #Dont set this too high because you can always increase later
# 		minFrags = 500,
# 		excludeChr = c("chrM") ,
# 		addTileMat = TRUE,
# 		addGeneScoreMat = TRUE ,
# 		verbose = TRUE ,
# 		subThreading = FALSE
# 	)
# 	# ArrowFiles <- file.path(output_dir,"TE013.arrow")
# }



## Loading ArrowFile  ----
# ArrowFiles <- "/Volumes/Data Transfer/data-dumps/TE013-multiome-data/TE013.arrow"
# ArrowFiles <- "/Users/av2729/Workspace/isc-project/experiments/2021-09-22-TE013-with-ArchR/reports/TE013.arrow"
# ArrowFiles <- "~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/arrows/TE013.arrow"
ArrowFiles <- "~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/ArchR-analysis/TE013.arrow"

TE013_ATAC.project <- ArchRProject( showLogo = FALSE ,
  ArrowFiles = ArrowFiles, 
  outputDirectory = reports.dir,
  # outputDirectory = output_dir,
  # outputDirectory = "~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/ArchR-analysis/"
  copyArrows = TRUE
)

object.size(TE013_ATAC.project) %>% format( units = "Mb")
getAvailableMatrices(TE013_ATAC.project)

## Inferring Doublets ----
doublet_score.obj <- addDoubletScores(
  input = ArrowFiles, 
  outDir = reports.dir , # This should be the Quality Control Directory
  k = 15, # Refers to how many cells near a "pseudo-doublet" to count.
  # knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search.
  knnMethod = "LSI",
  LSIMethod = 1 
)

print_msg_info(">>> Filtering Doublets")
{
	# log_warn("Doublets are detected but not yet removed !!!")
	TE013_ATAC.project <- ArchR::filterDoublets(TE013_ATAC.project)
}

print_msg_info(">>> Filtering Contaminant Cells")
{
	my_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-analysis-data.rds")
	gene_seurat <- readRDS(my_data_dir)
	
	CD45_cells <- colnames(gene_seurat)[ gene_seurat@assays$RNA@counts["Ptprc",] > 0 ]
	CD45_cells <- paste0("TE013#",CD45_cells)
	
	TE013_ATAC.project <- TE013_ATAC.project[ TE013_ATAC.project$cellNames[ !(TE013_ATAC.project$cellNames %in% CD45_cells) ] ]

	print_msg_info(">>> Removed " , length(CD45_cells), " cells")	

}

## addIterativeLSI ----
TE013_ATAC.project <- addIterativeLSI(
    ArchRProj = TE013_ATAC.project, 
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 5 , 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.5), 
        sampleCells = 10000, 
        n.start = 10 ,
        random.seed = my_seed
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30 ,
    force = TRUE ,
    seed = my_seed 
)

## addClusters ----
TE013_ATAC.project <- addClusters(
    input = TE013_ATAC.project,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "SeuratClusters",
    force = TRUE ,
    resolution = 0.75 ,
    seed = 666 
)
table(TE013_ATAC.project$SeuratClusters)

df <- tibble(cell_id=TE013_ATAC.project$cellNames,cluster_id=TE013_ATAC.project$SeuratClusters)
write_csv(df,"~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/processed/TE013-archr-clusters.csv")

# ## Integrating Viper Clusters with scATAC-Seq ----
# print_msg_info(">>> Loading Clusters from VIPER Analysis")
# {
# 	my_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-analysis-with-metacell-metadata.rds")
# 	my_sample_id <- "TE013"
# 	seurat_viper_analysis_data_filename <- file.path( my_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
# 	vp_data <- readRDS(seurat_viper_analysis_data_filename)
# 
# 	table(vp_data$seurat_clusters)
# 	vp_data$cell_id_new <- paste0("TE013#",vp_data$cell_id)
# 
# 	cells_in_common <- intersect(vp_data$cell_id_new,TE013_ATAC.project$cellNames)
# 
# 	TE013_ATAC.project <- TE013_ATAC.project[ cells_in_common ]
# 
# 	index <- match( TE013_ATAC.project$cellNames , vp_data$cell_id_new )
# 
# 	vp_data <- vp_data[,index]
# 	dim(vp_data)
# 
# 	stopifnot(identical( TE013_ATAC.project$cellNames , vp_data@meta.data$cell_id_new ))
# 
# 	TE013_ATAC.project@cellColData$ClustersVIPER <- as.character(vp_data@meta.data$seurat_clusters)
# 
# }

# print_msg_info(">>> Loading VIPER Clusters, run VIPER and assign them to new data")
# {
# 	filename <- "~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-viper-analysis.rds"
# 	vp_seurat <- readRDS(filename)
# 	vp_metadata <- vp_seurat@meta.data
# 	vp_metadata$cell_id <- rownames(vp_metadata)
# 	
# 	vp_metadata$cell_id_new <- paste0("TE013#",vp_metadata$cell_id)
# 	cells_in_common <- intersect(vp_metadata$cell_id_new,TE013_ATAC.project$cellNames)
# 
# 	TE013_ATAC.project <- TE013_ATAC.project[ cells_in_common ]
# 	# 	vp_data$cell_id_new <- paste0("TE013#",vp_data$cell_id)
# 
# 	index <- match( TE013_ATAC.project$cellNames , vp_metadata$cell_id_new )
# 
# 	vp_metadata <- vp_metadata[index,]
# 	dim(vp_metadata)
# 
# 	stopifnot(identical( TE013_ATAC.project$cellNames , vp_metadata$cell_id_new ))
# 
# 	TE013_ATAC.project@cellColData$ClustersVIPER <- as.character(vp_metadata$cluster_id_aREA_assigned)	
# 	table(TE013_ATAC.project@cellColData$ClustersVIPER)
# }

print_msg_info(">>> Loading VIPER Clusters from ingest Analysis")
{
	# my_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/multiome-metadata-ingest.csv")
	my_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE013/ingest-form-ermanno-atac-labels-paneth.csv")
	vp_metadata <- read_csv(my_data_dir)
	
	table(vp_metadata$cluster_id)
	table(vp_metadata$iter_cluster_id)
	vp_metadata$cell_id_new <- paste0("TE013#",vp_metadata$cell_id)
	cells_in_common <- intersect(vp_metadata$cell_id_new,TE013_ATAC.project$cellNames)

	TE013_ATAC.project <- TE013_ATAC.project[ cells_in_common ]
	# 	vp_data$cell_id_new <- paste0("TE013#",vp_data$cell_id)

	index <- match( TE013_ATAC.project$cellNames , vp_metadata$cell_id_new )

	vp_metadata <- vp_metadata[index,]
	dim(vp_metadata)

	stopifnot(identical( TE013_ATAC.project$cellNames , vp_metadata$cell_id_new ))

	# TE013_ATAC.project@cellColData$ClustersVIPER <- as.character(vp_metadata$cluster_id)
	TE013_ATAC.project@cellColData$ClustersVIPER <- as.character(vp_metadata$iter_cluster_id)

}

set.seed(666)
TE013_ATAC.project <- addUMAP(
    ArchRProj = TE013_ATAC.project, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 15, 
    minDist = 0.5, 
    force = TRUE ,
    metric = "euclidean" , 
    seed = 666 
)

p1 <- plotEmbedding(ArchRProj = TE013_ATAC.project, colorBy = "cellColData", name = "SeuratClusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = TE013_ATAC.project, colorBy = "cellColData", name = "ClustersVIPER", embedding = "UMAP")
# p3 <- plotEmbedding(ArchRProj = TE013_ATAC.project, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
# p4 <- ggAlignPlots(p1, p2, p3, type = "h",draw = FALSE)
p4 <- ggAlignPlots(p1, p2, type = "h",draw = FALSE)
plotPDF(p4, name = "Plot-UMAP-TE013-Clusters.pdf", ArchRProj = TE013_ATAC.project, addDOC = FALSE, width = 10, height = 4)
# plotPDF(p1, name = "Plot-UMAP-TE013-Clusters.pdf", ArchRProj = TE013_ATAC.project, addDOC = FALSE, width = 10, height = 4)

print_msg_info(">>> Computing accordance between doublets detected at single cell RNA and ATAC Seq")
{
	filename = "~/Clouds/Dropbox/Data/mouse-organoids/TE013/TE013-scanpy-obs.csv"	
	if ( file.exists(filename) )
	{
		scanpy.obs <- read_csv(filename)
		scanpy.obs <- scanpy.obs %>% dplyr::rename("cell_id"=X1)
		
		atac_doublets.tibble <- doublet_score.obj[[1]]@listData %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
		atac_doublets.tibble$cell_id <- gsub("TE013#","",atac_doublets.tibble$cell_id)
		
		# TE013_ATAC.project@embeddings@listData$UMAP@listData$df
		
		scanpy.obs <- left_join( scanpy.obs , atac_doublets.tibble , suffix = c("","ges_") )
		
		pas_umap.tibble <- vp_data@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
		
		scanpy.obs <- left_join( scanpy.obs , pas_umap.tibble , suffix = c("","pas_") )
		
		umap_plot <- ggplot( data = scanpy.obs %>% arrange(doubletScore), aes( UMAP_1 , UMAP_2 , color = doubletScore ) ) +
			geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
			geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
			viridis::scale_color_viridis(option = "C") +
			xlab( paste0( "UMAP 1") ) +
			ylab( paste0( "UMAP 2") ) +
			coord_fixed() +
			theme_minimal()
		
		pdf( file = file.path( reports.dir , "TE013-pas-umap-doublet-score-atac.pdf") , width = 6 , height = 6 )
			print(umap_plot)
		dev.off()	
		
		umap_plot <- ggplot( data = scanpy.obs %>% arrange(doublet_score), aes( UMAP_1 , UMAP_2 , color = doublet_score ) ) +
			geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
			geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
			viridis::scale_color_viridis(option = "C") +
			xlab( paste0( "UMAP 1") ) +
			ylab( paste0( "UMAP 2") ) +
			coord_fixed() +
			theme_minimal()
		
		pdf( file = file.path( reports.dir , "TE013-pas-umap-doublet-score-rna.pdf") , width = 6 , height = 6 )
			print(umap_plot)
		dev.off()			
		
		scanpy.obs$doubletEnrichFlag <- ifelse( scanpy.obs$doubletEnrich > 1.645 , "doublet" , "singlet" )
		table(scanpy.obs$doubletEnrichFlag)
		
		umap_plot <- ggplot( data = scanpy.obs %>% arrange(doubletEnrichFlag), aes( UMAP_1 , UMAP_2 , color = doubletEnrichFlag ) ) +
			geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
			geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
			scale_color_discrete("viridis") +
			xlab( paste0( "UMAP 1") ) +
			ylab( paste0( "UMAP 2") ) +
			coord_fixed() +
			theme_minimal()
		
		pdf( file = file.path( reports.dir , "TE013-pas-umap-doublet-score-atac-binary.pdf") , width = 6 , height = 6 )
			print(umap_plot)
		dev.off()	
		
		table(scanpy.obs$predicted_doublet)
		umap_plot <- ggplot( data = scanpy.obs %>% arrange(predicted_doublet), aes( UMAP_1 , UMAP_2 , color = predicted_doublet ) ) +
			geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
			geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
			scale_color_discrete("viridis") +
			xlab( paste0( "UMAP 1") ) +
			ylab( paste0( "UMAP 2") ) +
			coord_fixed() +
			theme_minimal()
		
		pdf( file = file.path( reports.dir , "TE013-pas-umap-doublet-score-rna-binary.pdf") , width = 6 , height = 6 )
			print(umap_plot)
		dev.off()			
		
		scanpy.obs %>% filter( predicted_doublet == "TRUE" & doubletEnrichFlag == "singlet" )
		table(scanpy.obs$doubletEnrichFlag)
		table(scanpy.obs$predicted_doublet)
		
		library("VennDiagram")
		universe_set = scanpy.obs$cell_id
		atac_set = scanpy.obs$cell_id[ scanpy.obs$doubletEnrichFlag == "doublet" ]
		rna_set = scanpy.obs$cell_id[ scanpy.obs$predicted_doublet == "TRUE" ]
		
		venn.diagram(
			x = list(universe_set, atac_set, rna_set),
			category.names = c("universe_set" , "atac_set" , "rna_set"),
			na = "remove" ,
			filename = 'doublet_detection_venn_diagramm.png',
			output=TRUE  , imagetype = "png" , height = 1750 , width = 1750 , 
		)
		
		scanpy.obs <- left_join( scanpy.obs , 
														 vp_data@meta.data %>% dplyr::select(cell_id,seurat_clusters) %>% as_tibble() , 
														 suffix = c("","_pas") )
		
		scanpy.obs %>% group_by(seurat_clusters) %>% summarize( doubletEnrichFlag == "doublet" )
		
		scanpy.obs %>%
			group_by(seurat_clusters) %>%
			mutate( countTot = n() ) %>%
			mutate( countC = sum(doubletEnrichFlag == "doublet") ) %>%
			group_by(seurat_clusters, add=TRUE) %>%
			mutate(percentage=paste0(round(100*countC/countTot,2),'%')) %>%
			dplyr::distinct(percentage)
		
		scanpy.obs %>%
			group_by(seurat_clusters) %>%
			mutate( countTot = n() ) %>%
			mutate( countC = sum(predicted_doublet == "TRUE") ) %>%
			group_by(seurat_clusters, add=TRUE) %>%
			mutate(percentage=paste0(round(100*countC/countTot,2),'%')) %>%
			dplyr::distinct(percentage)			
		
	} else {
		print_msg_warn("*** FILE NOT FOUND | SKIPPING ***")
		print_msg_warn(filename)
	}
	
	
}


## CHECK THIS ----
TE013_ATAC.project <- addImputeWeights(TE013_ATAC.project)

markerGenes  <- c("Lgr4",  "Lgr5" , "Atoh1", "Mki67" , "Ascl2" , "Sox4" , "Ung" , 
									"Fabp1" , "Chga" , "Dclk1","Dll1","Neurog3")

p <- plotEmbedding(
    ArchRProj = TE013_ATAC.project, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP" ,
    seed = 666 
    # imputeWeights = getImputeWeights(proj)
)

# Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
# do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
p3 <- cowplot::plot_grid( plotlist = p2 , ncol = 3)
pdf( file.path(reports.dir,"ATAC-Gene-Scores-Markers-on-UMAP.pdf") , width = 10, height = 10)
	print(p3)
dev.off()

# plotPDF(plotList = p2, 
#     name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
#     ArchRProj = TE013_ATAC.project, 
#     addDOC = FALSE, width = 5, height = 5)

## Motif Enrichment in Differential Peaks ----
TE013_ATAC.project <- addGroupCoverages(ArchRProj = TE013_ATAC.project, groupBy = "ClustersVIPER")
# TE013_ATAC.project <- addGroupCoverages(ArchRProj = TE013_ATAC.project, groupBy = "SeuratClusters")
my_pathToMacs2 <- findMacs2()
# my_pathToMacs2 <- "/Users/afpvax/opt/anaconda3/envs/single-cell-env/bin/macs2"
# my_pathToMacs2 <- "/Users/av2729/opt/anaconda3/bin/macs2"
# my_pathToMacs2 <- "/Applications/anaconda3/envs/single-cell-env/bin/macs2"
# my_pathToMacs2 <- "/Applications/anaconda3/bin/macs2"
TE013_ATAC.project <- addReproduciblePeakSet(
    ArchRProj = TE013_ATAC.project, 
    groupBy = "ClustersVIPER",
    # groupBy = "SeuratClusters",
    pathToMacs2 = my_pathToMacs2 ,
    seed = 666 
)
# getPeakSet(TE013_ATAC.project)

TE013_ATAC.project <- addPeakMatrix(TE013_ATAC.project)
TE013_ATAC.project@peakSet@metadata$PeakCallSummary %>% write_csv("~/Clouds/Dropbox/Data/isc/TE013/TE013-peakset.csv")
atac_clusters.tibble <- tibble( cell_id = gsub( "TE013#(.*)","\\1",rownames(TE013_ATAC.project@cellColData) ) ,
																SeuratClustersATAC = TE013_ATAC.project@cellColData$SeuratClusters )
atac_clusters.tibble %>% write_csv("~/Clouds/Dropbox/Data/isc/TE013/TE013-atac-clusters.csv")

##  Identifying Marker Peaks with ArchR ----
markersPeaks <- getMarkerFeatures(
	ArchRProj = TE013_ATAC.project, 
	useMatrix = "PeakMatrix", 
	groupBy = "ClustersVIPER",
	# groupBy = "SeuratClusters",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
)

my_cutoff_string <- "FDR <= 0.1 & Log2FC >= 0.5"
my_cutoff_string.down <- "FDR <= 0.1 & Log2FC <= -0.5"

# markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = my_cutoff_string )
markerList
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
# markerList

## Plotting Marker Peaks in ArchR ----
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks , 
  cutOff = my_cutoff_string,
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = TE013_ATAC.project, addDOC = FALSE)


selected_gene <- c("Lgr5","Ascl2","Dll1","Jag1","Notch1","Mki67")

p <- plotBrowserTrack(
    ArchRProj = TE013_ATAC.project, 
    groupBy = "ClustersVIPER",
    # groupBy = "SeuratClusters",
    geneSymbol = c(selected_gene),
    features = getMarkers(markersPeaks, cutOff = my_cutoff_string, returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)
# grid::grid.newpage()
# grid::grid.draw(p[[selected_gene]])
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 7, ArchRProj = TE013_ATAC.project, addDOC = FALSE)

## Printing Markers Plots
print_msg_info(">>> Printing Markers Plots")
{

	for (my_cluster in colnames(markersPeaks))
	{
		pma <- plotMarkers(seMarker = markersPeaks, name = my_cluster , cutOff = my_cutoff_string , plotAs = "MA")
		plotPDF(pma, name = paste0("Plot-Tracks-With-Features-C",my_cluster), width = 5, height = 5, ArchRProj = TE013_ATAC.project, addDOC = FALSE)			
	}
	
}

# ## Printing Markers Plots as differential between 2 clusters ----
# print_msg_info(">>> Printing Markers Plots as differential between 2 clusters")
# {
# 	
# 	markersPeaks
# 	markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# 	# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
# 	markerList
# 	# markerList$C7
# 	heatmapPeaks <- markerHeatmap(
# 		seMarker = markersPeaks, 
# 		cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
# 		transpose = TRUE
# 	)
# 	draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# 	plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 	
# 	pv <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
# 	# pv
# 	plotPDF(pv, pv, name = "C1-Volcano", width = 5, height = 5, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 	
# 	p <- plotBrowserTrack(
# 		ArchRProj = TE013_ATAC.project, 
# 		groupBy = "ClustersVIPER", 
# 		geneSymbol = c("AR"),
# 		features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C7"],
# 		upstream = 50000,
# 		downstream = 50000
# 	)
# 	plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 	
# 	
# 	table(TE013_ATAC.project$Clusters)
# 	table(TE013_ATAC.project$ClustersVIPER)
# 	
# 	# markersPeaks
# 	markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# 	sapply( markerList@listData , nrow )
# 	plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 	
# 	
# 	## Pairwise Testing Between Groups ----
# 	markersOf_2vsAll <- getMarkerFeatures(
# 		ArchRProj = TE013_ATAC.project, 
# 		useMatrix = "PeakMatrix",
# 		groupBy = "ClustersVIPER",
# 		testMethod = "wilcoxon",
# 		bias = c("TSSEnrichment", "log10(nFrags)"),
# 		useGroups = "2"
# 		# bgdGroups = "Progenitor"
# 	)
# 	pma <- markerPlot(seMarker = markersOf_2vsAll, name = "2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
# 	pma
# 	pv <- markerPlot(seMarker = markersOf_2vsAll, name = "2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
# 	pv
# 	plotPDF(pma, pv, name = "Vim-vs-all-MA-Volcano", width = 5, height = 5, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 	
# }


ggplot( TE013_ATAC.project@cellColData %>% as_tibble(), aes(ClustersVIPER,ReadsInPromoter) ) +
	geom_boxplot() + 
	theme_light()

ggplot( TE013_ATAC.project@cellColData %>% as_tibble(), aes(ClustersVIPER,nFrags) ) +
	geom_boxplot() + 
	theme_light()

# table( TE013_ATAC.project$ClustersVIPER , TE013_ATAC.project$ReadsInPromoter )

library(JASPAR2020)
TE013_ATAC.project <- addMotifAnnotations(ArchRProj = TE013_ATAC.project, 
																					# motifSet = "cisbp", species = "mus musculus" ,version = 2,
																					motifSet = "homer",species = "mus musculus" ,
																					# motifSet = "JASPAR2020", species = "mus musculus" ,
																					# motifSet = "encode",
																					force = TRUE ,
																					name = "Motif")

## Motif Enrichment in Differential Peaks ----
print_msg_info("Motif Enrichment in Differential Peaks")
{
	motifsUp <- peakAnnoEnrichment(
		seMarker = markersPeaks,
		# seMarker = markersOf_2vsAll, 
		ArchRProj = TE013_ATAC.project,
		peakAnnotation = "Motif",
		cutOff = my_cutoff_string 
	)
	
	motifsDo <- peakAnnoEnrichment(
		seMarker = markersPeaks,
		# seMarker = markersOf_2vsAll,
		ArchRProj = TE013_ATAC.project,
		peakAnnotation = "Motif",
		cutOff = my_cutoff_string.down
	)
	
	enrichMotifs <- peakAnnoEnrichment(
		seMarker = markersPeaks,
		# seMarker = markerTest,
		ArchRProj = TE013_ATAC.project,
		peakAnnotation = "Motif",
		cutOff = my_cutoff_string
	)
	
	print_msg_info(">>> Cycling on all clusters to get Up Motifs")
	for ( a_cluster in markersPeaks@colData@rownames )
	{
		print_msg_info(">>> >> Cluster: " ,  a_cluster )
		
		df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,a_cluster])
		df <- df[order(df$mlog10Padj, decreasing = TRUE),]
		df$rank <- seq_len(nrow(df))
		# print(head(df))
		
		ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
			geom_point(size = 1) +
			ggrepel::geom_label_repel(
				data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
				size = 3,
				nudge_x = 5,
				max.overlaps = 20 ,
				color = "black"
			) + theme_ArchR() + 
			ylab("-log10(P-adj) Motif Enrichment") + 
			xlab( paste0("Rank Sorted TFs Enriched for cluster " , a_cluster) ) +
			scale_color_gradientn(colors = paletteContinuous(set = "comet"))
		
		# ggUp
		
		df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,a_cluster])
		df <- df[order(df$mlog10Padj, decreasing = TRUE),]
		df$rank <- seq_len(nrow(df))
		
		ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
			geom_point(size = 1) +
			ggrepel::geom_label_repel(
				data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
				size = 3,
				nudge_x = 5,
				max.overlaps = 20 ,
				color = "black"
			) + theme_ArchR() + 
			ylab("-log10(FDR) Motif Enrichment") +
			xlab( paste0("Rank Sorted TFs Enriched for cluster " , a_cluster) ) +
			scale_color_gradientn(colors = paletteContinuous(set = "comet"))
		
		pdf( file.path(reports.dir,paste0(a_cluster,"-Markers-Motifs-Enriched-down.pdf") ) , width = 5 , height = 8 )
		print(ggDo)
		dev.off()
		pdf( file.path(reports.dir,paste0(a_cluster,"-Markers-Motifs-Enriched-up.pdf") ) , width = 5 , height = 8 )
		print(ggUp)
		dev.off()	
		
		plotPDF(ggUp, ggDo, name = paste0(a_cluster,"-Markers-Motifs-Enriched-both-first-up-second-down"), width = 5, height = 10, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
		
		enrichMotifs.sliced <- enrichMotifs
		N <- 20
		pos <- order(enrichMotifs.sliced[,a_cluster]@assays@data$mlog10Padj,decreasing = T)[1:N]
		enrichMotifs.sliced <- enrichMotifs.sliced[pos,]
		if ( sum(enrichMotifs.sliced[,a_cluster]@assays@data$mlog10Padj == 0) != N )
		{
			heatmapEM <- plotEnrichHeatmap(enrichMotifs.sliced, n = N, transpose = FALSE,
																		 pal = paletteContinuous(set = "comet", n = 10),
																		 labelRows=FALSE,
																		 returnMatrix = FALSE ,
																		 cutOff = 0,
																		 clusterCols = FALSE,
																		 binaryClusterRows=FALSE)
			# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
			plotPDF(heatmapEM, name = paste0(a_cluster,"-Motifs-Enriched-Marker-Heatmap"), width = 6, height = 12, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
			
		}
		
	}
	
}


## ChromVar Analysis ----
# print_msg_info("ChromVar Analysis")
# {
# 	if("Motif" %ni% names(TE013_ATAC.project@peakAnnotation)){
# 		TE013_ATAC.project <- addMotifAnnotations(ArchRProj = TE013_ATAC.project, 
# 																							motifSet = "cisbp", 
# 																							name = "Motif")
# 	}	
# 	TE013_ATAC.project <- addBgdPeaks(TE013_ATAC.project)
# 	TE013_ATAC.project <- addDeviationsMatrix(
# 		ArchRProj = TE013_ATAC.project, 
# 		peakAnnotation = "Motif",
# 		force = TRUE
# 	)	
# 	plotVarDev <- getVarDeviations(TE013_ATAC.project, name = "MotifMatrix", plot = TRUE)
# 	plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 15, height = 15, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 
# 	motifs <- c("AR", "ZEB1", "ZEB2", "MYC", "MYCN", "TWIST1" , "TWIST2")
# 	markerMotifs <- getFeatures(TE013_ATAC.project, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
# 	markerMotifs	
# 	
# 	markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# 	# markerMotifs <- markerMotifs[markerMotifs %in% "z:SREBF1_22"]
# 	markerMotifs
# 	
# 	# p <- plotGroups(ArchRProj = TE013_ATAC.project, 
# 	# 								groupBy = "ClusterVIPER", 
# 	# 								colorBy = "MotifMatrix", 
# 	# 								name = markerMotifs,
# 	# 								imputeWeights = getImputeWeights(TE013_ATAC.project)
# 	# )	
# 	p1 <- plotEmbedding(
# 		ArchRProj = TE013_ATAC.project, 
# 		colorBy = "MotifMatrix", 
# 		name = sort(markerMotifs), 
# 		embedding = "UMAP",
# 		imputeWeights = getImputeWeights(TE013_ATAC.project)
# 	)	
# 	p2 <- lapply(p1, function(x){
# 		x + guides(color = FALSE, fill = FALSE) + 
# 			theme_ArchR(baseSize = 2) +
# 			theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
# 			theme(
# 				axis.text.x=element_blank(), 
# 				axis.ticks.x=element_blank(), 
# 				axis.text.y=element_blank(), 
# 				axis.ticks.y=element_blank()
# 			)
# 	})
# 	
# 	p <- do.call(cowplot::plot_grid, c(list(ncol = 4),p1))
# 	plotPDF(p, name = "Variable-Motif-Deviation-Scores", width = 36, height = 36, ArchRProj = TE013_ATAC.project, addDOC = FALSE)
# 	
# }

## Saving ArchR Project ----
# print_msg_info(">>> Saving ArchR Project")
# {
# 	saveArchRProject(ArchRProj = TE013_ATAC.project, 
# 									 outputDirectory = reports.dir,
# 									 load = FALSE)
# 	
# }

## Integrate ChromVar and VIPER Analyses ----
# print_msg_info(">>> Integrate ChromVar and VIPER Analyses")
# {
# 	getAvailableMatrices(TE013_ATAC.project)
# 	motif_matrix <- getMatrixFromProject(TE013_ATAC.project,"MotifMatrix")
# 	str(motif_matrix,1)
# 	
# 	vpmat <- vp_data@assays$VIPER@scale.data %>% as.matrix()
# 	pas_clusters <- vp_data@meta.data$seurat_clusters
# 	names(pas_clusters) <- colnames(vpmat)
# 	vpmat.stouffered <- doStoufferFromClusters(vpmat,pas_clusters)	 
# 	
# 	motif_mat <- motif_matrix@assays@data$z %>% as.matrix()
# 	pas_clusters <- factor(motif_matrix$ClustersVIPER)
# 	names(pas_clusters) <- colnames(motif_mat)
# 	motifmat.stouffered <- doStoufferFromClusters(motif_mat,pas_clusters)	 	
# 	
# 	rownames(motifmat.stouffered) <- gsub("(.*)_(.*)","\\1",rownames(motifmat.stouffered) )
# 
# 	motifs_and_tfs_in_common <- intersect( rownames(motifmat.stouffered) , rownames(vpmat.stouffered) )
# 	length(motifs_and_tfs_in_common)
# 	
# 	vpmat.stouffered <- vpmat.stouffered[ motifs_and_tfs_in_common , ]
# 	dim(vpmat.stouffered)
# 	colnames(vpmat.stouffered)
# 	motifmat.stouffered <- motifmat.stouffered[ motifs_and_tfs_in_common , ]
# 	dim(motifmat.stouffered)
# 	colnames(motifmat.stouffered)
# 	
# 	print_msg_info(">>> VIPER and MOTIFs integration")
# 	{
# 		cluster_names_list <- colnames(vpmat.stouffered)
# 		
# 		stopifnot(identical(rownames(vpmat.stouffered),rownames(motifmat.stouffered)))
# 		
# 		for (my_cluster_id in cluster_names_list)
# 		{
# 			my_integration.tibble <- tibble( 
# 				tf = rownames(vpmat.stouffered) ,
# 				viper = vpmat.stouffered[,my_cluster_id] , 
# 				atac = motifmat.stouffered[,my_cluster_id] )
# 			
# 			p <- ggplot( my_integration.tibble , aes(x=viper,y=atac,label=tf)) +
# 				geom_point() +
# 				ggrepel::geom_label_repel() +
# 				theme_light()
# 			
# 			pdf( file.path( reports.dir , paste0("Cluster-",my_cluster_id,"-viper-atac-scatter.pdf") ) )
# 				print(p)				
# 			dev.off()
# 			
# 			## Priting Cluster Specific VIPER+ATAC Integration ----
# 			x <- doStouffer(my_integration.tibble %>% dplyr::select(viper,atac) %>% as.matrix())	
# 			names(x) <- my_integration.tibble$tf
# 			x <- sort(x,decreasing = T)
# 			my_tibble <- tibble( rank = 1:length(x) , "TF" = names(x) , z = x)
# 			# View(my_tibble)
# 			
# 			N <- 10
# 			my_tibble.top <- my_tibble %>% top_n(10)
# 			my_tibble.bottom <- my_tibble %>% top_n(10,wt = -z)
#       p <- ggplot( my_tibble , aes( x = rank , y = z , label = TF ) ) +
#         # xlim(-2000,40000) +
#         geom_text_repel( data = my_tibble.top , segment.color = "gray80" ,
#                          color = "darkred" ,
#                          min.segment.length = unit(0, 'lines') , point.padding = 1 , nudge_x = 200 , direction = "y" ,
#                          size = 5 , force = 1 , force_pull = 0.1 , max.overlaps = 40 , seed = 42) +
#         geom_text_repel( data = my_tibble.bottom , segment.color = "gray80" ,
#                          color = "darkblue" ,
#                          min.segment.length = unit(0, 'lines') , point.padding = 1 , nudge_x = -200 , direction = "y" ,
#                          size = 5 , force = 1 , force_pull = 5 , max.overlaps = 40 , seed = 42) +
#         geom_point( alpha = 1 , size = 0.25 ) +
#         # geom_hline( yintercept = high_threshold , alpha = 0.5 , color = "red" ) +
#         # geom_hline( yintercept = low_threshold , alpha = 0.5 , color = "red" ) +
#         ggtitle( "VIPER-Inferred TFs and ATAC-Seq Inferred Motifs (Integration)" ) +
#         xlab("TF Rank") +
#         ylab("Integration Stouffer's Score") +
#         theme_minimal()
# 			
# 			pdf( file.path( reports.dir , paste0("Cluster-",my_cluster_id,"-viper-atac-integration.pdf") ) , width = 4 , height = 8 )
# 				print(p)				
# 			dev.off()			
# 		}
# 	
# 	}
# 	
# 	print_msg_info(">>> UMAP plot of VIPER coordinates and TFs motif enrichment")
# 	{
# 		umap_vp_matrix <- vp_data@reductions$umap@cell.embeddings
# 		rownames(umap_vp_matrix)
# 		
# 		# table(vp_data@meta.data$seurat_clusters)
# 
# 		rownames(motif_mat) <- gsub("(.*)_(.*)","\\1",rownames(motif_mat) )
# 		colnames(motif_mat) <- gsub("(.*)#(.*)","\\2",colnames(motif_mat) )
# 		motif_mat <- t(motif_mat)
# 		
# 		motif_mat.tibble <- motif_mat %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
# 		umap_vp_matrix.tibble <- umap_vp_matrix %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
# 		
# 		my_tibble <- left_join(umap_vp_matrix.tibble,motif_mat.tibble)
# 
# 		my_motif <- "Ascl1"
# 		my_data <- my_tibble %>% dplyr::select(cell_id,UMAP_1,UMAP_2,"Ascl1")
# 		
# 		p <- ggplot( my_data %>% arrange(Ascl1), aes(x=UMAP_1,y=UMAP_2,color=Ascl1)) +
# 			geom_point( alpha = 0.75 , size = 0.15 ) +
# 			# scale_color_gradient2( low = "darkblue" , mid = "white" , high = "darkred" , midpoint = 0 ) +
# 			# scale_color_gradient2( low = "darkblue" , mid = "white" , high = "darkred" , midpoint = 0 ) +
# 			viridis::scale_color_viridis(option = "D") +
# 			theme_ArchR()
# 			# theme_light()
# 		
# 		pdf( file.path( reports.dir , paste0("Motif-",my_motif,"-viper-UMAP-scatter.pdf") ) , width = 4 , height = 4 )
# 			print(p)				
# 		dev.off()
# 			
# 		# print_msg_info(">>> Make 3D Scatter Plot")
# 		# {
# 		# 	# remotes::install_github("rstudio/webshot2")
# 		# 	library(webshot2)
# 		# 	library(car)
# 		# 	
# 		# 	df <- my_tibble	
# 		# 	
# 		# 	library(rgl)
# 		# 	# my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(my_tibble$) )
# 		# 	# df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]			
# 		# 	
# 		# 	my_colors <- viridis::viridis(100)
# 		# 	
# 		# 	rankimize <- function(arg0,n=100) {
# 		# 		arg0 <- round(rank(arg0)/length(arg0)*(n-1) ) + 1
# 		# 		print(summary(arg0))
# 		# 		return(arg0)
# 		# 	}
# 		# 	
# 		# 	colnames(df)
# 		# 	my_tf <- 'Nr3c1'
# 		# 	
# 		# 	require(viridis)
# 		# 	my_colors <- magma(100)
# 		# 	df$my_color <- my_colors[ rankimize( df[,my_tf] %>% pull() ) ]
# 		# 	
# 		# 	plot3d( 
# 		# 		x = df$UMAP_1, y = df$UMAP_2, z = df$UMAP_3,
# 		# 		# col = df$cluster_color , 
# 		# 		col = df$my_color ,
# 		# 		# type = 's', 
# 		# 		# radius = 10,
# 		# 		type = 's', size = 0.5 , 
# 		# 		xlab="UMAP 1", ylab="UMAP 2", zlab="UMAP 3" )
# 		# 	decorate3d(box = FALSE)
# 		# 	
# 		# 	my_filename_3d <- file.path(reports.dir,"TE013-3dscatter.html")
# 		# 	htmlwidgets::saveWidget( rglwidget(width = 400 , height = 400) , my_filename_3d , selfcontained = TRUE )
# 		# 	# writeWebGL( snapshot = FALSE , filename = my_filename_3d ) # ,  width=1000, height=1000)
# 		# 	
# 		# 	library(magick)
# 		# 	# Save like gif
# 		# 	movie3d(
# 		# 		movie="TE013-3dscatter-movie", 
# 		# 		fps = 10 ,
# 		# 		spin3d( axis = c(0, 0, 1), rpm = 3 ) ,
# 		# 		duration = 10 , 
# 		# 		convert = "convert -loop 0 -delay 1x%d %s*.png %s.%s" , 
# 		# 		dir = reports.dir ,
# 		# 		type = "gif", 
# 		# 		clean = TRUE
# 		# 	)						
# 		# 	
# 		# }
# 		
# 	}
# 	
# 		
# }


## Filtering Fragments ----
# atac_fragments.filename <- file.path("/Volumes/ac_lab_scratch/av2729/jia-nepc-organoids/TE013/analysis/TE013/cellranger_arc_count_outs/atac_fragments.tsv")
# atac_fragments.tibble <- read_tsv(atac_fragments.filename,col_names = c("chromosome","start","end","cell_id","counts") )
# table(atac_fragments.tibble$chromosome)
# print_msg_warn(">>> Removing weird chromosomes")
# {
# 	atac_fragments.tibble <- atac_fragments.tibble %>% filter( !grepl( "JH*|GL*",chromosome) )	
# }
# filtered_fragments.filename <- "/Volumes/ac_lab_scratch/av2729/jia-nepc-organoids/TE013/analysis/TE013/cellranger_arc_count_outs/atac_fragments_filtered.tsv"
# write_tsv(atac_fragments.tibble,filtered_fragments.filename,col_names = FALSE)

# print_msg_info(">>> Integrate with Clustering Info from VIPER Analysis (filtering not-matching cells from scRNA-Seq")
# {
#   require(tidyverse)
#   fragment_counts <- readRDS("/Volumes/ac_lab_scratch/av2729/jia-nepc-organoids/TE013/analysis/TE013/cellranger_arc_count_outs/atac_fragments.rds")
#   colnames(fragment_counts) <- c("chromosome","start","end","cell_id","counts")
# 	rna_meta <- readRDS("~/Clouds/Dropbox/Data/mouse-organoids/TE013/TE013-viper-analysis-with-metacell-metadata.rds")
# 	fragment_counts <- fragment_counts %>% filter( cell_id %in% rna_meta$sample_id )
# 	fragment_counts <- left_join( fragment_counts ,
# 																rna_meta %>% dplyr::select(sample_id,cluster_id,cytotrace_score_ges) ,
# 																by = c("cell_id"="sample_id") )
# 	
# 	x <- fragment_counts %>% group_by(cluster_id.x) %>% summarize(openness=sum(counts)/length(unique(cell_id)))
# 	x$cluster_id.x <- x$cluster_id.x %>% as.numeric()-1
# 	
# 	ggplot(x,aes(cluster_id.x,openness)) + 
# 	  geom_histogram(stat = "identity") +
# 	  theme_light(base_size = 20) +
# 	  xlab("Clusters") + 
# 	  ylab("Chromatin Openness")
# 	
# 	# x <- read_csv("~/Clouds/Dropbox/Data/mouse-organoids/TE013/TE013-adata-obs.csv")
# 	# viper_meta <- readRDS("~/Clouds/Dropbox/Data/mouse-organoids/TE013/")
# 	# viper_meta$
# 	
# 	
# 	# "/Users/av2729/Clouds/Dropbox/Data/mouse-organoids/TE013/TE013-adata-obs.csv"
# }




