## 
# Single cell Analysis of TE001 (VIPER Analysis) | Using 4 lists from network made in April 29, 2022
# --------------------------------------------------------------------------------------------------
# system.time({ source("sources/ermanno-data-analysis/TE001/TE001-viper-4-lists.R") })

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

create_workspace("TE001-4-lists-network-clustering")

library(logger)
log_threshold(DEBUG)
log_threshold(DEBUG,index=2)

log_filename <- file.path(reports.dir,"log.txt")
log_appender(appender = appender_file(log_filename),index=2)
# log_formatter(formatter_glue)
log_layout(layout_glue_colors)
log_warn('--- Start of Analysis ----')

require(tidyverse)
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-cpm.rds"
log_info('Readiimg CPM matrix: {filename}')
fs::file_info(filename)
cpm_matrix <- readRDS(filename)
dim(cpm_matrix)

filename <- "~/Clouds/Dropbox/Data/isc/TE001/networks/4-lists/TE001-metacells_unPruned.rds"
log_info('Readiimg network: {filename}')
isc_network <- readRDS(filename)

library(viper)
viper_method <- "none"
log_info("Signature for VIPER: {viper_method}")
log_info("Runnig VIPER analysis")
vpmat <- viper( cpm_matrix ,
								regulon = pruneRegulon(isc_network,50,adaptive = F,eliminate = T) , method = viper_method )
log_info("Runnig VIPER analysis [completed]")

log_info("Reading CytoTRACE analysis data")
x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-cytotrace-data.csv")

index <- match( colnames(vpmat) , x$cell_id )
x <- x[index,]
stopifnot(identical(x$cell_id,colnames(vpmat)))
x <- as.data.frame(x)
rownames(x) <- x$cell_id

library(Seurat)
vp_seurat <- CreateSeuratObject( counts = vpmat , 
																 assay = "VIPER" , 
																 meta.data = x )

# This is needed to fix this Seurat bug: "Error: Must request at least one colour from a hue palette." when you run DimPlot
rownames(vp_seurat@meta.data) <- colnames(vp_seurat) 

vp_seurat@assays$VIPER@scale.data=as.matrix(vp_seurat@assays$VIPER@data) # TODO: Scale vpmat
vp_seurat@meta.data$nCount_VIPER <- 1000
vp_seurat@meta.data$nFeature_VIPER <- nrow(vpmat)

vp_seurat <- RunPCA( vp_seurat , features = rownames(vp_seurat) , npcs = 30 , seed.use = 666 )

library(silhClust)
fc_params <- 	list(max.time=500,temperature=10^9)
log_info("Running FastClust with max.time={fc_params$max.time}, temperature={fc_params$temperature}")
log_info("Running FastClust with res.range=c(0.1,1), NN.range=c(15,51)")
system.time({
	res <- silhClust::SAClustering(vp_seurat,res.range = c(0.1,1),
																 reduction = TRUE ,
																 NN.range = c(15,51) , assay = "VIPER",
																 control = fc_params)
})
log_info("Running FastClust [completed]")

my_tibble <- tibble( cell_id = colnames(vp_seurat) , 
										 clusters = res$seurat_clusters ,
										 cytotrace = vp_seurat$cytotrace_score ,
										 gcs = vp_seurat$cytotrace_gcs )

table(res$seurat_clusters)
table(my_tibble$clusters)

View(res@assays$VIPER@misc$SA.history$par.history)

log_warn('--- END of Analysis ----')
