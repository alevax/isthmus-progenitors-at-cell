##
# CytoTrace Analysis
# ------------------
# system.time({source("sources/ermanno-data-analysis/TE001/TE001-cytotrace-analysis.R")})

# devtools::install_local("~/Downloads/CytoTRACE_0.3.3.tar.gz")
# library(CytoTRACE)

# results <- CytoTRACE(marrow_10x_expr, ncores = 8, subsamplesize = 1000)

# devtools::install_local("~//Downloads/OncoSig-master/",force = T)
# devtools::install_github("califano-lab/OncoSig",force = T)
# devtools::install_local("~/Downloads/OncoSig-master.zip")

library(dplyr)
library(Seurat)

source("~/Downloads/CytoTRACE/R/zzz.R")
source("~/Downloads/CytoTRACE/R/CytoTRACE.R")
source("~/Downloads/CytoTRACE/R/plotCytoGenes.R")
source("~/Downloads/CytoTRACE/R/plotCytoTRACE.R")

source("../vaxtools/R/utils.R")

my_sample_id <- "TE001"
create_workspace( paste0(my_sample_id,"-cytotrace-run") )

isc_data_dir <- file.path( "~/Clouds/Dropbox/Data/isc/" , my_sample_id )

print_msg_info(">>> Running CytoTrace on " , my_sample_id , " sample gene expression")
{
	print_msg_info(">>> >> Opening connection with GEX Data File ...")
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
	
	my_isc.sdata <- readRDS(seurat_analysis_data_filename)
	
	mat <- as.matrix(my_isc.sdata@assays$RNA@counts)
	print_msg_info(">>> >> Running CytoTrace  ...")
	cytotrace.data <- CytoTRACE(mat , enableFast = FALSE , ncores = 2 )
	# require(pryr)
	format(object.size(cytotrace.data),"MB")
	plotCytoGenes( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-ges-data")) ,numOfGenes = 10)
	
	df <- cbind(my_isc.sdata$UMAP_1.GES,my_isc.sdata$UMAP_2.GES)
	my_clusters <- as.character(my_isc.sdata$seurat_clusters)
	names(my_clusters) <- names(my_isc.sdata$seurat_clusters)
	plotCytoTRACE( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-ges-clusters")) , emb = df , phenotype = my_clusters )
	
	print_msg_info(">>> >> Opening connection with PAS Data File ...")
	seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
	
	my_isc.sdata <- readRDS(seurat_viper_analysis_data_filename)
	stopifnot( identical( names(my_isc.sdata$seurat_clusters) , colnames(my_isc.sdata) ))
	stopifnot( identical( names(cytotrace.data$CytoTRACE) , names(my_isc.sdata$seurat_clusters) ) )
	
	index <- match( names(my_isc.sdata$seurat_clusters) , names(cytotrace.data$CytoTRACE) )
	my_isc.sdata$cytotrace_score_on_ges <- cytotrace.data$CytoTRACE[index]

	print_msg_info(">>> >> Updating  PAS Data File with CytoTrace Scores ...")
	saveRDS(my_isc.sdata,seurat_viper_analysis_data_filename)
	
	stopifnot( identical( names(cytotrace.data$CytoTRACE) , names(cytotrace.data$GCS) ) )
	stopifnot( identical( names(cytotrace.data$CytoTRACE) , names(cytotrace.data$CytoTRACErank) ) )
	
	ct_tibble <- tibble(
		cell_id = names(cytotrace.data$CytoTRACE) ,
		cytotrace_score = cytotrace.data$CytoTRACE , 
		cytotrace_gcs = cytotrace.data$GCS ,
		cytotrace_rank = cytotrace.data$CytoTRACErank
	)
	
	write_csv(ct_tibble,file.path(isc_data_dir,"TE001-cytotrace-data.csv"))
	
	df <- cbind( my_isc.sdata@meta.data$PHATE_1.PAS , my_isc.sdata@meta.data$PHATE_2.PAS )
	rownames(df) <- rownames(my_isc.sdata@meta.data)
	my_clusters <- as.character(as.numeric(my_isc.sdata$seurat_clusters))
	names(my_clusters) <- names(my_isc.sdata$seurat_clusters)
	plotCytoTRACE( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-pas-clusters")) , emb = df , phenotype = my_clusters )
	
	print_msg_info(">>> >> Plotting PHATE 3D Axis with CytoTrace Scores ...")
	{
		
		df <- tibble( cell_id = rownames(my_isc.sdata@meta.data) , 
									PHATE1 = my_isc.sdata@meta.data$PHATE_1.PAS , 
									PHATE2 = my_isc.sdata@meta.data$PHATE_2.PAS ,
									PHATE3 = my_isc.sdata@meta.data$PHATE_3.PAS,
									cluster_id = my_isc.sdata@meta.data$pas_cluster_id ,
									ct_score = cytotrace.data$CytoTRACE )
		
		write_csv( df , file.path( isc_data_dir , "TE001-cytotrace-table.csv" ) )
		
		
		library(rgl)
		library(viridis)
		my_colors <- inferno(100)
		# my_colors <- viridis(100)
		df$ct_score_color <- my_colors[ ceiling(as.numeric(df$ct_score)*99)+1 ]			
		
		my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) )
		df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]			
		
		plot3d( 
			x = df$PHATE1, y = df$PHATE2, z = df$PHATE3,
			col = df$ct_score_color , 
			type = 's', size = 0.5 , 
			xlab="PHATE 1", ylab="PHATE 2", zlab="PHATE 3" )
		
		library(magick)
		# Save like gif
		movie3d(
			movie="PHATE-CytoTRACE-3dscatter", 
			fps = 10 ,
			spin3d( axis = c(0, 0, 1), rpm = 3 ) ,
			duration = 10 , 
			convert = "convert -loop 0 -delay 1x%d %s*.png %s.%s" , 
			dir = reports.dir ,
			type = "gif", 
			clean = TRUE
		)			
		
		writeWebGL( filename = file.path(reports.dir, paste0(my_sample_id,"-PHATE-CytoTrace-3dscatter.html") ) ,  width=1000, height=1000)		
		
		plot3d( 
			x = df$PHATE1, y = df$PHATE2, z = df$PHATE3,
			col = df$cluster_color , 
			type = 's', size = 0.5 , 
			xlab="PHATE 1", ylab="PHATE 2", zlab="PHATE 3" )
		decorate3d(box = FALSE)
		# decorate3d(sub = "Small Intestine Epithelial Cells")
		
		writeWebGL( filename = file.path(reports.dir, paste0(my_sample_id,"-PHATE-Clustering-3dscatter.html") ) ,  width=1000, height=1000)				
		
		# We can indicate the axis and the rotation velocity
		# play3d( spin3d( axis = c(0, 0, 1), rpm = 10), duration = 10 )
		
		library(magick)
		# Save like gif
		movie3d(
			movie = "PHATE-Clustering-3dscatter", 
			fps = 10 ,
			spin3d( axis = c(0, 0, 1), rpm = 3 ) ,
			duration = 10 , 
			convert = "convert -loop 0 -delay 1x%d %s*.png %s.%s" , 
			dir = reports.dir ,
			type = "gif", 
			clean = TRUE
		)	
	
		
	}

			
}
