## 
# Single cell Analysis of TE001 (PHATE Analysis)
# ----------------------------------------------
# source("sources/ermanno-data-analysis/TE001/TE001-phate-analysis.R")

source("../vaxtools/R/utils.R")
isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
create_workspace(run_dir = "TE001-phate-analysis-subnetworks-one-signature")
my_sample_id <- "TE001"
my_seed <- 42
knn_n_neighbor <- 33

library(Seurat)
library(viper)
# seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
# vp_seurat <- readRDS(seurat_viper_analysis_data_filename)
# x <- as.matrix(vp_seurat@assays$VIPER@data)
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data.rds"
vpmat_seurat <- readRDS(filename)
vpmat <- vpmat_seurat@assays$VIPER@scale.data
vpmat_metadata <- vpmat_seurat@meta.data
x <- vpmat
dim(x)
## Saving Metadata Preparation  ----
print_msg_info(">>> Saving Processed Data in CSV (for Scanpy)")
{
	stopifnot( identical( rownames(vpmat_metadata) , colnames(vpmat) ) )

	write_vpmat_to_3_csv( vpmat , "~/Clouds/Dropbox/Data/isc/TE001/" , filetag = paste0(my_sample_id,"-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data") , is_ensembl_to_sym = FALSE )
	write_csv( vpmat_metadata , "~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-metadata.csv" )
}

library(reticulate)
# use_python('/Applications/anaconda3/bin/python')
use_python('/Users/afpvax/opt/anaconda3/envs/single-cell-env/bin/python')
# py_config()
reticulate::py_discover_config(required_module = "phate")
reticulate::import("phate")

## PHATE Analysis ----
print_msg_info(">>> Running PHATE Analysis on Sham Epithelium Sample")
{
	# devtools::install_github("KrishnaswamyLab/phateR")
	library(phateR)
	
	# vpdist <- as.dist( 1-cor(x,method = "spe"))
	vpdist <- as.dist( viperSimilarity(x,50) )
	
	# run PHATE
	x_PHATE <- phate( 
		data = as.matrix(vpdist) ,
		knn.dist.method = "precomputed" ,
		# data = t(x) ,
		# knn.dist.method = "euclidean" ,
		knn = 33 ,
		ndim = 3 ,
		t = 100 ,
		# t = "auto" ,
		npca = 30 ,
		# knn.dist.method = "cosine" ,
		n.jobs = 5 , verbose = TRUE , seed = 666 )
	
	# str(x_PHATE,1)
	# x_PHATE$embedding
	df <- x_PHATE$embedding %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
	stopifnot( identical( df$cell_id , colnames(x) ) )
	df$cluster_id <- vp_seurat$seurat_clusters
	df$sample_id <- vp_seurat$cell_id
	
	stopifnot( identical( rownames(vp_seurat@meta.data) , names(df$sample_id) ) )
	
	vp_seurat@meta.data$PHATE_1.PAS = df$PHATE1
	vp_seurat@meta.data$PHATE_2.PAS = df$PHATE2
	vp_seurat@meta.data$PHATE_3.PAS = df$PHATE3
	
	write_csv(df,"~/Downloads/TE001-phate-pas.csv")
	
	# stopifnot(identical(df$cell_id,my_metadata.obj$sample_id))
	# df$sc_entropy <- 
	# df$stemness_index <- my_metadata.obj$stemness_index
	
	# df$Lgr5 <- x["Lgr5",]
	# df$Top2a <- x["Top2a",]
	# df$Krt19 <- x["Krt19",]
	
	# l <- levels(df$cluster_id)
	# cluster_colors <- c(brewer.pal(9,"Set1"),brewer.pal(6,"Set2"))[1:length(l)]
	# names(cluster_colors) <- l
	
	# df_3 <- df %>% pivot_longer( cols = contains("PHATE") , names_to = "phate_coordinates" , values_to = "phate_values" )
	
	# p <- ggplot( df , aes( x = PHATE2 , y = PHATE3 ,
	# 											# color = Lgr5
	# 											# fill = Top2a
	# 											fill = cluster_id
	# ) ) +
	# 	# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
	# 	# labs(color="Mpo") +
	# 	geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
	# 	# scale_fill_viridis() +
	# 	# scale_fill_discrete() +
	# 	scale_fill_brewer(palette = "Set1") +
	# 	# coord_fixed(ratio = 1) +
	# 	# facet_wrap(~phate_coordinates) +
	# 	theme_light()
	#
	# pdf( file.path(reports.dir,"TE001-pas-phate-2-3.pdf"))
	# 	print(p)
	# dev.off()
	
	p <- ggplot( df , aes(PHATE1, PHATE2 ,
												# color = Lgr5
												# fill = Top2a
												fill = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		theme_light()
	
	pdf( file.path(reports.dir,"TE001-pas-phate.pdf"))
		print(p)
	dev.off()
	
	p <- ggplot( df , aes(PHATE2, PHATE3 ,
												# color = Lgr5
												# fill = Top2a
												fill = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		theme_light()
	
	pdf( file.path(reports.dir,"TE001-pas-phate-2and3axis.pdf"))
		print(p)
	dev.off()		
	
	p <- ggplot( df , aes(PHATE1, PHATE2 ,
												# color = Lgr5
												# fill = Top2a
												fill = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "white",size = 1.5 , stroke = 0.25) +
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		facet_wrap(~cluster_id) +
		theme_light()
	
	pdf( file.path(reports.dir,"TE001-pas-phate-facet.pdf"))
		print(p)
	dev.off()	
	
	print_msg_info(">>> Make 3D Scatter Plot")
	{
		library(car)
		# scatter3d(formula, data)
		# scatter3d( x = df$PHATE1, 
		# 					 y = df$PHATE2, 
		# 					 z = df$PHATE3, 
		# 					 groups = df$cluster_id , 
		# 					 surface = FALSE ,
		# 					 surface.col = RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) ) )
		
		
		library(rgl)
		my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) )
		df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]			
		
		plot3d( 
			x = df$PHATE1, y = df$PHATE2, z = df$PHATE3,
			col = df$cluster_color , 
			# type = 's', 
			radius = 10,
			xlab="PHATE 1", ylab="PHATE 2", zlab="PHATE 3" )
		
		writeWebGL( filename = file.path(reports.dir,"TE001-3dscatter.html") ,  width=1000, height=1000)
		
	}
	
}