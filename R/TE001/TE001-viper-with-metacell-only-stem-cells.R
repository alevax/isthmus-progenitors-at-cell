## 
# Single cell Analysis of TE001 (VIPER Analysis)
# ----------------------------------------------
# New approach, metacells were created based on gene expression clustering
# ------------------------------------------------------------------------
# system.time({ source("sources/ermanno-data-analysis/TE001/TE001-viper-with-metacell-only-stem-cells.R") })

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

library(tidyverse)
library(readr)
library(pryr)
require(yaml)

require(viper)
require(ggplot2)
require(viridis)

library(reticulate)
use_python('/Applications/anaconda3/bin/python')
# py_config()

isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
create_workspace(run_dir = "TE001-viper-metacell-no-metaviper-sct-knn33-vpdistCorrFULL-all-regs-stem-cells-ct-score-0.75")
my_sample_id <- "TE001"
my_seed <- 42
knn_n_neighbor <- 33
regulon_n_targets <- 50
cytotrace_score_threshold <- 0.75
isComputingOptimalNumberOfClusters <- FALSE


## Loading counts matrix ----
print_msg_info(">>> Setting constants and paramters [in YAML]")
{
	my_params.list <- list()
	my_params.list$results_dir = reports.dir
	my_params.list$sample_id = my_sample_id
	# my_params.list$cpm_filename <- cpm_filename
	# my_params.list$seurat_original_data_filename <- seurat_original_data_filename
	# my_params.list$seurat_analysis_data_filename <- seurat_analysis_data_filename
	
	my_params.list$seed = my_seed
	my_params.list$knn_n = knn_n_neighbor
	
	# my_params.list$sct_variable_n_feats <- sct_variable_n_feats
	# my_params.list$sct_is_to_scale <- FALSE
	# my_params.list$n_pcs <- n_pcs
	# my_params.list$pcs_to_use <- pcs_to_use
	# my_params.list$clustering_algorithm <- clustering_algorithm
	
	my_yaml_file <- file.path(reports.dir,"run-parameters.yaml")
	yaml::write_yaml(my_params.list %>% as.list(),my_yaml_file)
	
}

## VIPER Run ----
print_msg_info(">>> VIPER Analysis")
{
	require(Seurat)
	require(viper)
	
	seurat_original_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-original-data.rds") )
	metadata_tibble_analysis_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-viper-analysis-with-metacell-metadata.rds") )
	
	print_msg_info(">>> Gene Expression computerd with SCT@scale.data")
	{
		my_gex_seurat <- readRDS(seurat_original_data_filename)
		my_metadata_from_viper_analysis <- readRDS(metadata_tibble_analysis_filename)
		
		selected_cells <- my_metadata_from_viper_analysis$sample_id[ my_metadata_from_viper_analysis$cytotrace_score_ges > cytotrace_score_threshold ]
		
		my_gex_seurat <- my_gex_seurat %>% 
			subset(cells=selected_cells) %>%
			SCTransform( do.scale = T , do.center = T , return.only.var.genes = FALSE ) %>% 			
			# SCTransform( variable.features.n = 3000 , 
			# 						 return.only.var.genes = FALSE ,
			# 						 do.scale = TRUE , 
			# 						 do.center = TRUE ,
			# 						 seed.use = my_seed ) %>% 
			RunPCA(nfeatures.print = 10) 
		my_ges.rankNormed <- as.matrix(my_gex_seurat@assays$SCT@scale.data)
		# print_msg_warn(">>> Using RankNorm")
		# my_ges.rankNormed <- rankNorm(as.matrix(my_gex_seurat@assays$RNA@counts))
	}	

	x <- read_csv( file.path(isc_data_dir,"lists/co-tfs-list-v3.csv") , col_names = FALSE )
	y <- read_csv( file.path(isc_data_dir,"lists/tfs-list-v3.csv") , col_names = FALSE )
	z1 <- read_csv( file.path(isc_data_dir,"lists/reg-extra-list-v3.csv") , col_names = FALSE )
	z2 <- read_csv( file.path(isc_data_dir,"lists/reg-membr-list-v3.csv") , col_names = FALSE )
	# tfs_and_cotfs <- unique(c( x$X1 , y$X1 , c("Rspo","Olfm4","Krt19")))
	tfs_and_cotfs <- unique(c( x$X1 , y$X1 , z1$X1 , z2$X1 ))
	length(tfs_and_cotfs)

	TE001_network_c1 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
	TE001_network_c1 <- pruneRegulon(TE001_network_c1,regulon_n_targets)
	length(TE001_network_c1)
	TE001_network_c1 <- TE001_network_c1[ names(TE001_network_c1) %in% tfs_and_cotfs ]
	length(TE001_network_c1)

	TE001_network_c2 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c2_unPruned.rds"))
	TE001_network_c2 <- pruneRegulon(TE001_network_c2,regulon_n_targets)
	length(TE001_network_c2)
	TE001_network_c2 <- TE001_network_c2[ names(TE001_network_c2) %in% tfs_and_cotfs ]
	length(TE001_network_c2)

	TE001_network_c3 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c3_unPruned.rds"))
	TE001_network_c3 <- pruneRegulon(TE001_network_c3,regulon_n_targets)
	length(TE001_network_c3)
	TE001_network_c3 <- TE001_network_c3[ names(TE001_network_c3) %in% tfs_and_cotfs ]
	length(TE001_network_c3)

	TE001_network_c3and4 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c3and4_unPruned.rds"))
	TE001_network_c3and4 <- pruneRegulon(TE001_network_c3and4,regulon_n_targets)
	length(TE001_network_c3and4)
	TE001_network_c3and4 <- TE001_network_c3and4[ names(TE001_network_c3and4) %in% tfs_and_cotfs ]
	length(TE001_network_c3and4)

	# my_ges.rankNormed <- rankNorm(my_counts)

	system.time({
		vpmat <- viper( eset = my_ges.rankNormed ,
										# regulon = list(TE001_network,
										regulon = list(TE001_network_c1,TE001_network_c2,TE001_network_c3and4) ,
										mvws = 5 , cores = 2 ,
										# regulon = TE001_network ,
										method = "none" , minsize = 30 , verbose = TRUE )
	})

	# vpmat_c1 <- viper( eset = my_ges.rankNormed , regulon = TE001_network_c1 , method = "none" , minsize = 30 , verbose = TRUE )
	# vpmat_c2 <- viper( eset = my_ges.rankNormed , regulon = TE001_network_c2 , method = "none" , minsize = 30 , verbose = TRUE )
	# vpmat_c3and4 <- viper( eset = my_ges.rankNormed , regulon = TE001_network_c3and4 , method = "none" , minsize = 30 , verbose = TRUE )
	# 
	# dim(vpmat_c1)
	# dim(vpmat_c2)
	# dim(vpmat_c3and4)
	# 
	# myMetaViper <- function( mat_list ) {
	# 
	# 	selected_net.matrix <- list()
	# 	i <- 1
	# 	for ( my_mat in mat_list ){
	# 		selected_net.matrix[[i]] <- colSums( abs(my_mat) ) / nrow(my_mat)
	# 		i <- i + 1
	# 	}
	# 	# return(selected_net.matrix)
	# 	x <- do.call( "cbind" , selected_net.matrix )
	# 	ret <- apply( x , 1 , which.max )
	# 	return(ret)
	# }
	# 
	# my_list <- list(vpmat_c1,vpmat_c2,vpmat_c3and4)
	# # apply( sapply( my_list , colnames , `[[` ) , 1 , function(x) all.equal.character(x) )
	# x <- myMetaViper(my_list)
	
	# # regs_in_common <- Reduce( "intersect" , list( rownames(c1.viper) , rownames(c2.viper) , rownames(c3.viper) , rownames(c4.viper) , rownames(c5_c6.viper) ) )
	# regs_in_common <- Reduce( "intersect" , list( rownames(c1.viper) , rownames(c2.viper) , rownames(c3.viper) , rownames(c3and4.viper) ) )
	# length(regs_in_common)
	# 
	# vpmat <- do.call(cbind, list( c1.viper[regs_in_common,] , c2.viper[regs_in_common,] , c3.viper[regs_in_common,] , c3and4.viper[regs_in_common,] ) )
	# str(vpmat,1)
	# 
	# vpmat <- vpmat[ , match( colnames(my_counts) , colnames(vpmat) ) ]
	# stopifnot( identical( colnames(my_counts) , colnames(vpmat) ) )
	
	vpsim <- viperSimilarity(vpmat)
	# vpdist <- as.dist( vpsim )
	vpdist <- as.dist( 1-cor(
		# vpmat[order(rowVars(vpmat),decreasing = T)[1:2000],] ,
		vpmat ,
		method = "pea")
	)

	# saveRDS(object = list(gesmat=my_ges.rankNormed,vpmat=vpmat,vpsim=vpsim,vpdist=vpdist) ,
	# 				file.path(isc_data_dir,"/TE001-viper-analysis-with-metacell.rds") )

}

message(">>> Humanizing data: Mouse --> Human Transformation")
{
	x <- vpmat
	feature_names <- getGeneSymbolsFromEnsemblId(getHumanOrthologousFromMouseEnsemblIds(getEnsemblIdsFromGeneSymbols(rownames(x),"mouse")))
	rownames(x) <- feature_names
	dim(x)
	x <- x[ !is.na(rownames(x)) , ]
	dim(x)
	vpmat_human <- x
}

message(">>> Stemness analysis")
{
	stemness_index <- getStemnessIndex(vpmat_human)	
}

## Patway Analysis ----
message(">>> Pathway Analysis on Hallmarks of Cancer")
{
	source("../vaxtools/R/pathway-analysis.R")
	
	samples.list <- list(vpmat_human)
	names(samples.list) <- "vpmat_human"
	h <- readRDS("../vaxtools/data/MSigDB-gene-sets/msigdb-h-as-regulon.rds")
	
	msigdb.list <- list(h)
	names(msigdb.list) <- c("h")
	
	system.time({
		pathway.results <- list()
		for( msigdb.index in seq_along(msigdb.list) )
		{
			pathway.results[[names(msigdb.list)[[msigdb.index]]]]
			for( i in seq_along(samples.list) )
			{
				filename <- paste0( names(samples.list)[[i]] , "-" , names(msigdb.list)[[msigdb.index]] )
				filename.pdf <- paste0(filename,".pdf")
				filename.xls <- paste0(filename,".xls")
				pathway.results[[ names(msigdb.list)[[msigdb.index]] ]][[i]] <- doPathwayAnalysis( samples = samples.list[[i]] , filename = filename.pdf , 
																																													 msigdb = NULL ,
																																													 regulon = msigdb.list[[msigdb.index]] ,
																																													 save2excel = FALSE , filename.xls = filename.xls , .DEBUG = FALSE )
				# doPathwayAnalysis( samples = samples.list[[i]] , filename = names(samples.list)[[i]] , msigdb = c5 , regulon = genesSetList.as.regulon , threshold = 5 )
			}
		}
	})        
	
	filetag <- "h"
	vpmat_pathways <- unlist(pathway.results,recursive = F)[[filetag]]
	# filename <- file.path( processed_data.dir , "vpmat-cancer-hallmarks.rds")
	# message("- Saving Clones Hallmark Of Cancer Analysis in " , filename )
	# saveRDS( vpmat_pathways , filename )
}

## Metadata Preparation (e.g. for scanpy) ----
message(">>> Metadata Preparation")
{
	vpmat_metadata <- tibble( cell_id = colnames(vpmat) ,
														stemness_index = stemness_index )
	
	pathways <- as_tibble(vpmat_pathways)
	pathways_names <- rownames(vpmat_pathways)
	cell_id <- colnames(vpmat_pathways)
	pathways <- as_tibble(t(pathways))
	colnames(pathways) <- pathways_names	
	pathways <- bind_cols( cell_id = cell_id , pathways )
	
	vpmat_metadata <- left_join( vpmat_metadata , pathways )		
	
	# vpmat_pathways <- t( vpmat_pathways )
	# index <- match( vpmat_metadata$cell_id , rownames(vpmat_pathways) )
	# vpmat_metadata <- as_tibble( do.call( cbind , list( vpmat_metadata , vpmat_pathways[ index , ] ) ) )
}

## Clustering with PAM ----
print_msg_info(">>> Clustering")
{
	require(fpc)
	require(factoextra)
	require(cluster)
	
	set.seed(my_seed)
	
	if (isComputingOptimalNumberOfClusters)
	{
		# Elbow method
		message(">>> Computing Optimal Number of Clusters using the Elbow method")
		pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-kmeans-elbow-optimal-n-clust.pdf") ))
		fviz_nbclust( as.matrix(vpmat) , diss = vpdist , kmeans , method = "wss" , verbose = TRUE ) +
			# geom_vline(xintercept = 4, linetype = 2)+
			labs(subtitle = "Elbow method using KMEANS")
		dev.off()
		pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-pam-elbow-optimal-n-clust.pdf") ))
		fviz_nbclust( as.matrix(vpmat) , diss = vpdist , pam , method = "wss" , verbose = TRUE ) +
			# geom_vline(xintercept = 4, linetype = 2)+
			labs(subtitle = "Elbow method using PAM")
		dev.off()
		
		# Silhouette method
		message(">>> Computing Optimal Number of Clusters using the Silhouette method")
		pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-kmeans-silhouette-optimal-n-clust.pdf") ))
		fviz_nbclust( x = as.matrix(vpdist) , diss = vpdist , kmeans, method = "silhouette" , verbose = TRUE ) +
			labs(subtitle = "Silhouette method using KMEANS")
		dev.off()
		pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-pam-silhouette-optimal-n-clust.pdf") ))
		fviz_nbclust( as.matrix(vpdist) , diss = vpdist , pam, method = "silhouette" , verbose = TRUE ) +
			labs(subtitle = "Silhouette method using PAM")
		dev.off()
		
		message(">>> Computing Optimal Number of Clusters using the Gap statistic")
		# Gap statistic
		# nboot = 50 to keep the function speedy.
		# recommended value: nboot= 500 for your analysis.
		# Use verbose = FALSE to hide computing progression.
		pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-pam-gap-stat-optimal-n-clust.pdf") ))
		fviz_nbclust( as.matrix(vpmat) , diss = vpdist , pam,  method = "gap_stat", nboot = 50 , verbose = TRUE ) +
			labs(subtitle = "Gap statistic method using PAM")
		dev.off()
		pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-kmeans-gap-stat-optimal-n-clust.pdf") ))
		fviz_nbclust( as.matrix(vpmat) , diss = vpdist , kmeans,  method = "gap_stat", nboot = 50 , verbose = TRUE) +
			labs(subtitle = "Gap statistic method using KMEANS")
		dev.off()
	}
	
	pam_k_obj <- pamk( vpdist, k = 2:10 , scaling = FALSE , stand = FALSE , diss = TRUE , metric = "euclidean" , seed = my_seed )
	# pam_k_obj <- pamk( t(vpmat[order(rowVars(vpmat),decreasing = T)[1:500],]) , k = 2:10 , scaling = FALSE , stand = FALSE , diss = FALSE , metric = "euclidean" , seed = my_seed )
	
	pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-pam-sil.pdf") ))
	print( fviz_silhouette( silhouette(pam_k_obj$pamobject) ) )
	dev.off()
	
	cluster_id <- pam_k_obj$pamobject$clustering
	index <- match( vpmat_metadata$cell_id , names(cluster_id) )
	cluster_id <- cluster_id[ index ]
	vpmat_metadata$cluster_id_pam <- as.factor( cluster_id )
	
	sil_score <- pam_k_obj$pamobject$silinfo$widths[,"sil_width"]
	index <- match( vpmat_metadata$cell_id , names(sil_score) )
	sil_score <- sil_score[ index ]
	vpmat_metadata$sil_score_pam <- sil_score
	
	stopifnot( identical( names(sil_score) , names(cluster_id) ) )
	stopifnot( identical( colnames(vpmat) , names(cluster_id) ) )
	
	# message(">>> Computing iter-PAM")
	# {
	#
	# 	vpmat_c1 <- viper( eset = rankNorm(my_counts[,cluster_id == 1]) , regulon = TE001_network , method = "none" , minsize = 30 , verbose = TRUE )
	# 	str(vpmat_c1)
	#
	# 	vpsim_c1 <- viperSimilarity(vpmat_c1,50)
	# 	vpdist_c1 <- as.dist( vpsim_c1 )
	# 	pam_k_obj_c1 <- pamk( vpdist_c1, k = 2:10 , scaling = FALSE , stand = FALSE , diss = TRUE , metric = "euclidean" , seed = my_seed )
	#
	# 	pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-c1-pam-sil.pdf") ))
	# 		print( fviz_silhouette( silhouette(pam_k_obj_c1$pamobject) ) )
	# 	dev.off()
	#
	# 	selected_cells_c1 <- names( ( pam_k_obj_c1$pamobject$silinfo$widths[,"sil_width"] > 0.5 ) == TRUE )
	#
	# 	vpmat_c2 <- viper( eset = rankNorm(my_counts[,cluster_id == 2]) , regulon = TE001_network , method = "none" , minsize = 30 , verbose = TRUE )
	# 	str(vpmat_c2)
	#
	# 	vpsim_c2 <- viperSimilarity(vpmat_c2,50)
	# 	vpdist_c2 <- as.dist( vpsim_c2 )
	# 	pam_k_obj_c2 <- pamk( vpdist_c2, k = 2:10 , scaling = FALSE , stand = FALSE , diss = TRUE , metric = "euclidean" , seed = my_seed )
	#
	# 	pdf( file.path(reports.dir,paste0(my_sample_id,"-vpmat-c2-pam-sil.pdf") ))
	# 		print( fviz_silhouette( silhouette(pam_k_obj_c2$pamobject) ) )
	# 	dev.off()
	#
	# 	selected_cells_c2 <- names( ( pam_k_obj_c2$pamobject$silinfo$widths[,"sil_width"] > 0.5 ) == TRUE )
	#
	# }
	
}

## UMAP ----
print_msg_info(">>> UMAP over PAS")
{
	library(ggplot2)
	library(ggrepel)
	library(umap)
	# library(tidyverse)
	
	stopifnot( identical( labels(vpdist),vpmat_metadata$cell_id ) )
	
	x <- vpdist
	# x <- t(vpmat[order(rowVars(vpmat),decreasing = T)[1:500],])
	
	# create a new settings object with n_neighbors set to 5
	custom.settings <- umap.defaults
	custom.settings$random_state <- 666
	custom.settings$input <- "dist"
	# custom.settings$input <- "data"
	custom.settings$n_neighbors <- 15
	custom.settings$n_components <- 2
	# custom.settings$metric <- "correlation"
	# custom.settings$metric <- "euclidean"
	# custom.settings$metric <- "manhattan"
	custom.settings$verbose <- TRUE
	
	set.seed(666)
	# y.umap <- umap( t(x) , config = custom.settings )
	y.umap <- umap( as.matrix(x) , config = custom.settings )
	
	umap.df <- as.data.frame( y.umap$layout )
	# umap.df$cell_line = ifelse( grepl( "LNCaP" , colnames(vpmat) ) , "clones" ," single_cells" )
	# 
	# umap.df$net_enrichment <- getNETenrichment(vpmat,net.regulon.obj)
	# umap.df$beltran_net_enrichment <- getNETenrichment(vpmat,beltran.net.regulon.obj)
	# umap.df$FOXM1 <- vpmat["Foxm1",]
	umap.df$stemness_index <- vpmat_metadata$stemness_index
	umap.df$cluster_id <- vpmat_metadata$cluster_id_pam
	umap.df$sil_score <- vpmat_metadata$sil_score
	
	# # umap.df <- cbind( umap.df , t(vpmat_pathways ))
	# umap.df$Lgr4 <- vpmat["Lgr4",]
	# umap.df$Lgr5 <- vpmat["Lgr5",]
	# umap.df$Olfm4 <- vpmat["Olfm4",]
	umap.df$cell_id <- vpmat_metadata$cell_id
	# umap.df$Atoh1 <- vpmat["Atoh1",]
	# umap.df$Ascl2 <- vpmat["Ascl2",]
	
	# library(mltools)
	
	# tmp <- model.matrix(~umap.df$model-1)
	# colnames(tmp) <- gsub("(umap.df\\$model)(.*)","\\2",colnames(tmp))
	# umap.df <- cbind(umap.df,tmp)
	# 
	# umap.df$NPK <- as.factor( umap.df$NPK )
	
	umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = stemness_index , shape = cluster_id ) ) +
		# umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = NPK , shape = castration_status ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 0.5 , size = 5 ) +
		scale_color_viridis() +
		geom_density_2d(lineend = "butt", color = "black",size=0.5,alpha=0.5) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "survival_group" , values = c("green","blue","red")) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-analysis.pdf") , width = 16 , height = 12 )
	print(umap_plot)
	dev.off()
}

# ## Saving Metadata Preparation  ----
# print_msg_info(">>> Saving Processed Data in CSV (for Scanpy)")
# {
# 	stopifnot( identical( rownames(umap.df) , colnames(vpmat) ) )
# 	stopifnot( identical( rownames(umap.df) , names(stemness_index) ) )
# 	
# 	write_vpmat_to_3_csv( vpmat , "~/Clouds/Dropbox/Data/isc/TE001/" , filetag = paste0(my_sample_id,"-with-metacell") , is_ensembl_to_sym = FALSE )
# 	write_csv( vpmat_metadata , "~/Clouds/Dropbox/Data/isc/TE001/TE001-with-metacell-metadata.csv" )
# }

## VIPER Seurat ----
message(">> VIPER Analysis with Seurat")
{
	require(Seurat)
	stopifnot( identical( colnames(vpmat) , vpmat_metadata$cell_id ) )
	vpmat_metadata.df <- vpmat_metadata %>% as.data.frame()
	rownames(vpmat_metadata.df) <- vpmat_metadata.df$cell_id
	# vpsim <- viperSimilarity(vpmat,50)
	# vpdist <- as.dist( vpsim )		
	vp_seurat <- CreateSeuratObject( counts = vpmat , 
																	 assay = "VIPER"  , 
																	 meta.data = vpmat_metadata.df )
	
	# This is needed to fix this Seurat bug: "Error: Must request at least one colour from a hue palette." when you run DimPlot
	rownames(vp_seurat@meta.data) <- colnames(vp_seurat) 
	
	# vp_seurat <- vp_seurat %>% ScaleData(  do.scale = TRUE,
	# 																			 do.center = TRUE,
	# 																			 scale.max = 10)
	vp_seurat@assays$VIPER@scale.data=as.matrix(vp_seurat@assays$VIPER@data) # TODO: Scale vpmat
	vp_seurat@meta.data$nCount_VIPER <- 1000
	vp_seurat@meta.data$nFeature_VIPER <- nrow(vpmat)
	
	# vp_seurat@meta.data$clusters_ges <- seurat_analysis$seurat_clusters
	
	vp_seurat <- RunPCA( vp_seurat , features = rownames(vp_seurat) , npcs = 30 , seed.use = my_seed )
	# vp_seurat <- RunICA(vp_seurat,features=rownames(vp_seurat))
	
	n_pcs <- 30	
	pdf( file.path(reports.dir,"pas-elbow-plot.pdf") )
	print(ElbowPlot(vp_seurat))
	dev.off()
	# vp_seurat <- JackStraw( vp_seurat , assay = "VIPER" , reduction = "pca",dims = n_pcs,verbose = TRUE)
	# # head(JS(object = sc_data[['pca']], slot = 'empirical'))	
	# vp_seurat <- ScoreJackStraw( vp_seurat , assay = "VIPER" , reduction = "pca",dims = 1:n_pcs)
	# pdf( file.path(reports.dir,"pas-jackstraw-plot.pdf") )
	# print( JackStrawPlot( vp_seurat , dims = 1:n_pcs , reduction = "pca", xmax = 0.1, ymax = 0.3) )
	# dev.off()
	# 
	# message(">>> Selecting only PC that are statistically significant using the JackStraw algorithm")
	# pca_significant_compounents <- vp_seurat@reductions$pca@jackstraw@overall.p.values
	# pcs_to_use <- pca_significant_compounents[,1][ pca_significant_compounents[,2] < 0.05 ]
	# 
	# x <- max(pcs_to_use)
	# x <- ifelse( x %% 2 == 0 , x , x+1 )
	# pcs_to_use <- 1:x
	# pcs_to_use
	# pcs_to_use <- 1:4
	pcs_to_use <- 1:30
	
	my_ssn_graph <- FindNeighbors( as.matrix(vpdist) , 
																 # dims = pcs_to_use , 
																 assay = "VIPER" , 
																 distance.matrix = TRUE , 
																 verbose = TRUE , 
																 k.param = knn_n_neighbor ,
																 # annoy.metric = "cosine" ,
																 annoy.metric = "euclidean" ,
																 compute.SNN = TRUE )
	vp_seurat@graphs <- my_ssn_graph		
	
	# vp_seurat <- FindNeighbors( vp_seurat , 
	# 													# dims = pcs_to_use , 
	# 													assay = "VIPER" , 
	# 													verbose = TRUE , 
	# 													k.param = knn_n_neighbor ,
	# 													compute.SNN = TRUE ,
	# 													annoy.metric = "cosine" ,
	# 													# distance.matrix = as.matrix(vpdist) ,
	# 													do.plot = FALSE )		
	
	## Silhouette Analysis ----
	vp_data_dist.matrix <- vpdist
	my_sil.df <- tibble( resolution = seq(0.01,1.5,by = 0.01) , sil_avg = 0 , sil_median = 0 , n_clust = 0 )
	for ( a_res in my_sil.df$resolution )
	{
		message("--- Silhouette score computation for resolution :" , a_res )
		vp_seurat <- FindClusters( vp_seurat , 
															 # graph.name = "VIPER_snn" ,
															 graph.name = "snn" ,
															 resolution = a_res , 
															 verbose = FALSE ,
															 modularity.fxn = 1 ,
															 # algorithm = 4 , # Leiden
															 algorithm = 1 , # Louvain
															 seed.use = my_seed
		)
		
		if ( nlevels(vp_seurat$seurat_clusters) == 1 ) next ;
		
		require(cluster)
		s <- silhouette( as.integer(vp_seurat$seurat_clusters) , vp_data_dist.matrix )
		# pdf( file.path(reports.dir,paste0("sil-res-",a_res,".pdf")) )
		require(factoextra)
		x <- fviz_silhouette(s,print.summary = FALSE)			
		y <- sapply( levels(x$data$cluster) , function(i) mean( x$data$sil_width[ x$data$cluster == i ] ) )
		z <- sapply( levels(x$data$cluster) , function(i) median( x$data$sil_width[ x$data$cluster == i ] ) )
		my_sil.df$sil_avg[ my_sil.df$resolution == a_res ] = mean(y)
		my_sil.df$sil_median[ my_sil.df$resolution == a_res ] = median(z)
		my_sil.df$n_clust[ my_sil.df$resolution == a_res ] = nlevels(x$data$cluster)
		# dev.off()
		
		# View(my_sil.df)    
	}
	
	message(">>> >> Printing Optimal Clustering Plot")
	{
		library(ggrepel)
		# p <- ggplot( my_sil.df , aes( x = resolution , y = sil_median ) ) +
		p <- ggplot( my_sil.df , aes( x = resolution , y = sil_avg ) ) +
			geom_line() +
			# geom_text( aes(x = resolution , y = 0.1 , label = n_clust) , label.size = 3 ) +
			# geom_text( aes(x = resolution , y = sil_median+0.01 , label = n_clust) , label.size = 3 ) +
			geom_text_repel( aes( label = n_clust ) , segment.color = "gray" ) +
			theme_light()
		
		pdf( file.path(reports.dir,"pas-optimal-clustering.pdf"))
		print(p)
		dev.off()
	}
	
	my_res <- my_sil.df$resolution[ which.max(my_sil.df$sil_avg) ]
	print_msg_warn(">>> Using MANUAL resolution")
	my_res <- 0.27
	# my_res <- 0.7
	# my_res <- 0.65
	
	# vp_seurat <- FindClusters(vp_seurat, resolution = 2 , verbose = TRUE,modularity.fxn = 2,algorithm = 4,random.seed = 666 ) 
	vp_seurat <- FindClusters( vp_seurat , 
														 assay = "VIPER"  ,
														 # graph.name = "VIPER_snn" ,
														 graph.name = "snn" ,
														 # graph.name = "VIPER_snn" ,
														 resolution = my_res , 
														 verbose = TRUE , 
														 modularity.fxn = 1 , 
														 # algorithm = 4 , # Leiden
														 group.singletons = TRUE ,
														 algorithm = 1 , # Louvain
														 random.seed = my_seed ) 
	
	print(table(vp_seurat$seurat_clusters))
	
	# print_msg_warn("\t *** Not Running RunUMAP from Seurat because continuously issues on umap-learn not recognized ***")
	vp_seurat <- RunUMAP( vp_seurat ,
												assay = "VIPER" ,
												n.components = 3L ,
												# features = prolif_genes ,
												# assay = "RNA" ,
												graph = "snn" ,
												# graph = "VIPER_snn" ,
												# graph.name = "snn" ,
												# n.neighbors = knn_n_neighbor ,
												# dims = pcs_to_use ,
												reduction = "pca" ,
												umap.method = "umap-learn" ,
												# umap.method = "uwot" ,
												# metric = "cosine" ,
												metric = "euclidean" ,
												verbose = FALSE , seed.use = my_seed )
	
	print_msg_info(">>> Creating metadata (clustering.tibble)")
	{
		my_clusters <- as.integer(vp_seurat@meta.data$seurat_clusters)
		names(my_clusters) <- colnames(vp_seurat)
		my_clusters <- factor(my_clusters)
		clustering.tibble <- tibble( sample_id = names(my_clusters) ,
																 cluster_id = my_clusters , 
																 sample_ordering = 1:length(my_clusters) ,
																 cluster_reliability = n1platform::clusterReliability(my_clusters,vpsim) ,
																 sil_score = silhouette(as.integer(my_clusters),vpsim)[,3] )
		clustering.tibble %>% group_by(cluster_id) %>% summarise(avg_sil=mean(sil_score))
		
		clustering.tibble <- clustering.tibble %>% 
			arrange(cluster_id,desc(sil_score)) %>%
			mutate(sil_score_order = 1:length(sil_score) )
		# arrange(cluster_id,desc(cluster_reliability)) %>%
		# mutate(cluster_reliability_order = 1:length(cluster_reliability) )
		clustering.tibble <- left_join( clustering.tibble , vpmat_metadata , by = c("sample_id"="cell_id") , suffix = c("","_vpmat_metadata") )	
		
		clustering.tibble <- clustering.tibble %>% arrange(sample_ordering)
		
		print_msg_info(">>> Integration with gene expression scores")
		{
			genes_to_add <- c("Lgr4","Lgr5","Tnfrsf19","Olfm4","Ascl2","Atoh1","Krt19","Top2a", 
												"Ung","Clu","Ly6a","Rad51","Smarcal1","Ihh",
												"Mki67","Fzd5","Fzd7","Dclk1","Rspo3","Wnt4","Wnt5a","Uri1")
			
			# x <- readRDS( file.path(isc_data_dir,"/TE001-seurat-analysis-data.rds") )
			x <- my_gex_seurat
			
			# tmp <- x@meta.data %>% dplyr::select(cell_id,UMAP_1.GES,UMAP_2.GES,clusters_ges=seurat_clusters,stemness_index.ges,singleR_labels)
			# vp_seurat@meta.data <- left_join(vp_seurat@meta.data,tmp,by=c("cell_id"="cell_id"))	%>% as.data.frame()
			
			rownames(vp_seurat@meta.data) <- colnames(vp_seurat)
			
			# x <- x@assays$RNA@counts %>% RelativeCounts(1e4) %>% as.matrix()
			# x <- rankNorm(x)
			
			x <- as.matrix(x@assays$SCT@scale.data)
			# x <- my_ges.rankNormed
			
			mat <- as.matrix(vp_seurat@assays$VIPER@data)
			stopifnot( identical( colnames(x) , colnames(mat) ) )
			x <- x[rownames(x) %in% genes_to_add,]
			mat_genes <- x
			rownames(x) <- paste0( "gene_", rownames(x) )
			x <- t(x) %>% as.data.frame() %>% rownames_to_column("sample_id") %>% as_tibble()
			clustering.tibble <- left_join( clustering.tibble , x , by = c("sample_id"="sample_id") )
			dim(clustering.tibble)
		}
		
		print_msg_info(">>> >> Integration with gene expression clustering")
		{
			x <- readRDS( file.path(isc_data_dir,"/TE001-seurat-analysis-data.rds") )
			gene_expression_clusters <- x$seurat_clusters
			my_mt_percent <- x@meta.data$mt_percent
			names(my_mt_percent) <- rownames(x@meta.data)
			clustering.tibble$gene_expression_clusters <- gene_expression_clusters[ match( clustering.tibble$sample_id , names(gene_expression_clusters) ) ]
			clustering.tibble$mt_percent <- my_mt_percent[ match( clustering.tibble$sample_id , names(my_mt_percent) ) ]
		}	
		
		print_msg_info(">>> >> Integration with protein avtivity clustering")
		{
			metadata_tibble_analysis_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-viper-analysis-with-metacell-metadata.rds") )
			pas_metadata <- readRDS(metadata_tibble_analysis_filename)
			
			print(table(pas_metadata$cluster_id))
			
			pas_clusters <- pas_metadata$cluster_id
			names(pas_clusters) <- pas_metadata$sample_id
			clustering.tibble$pas_clusters <- pas_clusters[ match( clustering.tibble$sample_id , names(pas_clusters) ) ]
		}	
		
		print_msg_info(">>> >> Integration with CytoTRACE stemness inference")
		{
			df <- read_csv(file.path(isc_data_dir, paste0(my_sample_id,"-cytotrace-table.csv")) )
			
			index <- match( clustering.tibble$sample_id , df$cell_id )
			df <- df[ index , ]
			stopifnot( identical( df$cell_id , clustering.tibble$sample_id ) )
			
			clustering.tibble$cytotrace_score_ges <- df$ct_score		
			
			stopifnot( identical( rownames(vp_seurat@meta.data) , clustering.tibble$sample_id ) )
			vp_seurat@meta.data$cytotrace_score_ges <- df$ct_score		
		}
	}
	
	## Stouffer's Integration ----
	print_msg_info(">>> Selecting markers with Stouffer's integrations")
	{
		source("../single-cell-pipeline/functions/cluster-functions.R")
		source("../single-cell-pipeline/functions/process-utils.R")
		source("../single-cell-pipeline/functions/viper-utils.R")
		
		mat <- as.matrix(vp_seurat@assays$VIPER@data)
		print_msg_info(">>> >> Applying Stouffer's Integration")
		index <-  match( colnames(mat) , clustering.tibble$sample_id )
		stopifnot( identical( colnames(mat) , clustering.tibble$sample_id[index] ) )
		clusters <- clustering.tibble$cluster_id[index]
		weigths <- clustering.tibble$sil_score[index]
		names(clusters) <- clustering.tibble$sample_id[ index ]
		names(weigths) <- clustering.tibble$sample_id[ index ]
		
		print_msg_info(">>> >> >> Selecting only samples with weights > 0.25 for integrations")
		{
			index <- weigths > 0.25
			mat <- mat[,index]
			clusters <- clusters[index]
			weigths <- weigths[index]
		}
		
		vpmat.stouffered <- doStoufferFromClusters( mat , clusters = clusters , weights = weigths )
		# vpmat.stouffered <- doStoufferFromClusters( mat , clusters = clusters )
		top_markers.pas <- CBCMRs(vpmat.stouffered,5)
		
		colnames(vpmat.stouffered) <- paste0("Cluster_", levels(clusters) )
		# write_csv( vpmat.stouffered %>% as.data.frame() %>% rownames_to_column("proteins") ,
		# 					 file.path(isc_data_dir , paste0(my_sample_id,"-viper-stouffer-with-metacell.csv") )
		# )
	}
	
	# # markers <- FindAllMarkers(vp_seurat, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.25 , test.use = 't' )
	# markers <- FindAllMarkers( vp_seurat , only.pos = TRUE )#, min.pct = 0.25, logfc.threshold = 0.25,test.use = 't')
	# # markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
	# # top_markers <- markers %>% group_by(cluster) %>% filter(myAUC > 0.99 & avg_diff > 4) %>% top_n(10,wt = myAUC )
	# # top_markers <- markers %>% group_by(cluster) %>% top_n(10,wt = avg_logFC )
	# n_top <- 10
	# top_markers <- markers %>% group_by(cluster) %>% top_n( n_top , wt = c(-p_val_adj) ) %>% top_n( n_top , wt = avg_logFC )
	# dim(top_markers)
	# top_markers.pas <- top_markers$gene
	
	rownames(vp_seurat@meta.data) <- colnames(vp_seurat) 
	
	p1 <- DimPlot( vp_seurat , reduction = "umap", pt.size = 0.5 , label = TRUE )
	
	require(scales)
	# p2 <- DoHeatmap(vp_seurat,features = features_to_show , raster = FALSE ,assay = "VIPER",draw.lines = TRUE,angle = 0,group.bar = TRUE ) + 
	p2 <- DoHeatmap(vp_seurat,features = top_markers.pas , raster = FALSE , assay = "VIPER" , draw.lines = TRUE , angle = 0 , group.bar = TRUE ) +
		scale_fill_distiller(palette = "RdBu") + 
		theme(axis.text=element_text(size=5) )
	
	pdf( file = file.path(reports.dir, paste0(my_sample_id,"-pas-clusters-umap.pdf")) )
	print(p1)
	dev.off()
	
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-pas-clusters-heatmap.pdf")) )
	print(p2)
	dev.off()
	
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-pas-clusters-combo.pdf")) , width = 12 , height = 6)
	print( cowplot::plot_grid(p1,p2,ncol = 2 ) )#, labels = c('A', 'B') ) )
	dev.off()
	
	g <- vp_seurat %>% FeaturePlot(features = c("Lgr4","Lgr5","Atoh1","Olfm4","Krt19") , pt.size = 0.5 , order = TRUE , cols = c("blue","yellow"), blend = FALSE )
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-pas-feat-plot.pdf")) )
	print(g)	
	dev.off()			
	
}

# Heatmap Printing on VIPER Data ----
print_msg_info(">>> Rendering clustering solution with ComplexHeatmap")
{
	require(ComplexHeatmap)
	require(circlize)
	require(viridis)
	require(dplyr)
	
	isUsingRaster <- TRUE
	
	n_top <- 20
	# selected_features <- CBCMRs(vpmat.stouffered,5)
	selected_features <- Reduce( union , apply( vpmat.stouffered , 2 , function(x) { x <- sort(x,decreasing = T) ; names(c(head(x,n_top))) } ) )
	length(selected_features)
	
	df_annot.df <- clustering.tibble 
	# index <- order( df_annot.df$sil_score_order )
	# df_annot.df <- df_annot.df[index,]
	# colnames(df_annot.df)
	
	sample_ordering_to_print <- order( df_annot.df$sil_score_order )
	
	l <- levels(df_annot.df$cluster_id)
	# cluster_colors <- brewer.pal(length(l),"Set1")
	cluster_colors <- c(brewer.pal(9,"Set1"),brewer.pal(6,"Set2"))[1:length(l)]
	names(cluster_colors) <- l
	
	l <- levels(df_annot.df$cluster_id_pam)
	# cluster_colors <- brewer.pal(length(l),"Set1")
	cluster_2_colors <- c(brewer.pal(6,"Set2"),brewer.pal(9,"Set1"))[1:length(l)]
	names(cluster_2_colors) <- l			
	
	l <- levels(df_annot.df$gene_expression_clusters)
	cluster_colors_exp <- hue_pal()(length(l))
	names(cluster_colors_exp) <- l			
	
	df <- df_annot.df %>% 
		arrange(sil_score_order) %>%
		dplyr::select(cluster_id,cluster_id_pam,gene_expression_clusters) %>% 
		as.data.frame()
	rownames(df) <- df_annot.df$sample_id
	df_annot_cols = HeatmapAnnotation( df = df ,
																		 col = list( 
																		 	cluster_id = cluster_colors ,
																		 	cluster_id_pam = cluster_2_colors ,
																		 	gene_expression_clusters = cluster_colors_exp
																		 	# pas_clusters = cluster_colors
																		 )
	)	
	
	my_color_func <- circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
	my_color_func_viridis <- circlize::colorRamp2(seq(0,1,length.out = 10), viridis(10))
	my_color_func_inferno <- circlize::colorRamp2(seq(0,1,length.out = 10), inferno(10))
	my_color_func_viridis_mt_percent <- circlize::colorRamp2(seq(0,20,length.out = 10), viridis(10))
	df <- df_annot.df %>% 
		arrange(sil_score_order) %>%
		dplyr::select(cytotrace_score_ges,stemness_index,HALLMARK_MITOTIC_SPINDLE,HALLMARK_G2M_CHECKPOINT,HALLMARK_WNT_BETA_CATENIN_SIGNALING,HALLMARK_DNA_REPAIR,mt_percent) %>% 
		as.data.frame() %>%
		dplyr::rename(MITOTIC_SPINDLE=HALLMARK_MITOTIC_SPINDLE,
									G2M_CHECKPOINT=HALLMARK_G2M_CHECKPOINT,
									WNT_BETA_CATENIN_SIGNALING=HALLMARK_WNT_BETA_CATENIN_SIGNALING,
									DNA_REPAIR=HALLMARK_DNA_REPAIR)
	rownames(df) <- df_annot.df$sample_id
	df_annot_cols_suppl = HeatmapAnnotation( df = df ,
																					 col = list( 
																					 	cytotrace_score_ges = my_color_func_inferno ,
																					 	stemness_index = my_color_func_viridis ,
																					 	MITOTIC_SPINDLE = my_color_func ,
																					 	G2M_CHECKPOINT = my_color_func ,
																					 	WNT_BETA_CATENIN_SIGNALING = my_color_func ,
																					 	DNA_REPAIR = my_color_func ,
																					 	mt_percent = my_color_func_viridis_mt_percent ,
																					 	heatmap_legend_param = list(title_position = "topcenter", nrow = 1)
																					 )
	)			
	
	# mat_genes <- mat_genes
	index <- match( df_annot.df$sample_id[ sample_ordering_to_print ] , colnames(mat_genes) )
	mat_genes <- mat_genes[,index]			
	
	my_colors <- viridis(10)
	genes_hm <- Heatmap( mat_genes , 
											 # col = x[3:length(x)] ,
											 # col = viridis(10) ,
											 # col = rev( brewer.pal(10,"RdBu") ) ,
											 col = circlize::colorRamp2(c(0, 2, 4), c(my_colors[1],my_colors[6],my_colors[10]) ),
											 # col = viridis(10),
											 # heatmap_legend_param = list(color_bar = "continuous") ,
											 heatmap_legend_param = list(title_position = "topcenter", nrow = 1) ,
											 cluster_rows = FALSE ,
											 cluster_columns = FALSE ,
											 column_split = df_annot.df$cluster_id[ sample_ordering_to_print ] ,
											 column_gap = unit(2, "mm") ,
											 # row_order = order(n_gemm2pats,decreasing = TRUE) ,
											 # column_order = order(n_pats2gemm,decreasing = TRUE) ,
											 # col = inferno(10,direction = -1) ,
											 use_raster = isUsingRaster ,
											 name = "Protein Activity", 
											 # rect_gp = gpar(col = "white", lwd = 0.5) ,
											 # top_annotation = df_annot_cols , 
											 # left_annotation = df_annot_rows ,
											 # bottom_annotation = df_annot_cols_bottom ,
											 # right_annotation = df_annot_rows_right ,
											 show_row_names = TRUE , show_column_names = FALSE ,
											 row_names_gp = gpar(fontsize = 6) 
	)
	
	
	# sample_ordering_to_print <- order( clustering.tibble$cluster_reliability_order )
	# sample_ordering_to_print <- order( df_annot.df$sil_score_order )
	
	df_annot_cols_cluster_score <- HeatmapAnnotation( 
		silhouette_score = anno_barplot(
			# cluster_rel_score = anno_barplot(
			# clustering.tibble$cluster_reliability[sample_ordering_to_print] ,
			df_annot.df$sil_score[sample_ordering_to_print] ,
			bar_width = 1,
			gp = gpar(col = NA, 
								fill = rep( cluster_colors , table(df_annot.df$cluster_id) ) ),
			border = FALSE,
			# axis_param = list(at = c(0, 5e5, 1e6, 1.5e6),
			#     labels = c("0", "500k", "1m", "1.5m")),
			width = unit(2, "cm") )
	)	
	
	mat <- vpmat[ selected_features , ]
	index <- match( df_annot.df$sample_id[ sample_ordering_to_print ] , colnames(mat) )
	mat <- mat[,index]
	
	# x <- viridis(10)
	main_hm <- Heatmap( mat , 
											# col = x[3:length(x)] ,
											# col = viridis(10) ,
											# col = rev( brewer.pal(10,"RdBu") ) ,
											# col = circlize::colorRamp2(c(-3, 0, 3), c("Darkblue", "white", "red")),
											col = circlize::colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3")),
											# heatmap_legend_param = list(color_bar = "continuous") ,                       
											cluster_rows = FALSE ,
											cluster_columns = FALSE ,
											column_split = df_annot.df$cluster_id[ sample_ordering_to_print ] ,
											column_gap = unit(2, "mm") ,
											# row_order = order(n_gemm2pats,decreasing = TRUE) ,
											# column_order = order(n_pats2gemm,decreasing = TRUE) ,
											# col = inferno(10,direction = -1) ,
											use_raster = isUsingRaster ,
											name = "Protein Activity", 
											# rect_gp = gpar(col = "white", lwd = 0.5) ,
											# top_annotation = df_annot_cols , 
											# left_annotation = df_annot_rows ,
											# bottom_annotation = df_annot_cols_bottom ,
											# right_annotation = df_annot_rows_right ,
											show_row_names = TRUE , show_column_names = FALSE ,
											row_names_gp = gpar(fontsize = 6) ,
											heatmap_legend_param = list(title_position = "topcenter", nrow = 1)
	)
	
	ht_list = df_annot_cols_cluster_score %v% df_annot_cols %v% df_annot_cols_suppl %v% genes_hm %v% main_hm
	# ht_list = df_annot_cols %v% df_annot_cols_bottom %v% main_hm
	
	pdf( file.path( reports.dir , "TE001-pas-clustering-with-metacell.pdf") , width = 8, height = 12 )
	draw( ht_list, column_km = 1 , column_gap = unit(2, "mm") , 
				heatmap_legend_side = "bottom" )
	dev.off()
	
}	

## UMAPs over PAS and GES and mixed ----
print_msg_info(">>> UMAP over PAS and GES and mixed")
{
	library(ggplot2)
	library(ggrepel)
	library(umap)
	
	my_meta <- vp_seurat@meta.data %>% as_tibble()
	
	x <- as.matrix(vpdist)
	stopifnot( identical( colnames(x) , my_meta$cell_id ) )
	stopifnot( identical( labels(vpdist) , colnames(vp_seurat) ) )
	
	# create a new settings object with n_neighbors set to 5
	custom.settings <- umap.defaults
	custom.settings$random_state <- 666
	custom.settings$input <- "dist"
	# custom.settings$input <- "data"
	custom.settings$n_neighbors <- knn_n_neighbor
	custom.settings$n_components <- 2
	# custom.settings$metric <- "correlation"
	custom.settings$metric <- "euclidean"
	# custom.settings$metric <- "manhattan"
	custom.settings$verbose <- TRUE
	
	set.seed(666)
	# y.umap <- umap( t(x) , config = custom.settings )
	y.umap <- umap( as.matrix(x) , config = custom.settings )
	
	my_umap_coords <- y.umap$layout
	colnames(my_umap_coords) <- c("PAS_UMAP_1","PAS_UMAP_2")
	stopifnot( identical( rownames(y.umap$layout) , my_meta$cell_id ) )
	umap.df <- cbind(my_meta,my_umap_coords)
	
	
	umap.df$pas_cluster_id <- my_meta$seurat_clusters
	# umap.df$ges_cluster_id <- my_meta$clusters_ges.x
	umap.df$ges_cluster_id <- my_meta$clusters_ges
	
	stopifnot( identical( rownames(umap.df) , rownames(vp_seurat@meta.data) ) )
	vp_seurat@meta.data <- umap.df %>% as.data.frame()
	
	umap_plot <- ggplot( data = umap.df , aes( PAS_UMAP_1 , PAS_UMAP_2 , color = cytotrace_score_ges ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 0.75 , size = 0.75 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.5 , contour_var = "ndensity") +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		# scale_color_brewer(palette = "Set1") +
		scale_color_viridis(option = "B") +
		xlim( c( min(umap.df$PAS_UMAP_1)-5 , max(umap.df$PAS_UMAP_1)+5) ) +
		ylim( c( min(umap.df$PAS_UMAP_2)-5 , max(umap.df$PAS_UMAP_2)+5) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-cytotrace-at-ges.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()	
	
	# umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = stemness_index , shape = cluster_id ) ) +
	umap_plot <- ggplot( data = umap.df , aes( PAS_UMAP_1 , PAS_UMAP_2 , color = pas_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 0.75 , size = 0.75 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.5 , contour_var = "ndensity") +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		scale_color_brewer(palette = "Set1") +
		xlim( c( min(umap.df$PAS_UMAP_1)-5 , max(umap.df$PAS_UMAP_1)+5) ) +
		ylim( c( min(umap.df$PAS_UMAP_2)-5 , max(umap.df$PAS_UMAP_2)+5) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-pas-clusters.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()
	
}

# ## 2D Density plots of Markers ----
print_msg_info(">>> 2D Density plots of Markers")
{
	my_df <- clustering.tibble %>% dplyr::select(gene_Lgr5,gene_Tnfrsf19,gene_Krt19,gene_Ung,cytotrace_score_ges,cluster_id)
	
	my_df <- my_df %>% pivot_longer(cols = starts_with("gene_") , names_to = "gene" , values_to = "expression")
	
	# library(hrbrthemes)
	# hrbrthemes::import_roboto_condensed()
	p <- ggplot(data=my_df, aes(x=expression, group=cluster_id, fill=cluster_id)) +
		# scale_fill_brewer(palette="Accent",direction = +1) +
		scale_fill_brewer(palette="Set1",direction = +1) +
		geom_density(adjust=1.5, alpha=.6) +
		theme_light() +
		facet_wrap( ~gene+cluster_id , scales = "free_y" , ncol = 3)
		# facet_wrap( ~gene , scales = "free_y" , ncol = 3)
	# theme_ipsum()
	
	pdf( file.path(reports.dir,"markers-by-gene-and-cluster-density-plot-facet.pdf") )
		print(p)
	dev.off()
	
	p <- ggplot(data=my_df, aes(x=expression, group=cluster_id, fill=cluster_id)) +
		# scale_fill_brewer(palette="Accent",direction = +1) +
		scale_fill_brewer(palette="Set1",direction = +1) +
		geom_density(adjust=1.5, alpha=.4) +
		theme_light() +
		facet_wrap( ~gene , scales = "free_y" , ncol = 3)
		# facet_wrap( ~gene , scales = "free_y" , ncol = 3)
	# theme_ipsum()
	
	pdf( file.path(reports.dir,"markers-by-cluster-density-plot.pdf") )
		print(p)
	dev.off()	
	
	# --
	
	my_df <- clustering.tibble %>% dplyr::select(gene_Lgr5,gene_Tnfrsf19,gene_Krt19,gene_Ung,cytotrace_score_ges,stemness_index,cluster_id)

	p <- ggplot(data=my_df, aes(x=cytotrace_score_ges, group=cluster_id, fill=cluster_id)) +
		# scale_fill_brewer(palette="Accent",direction = +1) +
		scale_fill_brewer(palette="Set1",direction = +1) +
		geom_density(adjust=1.5, alpha=.4) +
		theme_light() +
		facet_grid( ~cluster_id ) +
		theme(axis.text.x = element_text(size = 5 , angle = 90, vjust = 0.5, hjust=1))
		# facet_wrap( ~gene , scales = "free_y" , ncol = 3)
	# theme_ipsum()
	
	pdf( file.path(reports.dir,"cytotrace-score-by-cluster-density-plot.pdf") , width = 3, height = 2 )
		print(p)
	dev.off()
	
	p <- ggplot(data=my_df, aes(x=stemness_index, group=cluster_id, fill=cluster_id)) +
		# scale_fill_brewer(palette="Accent",direction = +1) +
		scale_fill_brewer(palette="Set1",direction = +1) +
		geom_density(adjust=1.5, alpha=.4) +
		theme_light() +
		facet_grid( ~cluster_id ) +
		theme(axis.text.x = element_text(size = 5 , angle = 90, vjust = 0.5, hjust=1))
	# facet_wrap( ~gene , scales = "free_y" , ncol = 3)
	# theme_ipsum()
	
	pdf( file.path(reports.dir,"stemness-index-score-by-cluster-density-plot.pdf") , width = 3, height = 2 )
		print(p)
	dev.off()	

	# # TODO : Check if the scores with high CytoTrace score are the same with high Krt19 expr	
	# x <- my_df %>% 
	# 	dplyr::filter(gene == "gene_Krt19") %>%
	# 	dplyr::filter(expression > 1) %>%
	# 	dplyr::filter(gene == "gene_cytotrace")
	
	# TODO : This !
	# per cluster do correlation
	# plot( clustering.tibble$gene_Krt19 , clustering.tibble$cytotrace_score_ges )
}

## Pseudo-Time Plots of Markers ----
print_msg_info(">>> Pseudo-Time Plots of Markers")
{
	my_matrix <- vp_seurat@assays$VIPER@scale.data %>% as.matrix()
	my_df_pas <- t(my_matrix) %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") 
	my_metadata <- vp_seurat@meta.data %>% 
		dplyr::select(!starts_with("snn")) %>% 
		dplyr::select(!orig.ident ) %>% 
		dplyr::rename(cluster_id=seurat_clusters)
	my_metadata$cluster_id <- factor(as.numeric(my_metadata$cluster_id))
	my_metadata$pas_cluster_id <- factor(as.numeric(my_metadata$pas_cluster_id))
	dim(my_df_pas)
	dim(my_metadata)
	my_df_pas <- my_df_pas %>% pivot_longer( cols = !starts_with("cell_id") , names_to = "protein" , values_to = "activity" )
	my_df_pas <- left_join(my_df_pas,my_metadata)
	
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
	gx_seurat <- readRDS(seurat_analysis_data_filename)
	object_size(gx_seurat)
	
	print_msg_info(">>> Extracting gene expression measurements")
	my_matrix <- gx_seurat@assays$SCT@scale.data %>% as.matrix()
	my_df_gex <- t(my_matrix) %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") 
	my_df_gex <- my_df_gex %>% pivot_longer( cols = !starts_with("cell_id") , names_to = "gene" , values_to = "expression" )
	my_df_gex <- right_join(my_df_gex, 
												 my_metadata %>% dplyr::select(cell_id,cytotrace_score_ges,pas_cluster_id))
	
	my_protein <- "Lgr5"
	# df_4_plot <- my_df %>% dplyr::filter(protein == "Olfm4")
	df_4_plot <- my_df_pas %>% dplyr::filter(protein == my_protein)
	
	p <- ggplot( data = df_4_plot , aes( x = cytotrace_score_ges , y = activity , group = cluster_id , 
																			 color = activity , fill = cluster_id )) +
		# scale_fill_brewer(palette="Accent",direction = +1) +
		scale_x_reverse() +
		scale_color_gradient2( low = "blue" , mid = "white" , high = "red" ) +
		geom_point( alpha = .75 , size = 0.75 ) +
		geom_line( alpha = .75 , size = 0.5 ) +
		geom_smooth( alpha = .95 , method = 'loess' , span = 0.75 ) +
		facet_grid( ~cluster_id) +
		theme_light() +
		theme(axis.text.x = element_text(size = 8 , angle = 90, vjust = 0.5, hjust=1)) +
		theme(legend.position="bottom" , legend.box = "horizontal",
					legend.title = element_text(color = "black", size = 8) ,
					legend.text = element_text(color = "black", size = 8)
		)
	
	pdf( file.path(reports.dir, paste0("TE001-isc-cells-cytotrace-ordered-for-",my_protein,"-protein.pdf")) , width = 7 , height = 3 )
		print(p)
	dev.off()		
	
	my_gene <- "Lgr5"
	df_4_plot <- my_df_gex %>% dplyr::filter(gene == my_gene)
	
	p <- ggplot( data = df_4_plot , aes( x = cytotrace_score_ges , y = expression , group = pas_cluster_id , 
																			 color = expression , fill = pas_cluster_id )) +
		scale_x_reverse() +
		scale_fill_brewer(palette="Set1",direction = +1) +
		scale_color_gradient2( low = "blue" , mid = "white" , high = "red" ) +
		geom_point( alpha = .75 , size = 0.75 ) +
		geom_line( alpha = .75 , size = 0.25 , ) +
		geom_smooth( alpha = .5 , method = 'loess' , span = 0.75 ) +
		# facet_wrap( ~cluster_id, nrow = 3) +
		facet_grid( ~pas_cluster_id) +
		theme_light() +
		theme(axis.text.x = element_text(size = 8 , angle = 90, vjust = 0.5, hjust=1)) +
		theme(legend.position="bottom" , legend.box = "horizontal",
					legend.title = element_text(color = "black", size = 8) ,
					legend.text = element_text(color = "black", size = 8)
		)
	
	pdf( file.path(reports.dir, paste0("TE001-isc-cells-cytotrace-ordered-for-",my_gene,"-gene.pdf")) , width = 8 , height = 3 )
		print(p)
	dev.off()		
}

## Pseudo-Time Plots of Markers ----
print_msg_info(">>> Pseudo-Time Plots of Markers")
{
	require(Seurat)
	require(pryr)
	
	metadata_tibble_analysis_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-viper-analysis-with-metacell-metadata.rds") )
	vp_metadata <- readRDS(metadata_tibble_analysis_filename)
	object_size(vp_metadata)
	vp_metadata <- vp_metadata %>% dplyr::select(cell_id=sample_id,sil_score)
	
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
	gx_seurat <- readRDS(seurat_analysis_data_filename)
	object_size(gx_seurat)
	
	seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
	vp_seurat_full <- readRDS(seurat_viper_analysis_data_filename)
	object_size(vp_seurat)

	print_msg_info(">>> Extracting protein activity measurements")
	my_matrix <- vp_seurat_full@assays$VIPER@scale.data %>% as.matrix()
	my_df <- t(my_matrix) %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") 
	my_metadata <- vp_seurat_full@meta.data %>% 
		dplyr::select(!starts_with("snn")) %>% 
		dplyr::select(!orig.ident ) %>% 
		dplyr::rename(cluster_id=seurat_clusters)
	my_metadata$cluster_id <- factor(as.numeric(my_metadata$cluster_id))
	dim(my_df)
	dim(my_metadata)
	my_df <- my_df %>% pivot_longer( cols = !starts_with("cell_id") , names_to = "protein" , values_to = "activity" )
	my_df_pas <- left_join(my_df,my_metadata)
	my_df_pas <- left_join(my_df_pas,vp_metadata)
	
	print_msg_info(">>> Extracting gene expression measurements")
	my_matrix <- gx_seurat@assays$SCT@scale.data %>% as.matrix()
	my_df_gex <- t(my_matrix) %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") 
	my_df_gex <- my_df_gex %>% pivot_longer( cols = !starts_with("cell_id") , names_to = "gene" , values_to = "expression" )
	my_df_gex <- left_join(my_df_gex,my_metadata)
	my_df_gex <- left_join(my_df_gex,vp_metadata)
	
	# print_msg_info(">>> Joining Data")
	# my_df <- left_join(my_df,my_df_gex)
	
	# df_4_plot <- my_df %>% dplyr::filter(protein == "Olfm4")
	
	my_protein <- "Lgr5"
	df_4_plot <- my_df_pas %>% dplyr::filter(protein == my_protein)
	
	p <- ggplot( data = df_4_plot , aes( x = cytotrace_score_ges , y = activity , group = cluster_id , 
																			 color = activity , fill = cluster_id )) +
		scale_x_reverse() +
		scale_fill_brewer(palette="Set1",direction = +1) +
		scale_color_gradient2( low = "blue" , mid = "white" , high = "red" ) +
		geom_point( alpha = .75 , size = 0.75 ) +
		geom_line( alpha = .75 , size = 0.25 ) +
		geom_smooth( alpha = .5 , method = 'loess' , span = 0.75 ) +
		# facet_wrap( ~cluster_id, nrow = 3) +
		facet_grid( ~cluster_id) +
		theme_light() +
		theme(axis.text.x = element_text(size = 8 , angle = 90, vjust = 0.5, hjust=1)) +
		theme(legend.position="bottom" , legend.box = "horizontal",
					legend.title = element_text(color = "black", size = 8) ,
					legend.text = element_text(color = "black", size = 8)
		)
	
	pdf( file.path(reports.dir, paste0("TE001-all-cells-cytotrace-ordered-for-",my_protein,"-protein.pdf")) , width = 8 , height = 3 )
		print(p)
	dev.off()	
	
	p <- ggplot( data = df_4_plot , aes( x = sil_score , y = activity , group = cluster_id , 
																			 color = activity , fill = cluster_id )) +
		scale_x_reverse() +
		scale_fill_brewer(palette="Set1",direction = +1) +
		scale_color_gradient2( low = "blue" , mid = "white" , high = "red" ) +
		geom_point( alpha = .75 , size = 0.75 ) +
		geom_line( alpha = .75 , size = 0.5 ) +
		geom_smooth( alpha = .5 , method = 'loess' , span = 0.75 ) +
		facet_wrap( ~cluster_id, nrow = 3) +
		theme_light()
	
	pdf( file.path(reports.dir, paste0("TE001-all-cells-silscore-ordered-silscore-for-",my_protein,"-protein.pdf")) )
		print(p)
	dev.off()		
	
	my_gene <- "Lgr5"
	df_4_plot <- my_df_gex %>% dplyr::filter(gene == my_gene)
	
	p <- ggplot( data = df_4_plot , aes( x = cytotrace_score_ges , y = expression , group = cluster_id , 
																			 color = expression , fill = cluster_id )) +
		scale_x_reverse() +
		scale_fill_brewer(palette="Set1",direction = +1) +
		scale_color_gradient2( low = "blue" , mid = "white" , high = "red" ) +
		geom_point( alpha = .75 , size = 0.75 ) +
		geom_line( alpha = .75 , size = 0.25 , ) +
		geom_smooth( alpha = .5 , method = 'loess' , span = 0.75 ) +
		# facet_wrap( ~cluster_id, nrow = 3) +
		facet_grid( ~cluster_id) +
		theme_light() +
		theme(axis.text.x = element_text(size = 8 , angle = 90, vjust = 0.5, hjust=1)) +
		theme(legend.position="bottom" , legend.box = "horizontal",
					legend.title = element_text(color = "black", size = 8) ,
					legend.text = element_text(color = "black", size = 8)
					)
	
	pdf( file.path(reports.dir, paste0("TE001-all-cells-cytotrace-ordered-for-",my_gene,"-gene.pdf")) , width = 8 , height = 3 )
		print(p)
	dev.off()		
	
	p <- ggplot( data = df_4_plot , aes( x = sil_score , y = expression , group = cluster_id , 
																			 color = expression , fill = cluster_id )) +
		scale_x_reverse() +
		scale_fill_brewer(palette="Set1",direction = +1) +
		scale_color_gradient2( low = "blue" , mid = "white" , high = "red" ) +
		geom_point( alpha = .75 , size = 0.75 ) +
		geom_line( alpha = .75 , size = 0.5 ) +
		geom_smooth( alpha = .5 , method = 'loess' , span = 0.75 ) +
		facet_wrap( ~cluster_id, nrow = 3) +
		theme_light()
	
	pdf( file.path(reports.dir, paste0("TE001-all-cells-silscores-ordered-for-",my_gene,"-gene.pdf")) )
		print(p)
	dev.off()			
	
}
