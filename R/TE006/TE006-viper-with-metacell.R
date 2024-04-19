## 
# Single cell Analysis of TE006 (VIPER Analysis)
# ----------------------------------------------
# New approach, metacells were created based on gene expression clustering
# ------------------------------------------------------------------------
# system.time({ source("sources/ermanno-data-analysis/TE006/TE006-viper-with-metacell.R") })

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

print_msg_info(">>> Setting Analysis Params")
{
	isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE006/")
	my_sample_id <- "TE006"
	my_seed <- 42
	knn_n_neighbor <- 15
	isComputingOptimalNumberOfClusters <- FALSE	
	regulon_n_targets <- 50
	ges_n_feature <- 15000
	metaviper_mvws <- "ClustSpecific"
}

## Loading counts matrix ----
print_msg_info(">>> Loading counts matrix")
{
	# x <- readRDS( file.path(isc_data_dir,"/TE006-seurat-data.rds") )
	# my_counts <- x@assays$SCT@data
	cpm_filename <- file.path(isc_data_dir,"TE006-cpm.rds")
	file.info(cpm_filename)["mtime"]
	x <- readRDS( cpm_filename )
	my_counts <- x
	dim(my_counts)
}

## Loading counts matrix ----
print_msg_info(">>> Setting constants and paramters [in YAML]")
{
	my_params.list <- list()
	my_params.list$sample_id = my_sample_id
	# my_params.list$seurat_original_data_filename <- seurat_original_data_filename
	# my_params.list$seurat_analysis_data_filename <- seurat_analysis_data_filename
	
	my_params.list$seed = my_seed
	my_params.list$knn_n = knn_n_neighbor
	
	my_params.list$is_SCT_do.center = TRUE 
	my_params.list$is_SCT_do.scale = TRUE
	my_params.list$ges_n_feature <- ges_n_feature
	my_params.list$viper_signature_method <- "none"
	# my_params.list$viper_signature_method <- "mad"
	my_params.list$regulon_minsize <- 50
	my_params.list$viper_cluster_specific <- TRUE
	my_params.list$viper_mvws_other_regs <- "allCellsNet"
	my_params.list$viper_dist_corr_method <- "pearson"
	# my_params.list$viper_dist_corr_method <- "spearman"
	
	# my_params.list$graph_type <- "nn"
	my_params.list$graph_type <- "snn"
	my_params.list$is_only_TFS <- TRUE
	# my_params.list$sct_variable_n_feats <- sct_variable_n_feats
	# my_params.list$sct_is_to_scale <- FALSE
	# my_params.list$n_pcs <- n_pcs
	my_params.list$pcs_to_use <- 1:10
	# my_params.list$clustering_algorithm <- clustering_algorithm
	
	tmp <- paste( my_params.list , collapse = "-" )
	create_workspace(tmp)
	
	my_params.list$results_dir = reports.dir
	my_params.list$cpm_filename <- cpm_filename	
	
	my_yaml_file <- file.path(reports.dir,"run-parameters.yaml")
	yaml::write_yaml(my_params.list %>% as.list(),my_yaml_file)
	
}

## VIPER Run ----
print_msg_info(">>> VIPER Analysis")
{
	require(Seurat)
	require(viper)

	seurat_original_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-original-data.rds") )
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
	
	seurat_analysis <- readRDS(seurat_analysis_data_filename)
	object_size(seurat_analysis)
	file.info(seurat_analysis_data_filename)

	seurat_analysis$cell_id <- colnames(seurat_analysis)
			
	# my_ges.rankNormed <- rankNorm(my_counts)
	
	# print_msg_info(">>> Gene Expression computerd with ranknorm")
	# {
	# my_counts <- as.matrix(seurat_analysis@assays$RNA@counts)
	# my_counts <- as.matrix(RelativeCounts(my_counts,1e6))
	# my_ges.rankNormed <- rankNorm(my_counts)
	# }
	
	print_msg_info(">>> Gene Expression computerd with SCT@scale.data")
	{
		my_gex_seurat <- readRDS(seurat_original_data_filename)
		my_gex_seurat <- my_gex_seurat %>% 
			SCTransform( do.scale = my_params.list$is_SCT_do.scale , 
									 do.center = my_params.list$is_SCT_do.center , 
									 return.only.var.genes = TRUE ,
									 variable.features.n = my_params.list$ges_n_feature
			) %>% 			
			# SCTransform( variable.features.n = 3000 , 
			# 						 return.only.var.genes = FALSE ,
			# 						 do.scale = TRUE , 
			# 						 do.center = TRUE ,
			# 						 seed.use = my_seed ) %>% 
			RunPCA() 
		my_ges.rankNormed <- as.matrix(my_gex_seurat@assays$SCT@scale.data)
		# rm(x)		
	}	

	x <- read_csv( file.path("~/Clouds/Dropbox/Data/isc/TE001/lists/co-tfs-list-v3.csv") , col_names = FALSE )
	y <- read_csv( file.path("~/Clouds/Dropbox/Data/isc/TE001/lists/tfs-list-v3.csv") , col_names = FALSE )
	z1 <- read_csv( file.path("~/Clouds/Dropbox/Data/isc/TE001/lists/reg-extra-list-v3.csv") , col_names = FALSE )
	z2 <- read_csv( file.path("~/Clouds/Dropbox/Data/isc/TE001/lists/reg-membr-list-v3.csv") , col_names = FALSE )
	if( my_params.list$is_only_TFS )
	{
		tfs_and_cotfs <- unique(c( x$X1 , y$X1 , c("Rspo","Olfm4","Krt19")))	
	} else {
		tfs_and_cotfs <- unique(c( x$X1 , y$X1 , z1$X1 , z2$X1 ))
	}
	print(length(tfs_and_cotfs))

	# TE006_network_full <- readRDS(file.path(isc_data_dir,"networks/TE006-full_unPruned.rds"))
	TE001_network_full <- readRDS(file.path("~/Clouds/Dropbox/Data/isc/TE001/networks/TE001-full_unPruned.rds"))
	TE001_network_full <- pruneRegulon(TE001_network_full,50)
	vpmat <- viper( my_ges.rankNormed , TE001_network_full , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
	
	## ----
	# TE006_network_c1 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE006_c1_unPruned.rds"))
	# TE006_network_c1 <- pruneRegulon(TE006_network_c1,50)
	# length(TE006_network_c1)
	# TE006_network_c1 <- TE006_network_c1[ names(TE006_network_c1) %in% tfs_and_cotfs ]
	# length(TE006_network_c1)
	# 
	# index <- colnames(my_ges.rankNormed) %in% seurat_analysis$cell_id[ (seurat_analysis$seurat_clusters == 0) ]
	# length(index)
	# sum(index)
	# c1.viper <- viper( my_ges.rankNormed[,index] , TE006_network_c1 , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
	# # c1.viper <- viper( rankNorm(my_counts[,index]) , TE006_network_c1 , method = "none" , minsize = 30 , verbose = TRUE)
	# 
	# TE006_network_c2 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE006_c2_unPruned.rds"))
	# TE006_network_c2 <- pruneRegulon(TE006_network_c2,50)
	# length(TE006_network_c2)
	# TE006_network_c2 <- TE006_network_c2[ names(TE006_network_c2) %in% tfs_and_cotfs ]
	# length(TE006_network_c2)
	# 
	# index <- colnames(my_ges.rankNormed) %in% seurat_analysis$cell_id[ (seurat_analysis$seurat_clusters == 1) ]
	# length(index)
	# sum(index)
	# c2.viper <- viper( my_ges.rankNormed[,index] , TE006_network_c2 , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
	# # c2.viper <- viper( rankNorm(my_counts[,index]) , TE006_network_c2 , method = "none" , minsize = 30 , verbose = TRUE)
	# 
	# TE006_network_c3 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE006_c3_unPruned.rds"))
	# TE006_network_c3 <- pruneRegulon(TE006_network_c3,50)
	# length(TE006_network_c3)
	# TE006_network_c3 <- TE006_network_c3[ names(TE006_network_c3) %in% tfs_and_cotfs ]
	# length(TE006_network_c3)
	# 
	# index <- colnames(my_ges.rankNormed) %in% seurat_analysis$cell_id[ (seurat_analysis$seurat_clusters == 2) ]
	# length(index)
	# sum(index)
	# c3.viper <- viper( my_ges.rankNormed[,index] , TE006_network_c3 , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
	# # c3.viper <- viper( rankNorm(my_counts[,index]) , TE006_network_c3 , method = "none" , minsize = 30 , verbose = TRUE)
	# 
	# 
	# TE006_network_c3and4 <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE006_c3and4_unPruned.rds"))
	# TE006_network_c3and4 <- pruneRegulon(TE006_network_c3and4,50)
	# length(TE006_network_c3and4)
	# TE006_network_c3and4 <- TE006_network_c3and4[ names(TE006_network_c3and4) %in% tfs_and_cotfs ]
	# length(TE006_network_c3and4)
	# 
	# index <- colnames(my_ges.rankNormed) %in% seurat_analysis$cell_id[ (seurat_analysis$seurat_clusters == 3) ]
	# length(index)
	# sum(index)
	# c3and4.viper <- viper( my_ges.rankNormed[,index] , TE006_network_c3and4 , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
	# # c4.viper <- viper( rankNorm(my_counts[,index]) , TE006_network_c4 , method = "none" , minsize = 30 , verbose = TRUE)
	# 
	# # index <- colnames(my_ges.rankNormed) %in% seurat_analysis$cell_id[ !(seurat_analysis$seurat_clusters %in% c(1,2,3,4)) ]
	# # length(index)
	# # sum(index)
	# # system.time({
	# # 	c5_c6.viper <- viper(
	# # 		eset = my_ges.rankNormed[,index] ,
	# # 		# eset = rankNorm(my_counts[,index]) ,
	# # 									# regulon = list(TE006_network,
	# # 									regulon = list(TE006_network_c1,TE006_network_c2,TE006_network_c3,TE006_network_c4) ,
	# # 									mvws = 3 , cores = 2 ,
	# # 									# regulon = TE006_network ,
	# # 									method = "none" , minsize = 30 , verbose = TRUE )
	# # })
	# 
	# # regs_in_common <- Reduce( "intersect" , list( rownames(c1.viper) , rownames(c2.viper) , rownames(c3.viper) , rownames(c4.viper) , rownames(c5_c6.viper) ) )
	# regs_in_common <- Reduce( "intersect" , list( rownames(c1.viper) , rownames(c2.viper) , rownames(c3.viper) , rownames(c3and4.viper) ) )
	# length(regs_in_common)
	# 
	# vpmat <- do.call(cbind, list( c1.viper[regs_in_common,] , c2.viper[regs_in_common,] , c3.viper[regs_in_common,] , c3and4.viper[regs_in_common,] ) )
	# str(vpmat,1)

	vpmat <- vpmat[ , match( colnames(my_counts) , colnames(vpmat) ) ]
	stopifnot( identical( colnames(my_counts) , colnames(vpmat) ) )
	
	# print_msg_info(">>> Running metaVIPER on the reminder of regulons")
	# {
	# 	extract_complement_regulons <- function(net,regs) {
	# 		index <- names(net) %in% regs
	# 		index <- !index
	# 		print_msg_info("*** Found " , sum(index) , " regulons" )
	# 		net <- net[index]
	# 		return(net)
	# 	}
	# 	TE006_network_c1.excl <- extract_complement_regulons(TE006_network_c1,regs_in_common)
	# 	length(TE006_network_c1.excl)
	# 	TE006_network_c2.excl <- extract_complement_regulons(TE006_network_c2,regs_in_common)
	# 	length(TE006_network_c2.excl)
	# 	TE006_network_c3.excl <- extract_complement_regulons(TE006_network_c3,regs_in_common)
	# 	length(TE006_network_c3.excl)
	# 	TE006_network_c3and4.excl <- extract_complement_regulons(TE006_network_c3and4,regs_in_common)
	# 	length(TE006_network_c3and4.excl)
	# 
	# 	system.time({
	# 		vpmat.excl <- viper( eset = my_ges.rankNormed ,
	# 										# regulon = list(TE006_network,
	# 										regulon = list(TE006_network_c1.excl,TE006_network_c2.excl,TE006_network_c3.excl,TE006_network_c3and4.excl) ,
	# 										mvws = my_params.list$viper_mvws_other_regs , cores = 2 ,
	# 										# regulon = TE006_network ,
	# 										method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE )
	# 	})
	# 
	# 	dim(vpmat.excl)
	# 	dim(vpmat)
	# 
	# 	# setdiff( rownames(vpmat.excl) , rownames(vpmat) )
	# 	intersect( rownames(vpmat.excl) , rownames(vpmat) )
	# 	stopifnot(identical(colnames(vpmat),colnames(vpmat.excl)))
	# 
	# 	print_msg_warn(">>> Original vpmat with cluster-specific VIPER run, feats: " , nrow(vpmat))
	# 	vpmat <- rbind(vpmat,vpmat.excl)
	# 	print_msg_warn(">>> Integrated with metaVIPER vpmat with cluster-specific VIPER run, feats: " , nrow(vpmat))
	# 
	# }
	
	vpsim <- viperSimilarity(vpmat)
	# vpdist <- as.dist( vpsim )
	vpdist <- as.dist( 1-cor(
		# vpmat[order(rowVars(vpmat),decreasing = T)[1:2000],] ,
		vpmat ,
		method = my_params.list$viper_dist_corr_method )
	)

	# saveRDS(object = list(gesmat=my_ges.rankNormed,vpmat=vpmat,vpsim=vpsim,vpdist=vpdist) ,
	# 				file.path(isc_data_dir,"/TE006-viper-analysis-with-metacell.rds") )

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

# ## UMAP over PAS ----
# print_msg_info(">>> UMAP over PAS")
# {
# 	library(ggplot2)
# 	library(ggrepel)
# 	library(umap)
# 	# library(tidyverse)
# 	
# 	stopifnot( identical( labels(vpdist),vpmat_metadata$cell_id ) )
# 	
# 	x <- vpdist
# 	# x <- t(vpmat[order(rowVars(vpmat),decreasing = T)[1:500],])
# 	
# 	# create a new settings object with n_neighbors set to 5
# 	custom.settings <- umap.defaults
# 	custom.settings$random_state <- 666
# 	custom.settings$input <- "dist"
# 	# custom.settings$input <- "data"
# 	custom.settings$n_neighbors <- 15
# 	custom.settings$n_components <- 2
# 	# custom.settings$metric <- "correlation"
# 	# custom.settings$metric <- "euclidean"
# 	# custom.settings$metric <- "manhattan"
# 	custom.settings$verbose <- TRUE
# 	
# 	set.seed(666)
# 	# y.umap <- umap( t(x) , config = custom.settings )
# 	y.umap <- umap( as.matrix(x) , config = custom.settings )
# 	
# 	umap.df <- as.data.frame( y.umap$layout )
# 	# umap.df$cell_line = ifelse( grepl( "LNCaP" , colnames(vpmat) ) , "clones" ," single_cells" )
# 	# 
# 	# umap.df$net_enrichment <- getNETenrichment(vpmat,net.regulon.obj)
# 	# umap.df$beltran_net_enrichment <- getNETenrichment(vpmat,beltran.net.regulon.obj)
# 	# umap.df$FOXM1 <- vpmat["Foxm1",]
# 	umap.df$stemness_index <- vpmat_metadata$stemness_index
# 	umap.df$cluster_id <- vpmat_metadata$cluster_id_pam
# 	umap.df$sil_score <- vpmat_metadata$sil_score
# 	
# 	# # umap.df <- cbind( umap.df , t(vpmat_pathways ))
# 	# umap.df$Lgr4 <- vpmat["Lgr4",]
# 	# umap.df$Lgr5 <- vpmat["Lgr5",]
# 	# umap.df$Olfm4 <- vpmat["Olfm4",]
# 	umap.df$cell_id <- vpmat_metadata$cell_id
# 	# umap.df$Atoh1 <- vpmat["Atoh1",]
# 	# umap.df$Ascl2 <- vpmat["Ascl2",]
# 	
# 	# library(mltools)
# 	
# 	# tmp <- model.matrix(~umap.df$model-1)
# 	# colnames(tmp) <- gsub("(umap.df\\$model)(.*)","\\2",colnames(tmp))
# 	# umap.df <- cbind(umap.df,tmp)
# 	# 
# 	# umap.df$NPK <- as.factor( umap.df$NPK )
# 	
# 	umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = stemness_index , shape = cluster_id ) ) +
# 		# umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = NPK , shape = castration_status ) ) +
# 		# geom_text_repel( aes( label = model ) , size = 5 ) +
# 		# geom_text_repel( aes( label = ID ) , size = 5 ) +
# 		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
# 		geom_point( alpha = 0.5 , stroke = 0.5 , size = 5 ) +
# 		scale_color_viridis() +
# 		geom_density_2d(lineend = "butt", color = "black",size=0.5,alpha=0.5) +
# 		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
# 		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
# 		# scale_colour_manual( name = "survival_group" , values = c("green","blue","red")) +
# 		xlab( paste0( "UMAP 1") ) +
# 		ylab( paste0( "UMAP 2") ) +
# 		theme_minimal()
# 	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
# 	# theme( legend.position = "none" )
# 	# coord_fixed()
# 	
# 	pdf( file = file.path( reports.dir , "TE006-pas-umap-analysis.pdf") , width = 16 , height = 12 )
# 		print(umap_plot)
# 	dev.off()
# }

# ## Saving Metadata Preparation  ----
# print_msg_info(">>> Saving Processed Data in CSV (for Scanpy)")
# {
# 	stopifnot( identical( rownames(umap.df) , colnames(vpmat) ) )
# 	stopifnot( identical( rownames(umap.df) , names(stemness_index) ) )
# 	
# 	write_vpmat_to_3_csv( vpmat , "~/Clouds/Dropbox/Data/isc/TE006/" , filetag = paste0(my_sample_id,"-with-metacell") , is_ensembl_to_sym = FALSE )
# 	write_csv( vpmat_metadata , "~/Clouds/Dropbox/Data/isc/TE006/TE006-with-metacell-metadata.csv" )
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
	pcs_to_use <- my_params.list$pcs_to_use
	
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
	my_sil.df <- tibble( resolution = seq(0.01,0.5,by = 0.01) , sil_avg = 0 , sil_median = 0 , n_clust = 0 )
	for ( a_res in my_sil.df$resolution )
	{
		message("--- Silhouette score computation for resolution :" , a_res )
		vp_seurat <- FindClusters( vp_seurat , 
															 # graph.name = "VIPER_snn" ,
															 graph.name = my_params.list$graph_type ,
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
	# print_msg_warn(">>> Using MANUAL resolution")
	# my_res <- 0.38
	# my_res <- 0.25
	
	# vp_seurat <- FindClusters(vp_seurat, resolution = 2 , verbose = TRUE,modularity.fxn = 2,algorithm = 4,random.seed = 666 ) 
		vp_seurat <- FindClusters( vp_seurat , 
															 # graph.name = "VIPER_snn" ,
															 graph.name = my_params.list$graph_type ,
															 resolution = my_res , 
															 verbose = FALSE ,
															 modularity.fxn = 1 ,
															 # algorithm = 4 , # Leiden
															 algorithm = 1 , # Louvain
															 seed.use = my_seed
		)
	
	print(table(vp_seurat$seurat_clusters))
	
	# print_msg_warn("\t *** Not Running RunUMAP from Seurat because continuously issues on umap-learn not recognized ***")
	vp_seurat <- RunUMAP( vp_seurat ,
												assay = "VIPER" ,
												n.components = 3L ,
												# features = prolif_genes ,
												# assay = "RNA" ,
												graph = my_params.list$graph_type ,
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
	
	print_msg_info(">>> Plotting UMAP at protein activity with sample_id clusters")
	{
		x <- vp_seurat@reductions$umap@cell.embeddings
		vp_seurat@meta.data$UMAP_1.PAS = x[,1]
		vp_seurat@meta.data$UMAP_2.PAS = x[,2]
		vp_seurat@meta.data$UMAP_3.PAS = x[,3]
		
		my_df <- vp_seurat@meta.data
		my_df$sample_id <- "X"
		my_cells_ids <- gsub( "(.*)_(.*)" , "\\1" , rownames(my_df) )
		# my_df$sample_id <- ifelse( my_cells_ids %in% colnames(pre_rad) , "MJ005" , my_df$sample_id )
		# my_df$sample_id <- ifelse( my_cells_ids %in% colnames(post_rad) , "MJ008" , my_df$sample_id )
		vp_seurat@meta.data$sample_id <- my_df$sample_id
		
		umap_plot <- ggplot( data = my_df , aes( UMAP_1.PAS , UMAP_2.PAS , color = sample_id , shape = sample_id ) ) +
			geom_point( alpha = 0.75 , stroke = 1 , size = 1 ) +
			# scale_color_manual(values = my_palette ) +
			xlab( paste0( "UMAP 1") ) +
			ylab( paste0( "UMAP 2") ) +
			theme_minimal()
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-pas-umap.pdf")) , width = 10 , height = 8 )
			print(umap_plot)
		dev.off()					
	}			
	
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
												"Mki67","Fzd5","Fzd7","Dclk1","Rspo3","Wnt4","Wnt5a","Uri1" ,
												"GFP","DTR")
			
			x <- readRDS( file.path(isc_data_dir,"/TE006-seurat-analysis-data.rds") )
			
			tmp <- x@meta.data %>% dplyr::select(cell_id,UMAP_1.GES,UMAP_2.GES,clusters_ges=seurat_clusters,stemness_index.ges,singleR_labels)
			vp_seurat@meta.data <- left_join(vp_seurat@meta.data,tmp,by=c("cell_id"="cell_id"))	%>% as.data.frame()
			
			rownames(vp_seurat@meta.data) <- colnames(vp_seurat)
			
			# x <- x@assays$RNA@counts %>% RelativeCounts(1e4) %>% as.matrix()
			# x <- rankNorm(x)
			
			x <- my_gex_seurat@assays$SCT@scale.data %>% as.matrix()
			# x <- as.matrix(x@assays$SCT@scale.data)
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
			x <- readRDS( file.path(isc_data_dir,"/TE006-seurat-analysis-data.rds") )
			gene_expression_clusters <- x$seurat_clusters
			my_mt_percent <- x@meta.data$mt_percent
			names(my_mt_percent) <- rownames(x@meta.data)
			clustering.tibble$gene_expression_clusters <- gene_expression_clusters[ match( clustering.tibble$sample_id , names(gene_expression_clusters) ) ]
			clustering.tibble$mt_percent <- my_mt_percent[ match( clustering.tibble$sample_id , names(my_mt_percent) ) ]
		}	
		
		# print_msg_info(">>> >> Integration with CytoTRACE stemness inference")
		# {
		# 	df <- read_csv(file.path(isc_data_dir, paste0(my_sample_id,"-cytotrace-table.csv")) )
		# 
		# 	index <- match( clustering.tibble$sample_id , df$cell_id )
		# 	na_tot <- sum(is.na(index))
		# 	print_msg_warn(">>> >> >> We are matching the cell labels, but there are " , na_tot , " NAs values")	
		# 	# index <- index[ !is.na(index) ]
		# 	df <- df[ index , ]
		# 	# stopifnot( identical( df$cell_id , clustering.tibble$sample_id ) )
		# 	
		# 	clustering.tibble$cytotrace_score_ges <- df$ct_score		
		# 	
		# 	stopifnot( identical( rownames(vp_seurat@meta.data) , clustering.tibble$sample_id ) )
		# 	vp_seurat@meta.data$cytotrace_score_ges <- df$ct_score		
		# }
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
			index <- weigths > 0.1
			mat <- mat[,index]
			clusters <- clusters[index]
			weigths <- weigths[index]
		}
		
		vpmat.stouffered <- doStoufferFromClusters( mat , clusters = clusters , weights = weigths )
		# vpmat.stouffered <- doStoufferFromClusters( mat , clusters = clusters )
		top_markers.pas <- CBCMRs(vpmat.stouffered,5)
		
		colnames(vpmat.stouffered) <- paste0("Cluster_", levels(clusters) )
		write_csv( vpmat.stouffered %>% as.data.frame() %>% rownames_to_column("proteins") , 
							 file.path(isc_data_dir , paste0(my_sample_id,"-viper-stouffer-with-metacell.csv") ) 
							 	)
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

## CytoTrace Analysis ----
print_msg_info(">>> CytoTrace Analysis")
{
	source("~/Downloads/CytoTRACE/R/zzz.R")
	source("~/Downloads/CytoTRACE/R/CytoTRACE.R")
	source("~/Downloads/CytoTRACE/R/plotCytoGenes.R")
	source("~/Downloads/CytoTRACE/R/plotCytoTRACE.R")
	
	mat <- as.matrix(my_gex_seurat@assays$RNA@counts)
	dim(mat)
	print_msg_info(">>> >> Running CytoTrace  ...")
	cytotrace.data <- CytoTRACE(mat , enableFast = FALSE , ncores = 2 )
	require(pryr)
	object_size(cytotrace.data)
	plotCytoGenes( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-ges-data")) ,numOfGenes = 10)
	
	vp_seurat$cytotrace_score_ges <- cytotrace.data$CytoTRACE[ match( colnames(vp_seurat) , names(cytotrace.data$CytoTRACE) ) ]
	clustering.tibble$cytotrace_score_ges <- cytotrace.data$CytoTRACE[ match( clustering.tibble$sample_id , names(cytotrace.data$CytoTRACE) ) ]
	# vp_seurat$cytotrace_score_ges <- cytotrace.data$CytoTRACE[ match( vp_seurat$cell_id , names(cytotrace.data$CytoTRACE) ) ]
	
	
	
	# df <- cbind(clustering.tibble$,clustering.tibble$UMAP_2.GES)
	# my_clusters <- as.character(my_isc.sdata$seurat_clusters)
	# names(my_clusters) <- names(my_isc.sdata$seurat_clusters)
	# plotCytoTRACE( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-ges-clusters")) , emb = df , phenotype = my_clusters )
	
}				

# Heatmap Printing on VIPER Data ----
print_msg_info(">>> Rendering clustering solution with ComplexHeatmap")
{
	require(ComplexHeatmap)
	require(circlize)
	require(viridis)
	require(dplyr)
	
	isUsingRaster <- TRUE
	
	n_top <- 10
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
	
	# l <- levels(df_annot.df$cluster_id_pam)
	# # cluster_colors <- brewer.pal(length(l),"Set1")
	# cluster_2_colors <- c(brewer.pal(6,"Set2"),brewer.pal(9,"Set1"))[1:length(l)]
	# names(cluster_2_colors) <- l			
	
	l <- levels(df_annot.df$gene_expression_clusters)
	cluster_colors_exp <- hue_pal()(length(l))
	names(cluster_colors_exp) <- l			
	
	df <- df_annot.df %>% 
		arrange(sil_score_order) %>%
		# dplyr::select(cluster_id,cluster_id_pam,gene_expression_clusters) %>%
			dplyr::select(cluster_id,gene_expression_clusters) %>%
		as.data.frame()
	rownames(df) <- df_annot.df$sample_id
	df_annot_cols = HeatmapAnnotation( df = df ,
																		 col = list( 
																		 	cluster_id = cluster_colors ,
																		 	# cluster_id_pam = cluster_2_colors ,
																		 	gene_expression_clusters = cluster_colors_exp
																		 )
	)	
	
	my_color_func <- circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
	my_color_func_viridis <- circlize::colorRamp2(seq(0,1,length.out = 10), viridis(10))
	my_color_func_inferno <- circlize::colorRamp2(seq(0,1,length.out = 10), inferno(10))
	my_color_func_viridis_mt_percent <- circlize::colorRamp2(seq(0,20,length.out = 10), viridis(10))
	df <- df_annot.df %>% 
		arrange(sil_score_order) %>%
		dplyr::select(cytotrace_score_ges,stemness_index,HALLMARK_MITOTIC_SPINDLE,HALLMARK_G2M_CHECKPOINT,HALLMARK_WNT_BETA_CATENIN_SIGNALING,HALLMARK_DNA_REPAIR,mt_percent) %>% 
		as.data.frame()
	rownames(df) <- df_annot.df$sample_id
	df_annot_cols_suppl = HeatmapAnnotation( df = df ,
																					 col = list( 
																					 	cytotrace_score_ges = my_color_func_inferno ,
																					 	stemness_index = my_color_func_viridis ,
																					 	HALLMARK_MITOTIC_SPINDLE = my_color_func ,
																					 	HALLMARK_G2M_CHECKPOINT = my_color_func ,
																					 	HALLMARK_WNT_BETA_CATENIN_SIGNALING = my_color_func ,
																					 	HALLMARK_DNA_REPAIR = my_color_func ,
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
											 name = "Gene Expression", 
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
	
	pdf( file.path( reports.dir , "TE006-pas-clustering-with-metacell.pdf") , width = 8, height = 11 )
		draw( ht_list, column_km = 1 , column_gap = unit(2, "mm") , 
					heatmap_legend_side = "bottom" )
	dev.off()
	
}	

## UMAPs over PAS and GES and mixed ----
print_msg_info(">>> UMAP over PAS and GES and mixed")
{
	print_msg_warn("*** ReRunning UMAP ***")
	library(ggplot2)
	library(ggrepel)
	library(umap)
	# library(tidyverse)
	
	my_meta <- vp_seurat@meta.data %>% as_tibble()
	my_meta <- my_meta %>% dplyr::rename(UMAP_1.PAS.with_seurat=UMAP_1.PAS,
																			 UMAP_2.PAS.with_seurat=UMAP_2.PAS,
																			 UMAP_3.PAS.with_seurat=UMAP_3.PAS)
	
	x <- as.matrix(vpdist)
	stopifnot( identical( colnames(x) , my_meta$cell_id ) )
	stopifnot( identical( labels(vpdist) , colnames(vp_seurat) ) )
	
	# x <- t(vpmat[order(rowVars(vpmat),decreasing = T)[1:500],])
	
	# create a new settings object with n_neighbors set to 5
	custom.settings <- umap.defaults
	custom.settings$random_state <- 666
	custom.settings$input <- "dist"
	# custom.settings$input <- "data"
	custom.settings$n_neighbors <- knn_n_neighbor
	custom.settings$n_components <- 3
	# custom.settings$metric <- "correlation"
	custom.settings$metric <- "euclidean"
	# custom.settings$metric <- "manhattan"
	custom.settings$verbose <- TRUE
	
	set.seed(666)
	# y.umap <- umap( t(x) , config = custom.settings )
	y.umap <- umap( as.matrix(x) , config = custom.settings )
	
	my_umap_coords <- y.umap$layout
	# colnames(my_umap_coords) <- c("PAS_UMAP_1","PAS_UMAP_2","PAS_UMAP_3")
	colnames(my_umap_coords) <- c("UMAP_1.PAS","UMAP_2.PAS","UMAP_3.PAS")
	print_msg_warn("*** RENAMING UMAP coordinates ***")
	stopifnot( identical( rownames(y.umap$layout) , my_meta$cell_id ) )
	umap.df <- cbind(my_meta,my_umap_coords)
	
	# umap.df <- as.data.frame( y.umap$layout )

	# umap.df$stemness_index <- vpmat_metadata$stemness_index
	# umap.df$cluster_id <- vpmat_metadata$cluster_id_pam
	# umap.df$sil_score <- vpmat_metadata$sil_score
	
	umap.df$pas_cluster_id <- my_meta$seurat_clusters
	# umap.df$ges_cluster_id <- my_meta$clusters_ges.x
	umap.df$ges_cluster_id <- my_meta$clusters_ges
	
	# umap.df$cell_id <- vpmat_metadata$cell_id
	# umap.df$Rspo3 <- as.matrix(vp_seurat@assays$VIPER@data)["Rspo3",]
	# umap.df$Ascl2 <- vpmat["Ascl2",]
	
	# library(mltools)
	
	# tmp <- model.matrix(~umap.df$model-1)
	# colnames(tmp) <- gsub("(umap.df\\$model)(.*)","\\2",colnames(tmp))
	# umap.df <- cbind(umap.df,tmp)
	# 
	# umap.df$NPK <- as.factor( umap.df$NPK )
	
	stopifnot( identical( rownames(umap.df) , rownames(vp_seurat@meta.data) ) )
	vp_seurat@meta.data <- umap.df %>% as.data.frame()
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = cytotrace_score_ges ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.5 , contour_var = "ndensity") +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		# scale_color_brewer(palette = "Set1") +
		scale_color_viridis(option = "B") +
		xlim( c( min(umap.df$UMAP_1.PAS)-5 , max(umap.df$UMAP_1.PAS)+5) ) +
		ylim( c( min(umap.df$UMAP_2.PAS)-5 , max(umap.df$UMAP_2.PAS)+5) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE006-pas-umap-with-cytotrace-at-ges.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()	
	
	# umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = stemness_index , shape = cluster_id ) ) +
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = pas_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.5 , contour_var = "ndensity") +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		scale_color_brewer(palette = "Set1") +
		xlim( c( min(umap.df$UMAP_1.PAS)-5 , max(umap.df$UMAP_1.PAS)+5) ) +
		ylim( c( min(umap.df$UMAP_2.PAS)-5 , max(umap.df$UMAP_2.PAS)+5) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE006-pas-umap-with-pas-clusters.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = singleR_labels ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 0.5 , size = 0.5 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.5 , contour_var = "ndensity") +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		# scale_color_brewer(palette = "Set1") +
		# scale_color_viridis(option = "B") +
		xlim( c( min(umap.df$UMAP_1.PAS)-5 , max(umap.df$UMAP_1.PAS)+5) ) +
		ylim( c( min(umap.df$UMAP_2.PAS)-5 , max(umap.df$UMAP_2.PAS)+5) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE006-pas-umap-with-singleR-at-ges.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = ges_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	
	
	pdf( file = file.path( reports.dir , "TE006-pas-umap-with-ges-clusters.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()	

	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.GES , UMAP_2.GES , color = pas_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		scale_color_brewer(palette = "Set1") +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	# coord_fixed()
	
	pdf( file = file.path( reports.dir , "TE006-ges-umap-with-pas-clusters.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.GES , UMAP_2.GES , color = ges_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		theme_minimal()
	# facet_wrap( ~model , scales = "free" , as.table = TRUE ) +
	# theme( legend.position = "none" )
	
	
	pdf( file = file.path( reports.dir , "TE006-ges-umap-with-ges-clusters.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.GES , UMAP_2.GES , color = ges_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		# scale_color_brewer(palette = "Set1") +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		facet_wrap( ~ges_cluster_id ) +
		theme_minimal()
	
	pdf( file = file.path( reports.dir , "TE006-ges-umap-with-ges-clusters-facets.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.GES , UMAP_2.GES , color = pas_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		scale_color_brewer(palette = "Set1") +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		coord_fixed() +
		facet_wrap( ~pas_cluster_id ) +
		theme_minimal()
	
	pdf( file = file.path( reports.dir , "TE006-ges-umap-with-pas-clusters-facets.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()	
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = pas_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		scale_color_brewer(palette = "Set1") +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		facet_wrap( ~pas_cluster_id ) +
		coord_fixed() +
		theme_minimal()
	
	pdf( file = file.path( reports.dir , "TE006-pas-umap-with-pas-clusters-facets.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()			
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = ges_cluster_id ) ) +
		# geom_text_repel( aes( label = model ) , size = 5 ) +
		# geom_text_repel( aes( label = ID ) , size = 5 ) +
		# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
		geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.25) +
		# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
		# scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
		# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
		xlab( paste0( "UMAP 1") ) +
		ylab( paste0( "UMAP 2") ) +
		facet_wrap( ~ges_cluster_id ) +
		coord_fixed() +
		theme_minimal()
	
	pdf( file = file.path( reports.dir , "TE006-pas-umap-with-ges-clusters-facets.pdf") , width = 6 , height = 6 )
		print(umap_plot)
	dev.off()		
	
	# umap_plot <- ggplot( data = umap.df , aes( PAS_UMAP_1 , PAS_UMAP_2 , color = Rspo3 ) ) +
	# 	# geom_text_repel( aes( label = model ) , size = 5 ) +
	# 	# geom_text_repel( aes( label = ID ) , size = 5 ) +
	# 	# geom_point( alpha = 0.5 , stroke = 0.5 , size = 2 ) +
	# 	geom_point( alpha = 0.5 , stroke = 1 , size = 3 ) +
	# 	geom_density_2d(lineend = "butt", color = "black",size=0.5,alpha=0.5) +
	# 	# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
	# 	scale_color_gradient2( breaks = c(-24,-16,-5,-1,0,1,5,16,24) , midpoint=0, low="blue", mid="white",high="red" ) +
	# 	xlab( paste0( "UMAP 1") ) +
	# 	ylab( paste0( "UMAP 2") ) +
	# 	theme_minimal()
	# 
	# pdf( file = file.path( reports.dir , "TE006-pas-umap-analysis-with-seurat-clusters-Rspo3.pdf") , width = 6 , height = 6 )
	# 	print(umap_plot)
	# dev.off()		
	
}

print_msg_info(">>> Feature plot on Lgr5 and other markers on UMAP done on PAS")
{
	x <- vp_seurat@reductions$umap@cell.embeddings
	str(x,1)
	x <- x %>% as.data.frame() %>% rownames_to_column("sample_id")
	df <- left_join( df_annot.df , x )
	
	# df <- df %>% arrange(gene_Lgr5)
	# df <- df %>% arrange(cluster_id)
	p <- ggplot(df,mapping = aes(x=UMAP_1,y=UMAP_2,fill=cluster_id)  ) +
		geom_point(shape = 21, colour = "gainsboro",size = 1.5, stroke = 0.35) +
		# scale_fill_viridis() +
		scale_fill_brewer(palette = "Set1") +
		theme_void()
	
	pdf( file.path( reports.dir , "TE006-umap-with-metacell.pdf") , width = 3 , height = 3)
		print(p)
	dev.off()		
	
	p <- ggplot(df,mapping = aes(x=UMAP_1,y=UMAP_2,fill=cluster_id)  ) +
		geom_point(shape = 21, colour = "black",size = 1.5, stroke = 0.35) +
		# scale_fill_viridis() +
		scale_fill_brewer(palette = "Set1") +
		facet_wrap(~cluster_id,nrow = 2) +
		theme_light()
	
	pdf( file.path( reports.dir , "TE006-umap-with-metacell-facet.pdf") , width = 6 , height = 4)
		print(p)
	dev.off()		
	
}

## Saving object ----
print_msg_warn(">>> Saving Seurat object and metadata object")
{
	seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
	metadata_tibble_analysis_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-viper-analysis-with-metacell-metadata.rds") )
	
	attr(clustering.tibble,"timestamp") <- date()
	attr(vp_seurat,"timestamp") <- date()

	my_metadata <- left_join(vpmat_metadata,vp_seurat@meta.data)
	my_metadata$sample_name <- my_sample_id
		
	# attr(clustering.tibble,"run_parameters.yaml") <- yaml::as.yaml(my_params.list)
	attr(clustering.tibble,"run_parameters.list") <- my_params.list
	
	identical( colnames(vp_seurat) , clustering.tibble$sample_id )
	# saveRDS( vp_seurat , seurat_viper_analysis_data_filename )
	# saveRDS( clustering.tibble , metadata_tibble_analysis_filename )
}

## Nice Plot Printing ----
print_msg_info(">>> Printing Nice Plots")
{
	# tmp %>% group_by(sample_name,cluster_id) %>% tally()
	# tmp <- tmp %>%
	# 	group_by(sample_name) %>%
	# 	mutate(n_total= n() ) %>%
	# 	group_by(cluster_id, add=TRUE) %>%
	# 	mutate(cluster_percentage=paste0(round(100*n()/n_total,2),'%')) %>%
	# 	distinct(sample_name,cluster_id,cluster_percentage,n_total)
	tmp <- my_metadata %>% mutate(pippo=1)
	umap_plot <- ggplot(tmp, aes(fill=pas_cluster_id, x=sample_name,y=pippo)) + 
		geom_bar(position="fill",stat = "identity") +
		scale_fill_brewer(palette = "Set1") +
		scale_y_continuous(labels = scales::percent_format()) +
		theme_light() +
		xlab("Sample") +
		ylab("Percent of cells per cluster") +
		theme( axis.text.x = element_text(face = "plain",size = 9,angle = 45) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12)
					 )
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-cluster-percent-cells.pdf") ) , width = 2.5 , height = 6 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = my_metadata , aes( UMAP_1.PAS , UMAP_2.PAS , color = pas_cluster_id ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_brewer(palette = "Set1") +
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
		)		
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-cluster-id.pdf") ) , width = 3.5 , height = 2.5 )
		print(umap_plot)
	dev.off()		
	
	# umap_plot <- ggplot( data = my_metadata , aes( Ascl2 , Ar , color = pas_cluster_id ) ) +
	# 	# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
	# 	geom_point( alpha = 1 , size = 0.5 ) +
	# 	geom_density_2d(lineend = "butt", color = "black",size=0.5,alpha=0.5) +
	# 	# geom_density_2d_filled(contour_var = "ndensity") + 
	# 	# geom_point( data = umap.df %>% filter(grepl("clones",cell_line)) , alpha = 0.5 , stroke = 2 , size = 4 ) +
	# 	# scale_color_gradient2( midpoint=0, low="steelblue", mid="white",high="red4" ) +
	# 	# scale_color_gradient2( midpoint=0, low="dodgerblue4", mid="white",high="red4" ) +
	# 	# scale_colour_manual( name = "cluster_id" , values = rep( cluster_colors , table(df_annot.df$cluster_id) ) ) +
	# 	scale_color_brewer(palette = "Set1") +
	# 	xlab( paste0( "Ascl2") ) +
	# 	ylab( paste0( "Ar") ) +
	# 	# labs(color = a_gene) + 
	# 	xlim(min(my_metadata$Ascl1)-2,max(my_metadata$Ascl1)+2) +
	# 	ylim(min(my_metadata$Ar)-2,max(my_metadata$Ar)+2) + 
	# 	# coord_fixed() +
	# 	theme_minimal() +
	# 	facet_wrap( ~sample_name , scales = "free" , as.table = TRUE )
	# # theme( legend.position = "none" )
	# # coord_fixed()
	# 
	# pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-Ascl1-vs-Ar.pdf") ) , width = 3.5 , height = 3 )
	# 	print(umap_plot)
	# dev.off()	
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(cytotrace_score_ges) , aes( UMAP_1.PAS , UMAP_2.PAS , color = cytotrace_score_ges ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_viridis(option = "B") +
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
		)		
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-cytotrace-ges.pdf") ) , width = 4 , height = 2.5 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(desc(stemness_index)) , aes( UMAP_1.PAS , UMAP_2.PAS , color = stemness_index ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_viridis() +
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
		)		
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-stemness-pas.pdf") ) , width = 3.5 , height = 2.5 )
		print(umap_plot)
	dev.off()			
	
	# umap_plot <- ggplot( data = my_metadata %>% arrange(gene_Mki67) , aes( UMAP_1.PAS , UMAP_2.PAS , color = gene_Mki67 ) ) +
	# 	# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
	# 	geom_point( alpha = 1 , size = 0.1 ) +
	# 	geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
	# 	# scale_color_viridis() +
	# 	scale_color_gradient(low="white",high="darkorange")	 +
	# 	xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
	# 	ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
	# 	# facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
	# 	theme_minimal() +
	# 	xlab("UMAP 1") +
	# 	ylab("UMAP 2") +
	# 	theme( axis.text.x = element_text(face = "plain",size = 9) , 
	# 				 axis.text.y = element_text(face = "plain",size = 9) ,
	# 				 axis.title.x = element_text(face = "bold",size = 12) ,
	# 				 axis.title.y = element_text(face = "bold",size = 12) 
	# 	)		
	# 
	# pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-Mki67-gene.pdf") ) , width = 3.5 , height = 2.5 )
	# print(umap_plot)
	# dev.off()		
	
	# umap_plot <- ggplot( data = my_metadata %>% arrange(S_Phase) , aes( UMAP_1.PAS , UMAP_2.PAS , color = S_Phase ) ) +
	# 	# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
	# 	geom_point( alpha = 0.75 , size = 0.25 ) +
	# 	geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
	# 	scale_color_viridis() +
	# 	xlim(min(my_meta.filtered$UMAP_1.PAS)-2,max(my_meta.filtered$UMAP_1.PAS)+2) +
	# 	ylim(min(my_meta.filtered$UMAP_2.PAS)-2,max(my_meta.filtered$UMAP_2.PAS)+2) + 
	# 	facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
	# 	theme_minimal() +
	# 	xlab("UMAP 1") +
	# 	ylab("UMAP 2") +
	# 	theme( axis.text.x = element_text(face = "plain",size = 9) , 
	# 				 axis.text.y = element_text(face = "plain",size = 9) ,
	# 				 axis.title.x = element_text(face = "bold",size = 12) ,
	# 				 axis.title.y = element_text(face = "bold",size = 12) 
	# 	)		
	# 
	# pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-S-Phase.pdf") ) , width = 6 , height = 2.5 )
	# print(umap_plot)
	# dev.off()			
	
	df <- my_metadata %>% arrange(abs(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) %>% dplyr::rename(EMT=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
	umap_plot <- ggplot( data = df , 
											 aes( UMAP_1.PAS , UMAP_2.PAS , color = EMT ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_gradient2(low = "deepskyblue3",mid = "white",high="brown3",midpoint = 0) + 
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		# facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		# guides(fill = guide_legend(title = "EMT", title.position = "left"))
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
					 # legend.position="bottom", legend.box = "horizontal"
		)		
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-EMT.pdf") ) , width = 3.5 , height = 2.5 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(desc(HALLMARK_WNT_BETA_CATENIN_SIGNALING)) , aes( UMAP_1.PAS , UMAP_2.PAS , color = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_gradient2(low = "deepskyblue3",mid = "white",high="brown3",midpoint = 0) + 
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		# facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		guides(color=guide_legend(title="WNT")) +
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
		)		
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-WNT.pdf") ) , width = 3.5 , height = 2.5 )
		print(umap_plot)
	dev.off()		
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(desc(HALLMARK_WNT_BETA_CATENIN_SIGNALING)) , aes( UMAP_1.PAS , UMAP_2.PAS , color = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_gradient2(low = "deepskyblue3",mid = "white",high="brown3",midpoint = 0) + 
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		# facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		guides(color=guide_colorbar(title="WNT")) +
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
		)		
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-WNT.pdf") ) , width = 3.5 , height = 2.5 )
		print(umap_plot)
	dev.off()			
	
	require(viridisLite)
	umap_plot <- ggplot( data = my_metadata , aes( UMAP_1.PAS , UMAP_2.PAS ) ) +
		geom_density_2d_filled(contour_var = "ndensity") + 
		scale_fill_viridis_d(option = "B") +
		xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
		ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
		# coord_fixed() +
		facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
		theme_minimal() +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		theme( axis.text.x = element_text(face = "plain",size = 9) , 
					 axis.text.y = element_text(face = "plain",size = 9) ,
					 axis.title.x = element_text(face = "bold",size = 12) ,
					 axis.title.y = element_text(face = "bold",size = 12) 
		)		

	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-density.pdf" )) , width = 3.5 , height = 2.5 )
		print(umap_plot)
	dev.off()		
	
	p <- ggplot( my_metadata , aes(pas_cluster_id, cytotrace_score_ges , 
												# color = sc_entropy ,
												# fill = Top2a
												fill = pas_cluster_id
	) ) +
		geom_boxplot() +
		scale_fill_brewer(palette = "Set1") +
		scale_x_discrete( position = "top" ) +
		theme_light() +
		theme( axis.text.x = element_text(angle = 45, hjust = 0, face = "bold") )
	
	pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-boxplot-cytotrace.pdf" )) , width = 4 , height = 2.5 )
		print(p)
	dev.off()		
	
	print_msg_info(">>> >> Printing Gene Expression on top of PAS UMAP")
	{
		gex.mat <- my_gex_seurat@assays$SCT@scale.data %>% as.matrix()
		
		gene_list <- c("Lgr5","Mki67",'Olfm4',"Ung","Sox4",'Tnfrsf19')
		
		for (my_gene in gene_list)
		{
			# my_gene <- "Lgr5"	
			my_gene.df <- gex.mat[my_gene,] %>% as.data.frame()
			colnames(my_gene.df) <- "a_gene_expr"
			my_gene.tibble <- my_gene.df %>% rownames_to_column("cell_id") %>% as_tibble()
			tibble4plot <- left_join( my_metadata , my_gene.tibble , by = c("cell_id"="cell_id") )
			
			df <- tibble4plot %>% arrange(a_gene_expr)
			p <- ggplot( data = df , aes(UMAP_1.PAS, UMAP_2.PAS , color = a_gene_expr ) ) +
				geom_point( alpha = 0.75 , size = 0.25 ) +
				# geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
				# scale_color_gradient2(low = "purple",mid = "white",high="orange",midpoint = 0) + 
				scale_color_viridis() +
				xlim(min(my_metadata$UMAP_1.PAS)-2,max(my_metadata$UMAP_1.PAS)+2) +
				ylim(min(my_metadata$UMAP_2.PAS)-2,max(my_metadata$UMAP_2.PAS)+2) + 
				# facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
				theme_minimal() +
				xlab("UMAP 1") +
				ylab("UMAP 2") +
				guides(color=guide_colorbar(title=my_gene)) +
				coord_equal() +
				theme( axis.text.x = element_text(face = "plain",size = 9) , 
							 axis.text.y = element_text(face = "plain",size = 9) ,
							 axis.title.x = element_text(face = "bold",size = 12) ,
							 axis.title.y = element_text(face = "bold",size = 12) 
				)		
			
			pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-expr-",my_gene,".pdf") ) , width = 4 , height = 2.5 )
				print(p)
			dev.off()				
		}
	}
	
}

## Printing Heatmap VPMAT with top cMRs ----
print_msg_info(">>> >>> Printing Heatmap VPMAT with top cMRs")
{
	require(circlize)
	require(ComplexHeatmap)
	
	mat <- vp_seurat@assays$VIPER@scale.data %>% as.matrix()
	stopifnot( identical( clustering.tibble$sample_id , colnames(mat) ) )
	clusters <- clustering.tibble$cluster_id
	names(clusters) <- clustering.tibble$sample_id		
	weigths <- clustering.tibble$sil_score
	names(weigths) <- clustering.tibble$sample_id	
	
	print_msg_warn(">>> >> Selecting only TFs and coTFs ...")	
	{
		print(dim(mat))
		tmp_1 <- read_csv("~/Clouds/Dropbox/Data/oncoloop/gemms-data/n1lists-pisces/lists/tfs.csv",col_names = FALSE)
		tmp_2 <- read_csv("~/Clouds/Dropbox/Data/oncoloop/gemms-data/n1lists-pisces/lists/cotfs.csv",col_names = FALSE)
		index <- (rownames(mat) %in% tmp_1$X1) | (rownames(mat) %in% tmp_2$X1)
		mat <- mat[index,]			
		print(dim(mat))
	}				
	print_msg_info(">>> >> >> Selecting only samples with weights > 0.25 for integrations")
	{
		index <- weigths > 0.25
		mat <- mat[,index]
		clusters <- clusters[index]
		weigths <- weigths[index]
	}		
	
	mat.stouffered <- doStoufferFromClusters( mat , clusters , weights = weigths )
	n_top <- 5
	# index <- names(sort(mat.stouffered["Ar",],decreasing = T))
	# mat.stouffered <- mat.stouffered[ , index ]
	selected_proteins <- Reduce( union , apply( mat.stouffered , 2 , function(x) { x <- sort(x,decreasing = T) ; names(c(head(x,n_top))) } ) )
	sort(selected_proteins)
	
	selected_cells <- sample( colnames(mat) , size = 1000 , replace = FALSE )
	mat <- mat[,selected_cells]
	
	df_annot.df <- vp_seurat@meta.data %>% dplyr::select(cell_id,cluster_id=seurat_clusters)
	df_annot.df$cluster_id <- factor(as.numeric(df_annot.df$cluster_id))
	df_annot.df <- df_annot.df[match(selected_cells,df_annot.df$cell_id),]
	
	stopifnot(identical( df_annot.df$cell_id , colnames(mat) ))
	
	x <- clustering.tibble %>% dplyr::select(sample_id,sil_score_order)
	df_annot.df <- left_join(df_annot.df,x,by=c("cell_id"="sample_id"))
	df_annot.df <- df_annot.df %>% arrange(sil_score_order)
	mat <- mat[,match(df_annot.df$cell_id,colnames(mat))]
	df_annot.df$sil_score_order <- NULL
	
	# x <- df_annot.df$cluster_id
	# names(x) <- df_annot.df$cell_id
	# mat.stouffered <- doStoufferFromClusters( mat , x )
	# n_top <- 15
	# index <- names(sort(mat.stouffered["Ar",],decreasing = T))
	# mat.stouffered <- mat.stouffered[ , index ]
	# selected_proteins <- Reduce( union , apply( mat.stouffered , 2 , function(x) { x <- sort(x,decreasing = T) ; names(c(head(x,n_top))) } ) )
	# sort(selected_proteins)
	
	df_annot.df <- df_annot.df %>% dplyr::select(cluster_id)
	n_clust <- nlevels(df_annot.df$cluster_id)
	cluster_colors <- brewer.pal(n_clust, "Set1")
	names(cluster_colors) <- 1:n_clust
	df_annot_cols = HeatmapAnnotation( df = df_annot.df ,
																		 col = list(cluster_id=cluster_colors)  , 
																		 show_legend = TRUE , na_col = "gray95"
	)  		
	
	col_fun = colorRamp2(c(-3, 0, 3), c("deepskyblue3", "white", "brown3"))
	
	nes_hm <- Heatmap( mat[selected_proteins,] ,
										 col = col_fun ,
										 use_raster = TRUE ,
										 heatmap_legend_param = list(color_bar = "continuous") ,	
										 # show_heatmap_legend = FALSE , 
										 # clustering_method_columns = "ward.D2" ,			
										 cluster_columns = FALSE ,
										 cluster_rows = FALSE ,
										 show_column_names = FALSE ,
										 # cluster_columns = columns_dend ,
										 name = "TF Activity" ,  
										 row_names_gp = gpar(fontsize = 10) , column_names_gp = gpar(fontsize = 8) ,
										 row_title = "TF Activity", 
										 row_title_side = "left" , row_title_gp = gpar(fontsize = 10, fontface = "bold")
	)
	
	ht_list = df_annot_cols %v% nes_hm
	
	p <- draw(ht_list, column_km = 1 , 
						# annotation_legend_list = pd ,
						heatmap_legend_side = "right", 
						annotation_legend_side = "bottom" )
	
	pdf( file.path( reports.dir , "clusters-tf-heatmap.pdf") , width = 5 , height = 9 )
	print(p)
	dev.off()  			
}

# message(">>> Benchmarking the silhouette score using Seurat")
# {
# 	source("../vaxtools/R/silhouette-score-benchmarking-optimal-number-of-clusters.R")
# 	mat <- vpmat
# 	df <- generate_optimal_clustering_solutions_on_seurat( mat , min_res = 0 , max_res = 0.5 , n_sim = 10 , n_iter = 5 , n_cells = 3000 )
# }
