## 
# Single cell Analysis of TE001 (VIPER Analysis) | SUBLUSTER-specific network but whole dataset signature
# ---------------------------------------------------------------
# New approach, metacells were created based on gene expression clustering, 
# 	but network was created on cluster-specific metacells that were merged in one dataset for one ARACNe run
# ------------------------------------------------------------------------
# system.time({ source("sources/ermanno-data-analysis/TE001/TE001-viper-with-subnetworks-one-signature.R") })

retainOnlyTFs <- function(network) {
	
	print_msg_info(">>> Filtering using my lists")
	x <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat" , col_names = FALSE )
	y <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat" , col_names = FALSE )
	# z <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat" , col_names = FALSE )
	
	selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))	
	
	index <- names(network) %in% selected_regulators
	return( network[index] )
	
}

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
# use_python('/Applications/anaconda3/bin/python')
# use_python('/Users/av2729/opt/anaconda3/bin/python')
# use_condaenv('/Users/afpvax/opt/anaconda3')
# py_config()

print_msg_info(">>> Setting Analysis Params")
{
	isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
	my_sample_id <- "TE001"
	my_seed <- 42
	knn_n_neighbor <- 15
	isComputingOptimalNumberOfClusters <- FALSE	
	regulon_n_targets <- 50
	ges_n_feature <- 20000
	metaviper_mvws <- 10
}

# ## Loading counts matrix ----
# print_msg_info(">>> Loading counts matrix")
# {
# 	# x <- readRDS( file.path(isc_data_dir,"/TE001-seurat-data.rds") )
# 	# my_counts <- x@assays$SCT@data
# 	cpm_filename <- file.path(isc_data_dir,"TE001-cpm.rds")
# 	file.info(cpm_filename)["mtime"]
# 	x <- readRDS( cpm_filename )
# 	my_counts <- x
# 	dim(my_counts)
# }

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
	my_params.list$is_SCT_do.scale = FALSE
	my_params.list$ges_n_feature <- ges_n_feature
	my_params.list$viper_signature_method <- "none"
	# my_params.list$viper_signature_method <- "mad"
	my_params.list$regulon_minsize <- 50
	my_params.list$viper_cluster_specific <- TRUE
	my_params.list$viper_mvws_other_regs <- metaviper_mvws
	my_params.list$viper_dist_corr_method <- "pearson"
	# my_params.list$viper_dist_corr_method <- "spearman"
	
	# my_params.list$graph_type <- "nn"
	my_params.list$graph_type <- "snn"
	my_params.list$is_only_TFS <- TRUE
	# my_params.list$sct_variable_n_feats <- sct_variable_n_feats
	# my_params.list$sct_is_to_scale <- FALSE
	# my_params.list$n_pcs <- n_pcs
	# my_params.list$pcs_to_use <- 1:5
	my_params.list$pcs_to_use <- 1:10
	# my_params.list$clustering_algorithm <- clustering_algorithm
	
	my_params.list$comments <- "metacellClusterSpecOneSignature"
	
	# tmp <- paste( my_params.list , collapse = "-" )
	# create_workspace(tmp)
	
	tmp <- "TE001-viper-di-ale-subnetworks-one-signature"
	create_workspace(tmp)
	
	my_params.list$results_dir = reports.dir
	# my_params.list$cpm_filename <- cpm_filename	
	
	my_yaml_file <- file.path(reports.dir,"run-parameters.yaml")
	yaml::write_yaml(my_params.list %>% as.list(),my_yaml_file)
	
}

## VIPER Run di Ale ----
print_msg_info(">>> VIPER Analysis")
{
	require(Seurat)
	require(viper)

	seurat_original_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-original-data.rds") )
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
	seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )

	print_msg_info(">>> Extracting Clustering Information from processed Seurat Object")

	seurat_analysis <- readRDS(seurat_analysis_data_filename)
	format(object.size(seurat_analysis),"Mb")
	file.info(seurat_analysis_data_filename)

	vp_seurat_analysis <- readRDS(seurat_viper_analysis_data_filename)
	format(object.size(vp_seurat_analysis),"Mb")
	file.info(seurat_viper_analysis_data_filename)

	seurat_analysis_original <- readRDS(seurat_original_data_filename)
	my_gex_seurat <- seurat_analysis_original %>% subset( cells = colnames(seurat_analysis) ) %>%
		SCTransform( do.scale = my_params.list$is_SCT_do.scale ,
								 do.center = my_params.list$is_SCT_do.center ,
								 return.only.var.genes = TRUE ,
								 variable.features.n = my_params.list$ges_n_feature
		)
	print_msg_info(">>> Computing Gene Expression Signatures over " , ncol(my_gex_seurat) , " cells")
	print_msg_info(">>> Computing Gene Expression Signatures over " , nrow(my_gex_seurat) , " genes")
	my_ges <- as.matrix(my_gex_seurat@assays$SCT@scale.data)
	dim(my_ges)

	print(table(vp_seurat_analysis$seurat_clusters))

	print_msg_info(">>> >> Statically assigning networks to clusters")
	{
		# clusters_networks.list <- vector("list",3)
		clusters_networks.list <- list()
		clusters_networks.list[["0"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
		clusters_networks.list[["1"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c2_unPruned.rds"))
		clusters_networks.list[["2"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c3and4_unPruned.rds"))
		str(clusters_networks.list,1)
	}

	vpmat.list <- vector("list",3)
	names(vpmat.list) <- c("0","1","2")
	for ( a_cluster_id in names(table(vp_seurat_analysis$seurat_clusters)) )
	{
		# index <- vp_seurat_analysis$seurat_clusters == a_cluster_id
		# selected_cells <- colnames(vp_seurat_analysis)[index]
		# length(selected_cells)

		# index <- match( selected_cells , colnames(my_ges) )
		# my_ges.subset <- my_ges[,index]
		my_ges.subset <- my_ges
		dim(my_ges.subset)

		x <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat") , col_names = FALSE )
		y <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat") , col_names = FALSE )
		z <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat") , col_names = FALSE )
		if( my_params.list$is_only_TFS )
		{
			selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))
		} else {
			selected_regulators <- unique(c( x$X1 , y$X1 , z$X1 ))
		}
		print(length(selected_regulators))

		print_msg_warn("*** Using Cluster Specific Network ***")

		print_msg_info(">>> Total regulatory proteins in the actual network: " , length(clusters_networks.list[[a_cluster_id]]) )
		index <- names(clusters_networks.list[[a_cluster_id]]) %in% selected_regulators
		print_msg_info(">>> Total regulatory proteins to expect from selected ones: " , sum(index) )
		clusters_networks.list[[a_cluster_id]] <- clusters_networks.list[[a_cluster_id]][index]
		print_msg_info(">>> Total regulatory proteins recovered in the network: " , length(clusters_networks.list[[a_cluster_id]]) )
		clusters_networks.list[[a_cluster_id]] <- pruneRegulon(clusters_networks.list[[a_cluster_id]],50)

		# vpmat <- viper( my_ges.subset , clusters_networks.list[[a_cluster_id]] , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
		vpmat <- viper( my_ges.subset , clusters_networks.list[[a_cluster_id]] , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
		print_msg_info(">>> VPMAT matrix shape: " , dim(vpmat) )

		vpmat <- vpmat[ , match( colnames(my_ges.subset) , colnames(vpmat) ) ]
		stopifnot( identical( colnames(my_ges.subset) , colnames(vpmat) ) )

		vpmat.list[[a_cluster_id]] <- vpmat

	}

	str(vpmat.list,1)

	cohort_list = vpmat.list
	vpmat <- mergeCohorts(cohort_list = cohort_list)
	print(dim(vpmat))
	regs_in_common <- rownames(vpmat)

	print_msg_info(">>> Selecting best matching network for each cell")
	net_selection.list <- vector("list",3)
	names(net_selection.list) <- c("0","1","2")
	for (a_cohort in names(cohort_list))
	{
		x_vpmat <- cohort_list[[a_cohort]]
		dim(x_vpmat)
		n_mrs <- 50
		net_selection.list[[a_cohort]] <- apply( x_vpmat , 2 , function(x) {
			# sum(abs( c(head(sort(x,decreasing = T),n_mrs),tail(sort(x,decreasing = T),n_mrs)) ))
			sum( abs( head(sort(x,decreasing = T),n_mrs) ))
		})
	}

	str(net_selection.list,1)
	x <- do.call("rbind",net_selection.list)
	selected_network_per_cell <- apply( x , 2 , which.max )
	print(table(selected_network_per_cell))

	print_msg_info(">>> Filtering in cell specific viper run")
	vpmat_selection.list <- vector("list",3)
	names(vpmat_selection.list) <- c("0","1","2")
	for (a_cohort in names(cohort_list))
	{
		x_vpmat <- cohort_list[[a_cohort]]
		dim(x_vpmat)

		index <- selected_network_per_cell %in% (as.numeric(a_cohort)+1)
		vpmat_selection.list[[a_cohort]] <- x_vpmat[,index]
	}
	str(vpmat_selection.list,1)

	vpmat.metaviper <- viper( my_ges.subset , clusters_networks.list ,
														mvws = my_params.list$viper_mvws_other_regs ,
														method = my_params.list$viper_signature_method ,
														minsize = my_params.list$regulon_minsize , verbose = TRUE)
	dim(vpmat.metaviper)

	tot_features <- rownames(vpmat.metaviper)

	print_msg_info(">>> Merging specific vpmat with metavipered one (Cluster specific) ...")
	vpmat_merged.list <- vector("list",3)
	names(vpmat_merged.list) <- c("0","1","2")
	for (a_cohort in names(cohort_list))
	{
		x_vpmat <- cohort_list[[a_cohort]]
		dim(x_vpmat)
		missing_features <- setdiff(tot_features,rownames(x_vpmat))
		vpmat_merged.list[[a_cohort]] <- rbind(x_vpmat,vpmat.metaviper[missing_features,colnames(x_vpmat)])
	}
	str(vpmat_merged.list,1)

	cohort_list = vpmat_merged.list
	vpmat <- mergeCohorts(cohort_list = vpmat_merged.list)
	print(dim(vpmat))
	regs_in_common <- rownames(vpmat)

	vpmat <- vpmat[ , match(colnames(vpmat.metaviper),colnames(vpmat))]
	stopifnot( identical( colnames(vpmat) , colnames(vpmat.metaviper) ) )
	stopifnot( identical( colnames(vpmat) , colnames(my_ges.subset) ) )

	vpsim <- viperSimilarity(vpmat)
	vpdist <- as.dist( vpsim )
	# vpdist <- as.dist( 1-cor(
	# 	# vpmat[order(rowVars(vpmat),decreasing = T)[1:200],] ,
	# 	vpmat ,
	# 	method = my_params.list$viper_dist_corr_method )
	# )

	# saveRDS(object = list(gesmat=my_ges.rankNormed,vpmat=vpmat,vpsim=vpsim,vpdist=vpdist) ,
	# 				file.path(isc_data_dir,"/TE001-viper-analysis-with-metacell.rds") )

}


# ## metaVIPER Run ----
# print_msg_info(">>> VIPER Analysis")
# {
# 	require(Seurat)
# 	require(viper)
# 	
# 	seurat_original_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-original-data.rds") )
# 	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
# 	seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
# 	
# 	print_msg_info(">>> Extracting Clustering Information from processed Seurat Object")
# 	
# 	seurat_analysis <- readRDS(seurat_analysis_data_filename)
# 	format(object.size(seurat_analysis),"Mb")
# 	file.info(seurat_analysis_data_filename)
# 	
# 	vp_seurat_analysis <- readRDS(seurat_viper_analysis_data_filename)
# 	format(object.size(vp_seurat_analysis),"Mb")
# 	file.info(seurat_viper_analysis_data_filename)
# 	
# 	seurat_analysis_original <- readRDS(seurat_original_data_filename)
# 	my_gex_seurat <- seurat_analysis_original %>% subset( cells = colnames(seurat_analysis) ) %>% 
# 		SCTransform( do.scale = my_params.list$is_SCT_do.scale , 
# 								 do.center = my_params.list$is_SCT_do.center , 
# 								 return.only.var.genes = TRUE ,
# 								 variable.features.n = my_params.list$ges_n_feature
# 		)
# 	print_msg_info(">>> Computing Gene Expression Signatures over " , ncol(my_gex_seurat) , " cells")
# 	print_msg_info(">>> Computing Gene Expression Signatures over " , nrow(my_gex_seurat) , " genes")
# 	my_ges <- as.matrix(my_gex_seurat@assays$SCT@scale.data)
# 	dim(my_ges)	
# 	
# 	print(table(vp_seurat_analysis$seurat_clusters))
# 	
# 	# clusters_networks.list <- vector("list",3)
# 	clusters_networks.list <- list()
# 	clusters_networks.list[["0"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
# 	clusters_networks.list[["1"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c2_unPruned.rds"))
# 	clusters_networks.list[["2"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c3and4_unPruned.rds"))
# 	str(clusters_networks.list,1)
# 	
# 	clusters_networks.list <- lapply( clusters_networks.list , retainOnlyTFs )
# 	clusters_networks.list <- lapply( clusters_networks.list , pruneRegulon , 50 )
# 	
# 	# vpmat <- viper( my_ges.subset , clusters_networks.list[[a_cluster_id]] , method = my_params.list$viper_signature_method , minsize = my_params.list$regulon_minsize , verbose = TRUE)
# 	vpmat <- viper( my_ges , clusters_networks.list , 
# 									method = my_params.list$viper_signature_method , 
# 									mvws = my_params.list$viper_mvws_other_regs ,
# 									minsize = my_params.list$regulon_minsize , verbose = TRUE)
# 	print_msg_info(">>> VPMAT matrix shape: " , dim(vpmat) )
# 	
# 	stopifnot( identical( colnames(my_ges) , colnames(vpmat) ) )				
# 	
# 	
# 	vpsim <- viperSimilarity(vpmat)
# 	# vpdist <- as.dist( vpsim )
# 	vpdist <- as.dist( 1-cor(
# 		# vpmat[order(rowVars(vpmat),decreasing = T)[1:200],] ,
# 		vpmat ,
# 		method = my_params.list$viper_dist_corr_method )
# 	)
# 	
# 	# saveRDS(object = list(gesmat=my_ges.rankNormed,vpmat=vpmat,vpsim=vpsim,vpdist=vpdist) ,
# 	# 				file.path(isc_data_dir,"/TE001-viper-analysis-with-metacell.rds") )
# 	
# }

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
	stemness_index.pas <- getStemnessIndex(vpmat_human)	
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
														stemness_index.pas = stemness_index.pas )
	
	pathways <- as_tibble(vpmat_pathways)
	pathways_names <- rownames(vpmat_pathways)
	cell_id <- colnames(vpmat_pathways)
	pathways <- as_tibble(t(pathways))
	colnames(pathways) <- pathways_names	
	pathways <- bind_cols( cell_id = cell_id , pathways )
	
	vpmat_metadata <- left_join( vpmat_metadata , pathways )		
	
	# colnames(seurat_analysis@meta.data)
	x <- seurat_analysis@meta.data %>% dplyr::select(cell_id,cytotrace_score.ges,stemness_index.ges,
																									 mt_percent,seurat_clusters.ges=seurat_clusters,
																									 singleR_labels.ges=singleR_labels)
	
	vpmat_metadata <- left_join( vpmat_metadata , x , by = c("cell_id"="cell_id") )

	print_msg_info(">>> >> Loading iterative clustering solution ...")
	{
		iter_clustering_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-iter-clustering-viper-analysis-with-metacell.tsv") )	
		iter_cluster.tibble <- read_delim(iter_clustering_filename,delim = ",")
		
		vpmat_metadata <- left_join( vpmat_metadata , iter_cluster.tibble , by = c( "cell_id" = "sample_id" ) )
		vpmat_metadata <- vpmat_metadata %>% dplyr::rename(iter_cluster_id=cluster_id)
	}

}

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
	
	## Clustering Analysis ----
	print_msg_info(">>> Clustering Analysis ...")
	{
		
	# 	my_ssn_graph <- FindNeighbors( as.matrix(vpdist) , 
	# 																 # dims = pcs_to_use , 
	# 																 assay = "VIPER" , 
	# 																 distance.matrix = TRUE , 
	# 																 verbose = TRUE , 
	# 																 k.param = knn_n_neighbor ,
	# 																 # annoy.metric = "cosine" ,
	# 																 annoy.metric = "euclidean" ,
	# 																 compute.SNN = TRUE )
	# 	vp_seurat@graphs <- my_ssn_graph		
	# 	
	# 	# vp_seurat <- FindNeighbors( vp_seurat , 
	# 	# 													# dims = pcs_to_use , 
	# 	# 													assay = "VIPER" , 
	# 	# 													verbose = TRUE , 
	# 	# 													k.param = knn_n_neighbor ,
	# 	# 													compute.SNN = TRUE ,
	# 	# 													annoy.metric = "cosine" ,
	# 	# 													# distance.matrix = as.matrix(vpdist) ,
	# 	# 													do.plot = FALSE )		
	# 	
	# 	## Silhouette Analysis ----
	# 	print_msg_info(">>> Selection of parameters for optimal clustering solution")
	# 	{
	# 		library(foreach)
	# 		library(doParallel)
	# 		
	# 		require(cluster)
	# 		require(factoextra)
	# 		
	# 		n_cores <- detectCores()-2
	# 		myCluster <- makeCluster(n_cores,type = "FORK")
	# 		registerDoParallel(myCluster)
	# 		
	# 		n_cells_to_subsample <- round(ncol(vp_seurat)*70/100)
	# 		
	# 		.resolutions <- seq(0.01,0.25,by = 0.01)
	# 		.bootstraps <- 1:2
	# 		.knns <- seq(5,15,by = 2)		
	# 		
	# 		print(system.time({ 
	# 			result <- foreach::foreach( a_res = .resolutions , .combine = 'rbind' ) %dopar% {
	# 				
	# 				knn = rep( .knns , length(.bootstraps) )
	# 				bootstrap = rep( .bootstraps , length(.knns) )
	# 				index <- 1:length(knn)
	# 				my_sil.df <- tibble( index = index ,
	# 														 bootstrap = bootstrap , knn = knn , resolution = 0 ,
	# 														 tot_neg_sil_cells = 0 ,
	# 														 sil_avg = 0 , sil_mean_median = 0 , n_clust = 0 )
	# 				nrow(my_sil.df)
	# 				my_dist <- as.matrix(as.dist( 1-cor( vp_seurat@assays$VIPER@scale.data , method = my_params.list$viper_dist_corr_method )))
	# 				for ( a_index in my_sil.df$index )
	# 				{
	# 					set.seed(a_index)
	# 					print_msg_info("--- Silhouette score computation resolution/bootstrap :" , a_res , "|" , my_sil.df$bootstrap[a_index] )
	# 					selected_samples <- sample(colnames(vp_seurat),size=n_cells_to_subsample,replace = FALSE)
	# 					x <- vp_seurat[ , colnames(vp_seurat) %in% selected_samples ]
	# 					x.dist <- my_dist[ colnames(my_dist) %in% selected_samples , colnames(my_dist) %in% selected_samples ]
	# 					# # Using PCA instead of distance matrix
	# 					# my_ssn_graph <- FindNeighbors( x , 
	# 					# 															 dims = pcs_to_use ,
	# 					# 															 distance.matrix = FALSE , 
	# 					# 															 verbose = TRUE , 
	# 					# 															 # k.param = knn_n_neighbor ,
	# 					# 															 k.param = my_sil.df$knn[a_index] ,
	# 					# 															 # annoy.metric = "cosine" ,
	# 					# 															 annoy.metric = "euclidean" ,
	# 					# 															 compute.SNN = TRUE )					
	# 					
	# 					my_ssn_graph <- FindNeighbors( x.dist ,
	# 																				 # dims = pcs_to_use ,
	# 																				 # assay = "VIPER" ,
	# 																				 distance.matrix = TRUE ,
	# 																				 verbose = TRUE ,
	# 																				 # k.param = knn_n_neighbor ,
	# 																				 k.param = my_sil.df$knn[a_index] ,
	# 																				 # annoy.metric = "cosine" ,
	# 																				 annoy.metric = "euclidean" ,
	# 																				 compute.SNN = TRUE )
	# 					x@graphs <- my_ssn_graph				
	# 					x <- FindClusters( x , 
	# 														 # graph.name = "VIPER_snn" ,
	# 														 graph.name = my_params.list$graph_type ,
	# 														 resolution = a_res  , 
	# 														 verbose = FALSE ,
	# 														 modularity.fxn = 1 ,
	# 														 # algorithm = 4 , # Leiden
	# 														 algorithm = 1 , # Louvain
	# 														 random.seed = my_seed
	# 					)
	# 					
	# 					if ( nlevels(x$seurat_clusters) == 1 ) next ;
	# 					
	# 					s <- silhouette( as.integer(x$seurat_clusters) , x.dist )
	# 					# pdf( file.path(reports.dir,paste0("sil-res-",a_res,".pdf")) )
	# 					x <- fviz_silhouette(s,print.summary = FALSE)			
	# 					y <- sapply( levels(x$data$cluster) , function(i) mean( x$data$sil_width[ x$data$cluster == i ] ) )
	# 					z <- sapply( levels(x$data$cluster) , function(i) median( x$data$sil_width[ x$data$cluster == i ] ) )
	# 					sil_neg <- sapply( levels(x$data$cluster) , function(i) sum( x$data$sil_width[ x$data$cluster == i ] < 0 ) )
	# 					my_sil.df$sil_avg[ a_index ] = mean(y)
	# 					my_sil.df$sil_mean_median[ a_index ] = mean(z)
	# 					my_sil.df$tot_neg_sil_cells[ a_index ] <- sum(sil_neg)
	# 					my_sil.df$n_clust[ a_index ] = nlevels(x$data$cluster)
	# 					
	# 					# dev.off()
	# 					# View(my_sil.df)    
	# 				}
	# 				
	# 				my_sil.df$resolution = a_res
	# 				
	# 				return(my_sil.df)
	# 				
	# 			} # End of dopar
	# 		})) # End of print
	# 		stopCluster(myCluster)		
	# 		
	# 		my_sil.df <- result
	# 	}
	# 	
	# 	message(">>> >> Printing Optimal Clustering Plot")
	# 	{
	# 		# my_sil.df %>% group_by(resolution) %>% summarize(mean(sil_avg,na.rm=TRUE))
	# 		
	# 		library(ggrepel)
	# 		# p <- ggplot( my_sil.df , aes( x = resolution , y = sil_median ) ) +
	# 		p <- ggplot( my_sil.df , aes( x = resolution , y = sil_avg , group = resolution) ) +
	# 			# geom_line() +
	# 			# geom_boxplot() +
	# 			geom_point() +
	# 			# geom_text( aes(x = resolution , y = 0.1 , label = n_clust) , label.size = 3 ) +
	# 			# geom_text( aes(x = resolution , y = sil_median+0.01 , label = n_clust) , label.size = 3 ) +
	# 			geom_text_repel( aes( label = n_clust ) , segment.color = "gray" ) +
	# 			theme_light()
	# 		
	# 		pdf( file.path(reports.dir,"pas-optimal-clustering.pdf") )
	# 		print(p)
	# 		dev.off()
	# 	}
	# 	
	# 	my_res <- my_sil.df$resolution[ which.max(my_sil.df$sil_avg) ]
	# 	# print_msg_warn(">>> Using MANUAL resolution")
	# 	# my_res <- 0.38
	# 	# my_res <- 0.2
	# 	
	# 	# vp_seurat <- FindClusters(vp_seurat, resolution = 2 , verbose = TRUE,modularity.fxn = 2,algorithm = 4,random.seed = 666 ) 
	# 	vp_seurat <- FindClusters( vp_seurat , 
	# 														 # graph.name = "VIPER_snn" ,
	# 														 graph.name = my_params.list$graph_type ,
	# 														 resolution = my_res , 
	# 														 verbose = FALSE ,
	# 														 modularity.fxn = 1 ,
	# 														 # algorithm = 4 , # Leiden
	# 														 algorithm = 1 , # Louvain
	# 														 seed.use = my_seed
	# 	)
	# 	
	}
	
	print_msg_info(">>> Injecting Clustering Info ...")	
	{
		vp_seurat@meta.data$seurat_clusters <- vp_seurat@meta.data$iter_cluster_id
	}
	print(table(vp_seurat$seurat_clusters))
	
	class(vpdist)
	my_ssn_graph <- FindNeighbors( vpdist ,
																 # dims = pcs_to_use ,
																 # distance.matrix = TRUE ,
																 verbose = TRUE ,
																 k.param = knn_n_neighbor ,
																 # k.param = my_sil.df$knn[a_index] ,
																 # annoy.metric = "cosine" ,
																 annoy.metric = "euclidean" ,
																 compute.SNN = TRUE )
	vp_seurat@graphs <- my_ssn_graph
	# print_msg_warn("\t *** Not Running RunUMAP from Seurat because continuously issues on umap-learn not recognized ***")
	vp_seurat <- RunUMAP( vp_seurat ,
												assay = "VIPER" ,
												n.components = 3L ,
												# graph.name = "VIPER_snn" ,
												# densmap = TRUE ,
												# features = prolif_genes ,
												# assay = "RNA" ,
												# graph = my_params.list$graph_type ,
												graph = "snn" ,
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
	
	## Creating metadata (clustering.tibble) ----
	print_msg_info(">>> Creating metadata (clustering.tibble)")
	{
		library(cluster)
		my_clusters <- as.integer( vp_seurat@meta.data$seurat_clusters %>% as.factor() )
		names(my_clusters) <- colnames(vp_seurat)
		my_clusters <- factor(my_clusters)
		clustering.tibble <- tibble( sample_id = names(my_clusters) ,
																 cluster_id = my_clusters , 
																 sample_ordering = 1:length(my_clusters) ,
																 cluster_reliability = n1platform::clusterReliability(my_clusters,vpsim) ,
																 sil_score = silhouette(as.integer(my_clusters),vpdist)[,3] )
		clustering.tibble %>% group_by(cluster_id) %>% summarise(avg_sil=mean(sil_score))
		
		clustering.tibble <- clustering.tibble %>% 
			arrange(cluster_id,desc(sil_score)) %>%
			mutate(sil_score_order = 1:length(sil_score) )
		# arrange(cluster_id,desc(cluster_reliability)) %>%
		# mutate(cluster_reliability_order = 1:length(cluster_reliability) )
		clustering.tibble <- left_join( clustering.tibble , vpmat_metadata , by = c("sample_id"="cell_id") , suffix = c("","_vpmat_metadata") )	
		
		clustering.tibble <- clustering.tibble %>% arrange(sample_ordering)
		
		print_msg_info(">>> >> Reordering clusters based on total number of cells")
		{
			x <- table(clustering.tibble$cluster_id) %>% sort(decreasing = TRUE)
			
			old_lables <- names(table(clustering.tibble$cluster_id))
			target_labels <- names(x)
			new_labels <- old_lables[ match(clustering.tibble$cluster_id,target_labels) ]
			clustering.tibble$cluster_id <- factor(new_labels) #,levels = target_labels,ordered = TRUE)			
			table(clustering.tibble$cluster_id)
		}
		
		print_msg_info(">>> Integration with gene expression scores")
		{
			genes_to_add <- c("Lgr4","Lgr5","Tnfrsf19","Olfm4","Ascl2","Atoh1","Krt19","Top2a", "Chga" ,
												"Ung","Clu","Ly6a","Rad51","Smarcal1","Ihh", "Fgfbp1" , "Dll1" ,"Lyz1",'Nupr1', 'Dll4', 'Mmp7',
												"Mki67","Fzd5","Fzd7","Dclk1","Rspo3","Wnt4","Wnt5a","Uri1")
			
			# x <- readRDS( file.path(isc_data_dir,"/TE001-seurat-analysis-data.rds") )
			x <- seurat_analysis
			
			tmp <- x@meta.data %>% dplyr::select(cell_id,UMAP_1.GES,UMAP_2.GES,clusters_ges=seurat_clusters,stemness_index.ges,singleR_labels)
			vp_seurat@meta.data <- left_join(vp_seurat@meta.data,tmp,by=c("cell_id"="cell_id"))	%>% as.data.frame()
			
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
			# x <- readRDS( file.path(isc_data_dir,"/TE001-seurat-analysis-data.rds") )
			x <- seurat_analysis
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
		# 	clustering.tibble$cytotrace_score.ges <- df$ct_score		
		# 	
		# 	stopifnot( identical( rownames(vp_seurat@meta.data) , clustering.tibble$sample_id ) )
		# 	vp_seurat@meta.data$cytotrace_score.ges <- df$ct_score		
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
		# 					 	)
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
	
	g <- vp_seurat %>% FeaturePlot(features = c("Lgr4","Lgr5","Atoh1","Olfm4","Krt19",'Dclk1',"Chga") , pt.size = 0.5 , order = TRUE , cols = c("blue","yellow"), blend = FALSE )
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
		dplyr::select(cytotrace_score.ges,stemness_index.pas,HALLMARK_MITOTIC_SPINDLE,HALLMARK_G2M_CHECKPOINT,HALLMARK_WNT_BETA_CATENIN_SIGNALING,HALLMARK_DNA_REPAIR,mt_percent) %>% 
		as.data.frame()
	rownames(df) <- df_annot.df$sample_id
	df_annot_cols_suppl = HeatmapAnnotation( df = df ,
																					 col = list( 
																					 	cytotrace_score.ges = my_color_func_inferno ,
																					 	stemness_index.pas = my_color_func_viridis ,
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
											 cluster_rows = TRUE ,
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
	
	# my_color_vector = rep( cluster_colors , as.integer(table(df_annot.df$cluster_id)) )
	# names(my_color_vector) <- factor(names(my_color_vector),levels = levels(df_annot.df$cluster_id),ordered = TRUE)
	
	df_annot_cols_cluster_score <- HeatmapAnnotation( 
		silhouette_score = anno_barplot(
			# cluster_rel_score = anno_barplot(
			# clustering.tibble$cluster_reliability[sample_ordering_to_print] ,
			df_annot.df$sil_score[sample_ordering_to_print] ,
			bar_width = 1,
			gp = gpar(col = "gray50"),
			# gp = gpar(col = NA, 
			# 					# fill = rep( cluster_colors , as.integer(table(df_annot.df$cluster_id)) ) ),
			# 					fill = rep( brewer.pal(nlevels(df_annot.df$cluster_id),"Set1") , table(df_annot.df$cluster_id) ) ),
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
											col = circlize::colorRamp2(c(-5, 0, 5), c("deepskyblue3", "white", "brown3")),
											# heatmap_legend_param = list(color_bar = "continuous") ,                       
											cluster_rows = FALSE ,
											cluster_columns = FALSE ,
											cluster_row_slices = TRUE,
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
	
	pdf( file.path( reports.dir , "TE001-pas-clustering-heatmap-with-metacell.pdf") , width = 16, height = 20 )
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
	y.umap <- umap::umap( as.matrix(x) , config = custom.settings )
	
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
	
	umap_plot <- ggplot( data = umap.df , aes( UMAP_1.PAS , UMAP_2.PAS , color = cytotrace_score.ges ) ) +
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
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-cytotrace-at-ges.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-pas-clusters.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-singleR-at-ges.pdf") , width = 6 , height = 6 )
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
	
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-ges-clusters.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-ges-umap-with-pas-clusters.pdf") , width = 6 , height = 6 )
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
	
	
	pdf( file = file.path( reports.dir , "TE001-ges-umap-with-ges-clusters.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-ges-umap-with-ges-clusters-facets.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-ges-umap-with-pas-clusters-facets.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-pas-clusters-facets.pdf") , width = 6 , height = 6 )
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
	
	pdf( file = file.path( reports.dir , "TE001-pas-umap-with-ges-clusters-facets.pdf") , width = 6 , height = 6 )
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
	# pdf( file = file.path( reports.dir , "TE001-pas-umap-analysis-with-seurat-clusters-Rspo3.pdf") , width = 6 , height = 6 )
	# 	print(umap_plot)
	# dev.off()		
	
}

# ## PHATE Analysis ----
# print_msg_info(">>> Running PHATE Analysis")
# {
# 	# devtools::install_github("KrishnaswamyLab/phateR")
# 	library(phateR)
# 	
# 	x <- as.matrix(vp_seurat@assays$VIPER@scale.data)
# 	
# 	# library(reticulate)
# 	# reticulate::py_discover_config(required_module = "phate")
# 	# reticulate::import("phate")
# 	
# 	# vpdist <- as.dist( 1-cor(x,method = "spe"))
# 	vpdist <- as.dist( viperSimilarity(x) )
# 	
# 	# run PHATE
# 	x_PHATE <- phate( 
# 		data = as.matrix(vpdist) ,
# 		knn.dist.method = "precomputed" ,
# 		# data = t(x) ,
# 		# knn.dist.method = "euclidean" ,
# 		knn = 15 ,
# 		ndim = 3 ,
# 		t = 100 ,
# 		# t = "auto" ,
# 		npca = 30 ,
# 		# knn.dist.method = "cosine" ,
# 		n.jobs = 5 , verbose = TRUE , seed = 666 )		
# 	# # run PHATE
# 	# x_PHATE <- phate( t(x) ,
# 	# 									knn = knn_n_neighbor ,
# 	# 									ndim = 3 ,
# 	# 									t = 100 ,
# 	# 									npca = 30 ,
# 	# 									# knn.dist.method = "cosine" ,
# 	# 									knn.dist.method = "euclidean" ,
# 	# 									n.jobs = 5 , verbose = TRUE , seed = 666 )
# 	
# 	df <- x_PHATE$embedding %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
# 	stopifnot( identical( df$cell_id , colnames(x) ) )
# 	df$cluster_id <- vp_seurat$seurat_clusters
# 	df$sample_id <- vp_seurat$sample_id
# 	
# 	
# 	# stopifnot(identical(df$cell_id,my_metadata.obj$sample_id))
# 	# df$sc_entropy <- 
# 	# df$stemness_index <- my_metadata.obj$stemness_index
# 	
# 	# df$Lgr5 <- x["Lgr5",]
# 	# df$Top2a <- x["Top2a",]
# 	# df$Krt19 <- x["Krt19",]
# 	
# 	# l <- levels(df$cluster_id)
# 	# cluster_colors <- c(brewer.pal(9,"Set1"),brewer.pal(6,"Set2"))[1:length(l)]
# 	# names(cluster_colors) <- l
# 	
# 	# df_3 <- df %>% pivot_longer( cols = contains("PHATE") , names_to = "phate_coordinates" , values_to = "phate_values" )
# 	
# 	# p <- ggplot( df , aes( x = PHATE2 , y = PHATE3 ,
# 	# 											# color = Lgr5
# 	# 											# fill = Top2a
# 	# 											fill = cluster_id
# 	# ) ) +
# 	# 	# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 	# 	# labs(color="Mpo") +
# 	# 	geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
# 	# 	# scale_fill_viridis() +
# 	# 	# scale_fill_discrete() +
# 	# 	scale_fill_brewer(palette = "Set1") +
# 	# 	# coord_fixed(ratio = 1) +
# 	# 	# facet_wrap(~phate_coordinates) +
# 	# 	theme_light()
# 	#
# 	# pdf( file.path(reports.dir,"TE001-pas-phate-2-3.pdf"))
# 	# 	print(p)
# 	# dev.off()
# 	
# 	p <- ggplot( df , aes(PHATE1, PHATE2 ,
# 												# color = Lgr5
# 												# fill = Top2a
# 												fill = cluster_id
# 	) ) +
# 		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 		# labs(color="Mpo") +
# 		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
# 		# scale_fill_viridis() +
# 		# scale_fill_discrete() +
# 		scale_fill_brewer(palette = "Set1") +
# 		coord_fixed(ratio = 1) +
# 		theme_light()
# 	
# 	pdf( file.path(reports.dir,"epithelium-integrated-pas-phate.pdf"))
# 	print(p)
# 	dev.off()
# 	
# 	p <- ggplot( df , aes(PHATE2, PHATE3 ,
# 												# color = Lgr5
# 												# fill = Top2a
# 												fill = cluster_id
# 	) ) +
# 		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 		# labs(color="Mpo") +
# 		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
# 		# scale_fill_viridis() +
# 		# scale_fill_discrete() +
# 		scale_fill_brewer(palette = "Set1") +
# 		coord_fixed(ratio = 1) +
# 		theme_light()
# 	
# 	pdf( file.path(reports.dir,"epithelium-integrated-pas-phate-2and3axis.pdf"))
# 	print(p)
# 	dev.off()		
# 	
# 	p <- ggplot( df , aes(PHATE1, PHATE2 ,
# 												# color = Lgr5
# 												# fill = Top2a
# 												fill = cluster_id
# 	) ) +
# 		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 		# labs(color="Mpo") +
# 		geom_point(shape = 21, colour = "white",size = 1.5 , stroke = 0.25) +
# 		# scale_fill_viridis() +
# 		# scale_fill_discrete() +
# 		scale_fill_brewer(palette = "Set1") +
# 		coord_fixed(ratio = 1) +
# 		facet_wrap(~sample_id) +
# 		theme_light()
# 	
# 	pdf( file.path(reports.dir,"epithelium-integrated-pas-phate-facet.pdf"))
# 	print(p)
# 	dev.off()	
# 	
# 	print_msg_info(">>> Make 3D Scatter Plot")
# 	{
# 		library(car)
# 		# scatter3d(formula, data)
# 		# scatter3d( x = df$PHATE1, 
# 		# 					 y = df$PHATE2, 
# 		# 					 z = df$PHATE3, 
# 		# 					 groups = df$cluster_id , 
# 		# 					 surface = FALSE ,
# 		# 					 surface.col = RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) ) )
# 		
# 		
# 		library(rgl)
# 		my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) )
# 		df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]			
# 		
# 		plot3d( 
# 			x = df$PHATE1, y = df$PHATE2, z = df$PHATE3,
# 			col = df$cluster_color , 
# 			# type = 's', 
# 			radius = 10,
# 			xlab="PHATE 1", ylab="PHATE 2", zlab="PHATE 3" )
# 		
# 		writeWebGL( filename = file.path(reports.dir,"3dscatter.html") ,  width=1000, height=1000)
# 		
# 	}
# }

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
	
	pdf( file.path( reports.dir , "TE001-umap-with-metacell.pdf") , width = 3 , height = 3)
	print(p)
	dev.off()		
	
	p <- ggplot(df,mapping = aes(x=UMAP_1,y=UMAP_2,fill=cluster_id)  ) +
		geom_point(shape = 21, colour = "black",size = 1.5, stroke = 0.35) +
		# scale_fill_viridis() +
		scale_fill_brewer(palette = "Set1") +
		facet_wrap(~cluster_id,nrow = 2) +
		theme_light()
	
	pdf( file.path( reports.dir , "TE001-umap-with-metacell-facet.pdf") , width = 6 , height = 4)
	print(p)
	dev.off()		
	
}

## Saving object ----
print_msg_warn(">>> Saving Seurat object and metadata object")
{
	seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data.rds") )
	metadata_tibble_analysis_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-subnetworks-one-signature-viper-analysis-with-metacell-metadata.rds") )
	# clustering_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-subnetworks-one-signature-viper-analysis-with-metacell.tsv") )
	# iter_clustering_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-subnetworks-one-signature-iter-clustering-viper-analysis-with-metacell.tsv") )

	attr(clustering.tibble,"timestamp") <- date()
	attr(vp_seurat,"timestamp") <- date()

	my_metadata <- left_join(vpmat_metadata,vp_seurat@meta.data)
	my_metadata$sample_name <- paste0(my_sample_id,"-7CLUST")

	# attr(clustering.tibble,"run_parameters.yaml") <- yaml::as.yaml(my_params.list)
	attr(clustering.tibble,"run_parameters.list") <- my_params.list

	identical( colnames(vp_seurat) , clustering.tibble$sample_id )
	saveRDS( vp_seurat , seurat_viper_analysis_data_filename )
	saveRDS( clustering.tibble , metadata_tibble_analysis_filename )

	x <- clustering.tibble %>% dplyr::select(cell_id=sample_id,contains("gene_"))
	vp_seurat@meta.data = left_join( vp_seurat@meta.data , x , by=c("cell_id"="cell_id") )
	rownames(vp_seurat@meta.data) <- vp_seurat@meta.data$cell_id
	identical( rownames(vp_seurat@meta.data) , colnames(vp_seurat) )

	markers_4_excel.filename <- file.path( isc_data_dir , paste0(my_sample_id,"-subnetworks-one-signature-iter-clustering-viper-analysis-with-metacell.xlsx") )
	x <- vpmat.stouffered %>% round(2) %>% as.data.frame() %>% rownames_to_column("Regulator")
	writexl::write_xlsx( x , markers_4_excel.filename )
	
	## Renaming identities ----
	Idents(vp_seurat) <- vp_seurat@meta.data$iter_cluster_id
	table(Idents(vp_seurat))
	# new.cluster.ids <- c("stem-progenitors","absorptive","secretory")
	# names(new.cluster.ids) <- levels(vp_seurat)
	# vp_seurat <- RenameIdents(vp_seurat, new.cluster.ids)		
	my_gex_seurat@meta.data$viper_clusters <- vp_seurat@meta.data$pas_cluster_id[ match( colnames(vp_seurat) , colnames(my_gex_seurat) ) ]
	Idents(my_gex_seurat) <- my_gex_seurat@meta.data$viper_clusters
	table(Idents(my_gex_seurat))
	# my_gex_seurat <- RenameIdents(my_gex_seurat, new.cluster.ids)		
	# table(Idents(my_gex_seurat))
	
	markers_table <- FindAllMarkers(my_gex_seurat,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	markers_table <- markers_table %>% dplyr::select(gene,everything())
	writexl::write_xlsx( markers_table , 
											 file.path(isc_data_dir , paste0(my_sample_id,"-expression-viper-clusters-7-clusters.xlsx") ) 
	)

	# print_msg_warn(">>> >> Converting Seurat object for Scanpy")
	# {
	# 	library(SeuratDisk)
	# 	my_filename <- "~/Downloads/TE001-viper.h5Seurat"
	# 	SaveH5Seurat(vp_seurat, filename = my_filename , verbose = T , overwrite = TRUE )
	# 	Convert(my_filename, dest = "h5ad",overwrite = TRUE)
	# }
	
	## Saving Metadata Preparation  ----
	print_msg_info(">>> Saving Processed Data in CSV (for Scanpy)")
	{
		stopifnot( identical( rownames(umap.df) , colnames(vpmat) ) )
		stopifnot( identical( rownames(umap.df) , names(stemness_index.pas) ) )
		
		write_vpmat_to_3_csv( vpmat , "~/Clouds/Dropbox/Data/isc/TE001/" , filetag = paste0(my_sample_id,"-viper-with-subnetworks-one-signature") , is_ensembl_to_sym = FALSE )
		write_csv( vpmat_metadata , "~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-with-subnetworks-one-signature-metadata.csv" )
	}

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
	
	umap_plot <- ggplot( data = my_metadata , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = pas_cluster_id ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_brewer(palette = "Set1") +
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(cytotrace_score.ges) , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = cytotrace_score.ges ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_viridis(option = "B") +
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(desc(stemness_index.pas)) , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = stemness_index.pas ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_viridis() +
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	# umap_plot <- ggplot( data = my_metadata %>% arrange(gene_Mki67) , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = gene_Mki67 ) ) +
	# 	# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
	# 	geom_point( alpha = 1 , size = 0.1 ) +
	# 	geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
	# 	# scale_color_viridis() +
	# 	scale_color_gradient(low="white",high="darkorange")	 +
	# 	xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
	# 	ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	# umap_plot <- ggplot( data = my_metadata %>% arrange(S_Phase) , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = S_Phase ) ) +
	# 	# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
	# 	geom_point( alpha = 0.75 , size = 0.25 ) +
	# 	geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
	# 	scale_color_viridis() +
	# 	xlim(min(my_meta.filtered$UMAP_1.PAS.with_seurat)-2,max(my_meta.filtered$UMAP_1.PAS.with_seurat)+2) +
	# 	ylim(min(my_meta.filtered$UMAP_2.PAS.with_seurat)-2,max(my_meta.filtered$UMAP_2.PAS.with_seurat)+2) + 
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
											 aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = EMT ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_gradient2(low = "deepskyblue3",mid = "white",high="brown3",midpoint = 0) + 
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(desc(HALLMARK_WNT_BETA_CATENIN_SIGNALING)) , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_gradient2(low = "deepskyblue3",mid = "white",high="brown3",midpoint = 0) + 
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	umap_plot <- ggplot( data = my_metadata %>% arrange(desc(HALLMARK_WNT_BETA_CATENIN_SIGNALING)) , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat , color = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION ) ) +
		# geom_point( alpha = 0.5 , stroke = 1 , size = 1 ) +
		geom_point( alpha = 0.75 , size = 0.25 ) +
		geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
		scale_color_gradient2(low = "deepskyblue3",mid = "white",high="brown3",midpoint = 0) + 
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	umap_plot <- ggplot( data = my_metadata , aes( UMAP_1.PAS.with_seurat , UMAP_2.PAS.with_seurat ) ) +
		geom_density_2d_filled(contour_var = "ndensity") + 
		scale_fill_viridis_d(option = "B") +
		xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
		ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
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
	
	p <- ggplot( my_metadata , aes(pas_cluster_id, cytotrace_score.ges , 
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
		
		gene_list <- c("Lgr4","Lgr5","Mki67",'Olfm4',"Ung","Sox4",'Tnfrsf19','Dclk1',"Fgfbp1")
		
		for (my_gene in gene_list)
		{
			# my_gene <- "Lgr5"	
			my_gene.df <- gex.mat[my_gene,] %>% as.data.frame()
			colnames(my_gene.df) <- "a_gene_expr"
			my_gene.tibble <- my_gene.df %>% rownames_to_column("cell_id") %>% as_tibble()
			tibble4plot <- left_join( my_metadata , my_gene.tibble , by = c("cell_id"="cell_id") )
			
			df <- tibble4plot %>% arrange(a_gene_expr)
			p <- ggplot( data = df , aes(UMAP_1.PAS.with_seurat, UMAP_2.PAS.with_seurat , color = a_gene_expr ) ) +
				geom_point( alpha = 0.75 , size = 0.25 ) +
				# geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
				# scale_color_gradient2(low = "purple",mid = "white",high="orange",midpoint = 0) + 
				scale_color_viridis() +
				xlim(min(my_metadata$UMAP_1.PAS.with_seurat)-2,max(my_metadata$UMAP_1.PAS.with_seurat)+2) +
				ylim(min(my_metadata$UMAP_2.PAS.with_seurat)-2,max(my_metadata$UMAP_2.PAS.with_seurat)+2) + 
				# facet_wrap( ~sample_name , scales = "free" , as.table = TRUE ) +
				theme_minimal() +
				xlab("UMAP 1") +
				ylab("UMAP 2") +
				guides(color=guide_colorbar(title=my_gene)) +
				# coord_equal() +
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
	print_msg_info(">>> >> >> Selecting only samples with weights > 0.1 for integrations") # It should be 0.25 but we often don't have enough samples with that threshold
	{
		index <- weigths > 0.1
		mat <- mat[,index]
		clusters <- clusters[index]
		weigths <- weigths[index]
	}		
	
	mat.stouffered <- doStoufferFromClusters( mat , clusters , weights = weigths )
	
	writexl::write_xlsx(as.data.frame(mat.stouffered),file.path(isc_data_dir,"TE001-viper-stouffer-7-clusters.xlsx"))
	saveRDS(mat.stouffered,file.path(isc_data_dir,"TE001-viper-stouffer-7-clusters.rds"))
	
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
	cluster_colors <- cluster_colors[ !is.na(names(cluster_colors)) ]
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
