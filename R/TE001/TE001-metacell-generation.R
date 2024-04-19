##
# Metacell generation for TE001
# -----------------------------

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

library(pryr)
library(Seurat)

create_workspace( isWorkspaceToClean = FALSE , 
									experiments_dir = "experiments" , 
									run_dir = "TE001-metacell-generation")

isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")

print_msg_info(">>> Setting constants and paramters")
{
	my_sample_id <- "TE001"
	my_seed <- 42
	# knn_n_neighbor <- 15
	
	# my_path <- paste0("/Volumes/ac_lab_archive/av2729/isc/isc-ermanno-single-cell-data/191101_TIMOTHY_ERMANNO_4_MOUSE_10X/",my_sample_id,"/cellranger_analysis_outs/filtered_feature_bc_matrix/")
	# my_path <- paste0("/Volumes/av2729/isc/isc-ermanno-single-cell-data/191101_TIMOTHY_ERMANNO_4_MOUSE_10X/",my_sample_id,"/cellranger_analysis_outs/filtered_feature_bc_matrix/")
	my_path <- paste0("~/Clouds/Dropbox/Data/isc/TE001/filtered_feature_bc_matrix/")
	cpm_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-cpm.rds") )
	seurat_original_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-original-data.rds") )
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
}

print_msg_info(">>> Metacell dataset generation for ARACNe network inference")
{
	sc_data <- readRDS(seurat_analysis_data_filename)
	print(object_size(sc_data))
}

# print_msg_info(">>> Pruning BAD cells")
# {
# 	print(dim(sc_data))
# 	print_msg_warn(">>> >> Removing B Cells")	
# 	sc_data <- sc_data %>% subset( subset = singleR_labels != "B_cell" )
# 	print_msg_warn(">>> >> Removing Epithelial Cells")	
# 	index <- as.matrix(sc_data@assays$RNA@counts)["Epcam",] > 0 ; sum(index)
# 	sc_data <- sc_data[,!index]
# 	print(dim(sc_data))
# 	
# 	table( sc_data$seurat_clusters )
# 	
# }

## Metacell Generation ----
print_msg_info(">>> Metacell dataset generation for ARACNe network inference")
{
	source("../single-cell-pipeline/functions/cluster-functions.R")
	source("../single-cell-pipeline/functions/process-utils.R")
	source("../single-cell-pipeline/functions/viper-utils.R")
	
	# d <- sc_data@graphs$VIPER_snn
	mat <- as.matrix(sc_data@assays$SCT@scale.data)
	print(dim(mat))
	
	set.seed(my_seed)
	my_list <- list()
	my_list.counts <- list()
	for ( a_cluster in as.numeric(levels(sc_data$seurat_clusters)) )
	{
		message(">>> Generating metacells for Cluster " , a_cluster )
		# head( summary(sc_data@graphs$VIPER_snn) )
		c_id <- a_cluster
		n_cells <- 20
		cell_ids <- colnames(sc_data)[ colnames(sc_data) %in% names(sc_data$seurat_clusters)[ sc_data$seurat_clusters == c_id] ]
		tot_cells <- length(cell_ids)
		tot_cells
		n_metacells <- round(tot_cells/3)
		n_metacells
		
		# d <- as.matrix(sc_data@graphs$SCT_nn)[,cell_ids]
		d <- as.matrix(sc_data@graphs$SCT_snn)[,cell_ids]
		# d <- cor( mat[,cell_ids] )
		dim(d)
		
		my_metacell.tibble <- tibble( metacell_id = 1:n_metacells , cell_id_list = vector(length(n_metacells),mode = "list") )
		medioid_cells.vector <- sample( cell_ids , size = n_metacells , replace = FALSE )
		for ( i in 1:n_metacells)
		{
			medioid_cell <- medioid_cells.vector[i]
			my_metacell.tibble$cell_id_list[[i]] <- head( sort( d[,medioid_cell] , decreasing = T ) , n_cells )
		}
		
		tmp <- my_metacell.tibble %>% unnest(cell_id_list) %>% group_by(metacell_id) %>% summarize( min_cor_value=min(cell_id_list) , 
																																																max_m1_cor_value=max(cell_id_list[-which.max(cell_id_list)]) ,
																																																mean_m1_cor_value=mean(cell_id_list[-which.max(cell_id_list)]) )
		head(tmp)
		print_msg_info("[ MIN Corr Score Among the lowest per cell | " , round( min(tmp$min_cor_value) , 3 )  , "]")
		print_msg_info("[ MIN Corr Score Among the lowest per cell | " , round( max(tmp$min_cor_value) , 3 )  , "]")
		print_msg_info("[ Quantile | " , round( quantile(tmp$min_cor_value) , 3 )  , "]")
		print_msg_info("[ MIN Corr Score Among the highest per cell | " , round( min(tmp$max_m1_cor_value) , 3 )  , "]")
		print_msg_info("[ MAX Corr Score Among the highest per cell | " , round( max(tmp$max_m1_cor_value) , 3 )  , "]")
		print_msg_info("[ Quantile | " , round( quantile(tmp$max_m1_cor_value) , 3 )  , "]")
		
		quantile(tmp$mean_m1_cor_value)
		
		my_t <- table( unlist( lapply( my_metacell.tibble$cell_id_list , names ) ) )
		print_msg_info("[ Total Cells Used | " , sum(my_t) , "]")
		print_msg_info("[ Unique Cells Used | " , length(unique(names(my_t))) , "]")
		
		head( sort( my_t , decreasing = T ) , 20 )
		
		p <- ggplot( as.data.frame(my_t) , aes(x = Freq) ) +
			geom_histogram( binwidth = 1 , stat = "count" , col="orange", fill="green", alpha=.2 ) +
			# stat_bin( aes(label=..count..), geom="text", vjust=-.5) +
			theme_light()
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-metacell-histogram-cluster-",a_cluster,".pdf")) , height = 4 , width = 2)
			print(p)
		dev.off()
		
		umi_counts <- as.matrix(sc_data@assays$RNA@counts)
		x.counts <- do.call( cbind , lapply( lapply( my_metacell.tibble$cell_id_list , names ) , function(i) rowSums(umi_counts[,i]) ) )
		c_names <- sapply( lapply( my_metacell.tibble$cell_id_list , names ) , '[[' , 1 )
		x.cpm <- as.matrix(RelativeCounts(x.counts,1e6))
		colnames(x.cpm) <- make.unique(c_names)
		colnames(x.counts) <- make.unique(c_names)
		# head(colSums(x.cpm))
		print_msg_warn("*** Adding 1 to clusters to let them start C=1")
		my_list[[a_cluster+1]] <- x.cpm
		my_list.counts[[a_cluster+1]] <- x.counts
	}
	
	metacell_cpm <- do.call( cbind , my_list )
	dim(metacell_cpm)
	sum(duplicated(colnames(metacell_cpm)))
	
	metacell_counts <- do.call( cbind , my_list.counts )
	dim(metacell_counts)
	sum(duplicated(colnames(metacell_counts)))	
	
	x <- apply( metacell_cpm , 2 , function(x) { sum(x > 0) } )
	print_msg_info("[ Metacells Quantile | " , round( quantile(x) , 3 )  , "]")
	y <- apply( as.matrix(sc_data@assays$RNA@counts) , 2 , function(x) { sum(x > 0) } )
	print_msg_info("[ Single Cells Quantile | " , round( quantile(y) , 3 )  , "]")
	
	sapply( my_list , function(z) dim(z) )
	
	print_msg_warn("*** Adding 1 to clusters to let them start C=1")
	lapply( as.numeric(levels(sc_data$seurat_clusters))+1 ,
					function(z) saveRDS( my_list[[z]] ,
															 file.path( processed_data.dir , paste0(my_sample_id,"-cluster-",z,"-cpm.rds" ) ) ) )
															 # paste0("~/Clouds/OneDrive - cumc.columbia.edu/cumc.columbia.edu/Malagola, Ermanno - wang-califano-isc-collaboration/data/aracne-run-TE001-metacell-4-clusters-04-02-2020/cluster-",z,"-cpm.rds" ) ) )

	saveRDS( metacell_cpm , file.path( processed_data.dir , paste0(my_sample_id,"-all-metacells-cpm.rds" ) ) )
	
}

# ## Metacell Generation PISCES ----
# print_msg_info(">>> Metacell dataset generation for ARACNe network inference")
# {
# 	# source("../single-cell-pipeline/functions/cluster-functions.R")
# 	# source("../single-cell-pipeline/functions/process-utils.R")
# 	# source("../single-cell-pipeline/functions/viper-utils.R")
# 	source("../PISCES/R/metacell_funcs.R")
# 
# 	table(sc_data$seurat_clusters)
# 
# 	metacell_cpm <- MetaCells( counts.mat = as.matrix(sc_data@assays$RNA@counts) ,
# 														 wReplace = TRUE ,
# 														 dist.mat = as.dist(as.matrix(sc_data@graphs$SCT_snn)) ,
# 														 clust.vect = sc_data$seurat_clusters ,
# 														 subset = 100 ,
# 														 min.samps = 20 ,
# 														 num.neighbors = 5 )
# 	str(metacell_cpm,1)
# 	metacell_cpm <- do.call( cbind , metacell_cpm )
# 	str(metacell_cpm,1)
# 	colnames(metacell_cpm) <- make.names(colnames(metacell_cpm))
# 	colnames(metacell_cpm) <- make.unique(colnames(metacell_cpm))
# 	colnames(metacell_cpm) <- gsub("\\.","-",colnames(metacell_cpm))
# 	identical( rownames(metacell_cpm[[1]]) , rownames(metacell_cpm[[2]]) )	
# }

## Analysis with Seurat of Metacell Dataset ----
print_msg_info(">>> Analysis with Seurat of Metacell Dataset")
{
	sc_data.metacell <- CreateSeuratObject(metacell_counts , min.cells = param_min_cells , min.features = param_min_feats ) %>%
		PercentageFeatureSet(pattern = "^mt-", col.name = "mt_percent") %>% 
		subset( mt_percent < param_mt_percent_threshold & nFeature_RNA > param_min_nFeature_RNA ) # & nFeature_RNA < param_max_nFeature_RNA )
	
	dim(sc_data.metacell)
	
	sc_data.metacell <- sc_data.metacell %>% 
		SCTransform( variable.features.n = 2000 , 
								 do.scale = FALSE , 
								 do.center = TRUE , 
								 seed.use = my_seed ) %>% 
		RunPCA(nfeatures.print = 10) 
	
	# sc_data$CC.Difference <- sc_data$S.Score - sc_data$G2M.Score
	# sc_data <- ScaleData(sc_data, vars.to.regress = "CC.Difference", features = rownames(sc_data))
	
	# n_pcs <- 30	
	pdf( file.path(reports.dir,"elbow-plot-metacell.pdf") )
		print(ElbowPlot(sc_data.metacell))
	dev.off()
	
	# n_pcs <- 30	
	# message(">>> Run JackStraw Algorithm")
	# {
	# 	pdf( file.path(reports.dir,"elbow-plot.pdf") )
	# 	print(ElbowPlot(sc_data))
	# 	dev.off()
	# 	sc_data <- JackStraw( sc_data , assay = "RNA" , reduction = "pca",dims = n_pcs,verbose = TRUE)
	# 	# head(JS(object = sc_data[['pca']], slot = 'empirical'))	
	# 	sc_data <- ScoreJackStraw( sc_data , assay = "RNA" , reduction = "pca",dims = 1:n_pcs)
	# 	pdf( file.path(reports.dir,"jackstraw-plot.pdf") )
	# 	print( JackStrawPlot( sc_data , dims = 1:n_pcs , reduction = "pca", xmax = 0.1, ymax = 0.3) )
	# 	dev.off()
	# 	
	# 	message(">>> Selecting only PC that are statistically significant using the JackStraw algorithm")
	# 	pca_significant_compounents <- sc_data@reductions$pca@jackstraw@overall.p.values
	# 	pcs_to_use <- pca_significant_compounents[,1][ pca_significant_compounents[,2] < 0.05 ]
	# 	
	# 	x <- max(pcs_to_use)
	# 	x <- ifelse( x %% 2 == 0 , x , x+1 )
	# 	pcs_to_use <- 1:x
	# 	# pcs_to_use <- 1:6
	# }
	# pcs_to_use <- 1:6
	
	fs <- sc_data.metacell %>% FeatureScatter("nFeature_RNA","mt_percent")
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-mt-vs-detected-genes-scatter-plot-metacell.pdf")) )
		print(fs)	
	dev.off()
	
	g <- VlnPlot(object = sc_data.metacell, features = c("mt_percent","Lgr4","Lgr5") ) 
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-mt-violin-plot-metacell.pdf")) )
		print(g)	
	dev.off()
	
	g <- sc_data.metacell %>% FeaturePlot(features = c("Epcam","Sox6","Jam3","Jam2","Cd81","Cd38","Dll1") , pt.size = 0.5 , order = TRUE , cols = c("blue","yellow"), blend = FALSE )
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-feat-plot-metacell.pdf")) ,width = 12)
		print(g)	
	dev.off()	
	
	# prolif_genes <- c(s.genes %>% human_to_mouse() , g2m.genes %>% human_to_mouse())
	# prolif_genes <- prolif_genes[ prolif_genes %in% rownames(sc_data@assays$SCT) ]
	# sc_data <- CellCycleScoring(sc_data, s.features = s.genes %>% human_to_mouse() , g2m.features = g2m.genes %>% human_to_mouse(), set.ident = TRUE)
	# RidgePlot(sc_data, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)		
	
	# sc_data <- RunUMAP( sc_data , 
	# 	assay = "SCT" ,
	# 	features = prolif_genes ,
	# 	# assay = "RNA" , 
	# 	# graph = "nn" ,
	# 	# graph = "SCT_snn" ,
	# 	# n.neighbors = knn_n_neighbor , 
	# 	# dims = pcs_to_use ,
	# 	reduction = "pca" , umap.method = "umap-learn" ,
	# 	metric = "cosine" , 
	# 	verbose = FALSE , seed.use = my_seed )	
	# 
	# p1 <- DimPlot( sc_data , reduction = "umap", pt.size = 0.5 , label = TRUE )
	# pdf( file = file.path(reports.dir, paste0(my_sample_id,"-ges-cell-cycle-umap.pdf")) )
	# 	print(p1)
	# dev.off()		
	
	sc_data.metacell <- FindNeighbors( sc_data.metacell , 
														dims = pcs_to_use ,
														assay = "SCT" , 
														verbose = TRUE , 
														k.param = knn_n_neighbor ,
														compute.SNN = TRUE ,
														# annoy.metric = "cosine" ,
														annoy.metric = "euclidean" ,
														force.recalc = TRUE ,
														# distance.matrix = as.matrix(vpdist) ,
														do.plot = FALSE )		
	
	## Silhouette Analysis ----
	sc_data_dist.matrix <- as.dist(1-cor(as.matrix(sc_data.metacell@assays$SCT@scale.data),method = "pea"))
	my_sil.df <- tibble( resolution = seq(0.01,0.5,by = 0.01) , sil_avg = 0 , sil_median = 0 , n_clust = 0 )
	for ( a_res in my_sil.df$resolution )
	{
		message("--- Silhouette score computation for resolution :" , a_res )
		sc_data.metacell <- FindClusters( sc_data.metacell , 
														 graph.name = "SCT_snn" ,
														 # graph.name = "nn" ,
														 # graph.name = "VIPER_snn" ,
														 # graph.name = "VIPER_nn" ,
														 # resolution = resolution_value , 
														 resolution = a_res , 
														 # verbose = TRUE , 
														 modularity.fxn = 1 ,
														 # algorithm = 4 , # Leiden
														 algorithm = 1 , # Louvain
														 # seed.use = my_seed 
		)
		
		if ( nlevels(sc_data.metacell$seurat_clusters) == 1 ) next ;
		
		require(cluster)
		s <- silhouette( as.integer(sc_data.metacell$seurat_clusters) , sc_data_dist.matrix )
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
	
	print_msg_info(">>> >> Printing Optimal Clustering Plot")
	{
		library(ggrepel)
		# p <- ggplot( my_sil.df , aes( x = resolution , y = sil_median ) ) +
			p <- ggplot( my_sil.df , aes( x = resolution , y = sil_avg ) ) +
			geom_line() +
			# geom_text( aes(x = resolution , y = 0.1 , label = n_clust) , label.size = 3 ) +
			# geom_text( aes(x = resolution , y = sil_median+0.01 , label = n_clust) , label.size = 3 ) +
			geom_text_repel( aes( label = n_clust ) , segment.color = "gray" ) +
			theme_light()
		
		pdf( file.path(reports.dir,"optimal-clustering-metacell.pdf"))
			print(p)
		dev.off()
	}
	
	my_res <- my_sil.df$resolution[ which.max(my_sil.df$sil_avg) ]
	# my_res <- my_sil.df$resolution[ which.max(my_sil.df$sil_median) ]
	# my_res <- 0.26
	print_msg_info(">>> Using Resolution: " , my_res )
	
	# vp_seurat <- FindClusters(vp_seurat, resolution = 2 , verbose = TRUE,modularity.fxn = 2,algorithm = 4,random.seed = 666 ) 
	sc_data.metacell <- FindClusters( sc_data.metacell , 
													 assay = "SCT"  ,
													 graph.name = "SCT_snn" ,
													 # graph.name = "VIPER_snn" ,
													 resolution = my_res , 
													 verbose = TRUE , 
													 modularity.fxn = 1 , 
													 algorithm = 4 , # Leiden
													 group.singletons = TRUE ,
													 # algorithm = 1 , # Louvain
													 random.seed = my_seed ) 
	
	table(sc_data.metacell$seurat_clusters)		
	
	## Run UMAP ----
	sc_data.metacell <- RunUMAP( sc_data.metacell , 
											assay = "SCT" ,
											# features = prolif_genes ,
											# assay = "RNA" , 
											# graph = "nn" ,
											# graph = "SCT_snn" ,
											# n.neighbors = knn_n_neighbor , 
											dims = pcs_to_use ,
											reduction = "pca" , 
											# umap.method = "umap-learn" ,
											# metric = "cosine" , 
											metric = "euclidean" , 
											verbose = FALSE , seed.use = my_seed )		
	# sc_data <- RunUMAP( sc_data , 
	# 										assay = "SCT" ,
	# 										# features = prolif_genes ,
	# 										# assay = "RNA" , 
	# 										# graph = "nn" ,
	# 										graph = "SCT_snn" ,
	# 										# n.neighbors = knn_n_neighbor , 
	# 										# dims = pcs_to_use ,
	# 										reduction = "pca" , 
	# 										# umap.method = "umap-learn" ,
	# 										metric = "cosine" , 
	# 										verbose = FALSE , seed.use = my_seed )				
	
	p1 <- DimPlot( sc_data.metacell , reduction = "umap", pt.size = 0.5 , label = TRUE )
	pdf( file = file.path(reports.dir, paste0(my_sample_id,"-ges-clusters-umap-metacell.pdf")) )
		print(p1)
	dev.off()			
	
	# Identify the 10 most highly variable genes
	top10 <- head(VariableFeatures(sc_data.metacell), 10)
	
	# plot variable features with and without labels
	plot1 <- VariableFeaturePlot(sc_data.metacell)
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	p <- CombinePlots(plots = list(plot1, plot2))
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-top-variables-features-metacell.pdf")) , width = 16)
		print(p)
	dev.off()		
	
	# markers <- FindAllMarkers(vp_seurat, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.25 , test.use = 't' )
	# markers <- FindAllMarkers( sc_data.metacell , only.pos = FALSE , test.use = "negbinom" )#, min.pct = 0.25, logfc.threshold = 0.25,test.use = 't')
	markers <- FindAllMarkers( sc_data.metacell , only.pos = TRUE , assay = "RNA" )#, min.pct = 0.25, logfc.threshold = 0.25,test.use = 't')
	# markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
	# top_markers <- markers %>% group_by(cluster) %>% filter(myAUC > 0.99 & avg_diff > 4) %>% top_n(10,wt = myAUC )
	# top_markers <- markers %>% group_by(cluster) %>% top_n(10,wt = avg_logFC )
	n_top <- 20
	top_markers <- markers %>% group_by(cluster) %>% top_n( n_top , wt = c(-p_val_adj) ) %>% top_n( n_top , wt = avg_log2FC )
	dim(top_markers)
	
	# markers %>% as_tibble() %>% filter(cluster==5) %>% arrange(p_val_adj)
	
	require(scales)
	require(viridis)
	
	# p2 <- DoHeatmap(vp_seurat,features = features_to_show , raster = FALSE ,assay = "VIPER",draw.lines = TRUE,angle = 0,group.bar = TRUE ) + 
	p2 <- DoHeatmap( sc_data.metacell , features = top_markers$gene , 
									 raster = TRUE , 
									 assay = "SCT" ,
									 draw.lines = TRUE , angle = 0 , group.bar = TRUE ) +
		scale_fill_viridis() +
		theme( text = element_text(size=7) ) +
		theme(axis.text=element_text(size=5) )
	# scale_fill_distiller(palette = "RdBu") + 
	

	
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-clusters-heatmap-metacell.pdf")) )
		print(p2)
	dev.off()
	
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-clusters-combo-metacell.pdf")) , width = 14 , height = 6)
		print( cowplot::plot_grid(p1,p2,ncol = 2 ) )
	dev.off()		
	
	source("../vaxtools/R/cross-species-utils.R")
	stm_idx <- getStemnessIndex( as.matrix( sc_data.metacell@assays$SCT@data ) %>% mouse_to_human() )
	
	df <- as_tibble( sc_data.metacell@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_id") )
	stopifnot( identical ( df$cell_id , names(stm_idx) ) )
	df$stemness_index <- stm_idx
	
	p <- ggplot( df , aes( UMAP_1 , UMAP_2 , color = stemness_index ) ) +
		geom_point(size=1) +
		scale_color_viridis() +
		theme_light()
	
	pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-stemmness-index-metacell.pdf")) ) #, width = 12 , height = 6)
		print(p)
	dev.off()
	
	## Running SingleR ----
	print_msg_info(">>> Running SingleR")
	{
		# devtools::install_github('dviraran/SingleR')
		require(SingleR)
		
		hpca.se <- HumanPrimaryCellAtlasData()
		# mrsd <- MouseRNAseqData()
		
		pred.isc <- SingleR( test = as.matrix(sc_data.metacell@assays$RNA@data) %>% mouse_to_human() , 
												 # ref = mrsd , labels = mrsd$label.main
												 ref = hpca.se , labels = hpca.se$label.main ,
												 quantile = 0.9 , tune.thresh = 0.01
		)
		pred.isc
		
		table(pred.isc$labels)
		
		p <- plotScoreHeatmap(pred.isc)
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-singleR-metacell.pdf") ) , width = 64 , height = 6)
			print(p)
		dev.off()
		
		umap.df <- tibble(
			cell_id = colnames(sc_data.metacell) ,
			V1 = sc_data.metacell@reductions$umap@cell.embeddings[,1] ,
			V2 = sc_data.metacell@reductions$umap@cell.embeddings[,2] ,
			labels = pred.isc$labels ,
			cluster_id = sc_data.metacell$seurat_clusters
		)
		
		library(wesanderson)
		my_palette <- c( wes_palette("Darjeeling1", n = 5) , 
										 wes_palette("Darjeeling2", n = 5) , 
										 wes_palette("FantasticFox1", n = 5) )
		
		n <- length(table(umap.df$labels))
		umap_plot <- ggplot( data = umap.df , aes( V1 , V2 , color = labels , shape = cluster_id ) ) +
			geom_point( alpha = 0.75 , stroke = 0.5 , size = 1 ) +
			scale_color_manual(values = my_palette ) +
			xlab( paste0( "UMAP 1") ) +
			ylab( paste0( "UMAP 2") ) +
			theme_minimal()
		
		stopifnot( identical( umap.df$cell_id , colnames(sc_data.metacell) ) ) # This is really redoundant
		sc_data.metacell@meta.data <- cbind( sc_data.metacell@meta.data , umap.df )
		
		# pdf( file = file.path( reports.dir , "sc-untr+clones-pas-umap-analysis.pdf") , width = 16 , height = 12 )
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-singleR-umap-with-labels-metacell.pdf") ) , width = 10 , height = 8 )
			print(umap_plot)
		dev.off()		
		
	}	
	
	# saveRDS( sc_data , seurat_analysis_data_filename )
	
}





