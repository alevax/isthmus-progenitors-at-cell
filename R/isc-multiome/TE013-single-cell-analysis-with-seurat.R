##
# Analysis of TExxx samples with SEURAT
# -------------------------------------
# system.time({source("sources/isc-multiome/TE013-single-cell-analysis-with-seurat.R")})

	source("../vaxtools/R/utils.R")
	source("../vaxtools/R/cross-species-utils.R")

	# init_python_on_laptop()
	
	library(Seurat)
	require(dplyr)
	require(ggplot2)
	require(tidyverse)
	require(umap)
	require(yaml)

	print_msg_info(">>> Setting constants and paramters")
	{
		my_sample_id <- "TE013"
		my_seed <- 42
		knn_n_neighbor <- 29
		my_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/",my_sample_id)
		# my_path <- paste0("/Volumes/ac_lab_archive/av2729/isc/isc-ermanno-single-cell-data/191101_TIMOTHY_ERMANNO_4_MOUSE_10X/",my_sample_id,"/cellranger_analysis_outs/filtered_feature_bc_matrix/")
		# my_path <- paste0("/Volumes/av2729/isc/isc-ermanno-single-cell-data/191101_TIMOTHY_ERMANNO_4_MOUSE_10X/",my_sample_id,"/cellranger_analysis_outs/filtered_feature_bc_matrix/")
		my_path <- file.path(my_data_dir,"/filtered_feature_bc_matrix/")
		cpm_filename <- file.path( my_data_dir , paste0(my_sample_id,"-cpm.rds") )
		seurat_original_data_filename <- file.path( my_data_dir , paste0(my_sample_id,"-seurat-original-data.rds") )
		seurat_analysis_data_filename <- file.path( my_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
		
		isPerformingPathwayAnalysis <- FALSE
		isComputingOptimalNumberOfClusters <- FALSE
		
		param_min_nFeature_RNA <- 1000
		param_max_nFeature_RNA <- 20000
		param_max_nCount_RNA <- 100000
		param_mt_percent_threshold <- 20
		param_min_cells <- 5
		param_min_feats <- 500
		
		sct_variable_n_feats <- 2000
		scr_return_only_var_genes <- TRUE
		sct_do.scale <- FALSE
		sct_do.center <- TRUE
		# sct_graph_name <- "SCT_nn"
		sct_graph_name <- "SCT_snn"
		n_pcs <- 30
		pcs_to_use <- 1:30
		
		# sct_is_to_scale <- FALSE
		clustering_algorithm <- "louvain"
		# clustering_algorithm <- "leiden"
	}
	
	create_workspace( isWorkspaceToClean = FALSE , 
										experiments_dir = "experiments" , 
										run_dir = paste0(my_sample_id,"-gene-expression-analysis"))

	print_msg_info(">>> Setting constants and paramters [in YAML]")
	{
		my_params.list <- list()
		my_params.list$results_dir = reports.dir
		my_params.list$sample_id = my_sample_id
		my_params.list$cpm_filename <- cpm_filename
		my_params.list$seurat_original_data_filename <- seurat_original_data_filename
		my_params.list$seurat_analysis_data_filename <- seurat_analysis_data_filename
		
		my_params.list$seed = my_seed
		my_params.list$knn_n = knn_n_neighbor
		my_params.list$param_min_nFeature_RNA <- param_min_nFeature_RNA
		my_params.list$param_max_nFeature_RNA <- param_max_nFeature_RNA
		my_params.list$param_max_nCount_RNA <- param_max_nCount_RNA
		my_params.list$param_mt_percent_threshold <- param_mt_percent_threshold
		my_params.list$param_min_cells <- param_min_cells
		my_params.list$param_min_feats <- param_min_feats
		
		my_params.list$sct_variable_n_feats <- sct_variable_n_feats
		my_params.list$scr_return_only_var_genes <- scr_return_only_var_genes
		my_params.list$sct_do.scale <- sct_do.scale
		my_params.list$sct_do.center <- sct_do.center
		my_params.list$sct_graph_name <- sct_graph_name
		my_params.list$n_pcs <- n_pcs
		my_params.list$pcs_to_use <- pcs_to_use
		my_params.list$clustering_algorithm <- clustering_algorithm

		my_yaml_file <- file.path(reports.dir,"run-parameters.yaml")
		yaml::write_yaml(my_params.list %>% as.list(),my_yaml_file)
		
	}
	
	print_msg_info(">>> Loading umi counts matrix")
	{
		if ( !file.exists(cpm_filename) )
		{
			print_msg_warn(">>> >>> CPM File NOT Found, reading data ... ")
			
			data <- Read10X( my_path )
			
			print_msg_warn(">>> >>> >> Extracting only gene expression data ... ")
			sc_data <- CreateSeuratObject(data$`Gene Expression` , min.cells = param_min_cells , min.features = param_min_feats ) %>%
				PercentageFeatureSet(pattern = "^mt-", col.name = "mt_percent") %>% 
				subset( mt_percent < param_mt_percent_threshold & nFeature_RNA > param_min_nFeature_RNA & nFeature_RNA < param_max_nFeature_RNA & nCount_RNA < param_max_nCount_RNA)		
			
			# x <- read_csv("~/Clouds/Dropbox/Data/isc/TE013/TE013-doublets-df.csv")
			# sum(x$predicted_doublet == TRUE)
			
			# tmp <- as.matrix(sc_data@assays$RNA@counts)
			# print_msg_warn("*** REMOVING IMMUNE CELLS")
			# {
			# 	print_msg_warn("Removed " , sum( tmp["Ptprc",] > 0 ) , " cells")	
			# 	sc_data <- sc_data %>% subset( subset = Ptprc == 0)			
			# }
			
			require(pryr)
			my_counts <- as.matrix(sc_data@assays$RNA@counts)
			my_counts <- as.matrix(RelativeCounts(my_counts,1e6))
			dim(my_counts)
			index <- rowSums(my_counts) > 0
			my_counts <- my_counts[ index , ]
			print(dim(my_counts))
			# set.seed(my_seed)
			# selected_cells <- sample(colnames(my_counts),size = 1500,replace = FALSE)
			# my_counts <- my_counts[ ,selected_cells ]
			# dim(my_counts)
			object_size(my_counts)
			attr(my_counts,"nFeature_RNA") <- paste0("nFeature_RNA > ", param_min_nFeature_RNA ," & nFeature_RNA < ", param_max_nFeature_RNA )
			attr(my_counts,"mt_percent") <- paste0("mt_percent < ", param_mt_percent_threshold )
			attr(my_counts,"min.features") <- param_min_feats
			attr(my_counts,"min.cells") <- param_min_cells
			
			saveRDS( my_counts , cpm_filename )
			saveRDS( sc_data , seurat_original_data_filename )
			
		} else {
			print_msg_warn(">>> >>> CPM File Found, loading ... ")
			my_counts <- readRDS(cpm_filename)
			sc_data <- readRDS(seurat_original_data_filename)
			
			print(dim(my_counts))
		}
	}
	
	## CytoTrace Analysis ----
	print_msg_info(">>> CytoTrace Analysis")
	{
		source("~/Downloads/CytoTRACE/R/zzz.R")
		source("~/Downloads/CytoTRACE/R/CytoTRACE.R")
		source("~/Downloads/CytoTRACE/R/plotCytoGenes.R")
		source("~/Downloads/CytoTRACE/R/plotCytoTRACE.R")
		
		mat <- as.matrix(sc_data@assays$RNA@counts)
		dim(mat)
		print_msg_info(">>> >> Running CytoTrace  ...")
		cytotrace.data <- CytoTRACE(mat , enableFast = FALSE , ncores = 2 )
		require(pryr)
		object_size(cytotrace.data)
		plotCytoGenes( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-ges-data")) ,numOfGenes = 10)
		
		sc_data$cytotrace_score.ges <- cytotrace.data$CytoTRACE[ match( colnames(sc_data) , names(cytotrace.data$CytoTRACE) ) ]
		
		tmp <- tibble( cell_id = colnames(sc_data) , cytotrace_score = sc_data$cytotrace_score.ges %>% round(4) )
		write_csv(tmp, file.path(my_data_dir,paste0(my_sample_id,"-cytotrace.csv")) )
		
		# df <- cbind(clustering.tibble$,clustering.tibble$UMAP_2.GES)
		# my_clusters <- as.character(my_isc.sdata$seurat_clusters)
		# names(my_clusters) <- names(my_isc.sdata$seurat_clusters)
		# plotCytoTRACE( cytotrace.data , outputDir = file.path(reports.dir,paste0(my_sample_id,"cytotrace-ges-clusters")) , emb = df , phenotype = my_clusters )
		
	}	
	
	## Seurat Analysis ----
	print_msg_info(">>> Seurat Analysis")
	{
		# print_msg_info(">>> Data to Regress Out Cell Cycle")
		# {
		# 	# Read in the expression matrix The first row is a header row, the first column is rownames
		# 	exp.mat <- read.table(file = "~/Downloads/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
		# 	    as.is = TRUE, row.names = 1)
		# 	
		# 	# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
		# 	# segregate this list into markers of G2/M phase and markers of S phase
		# 	s.genes <- cc.genes$s.genes
		# 	g2m.genes <- cc.genes$g2m.genes			
		# }
		
		# sc_data <- sc_data %>% 
		# 	SCTransform( 
		# 		variable.features.n = sct_variable_n_feats ,
		# 							 do.scale = sct_is_to_scale ,
		# 							 do.center = TRUE ,
		# 							 return.only.var.genes = scr_return_only_var_genes ,
		# 							 seed.use = my_seed ) %>% 
		# 	RunPCA(nfeatures.print = 10) 
		
		sc_data <- sc_data %>% 
			SCTransform( variable.features.n = my_params.list$sct_variable_n_feats , 
									 return.only.var.genes = my_params.list$scr_return_only_var_genes ,
									 do.scale = my_params.list$sct_do.scale , do.center = my_params.list$sct_do.center , 
									 seed.use = my_params.list$seed ) %>% 
			RunPCA()
		
		# sc_data$CC.Difference <- sc_data$S.Score - sc_data$G2M.Score
		# sc_data <- ScaleData(sc_data, vars.to.regress = "CC.Difference", features = rownames(sc_data))
		
		## Cell Cycle Scoring Analysis ----
		print_msg_info(">>> Cell Cycle Scoring Analysis")
		{
			# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
			# segregate this list into markers of G2/M phase and markers of S phase
			s.genes <- cc.genes$s.genes %>% human_to_mouse()
			g2m.genes <- cc.genes$g2m.genes	%>% human_to_mouse()
			
			sc_data <- CellCycleScoring(sc_data, 
																	s.features = s.genes, 
																	g2m.features = g2m.genes, 
																	set.ident = FALSE)
		
			my_tibble <- tibble( cell_id = colnames(sc_data) ,
													 s_score = sc_data@meta.data$S.Score ,
													 g2m_score = sc_data@meta.data$G2M.Score ,
													 phase_score = sc_data@meta.data$Phase ,
			)
			
			write_tsv( my_tibble , file.path(my_data_dir,paste0(my_sample_id,"-cell-cycle-scoring-table.csv")) )						
				
		}
	
		# message(">>> Run JackStraw Algorithm")
		# {
			pdf( file.path(reports.dir,"elbow-plot.pdf") )
				print(ElbowPlot(sc_data))
			dev.off()
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
		# pcs_to_use <- 1:10
		
		fs <- sc_data %>% FeatureScatter("nFeature_RNA","mt_percent")
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-mt-vs-detected-genes-scatter-plot.pdf")) )
			print(fs)	
		dev.off()
		
		g <- VlnPlot(object = sc_data, features = c("mt_percent","Lgr4","Lgr5") ) 
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-mt-violin-plot.pdf")) )
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
		
		sc_data <- FindNeighbors( sc_data , 
															# dims = pcs_to_use , 
															assay = "SCT" , 
															verbose = TRUE , 
															k.param = knn_n_neighbor ,
															compute.SNN = TRUE ,
															# annoy.metric = "cosine" ,
															annoy.metric = "euclidean" ,
															# force.recalc = TRUE ,
															# distance.matrix = as.matrix(vpdist) ,
															do.plot = FALSE )		
		
		## Silhouette Analysis ----
		sc_data_dist.matrix <- as.dist(1-cor(as.matrix(sc_data@assays$SCT@scale.data),method = "pea"))
		my_sil.df <- tibble( resolution = seq(0.01,0.5,by = 0.01) , sil_avg = 0 , sil_median = 0 , n_clust = 0 )
		for ( a_res in my_sil.df$resolution )
		{
			message("--- Silhouette score computation for resolution :" , a_res )
			sc_data <- FindClusters( sc_data , 
															 graph.name = my_params.list$sct_graph_name ,
																 # graph.name = "nn" ,
																 # graph.name = "VIPER_snn" ,
																 # graph.name = "VIPER_nn" ,
																 # resolution = resolution_value , 
																 resolution = a_res , 
																 # verbose = TRUE , 
																 modularity.fxn = 1 ,
																algorithm = ifelse( clustering_algorithm == "louvain" , 1 , 4 )
																 # algorithm = 4 , # Leiden
																 # algorithm = 1 , # Louvain
																 # seed.use = my_seed 
			)
			
			if ( nlevels(sc_data$seurat_clusters) == 1 ) next ;
			
			require(cluster)
			s <- silhouette( as.integer(sc_data$seurat_clusters) , sc_data_dist.matrix )
			# pdf( file.path(reports.dir,paste0("sil-res-",a_res,".pdf")) )
			require(factoextra)
			x <- fviz_silhouette(s,print.summary = FALSE)			
			y <- sapply( levels(x$data$cluster) , function(i) mean( x$data$sil_width[ x$data$cluster == i ] ) )
			z <- sapply( levels(x$data$cluster) , function(i) median( x$data$sil_width[ x$data$cluster == i ] ) )
			my_sil.df$sil_avg[ my_sil.df$resolution == a_res ] = mean(y)
			my_sil.df$sil_median[ my_sil.df$resolution == a_res ] = mean(z)
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
			
			pdf( file.path(reports.dir,paste0(my_sample_id,"-optimal-clustering.pdf")) )
				print(p)
			dev.off()
		}
		
		# my_res <- my_sil.df$resolution[ which.max(my_sil.df$sil_median) ]
		my_res <- my_sil.df$resolution[ which.max(my_sil.df$sil_avg) ]
		# my_res <- 0.30
		
		# vp_seurat <- FindClusters(vp_seurat, resolution = 2 , verbose = TRUE,modularity.fxn = 2,algorithm = 4,random.seed = 666 ) 
		sc_data <- FindClusters( sc_data , 
														 assay = "SCT"  ,
														 graph.name = my_params.list$sct_graph_name ,
														 # graph.name = "VIPER_snn" ,
														 resolution = my_res , 
														 verbose = TRUE , 
														 modularity.fxn = 1 , 
														 # algorithm = 4 , # Leiden
														 group.singletons = FALSE ,
														 algorithm = ifelse( clustering_algorithm == "louvain" , 1 , 4 ) ,
														 random.seed = my_seed ) 
		
		table(sc_data$seurat_clusters)		
		
		sc_data <- RunUMAP( sc_data , 
												assay = "SCT" ,
												n.components = 3L ,
												# features = prolif_genes ,
												# graph = my_params.list$sct_graph_name ,
												# graph.name = my_params.list$sct_graph_name ,
												# n.neighbors = knn_n_neighbor , 
												dims = pcs_to_use ,
												reduction = "pca" , 
												# umap.method = "umap-learn" ,
												# annoy.metric = "cosine" ,
												annoy.metric = "euclidean" ,
												verbose = FALSE , seed.use = my_seed )				
		
		p1 <- DimPlot( sc_data , reduction = "umap", 
									 group.by = "seurat_clusters" , 
									 pt.size = 0.5 , label = TRUE )
		pdf( file = file.path(reports.dir, paste0(my_sample_id,"-ges-clusters-umap.pdf")) )
			print(p1)
		dev.off()

		p1 <- DimPlot( sc_data , reduction = "umap", 
									 group.by = "seurat_clusters" , 
									 pt.size = 0.5 , label = TRUE )		
		p2 <- DimPlot( sc_data , reduction = "umap", 
									 group.by = "Phase" , 
									 pt.size = 0.5 , label = TRUE )
		p <- CombinePlots(plots = list(p1, p2) , )
		pdf( file = file.path(reports.dir, paste0(my_sample_id,"-ges-cell-cycle-phase-and-clusters-umap.pdf")) , width = 8 ,  height = 3.5 )
			print(p)
		dev.off()
		
		p1 <- DimPlot( sc_data , reduction = "umap", 
									 group.by = "Phase" ,
									 pt.size = 0.5 , label = TRUE )		
		p2 <- DimPlot( sc_data , reduction = "umap", 
									 split.by = "Phase" ,
									 pt.size = 0.5 , label = TRUE )
		p <- CombinePlots(plots = list(p1, p2) , rel_widths = c(1,2))
		pdf( file = file.path(reports.dir, paste0(my_sample_id,"-ges-cell-cycle-phase-umap.pdf")) , width = 12 ,  height = 3.75 )
		print(p)
		dev.off()		
		
		g <- sc_data %>% FeaturePlot(features = c("Epcam","Sox6","Lgr4","Lgr5","Tnfrsf19","Krt19","Olfm4","Mki67","Ung","Dll1") , 
																 pt.size = 0.25 , order = TRUE , cols = c("blue","yellow"), blend = FALSE ,ncol = 2)
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-feat-plot-1.pdf")),height = 12)
			print(g)	
		dev.off()		
		
		g <- sc_data %>% FeaturePlot(features = c("Lgr4","Lgr5","Tnfrsf19","Ung","Krt19","Olfm4","Atoh1","Dclk1") , 
																 pt.size = 0.25 , order = TRUE , cols = c("blue","yellow"), blend = FALSE ,ncol = 2)
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-sct-feat-plot-2.pdf")))
		print(g)	
		dev.off()				
		
		# Identify the 10 most highly variable genes
		top10 <- head(VariableFeatures(sc_data), 10)
		
		# plot variable features with and without labels
		plot1 <- VariableFeaturePlot(sc_data)
		plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
		p <- CombinePlots(plots = list(plot1, plot2))
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-top-variables-features.pdf")) , width = 16)
		print(p)
		dev.off()		
		
		# markers <- FindAllMarkers(vp_seurat, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.25 , test.use = 't' )
		markers <- FindAllMarkers( sc_data , only.pos = TRUE )#, min.pct = 0.25, logfc.threshold = 0.25,test.use = 't')
		# markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
		# top_markers <- markers %>% group_by(cluster) %>% filter(myAUC > 0.99 & avg_diff > 4) %>% top_n(10,wt = myAUC )
		# top_markers <- markers %>% group_by(cluster) %>% top_n(10,wt = avg_logFC )
		n_top <- 10
		top_markers <- markers %>% group_by(cluster) %>% top_n( n_top , wt = c(-p_val_adj) ) %>% top_n( n_top , wt = avg_log2FC )
		dim(top_markers)
		
		# markers %>% as_tibble() %>% filter(cluster==5) %>% arrange(p_val_adj)
		
		require(scales)
		require(viridis)
		
		# p2 <- DoHeatmap(vp_seurat,features = features_to_show , raster = FALSE ,assay = "VIPER",draw.lines = TRUE,angle = 0,group.bar = TRUE ) + 
		p2 <- DoHeatmap( sc_data , features = top_markers$gene , 
			raster = TRUE , 
			assay = "SCT" , draw.lines = TRUE , angle = 0 , group.bar = TRUE ) +
			scale_fill_viridis() +
			theme( text = element_text(size=7) ) +
			# scale_fill_distiller(palette = "RdBu") + 
			theme(axis.text=element_text(size=5) )
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-clusters-heatmap.pdf")) )
			print(p2)
		dev.off()
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-clusters-combo.pdf")) , width = 14 , height = 6)
			print( cowplot::plot_grid(p1,p2,ncol = 2 ) )
		dev.off()		
		
		source("../vaxtools/R/cross-species-utils.R")
		stm_idx <- getStemnessIndex( as.matrix( sc_data@assays$SCT@scale.data ) %>% mouse_to_human() )
		stopifnot( identical( colnames(sc_data) , names(stm_idx) ) )
		sc_data@meta.data$stemness_index.ges <- stm_idx		
		
		stopifnot( identical( colnames(sc_data) , rownames(sc_data@meta.data) ) )
		df <- as_tibble( sc_data@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_id") )
		df <- df %>% dplyr::rename(UMAP_1.GES=UMAP_1,UMAP_2.GES=UMAP_2,UMAP_3.GES=UMAP_3)
		
		sc_data@meta.data <- sc_data@meta.data %>% rownames_to_column("cell_id")
		tmp <- left_join( sc_data@meta.data , df , by = c("cell_id"="cell_id") )
		sc_data@meta.data <- tmp[ match( sc_data@meta.data$cell_id , tmp$cell_id ) , ]

		sc_data@meta.data <- as.data.frame(sc_data@meta.data)
		rownames(sc_data@meta.data) <- colnames(sc_data)
		
		p <- ggplot( sc_data@meta.data , aes( UMAP_1.GES , UMAP_2.GES , color = stemness_index.ges ) ) +
			geom_point(size=1) +
			scale_color_viridis() +
			theme_light()
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-stemness-index.pdf")) ) #, width = 12 , height = 6)
			print(p)
		dev.off()
		
		p <- ggplot( sc_data@meta.data , aes( UMAP_1.GES , UMAP_2.GES , color = cytotrace_score.ges ) ) +
			geom_point(size=1) +
			scale_color_viridis() +
			theme_light()
		
		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-cytotrace-score.pdf")) ) #, width = 12 , height = 6)
			print(p)
		dev.off()		
		
		print_msg_info(">>> Plotting UMAP with clusters and mitocontrial genes content")
		{
			df <- as_tibble(sc_data@meta.data) %>% dplyr::select(cell_id,mt_percent,UMAP_1.GES,UMAP_2.GES)
			
			x <- sc_data@meta.data
			x %>% group_by(seurat_clusters) %>% summarise(mean_mt=mean(mt_percent),max_mt=max(mt_percent),count=n())
			
			p <- ggplot( df , aes( UMAP_1.GES , UMAP_2.GES , color = mt_percent ) ) +
				geom_point(size=0.5,alpha=0.5) +
				scale_color_viridis() +
				theme_light()
			
			pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-mt-percent-index.pdf")) ) #, width = 12 , height = 6)
			print(p)
			dev.off()
		}
		
		## Running SingleR ----
		print_msg_info(">>> Running SingleR")
		{
			library(SingleR)
			library(celldex)
			bped.se <- BlueprintEncodeData()
			bped.se
			
			my_mat <- as.matrix(sc_data@assays$RNA@counts)
			system.time({
				sR_BP_predicted_labels <- SingleR( test = my_mat %>% mouse_to_human(), 
																					 ref = bped.se, 
																					 # fine.tune = FALSE,prune=T,tune.thresh=0.1,
																					 # genes = "sd" ,
																					 # assay.type.test=1,
																					 labels = bped.se$label.main)
			})
			table(sR_BP_predicted_labels$labels)
			
			stopifnot( identical( colnames(sc_data) , rownames(sc_data@meta.data) ) )
			sc_data@meta.data$singleR_labels <- sR_BP_predicted_labels$labels 
			
			umap.df <- as_tibble(sc_data@meta.data)
			
			library(wesanderson)
			my_palette <- c( wes_palette("Darjeeling1", n = 5) , 
											 wes_palette("Darjeeling2", n = 5) , 
											 wes_palette("FantasticFox1", n = 5) )
			
			n <- length(table(umap.df$labels))
			umap_plot <- ggplot( data = umap.df , aes( UMAP_1.GES , UMAP_2.GES , color = singleR_labels , shape = seurat_clusters ) ) +
				geom_point( alpha = 0.75 , stroke = 0.5 , size = 1 ) +
				scale_color_manual(values = my_palette ) +
				xlab( paste0( "UMAP 1") ) +
				ylab( paste0( "UMAP 2") ) +
				theme_minimal()
			
			# pdf( file = file.path( reports.dir , "sc-untr+clones-pas-umap-analysis.pdf") , width = 16 , height = 12 )
			pdf( file = file.path(reports.dir,paste0(my_sample_id,"-ges-singleR-umap-with-labels.pdf") ) , width = 10 , height = 8 )
			print(umap_plot)
			dev.off()		
			
		}			
		
		saveRDS( sc_data , seurat_analysis_data_filename )
		
	}
		
	# ## Pathway Analysis ----
	# if (isPerformingPathwayAnalysis)
	# {
	# 	##
	# 	# Pathway Analysis
	# 	# ----------------
	# 	message(">>> Pathway Analysis on Hallmarks of Cancer")
	# 	{
	# 		source("../vaxtools/R/pathway-analysis.R")
	# 		
	# 		# x <- readRDS("~/Clouds/Dropbox/Data/isc/TE013-cpm.rds")
	# 		x <- my_counts
	# 		feature_names <- getGeneSymbolsFromEnsemblId(getHumanOrthologousFromMouseEnsemblIds(getEnsemblIdsFromGeneSymbols(rownames(x),"mouse")))
	# 		rownames(x) <- feature_names
	# 		dim(x)
	# 		x <- x[ !is.na(rownames(x)) , ]
	# 		dim(x)
	# 		samples.list <- list(x)
	# 		names(samples.list) <- "vpmat"
	# 		h <- readRDS("../vaxtools/data/MSigDB-gene-sets/msigdb-h-as-regulon.rds")
	# 		
	# 		msigdb.list <- list(h)
	# 		names(msigdb.list) <- c("h")
	# 		
	# 		system.time({
	# 			pathway.results <- list()
	# 			for( msigdb.index in seq_along(msigdb.list) )
	# 			{
	# 				pathway.results[[names(msigdb.list)[[msigdb.index]]]]
	# 				for( i in seq_along(samples.list) )
	# 				{
	# 					filename <- paste0( names(samples.list)[[i]] , "-" , names(msigdb.list)[[msigdb.index]] )
	# 					filename.pdf <- paste0(filename,".pdf")
	# 					filename.xls <- paste0(filename,".xls")
	# 					pathway.results[[ names(msigdb.list)[[msigdb.index]] ]][[i]] <- doPathwayAnalysis( samples = samples.list[[i]] , filename = filename.pdf , 
	# 																																														 msigdb = NULL ,
	# 																																														 regulon = msigdb.list[[msigdb.index]] ,
	# 																																														 save2excel = FALSE , filename.xls = filename.xls , .DEBUG = FALSE )
	# 					# doPathwayAnalysis( samples = samples.list[[i]] , filename = names(samples.list)[[i]] , msigdb = c5 , regulon = genesSetList.as.regulon , threshold = 5 )
	# 				}
	# 			}
	# 		})        
	# 		
	# 		filetag <- "h"
	# 		mat_pathways <- unlist(pathway.results,recursive = F)[[filetag]]
	# 		# filename <- file.path( processed_data.dir , "vpmat-cancer-hallmarks.rds")
	# 		# message("- Saving Clones Hallmark Of Cancer Analysis in " , filename )
	# 		# saveRDS( vpmat_pathways , filename )
	# 		
	# 		
	# 		# g <- init_pheatmap_color_palette(z.limit = 6,z.satur = 3,z.signific = 1.6)							
	# 		# pheatmap(mat_pathways,scale = "row",
	# 		# 				 color = g$color.palette ,
	# 		# 				 breaks = g$palette.breaks )
	# 	}			
	# }
	
	# ## MetaCell Analysis ----
	# message(">>> Metacell dataset generation for ARACNe network inference")
	# {
	# 	source("../single-cell-pipeline/functions/cluster-functions.R")
	# 	source("../single-cell-pipeline/functions/process-utils.R")
	# 	source("../single-cell-pipeline/functions/viper-utils.R")
	# 
	# 	sc_data <- readRDS(  seurat_analysis_data_filename )
	# 
	# 	table(sc_data$seurat_clusters)
	# 
	# 	set.seed(my_seed)
	# 	my_list <- list()
	# 	for ( a_cluster in as.numeric(levels(sc_data$seurat_clusters)) )
	# 	{
	# 		# head( summary(vp_seurat@graphs$VIPER_snn) )
	# 		c_id <- a_cluster
	# 		n_cells <- 20
	# 		# cell_ids <- colnames(vp_seurat)[ colnames(vp_seurat) %in% names(vp_seurat$seurat_clusters)[ vp_seurat$seurat_clusters == c_id] ]
	# 		cell_ids <- colnames(sc_data)[ colnames(sc_data) %in% names(sc_data$seurat_clusters)[ sc_data$seurat_clusters == c_id] ]
	# 		tot_cells <- length(cell_ids)
	# 		tot_cells
	# 		n_metacells <- round(tot_cells/3)
	# 		n_metacells
	# 
	# 		d <- sc_data@graphs$SCT_snn
	# 		dim(d)
	# 		d <- d[cell_ids,cell_ids]
	# 		dim(d)
	# 
	# 		my_metacell.tibble <- tibble( metacell_id = 1:n_metacells , cell_id_list = vector(length(n_metacells),mode = "list") )
	# 		medioid_cells.vector <- sample( cell_ids , size = n_metacells , replace = FALSE )
	# 		for ( i in 1:n_metacells)
	# 		{
	# 			medioid_cell <- medioid_cells.vector[i]
	# 			my_metacell.tibble$cell_id_list[[i]] <- head( sort( d[,medioid_cell] , decreasing = T ) , n_cells )
	# 		}
	# 
	# 		my_metacell.tibble %>% unnest() %>% group_by(metacell_id) %>% summarize( min(cell_id_list) )
	# 
	# 		my_t <- table( unlist( lapply( my_metacell.tibble$cell_id_list , names ) ) )
	# 		sum(my_t)
	# 		head( sort( my_t , decreasing = T ) , 20 )
	# 
	# 		p <- ggplot( as.data.frame(table(my_t)) , aes(x = my_t,y = Freq) ) +
	# 			geom_histogram( stat = "identity" , col="orange", fill="green", alpha=.2 ) +
	# 			theme_light()
	# 
	# 		pdf( file = file.path(reports.dir,paste0(my_sample_id,"-metacell-histogram-cluster-",a_cluster,".pdf")) )
	# 			print(p)
	# 		dev.off()
	# 
	# 		umi_counts <- as.matrix(sc_data@assays$RNA@counts)
	# 		x <- do.call( cbind , lapply( lapply( my_metacell.tibble$cell_id_list , names ) , function(i) rowSums(umi_counts[,i]) ) )
	# 		c_names <- sapply( lapply( my_metacell.tibble$cell_id_list , names ) , '[[' , 1 )
	# 		x.cpm <- as.matrix(RelativeCounts(x,1e4))
	# 		colnames(x.cpm) <- make.unique(c_names)
	# 		# head(colSums(x.cpm))
	# 		my_list[[a_cluster+1]] <- x.cpm
	# 	}
	# 
	# 	metacell_cpm <- do.call( cbind , my_list )
	# 	dim(metacell_cpm)
	# 
	# 	# lapply( as.numeric(levels(vp_seurat$seurat_clusters)) ,
	# 	# 	function(z) saveRDS( my_list[[z+1]] ,
	# 	# 											 file.path( my_data_dir , paste0("cluster-",z,"-cpm.rds" ) ) )
	# 	# 	)
	# 
	# }
	# 
	# sapply( my_list , function(z) dim(z) )



