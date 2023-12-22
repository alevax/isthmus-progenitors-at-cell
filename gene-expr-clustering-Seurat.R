# Luca Zanella 
# 12.21.2023

# Gene expression based clustering
# INPUTS (MANDATORY): 
# DataFolder: path to data are stored (scratch)


rm(list = ls()) # cleans workspace
gc()
cat("\014")

# Set paths
cat("YOU MUST SET DataFolder and decide whether to run clustering or not (by using predetermined values)..\n")
DataFolder <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001/filtered_feature_bc_matrix/" # "M:/shares/hcmi/" #
metadata_path <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001-h5ad/"
figures_path = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/TE001_gExpr_clustering/"


running_clustering <- FALSE # do not re-run clustering
cat("'runnig_clustering' is set to: ", running_clustering)

# optimal clustering parameters previously determined
# step 1
opt_param <- c(0.03013859, 22) # res NN 
# step 2 (subclustering)
optim_par_list <- matrix(c(0.05655232, 7, 0.0876204, 8, 0.06659604, 7),byrow=TRUE,nrow=3,ncol=2)
colnames(optim_par_list) <- c("res", "NN")
rownames(optim_par_list) <- c("1", "3", "2")


cytotrace_markers <- c('Smarca5','Rbbp7','Tcerg1','Hnrnpd','Hmg20b','Nelfe','Ube2i','Etv5','Ubn1','Mbd3','Dek','Maz',
                       'Itgb3bp','Ilf2','Pa2g4', 'Id3','Hnf4g','Atoh1','Spdef','Neurod1')
  

# Filtering parameters
{
  min_cells <- 5
  min_feats <- 1000
  
  min_nFeature_RNA <- 1500
  max_nFeature_RNA <- 20000
  max_nCount_RNA <- 100000
  mt_percent_threshold <- 10
}


# Libraries and dependencies
{
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(matrixStats)))
  suppressWarnings(suppressMessages(library(Seurat)))
  suppressWarnings(suppressMessages(library(tibble)))
  suppressWarnings(suppressMessages(library(acdc)))
  suppressWarnings(suppressMessages(library(cluster)))
  suppressWarnings(suppressMessages(library(ggplot2)))
  suppressWarnings(suppressMessages(library(factoextra)))
  suppressWarnings(suppressMessages(library(ComplexHeatmap)))
  suppressWarnings(suppressMessages(library(RColorBrewer)))
  suppressWarnings(suppressMessages(library(circlize)))
  suppressWarnings(suppressMessages(library(gridExtra)))
  
}


# Loading data and metadata
{
  
  TE001_data <- Read10X(file.path(data.dir=file.path(DataFolder)))
  
  
  metadata_csv <- file.path(metadata_path, "TE001-metadata-umap-and-clusters-for-paper.csv")
  metadata <- read.csv(metadata_csv)
  
  TE001 <- CreateSeuratObject(counts=TE001_data,
                              project="SeuratObject",
                              assay="RNA",
                              min.cells = min_cells,
                              min.features = min_feats)
  
  TE001@meta.data$cell_id <- rownames(TE001@meta.data) 
  
  TE001@meta.data <- TE001@meta.data %>% dplyr::select(.,cell_id, everything())
  
  # merge metadata with stored metadata
  TE001@meta.data <- TE001@meta.data %>% dplyr::left_join(., metadata, 
                                                          by = c("cell_id" = "cell_id"))
  rownames(TE001@meta.data) <- TE001@meta.data$cell_id
  
}


# Filtering and Quality Control
{
  TE001[["percent.mt"]] <- PercentageFeatureSet(TE001, pattern="^mt")
  
  # Visualize QC metrics as a violin plot BEFORE filtering
  VlnPlot(TE001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
  
  TE001 <- subset(TE001, subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & nCount_RNA < max_nCount_RNA & percent.mt < mt_percent_threshold)
  
 
  # Visualize QC metrics as a violin plot AFTER filtering
  VlnPlot(TE001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(TE001, feature1 = "nCount_RNA", feature2 = "percent.mt", col="darkgray")
  plot2 <- FeatureScatter(TE001, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", col="darkgray")
  grid.arrange(plot1, plot2, ncol = 2)
  
  dev.off()
  
}

# Basic Data Normalization, Feature Selection and Clustering with default Seurat parameters 
{
  # Normalize Data
  TE001 <- NormalizeData(TE001, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Finding Variable features
  TE001 <- FindVariableFeatures(TE001, selection.method = "vst", nfeatures = 2000)
  top10_hvg <- head(VariableFeatures(TE001), 10)
  plot1 <- VariableFeaturePlot(TE001)
  plot2 <- LabelPoints(plot = plot1, points = top10_hvg, repel = TRUE, xnudge=0, ynudge=0)
  plot2
  dev.off()
  
  # Scaling the data
  all.genes <- rownames(TE001)
  TE001 <- ScaleData(TE001, features = all.genes)
  
  # Perform linear dimensionality reduction
  TE001 <- RunPCA(TE001, features = VariableFeatures(object = TE001))

  # Cluster the cells with default parameters
  TE001 <- FindNeighbors(TE001, dims = 1:10)
  TE001 <- FindClusters(TE001, resolution = 0.5)
  
  levels(TE001@meta.data$seurat_clusters) <- as.character(as.numeric(levels(TE001@meta.data$seurat_clusters)) + 1)

  TE001 <- RunUMAP(TE001, dims = 1:10)
  
  pdf(file.path(figures_path, "umap_clusters_default.pdf") ,width = 12, height = 6)
  DimPlot(TE001, reduction = "umap", pt.size=1.5, group.by = c("seurat_clusters","iter_cluster_id_with_paneth"))
  dev.off()
  
  pdf(file.path(figures_path, "umap_default_cytotrace.pdf") , width = 6 , height = 6)
  FeaturePlot(TE001, features = "cytotrace_te001", pt.size=1.5) & scale_color_viridis_c()
  dev.off()
  
  # evaluate clustering solution
  PCs<-as.data.frame(TE001$pca@cell.embeddings[,1:10])
  euc_distance<-dist(PCs, method="euclidean")
  
  
  s <- silhouette( as.integer(TE001$seurat_clusters) , euc_distance )
  pdf(file.path(figures_path, "silhouette_default.pdf") , width = 6 , height = 8)
  fviz_silhouette(s,print.summary = FALSE)	
  dev.off()
  
  # Find differentially expressed genes
  cat("Perform differential gene expression analysis on the identified clusters with default Seurat pipeline..\n")
  Idents(TE001) <- "seurat_clusters"
    
  # FindAllMarkers with default parameters
  TE001_markers <- FindAllMarkers(TE001, 
                          assay = "RNA",
                          slot="data",
                          test.use = "wilcox", 
                          only.pos = TRUE)
    
    
  top10_markers <- TE001_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  genes_default <- top10_markers$gene
  
  #######################################
  ####### Heatmap for DE analysis #######
  #######################################
  
  # Visualize clusters and top DE genes on a heatmap
  mat_to_show_default <- as.matrix(TE001@assays$RNA@data[genes_default,])
  Clusters <- TE001$seurat_clusters
  colors <- rainbow(nlevels(Clusters)) 
  names(colors) <- levels(Clusters)
  
  # prepare df_annot.df 
  per.cell.sil.width <- s[,"sil_width"]
  names(per.cell.sil.width) <- colnames(TE001)
  
  df_annot.df <- TE001@meta.data %>%
    dplyr::select(cluster_id=seurat_clusters, cytotrace_score = cytotrace_te001) %>%
    tibble::rownames_to_column(.,var="cell_id")
  df_annot.df$cluster_id <- factor(as.numeric(df_annot.df$cluster_id))
   
  x <- as.data.frame(per.cell.sil.width) %>%
    tibble::rownames_to_column(.,var="sample_id") 
  
  df_annot.df <- left_join(df_annot.df,x,by=c("cell_id"="sample_id"))
  
  sample_ordering_to_print <- rev( order( df_annot.df$per.cell.sil.width ) )
  df_annot.df <- df_annot.df[sample_ordering_to_print,]
  
  mat_to_show_default <- mat_to_show_default[,match(df_annot.df$cell_id,colnames(mat_to_show_default))]
  
  stopifnot(identical( df_annot.df$cell_id, colnames(mat_to_show_default) ))
  
  # the silhouette score
  sil.fill <- rep(NA, nrow(df_annot.df))
  unique_clusters <- unique(df_annot.df$cluster_id)
  n_colors <- length(colors)
  
  color_index <- 1
  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]
    sil.fill[df_annot.df$cluster_id == cluster] <- colors[names(colors) == cluster]
    color_index <- i + 1  # Use modulo to cycle through colors
  }
  
  df_annot_cols_cluster_score <- HeatmapAnnotation( 
    silhouette_score = anno_barplot(
      df_annot.df$per.cell.sil.width,
      bar_width = 1,
      gp = gpar(col = NA, 
                fontsize = 6 ,
                cex = 0.8 ,
                fill =  sil.fill 
      ),
      border = FALSE,
      height = unit(1.25, "cm") ,
      width = unit(0.5, "cm")
    )
  )
  
  
  # heatmap annotation
  df_annot.df <- df_annot.df %>% dplyr::select(cluster_id)
  cluster_colors <- colors[ !is.na(names(colors)) ]
  df_annot_cols = HeatmapAnnotation( df = df_annot.df,
                                     col = list(cluster_id=cluster_colors)  , 
                                     show_legend = TRUE , na_col = "gray95",
                                     annotation_name_gp = gpar(fontsize=0.1)
                                    )
  
  top10_default_hm <- Heatmap( mat_to_show_default,
                       col = circlize::colorRamp2(range(mat_to_show_default)/2, c("darkslateblue", "yellow")),
                       use_raster = FALSE ,
                       column_split = df_annot.df$cluster_id,
                       heatmap_legend_param = list(color_bar = "continuous") ,	
                       # show_heatmap_legend = FALSE , 
                       # clustering_method_columns = "ward.D2" ,			
                       cluster_columns = FALSE ,
                       cluster_rows = FALSE ,
                       show_column_names = FALSE ,
                       # cluster_columns = columns_dend ,
                       name = "Gene Expression (log1p)" ,  
                       row_names_gp = gpar(fontsize = 10) , column_names_gp = gpar(fontsize = 8) ,
                       row_title_side = "left" , row_title_gp = gpar(fontsize = 10, fontface = "bold")
                      )
  
  
  
  ht_top10_list = df_annot_cols_cluster_score %v% df_annot_cols %v% top10_default_hm
  
  pdf(file.path(figures_path,"top_10_default_heatmap.pdf"), width = 10, height = 15)  # Adjust width and height as needed
  p <- draw(ht_top10_list, column_km = 1 , 
            heatmap_legend_side = "right", 
            annotation_legend_side = "bottom" )
  dev.off()  

}

# DotPlot - default solution 
{
  top5_markers <- TE001_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  genes_dp_default <- top5_markers$gene %>% unique()
  
  # intialize lists containing dataframes 
  Avgs_default <- vector( mode="list", length=length(levels(TE001$seurat_clusters)) )
  names(Avgs_default) <- levels(TE001$seurat_clusters)
  
  Pct_default <- vector( mode="list", length=length(levels(TE001$seurat_clusters)) )
  names(Pct_default) <- levels(TE001$seurat_clusters)
  
  for (cluster in levels(TE001$seurat_clusters)){
    
    cells_in_cluster <- TE001$seurat_clusters == cluster
    
    # Avg_default
    Avgs_default[[cluster]] <- TE001[["RNA"]]@data[genes_dp_default, cells_in_cluster] %>%
      rowMeans() %>% 
      as.data.frame() %>% 
      rownames_to_column() %>%
      dplyr::rename(.,"Gene"=rowname,"mean"=.) %>% 
      dplyr::mutate(cluster=cluster)
    
    
    # Pct
    Pct_default[[cluster]] <- TE001[["RNA"]]@data[genes_dp_default, cells_in_cluster] %>%
      apply(.,1, function(x){ sum(x>0)/length(x)*100 }) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      dplyr::rename(.,"Gene"=rowname,"Pct"=.) %>%
      dplyr::mutate(cluster=cluster)
    
  }
  
  # Mean
  df_mean.GEx_default <- do.call(rbind, Avgs_default)
  rownames(df_mean.GEx_default) <- NULL
  # Pct
  df_pct.GEx_default <- do.call(rbind, Pct_default)
  rownames(df_pct.GEx_default) <- NULL
  
  
  # Summary
  summary.GEx_default <- cbind(df_mean.GEx_default[,1:2],df_pct.GEx_default[,2:3]) 
  summary.GEx_default$Gene <- factor(summary.GEx_default$Gene, levels=genes_dp_default, ordered = TRUE) 
  
  
  
  # Printing dotplot for genes
  cat("Printing dotplots..\n")
  GEx.DGE_DotPlot_default <- ggplot(summary.GEx_default, aes(x=Gene,y=cluster)) +
    geom_point(aes(size = Pct, fill= mean), color="black", shape=21) +
    scale_size("% detected", range = c(0,15), breaks=c(0,25,50,75,100), limits=c(0,100)) +
    scale_fill_gradient2(low="steelblue1", high = "red", mid="white",midpoint=0,
                         limits=range(summary.GEx_default$mean),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Average Expression\nlog1p(counts)") +
    ylab("Cluster") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title = element_text(size=14) ,
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()
          ) 
  
  
  pdf(file.path(figures_path, "GEx-DGE-dotplot_default.pdf") , width = 25 , height = 9)
  print(GEx.DGE_DotPlot_default)
  dev.off()
  
}




# Compute SCTransform using all genes and run PCA
{
  TE001@reductions$pca <- NULL
  TE001@reductions$umap <- NULL
  
  
  TE001 <- SCTransform(object=TE001,
                       assay="RNA",
                       ncells=ncol(TE001),
                       variable.features.n = length(all.genes),
                       do.center = TRUE,
                       do.scale = TRUE,
                       conserve.memory = TRUE)
  
  TE001 <- RunPCA(TE001, features=all.genes)
  
}




# Compute clustering solution at gene expression
{
  
  if (running_clustering == TRUE){
    cat("Computing clustering solution at gene expression..\n")
    control <- list(maxit=1e4,
                    smooth=FALSE,
                    max.time=60*30,
                    temperature=1e9)
    
    TE001 <- SAClustering(S.obj = TE001,
                          assay="SCT",
                          control=control,
                          type.fun = "group.mean.silhouette")
    
    opt_param <- TE001@assays$SCT@misc$SA.history$optim.par # res, NN 
  } else if (running_clustering == FALSE){
    cat("Using previously identified optimal solution at gene expression, without recomputing it..\n")
    
    TE001 <- getFinal(S.obj = TE001,
                      assay="SCT",
                      res=opt_param[1],
                      NN=opt_param[2],
                      type.fun="group.mean.silhouette")
    
  }
  
  # rename clusters such that levels start from 1 instead of 0
  levels(TE001$seurat_clusters) <- as.character(as.numeric(levels(TE001$seurat_clusters)) + 1)
  
  npcs <- dim( TE001@reductions$pca@cell.embeddings )[2]
  
  PCs<-as.data.frame(TE001$pca@cell.embeddings[,1:npcs])
  euc_distance<-dist(PCs, method="euclidean")
  s <- silhouette( as.integer(TE001$seurat_clusters) , euc_distance )
  
  pdf(file.path(figures_path, "step_1_silhouette.pdf") , width = 6 , height = 8)
  fviz_silhouette(s,print.summary = FALSE)			
  dev.off()
  
  
  npcs <- dim( TE001@reductions$pca@cell.embeddings )[2]
  
  TE001 <-  RunUMAP(TE001, dims = 1:npcs)
  
  pdf(file.path(figures_path,"step_1_clustering_umap.pdf"), width = 12, height = 6)  # Adjust width and height as needed
  DimPlot(TE001, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), reduction="umap", pt.size=1.5)
  dev.off()
  
  TE001@reductions$umap@cell.embeddings <- as.matrix( TE001@meta.data[c("UMAP_1_scanpy","UMAP_2_scanpy")] )
  colnames(TE001@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  
  
  pdf(file.path(figures_path,"step_1_clustering_umap_pa.pdf"), width = 12, height = 6)  # Adjust width and height as needed
  DimPlot(TE001, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), reduction="umap", pt.size=1.5)
  dev.off()
  
}


# Split by cluster and subcluster
{
  TE001_list <- SplitObject(object=TE001, split.by="seurat_clusters")
  
  TE001_clusters <- TE001$seurat_clusters %>% unique()
  
  cat("Subclustering only clusters with at least 100 cells..\n")
  n_cells_per_cluster <- sapply(TE001_list, function(x){dim(x)[2]}) 
  TE001_clusters_to_subcluster <- TE001_clusters[n_cells_per_cluster >=100]
  
  
  if (running_clustering == TRUE){
    optim_par_list <- matrix(NA, nrow=length(TE001_clusters_to_subcluster), ncol=2)
    colnames(optim_par_list) <- c("res", "NN")
    rownames(optim_par_list) <- as.character(TE001_clusters_to_subcluster)
  
    cat("Subclustering on a cluster by cluster basis..\n")
    for (k in TE001_clusters_to_subcluster){
    
      TE001_list[[k]] <- SAClustering(S.obj = TE001_list[[k]],
                                 assay="SCT",
                                 control=control,
                                 type.fun = "group.mean.silhouette")
    
      # rename clusters such that levels start from 1 instead of 0
      levels(TE001_list[[k]]@meta.data$seurat_clusters) <- as.character(as.numeric(levels(TE001_list[[k]]@meta.data$seurat_clusters)) + 1)
      # store name for branching
      TE001_list[[k]]@meta.data$seurat_clusters <- paste0(k,"_",TE001_list[[k]]@meta.data$seurat_clusters)
    
      optim_par_list[k,] <-  TE001_list[[k]]@assays$SCT@misc$SA.history$optim.par # res, NN 
    }
  } else if (running_clustering == FALSE){
    cat("Subclustering on a cluster by cluster basis by using the previously optimized clustering parameters..\n")
    
    for (k in TE001_clusters_to_subcluster){
      
      opt_res <- optim_par_list[k,"res"]
      opt_NN <- optim_par_list[k,"NN"]
      TE001_list[[k]] <- getFinal(S.obj = TE001_list[[k]],
                                      assay="SCT",
                                      res=opt_res,
                                      NN=opt_NN,
                                      type.fun = "group.mean.silhouette")
      
      levels(TE001_list[[k]]@meta.data$seurat_clusters) <- as.character(as.numeric(levels(TE001_list[[k]]@meta.data$seurat_clusters)) + 1)
      TE001_list[[k]]@meta.data$seurat_clusters <- paste0(k,"_",TE001_list[[k]]@meta.data$seurat_clusters)
   }
  }
  
  TE001_list[[k]]@meta.data$seurat_clusters <- as.factor(TE001_list[[k]]@meta.data$seurat_clusters)
  
}

# Remerge all clusters into a Seurat object and display it
{
  TE001_subclustering <- merge(x=TE001_list[[1]],
                                y=TE001_list[2:length(TE001_list)])

  # re-run PCA to visualize it
  TE001_subclustering <- RunPCA(TE001_subclustering, features=all.genes)
  
  npcs <- dim( TE001_subclustering@reductions$pca@cell.embeddings )[2]
  
  
  PCs<-as.data.frame(TE001_subclustering$pca@cell.embeddings[,1:npcs])
  euc_distance<-dist(PCs, method="euclidean")
  
  TE001_subclustering@meta.data$seurat_clusters <- as.factor(TE001_subclustering@meta.data$seurat_clusters)
  
  
  s <- silhouette( as.integer(TE001_subclustering$seurat_clusters) , euc_distance )
  
  
  pdf(file.path(figures_path, "step_2_silhouette.pdf") , width = 6 , height = 8)
  fviz_silhouette(s,print.summary = FALSE)			
  dev.off()
  
  TE001_subclustering <-  RunUMAP(TE001_subclustering, dims = 1:npcs)
  
  
  
  pdf(file.path(figures_path,"step_2_clustering_umap.pdf"), width = 12, height = 6)  # Adjust width and height as needed
  DimPlot(TE001_subclustering, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), 
          reduction="umap", pt.size=1.5 )
  dev.off()
  
  
  TE001_subclustering@reductions$umap@cell.embeddings <- as.matrix( TE001_subclustering@meta.data[c("UMAP_1_scanpy","UMAP_2_scanpy")] )
    colnames(TE001_subclustering@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  
  pdf(file.path(figures_path,"step_2_clustering_umap_pa.pdf"), width = 12, height = 6)  # Adjust width and height as needed
  DimPlot(TE001_subclustering, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), 
          reduction="umap", pt.size=1.5 )
  dev.off()
  
  pdf(file.path(figures_path,"step_2_clustering_umap_pa_cytotrace.pdf"), width = 6, height = 6)  # Adjust width and height as needed
  FeaturePlot(TE001_subclustering, features = "cytotrace_te001", pt.size=1.5) & scale_color_viridis_c()
  dev.off()
  
}


# Perform differential gene expression analysis
{
  cat("Perform differential gene expression analysis on the identified clusters..\n")
  Idents(TE001_subclustering) <- "seurat_clusters"
  
  # Use RNA slot for calculation of DE genes
  TE001_subclustering <- NormalizeData(TE001_subclustering, assay="RNA", normalization.method = "LogNormalize")
  
  
  # FindAllMarkers
  TE001_subclustering_markers <- FindAllMarkers(TE001_subclustering, 
                                                assay = "RNA",
                                                slot="data",
                                                test.use = "wilcox", 
                                                only.pos = TRUE, 
                                                min.pct = 0.25)
  
  
  
  top10 <- TE001_subclustering_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

  
}


# Heatmap for DE analysis 
{
  genes <- top10$gene
  mat_to_show <- as.matrix(TE001_subclustering@assays$SCT@scale.data[genes,])
  
  Clusters <- TE001_subclustering$seurat_clusters
  colors <- rainbow(nlevels(Clusters)) 
  names(colors) <- levels(Clusters)
  
  # prepare df_annot.df 
  per.cell.sil.width <- s[,"sil_width"]
  names(per.cell.sil.width) <- colnames(TE001_subclustering)
  
  df_annot.df <- TE001_subclustering@meta.data %>%
    dplyr::select(cluster_id=seurat_clusters, cytotrace_score = cytotrace_te001) %>%
    tibble::rownames_to_column(.,var="cell_id")
  df_annot.df$cluster_id <- factor(as.numeric(df_annot.df$cluster_id))
  
  x <- as.data.frame(per.cell.sil.width) %>%
    tibble::rownames_to_column(.,var="sample_id")
   
  df_annot.df <- left_join(df_annot.df,x,by=c("cell_id"="sample_id"))
  
  df_annot.df <- merge(df_annot.df, TE001_subclustering@meta.data[, c("cell_id", "seurat_clusters")], by = "cell_id", all.x = TRUE) # rename clusters as a tree-like structure
  df_annot.df <- df_annot.df[,-2] %>% dplyr::rename(cluster_id=seurat_clusters)

  sample_ordering_to_print <- rev( order( df_annot.df$per.cell.sil.width ) )
  
  df_annot.df <- df_annot.df[sample_ordering_to_print,]
  
  mat_to_show <- mat_to_show[,match(df_annot.df$cell_id,colnames(mat_to_show))]
  
  
  stopifnot(identical( df_annot.df$cell_id, colnames(mat_to_show) ))
  
  
  
  # the silhouette score
  sil.fill <- rep(NA, nrow(df_annot.df))
  unique_clusters <- unique(df_annot.df$cluster_id)
  n_colors <- length(colors)
  
  color_index <- 1
  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]
    sil.fill[df_annot.df$cluster_id == cluster] <- colors[names(colors) == cluster]
    color_index <- i + 1  # Use modulo to cycle through colors
  }
  
  df_annot_cols_cluster_score <- HeatmapAnnotation( 
    silhouette_score = anno_barplot(
      df_annot.df$per.cell.sil.width,
      bar_width = 1,
      gp = gpar(col = NA, 
                fontsize = 6 ,
                cex = 0.8 ,
                fill =  sil.fill 
      ),
      border = FALSE,
      height = unit(1.25, "cm") ,
      width = unit(0.5, "cm")
    )
  )
  
  
  # heatmap annotation
  df_annot.df <- df_annot.df %>% dplyr::select(cluster_id)
  cluster_colors <- colors[ !is.na(names(colors)) ]
  df_annot_cols = HeatmapAnnotation( df = df_annot.df,
                                     col = list(cluster_id=cluster_colors)  , 
                                     show_legend = TRUE , na_col = "gray95",
                                     annotation_name_gp = gpar(fontsize=0.1)
  )
  
  
  top10_hm <- Heatmap( mat_to_show,
                     col = circlize::colorRamp2(range(mat_to_show_default)/2, c("darkslateblue", "yellow")),
                     use_raster = FALSE ,
                     column_split = df_annot.df$cluster_id,
                     heatmap_legend_param = list(color_bar = "continuous") ,	
                     # show_heatmap_legend = FALSE , 
                     # clustering_method_columns = "ward.D2" ,			
                     cluster_columns = FALSE ,
                     cluster_rows = FALSE ,
                     show_column_names = FALSE ,
                     # cluster_columns = columns_dend ,
                     name = "Gene Expression (log1p)" ,  
                     row_names_gp = gpar(fontsize = 10) , column_names_gp = gpar(fontsize = 8) ,
                     row_title_side = "left" , row_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  ht_list = df_annot_cols_cluster_score %v% df_annot_cols %v% top10_hm
  
  pdf(file.path(figures_path,"top_10_heatmap.pdf"), width = 10, height = 12)  # Adjust width and height as needed
  p <- draw(ht_list, column_km = 1 , 
            heatmap_legend_side = "right", 
            annotation_legend_side = "bottom" )
  dev.off()  
}


# DotPlot
{
  
  
  top5 <- TE001_subclustering_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  genes_dp <- top5$gene %>% unique()
  
  # intialize lists containing dataframes 
  Avgs <- vector( mode="list", length=length(levels(TE001_subclustering$seurat_clusters)) )
  names(Avgs) <- levels(TE001_subclustering$seurat_clusters)
  
  Pct <- vector( mode="list", length=length(levels(TE001_subclustering$seurat_clusters)) )
  names(Pct) <- levels(TE001_subclustering$seurat_clusters)
  
  for (cluster in levels(TE001_subclustering$seurat_clusters)){
  
    cells_in_subcluster <- TE001_subclustering$seurat_clusters == cluster
    
    # Avg 
    Avgs[[cluster]] <- TE001_subclustering[["RNA"]]@data[genes_dp, cells_in_subcluster] %>%
      rowMeans() %>% 
      as.data.frame() %>% 
      rownames_to_column() %>%
      dplyr::rename(.,"Gene"=rowname,"mean"=.) %>% 
      dplyr::mutate(cluster=cluster)
      
  
    # Pct
    Pct[[cluster]] <- TE001_subclustering[["RNA"]]@data[genes_dp, cells_in_subcluster] %>%
      apply(.,1, function(x){ sum(x>0)/length(x)*100 }) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      dplyr::rename(.,"Gene"=rowname,"Pct"=.) %>%
      dplyr::mutate(cluster=cluster)
    
  }
  
  # Mean
  df_mean.GEx <- do.call(rbind, Avgs)
    rownames(df_mean.GEx) <- NULL
  # Pct
  df_pct.GEx <- do.call(rbind, Pct)
  rownames(df_pct.GEx) <- NULL
  
  
  # Summary
  summary.GEx <- cbind(df_mean.GEx[,1:2],df_pct.GEx[,2:3]) 
  summary.GEx$Gene <- factor(summary.GEx$Gene, levels=genes_dp, ordered = TRUE) 
  
  
  
  # Printing dotplot for genes
  cat("Printing dotplots..\n")
  GEx.DGE_DotPlot <- ggplot(summary.GEx, aes(x=Gene,y=cluster)) +
    geom_point(aes(size = Pct, fill= mean), color="black", shape=21) +
    scale_size("% detected", range = c(0,15), breaks=c(0,25,50,75,100), limits=c(0,100)) +
    scale_fill_gradient2(low="steelblue1", high = "red", mid="white",midpoint=0,
                         limits=range(summary.GEx$mean),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Average Expression\nlog1p(counts)") +
    ylab("Cluster") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title = element_text(size=14) ,
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) 
  
  
  pdf(file.path(figures_path, "GEx-DGE-dotplot.pdf") , width = 22 , height = 9)
  print(GEx.DGE_DotPlot)
  dev.off()
  
}

# Save clustering result into Seurat object
{ 
  saveRDS(object=TE001_subclustering, file=file.path(DataFolder, "TE001_iter_subclustering_gene_expr.rds"))
}

