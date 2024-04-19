# Luca Zanella 
# 12.19.2023

# Gene expression based clustering
# INPUTS (MANDATORY): 
# DataFolder: path to data are stored (scratch)


rm(list = ls()) # cleans workspace
gc()
cat("\014")

# Set paths
cat("YOU MUST SET DataFolder and decide whether to run clustering or not (by using predetermined values)..\n")
DataFolder <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001/" # "M:/shares/hcmi/" #
h5ad_path <- "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001-h5ad/"
figures_path = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/TE001_gExpr_clustering/"


running_clustering <- FALSE # do not re-run clustering
cat("'runnig_clustering' is set to: ", running_clustering)

# optimal clustering parameters previously determined
# step 1
opt_param <- c(0.02213701, 7) # res NN 
# step 2 (subclustering)
optim_par_list <- matrix(c(0.02839557, 27, 0.05340340, 4, 0.04845794, 6),byrow=TRUE,nrow=3,ncol=2)
colnames(optim_par_list) <- c("res", "NN")
rownames(optim_par_list) <- c("1", "3", "2")


cytotrace_markers <- c('Smarca5','Rbbp7','Tcerg1','Hnrnpd','Hmg20b','Nelfe','Ube2i','Etv5','Ubn1','Mbd3','Dek','Maz',
                       'Itgb3bp','Ilf2','Pa2g4', 'Id3','Hnf4g','Atoh1','Spdef','Neurod1')
  
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
  
  
  
  
}



# Loading data and metadata
{
  TE001 <- readRDS(file.path(DataFolder, "TE001-seurat-analysis-data.rds"))
  
  metadata_csv <- file.path(h5ad_path, "TE001-metadata-umap-and-clusters-for-paper.csv")
  metadata <- read.csv(metadata_csv)
  
  TE001@meta.data <- TE001@meta.data %>% dplyr::left_join(., metadata, 
                                                          by = c("cell_id" = "cell_id"))
  rownames(TE001@meta.data) <- TE001@meta.data$cell_id
  
}

# Recompute SCTransform using all genes and PCA
{
  
  DefaultAssay(TE001) <- "RNA"
  TE001[["SCT"]] <- NULL
  TE001@reductions$umap <- NULL
  
  genes.list <- rownames(TE001@assays$RNA@counts)
  
  TE001 <- SCTransform(object=TE001,
                       assay="RNA",
                       ncells=ncol(TE001),
                       variable.features.n = length(genes.list),
                       do.center = TRUE,
                       do.scale = TRUE,
                       conserve.memory = TRUE)
  
  TE001 <- RunPCA(TE001, features=genes.list)
  
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
  }
  else if (running_clustering == FALSE){
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
  
  TE001@reductions$umap@cell.embeddings <- as.matrix( TE001@meta.data[c("UMAP_1_scanpy","UMAP_2_scanpy")] )
  colnames(TE001@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  
  
  pdf(file.path(figures_path,"step_1_clustering_umap.pdf"), width = 12, height = 6)  # Adjust width and height as needed
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
  TE001_subclustering <- RunPCA(TE001_subclustering, features=genes.list)
  
  npcs <- dim( TE001_subclustering@reductions$pca@cell.embeddings )[2]
  
  
  PCs<-as.data.frame(TE001_subclustering$pca@cell.embeddings[,1:npcs])
  euc_distance<-dist(PCs, method="euclidean")
  
  TE001_subclustering@meta.data$seurat_clusters <- as.factor(TE001_subclustering@meta.data$seurat_clusters)
  
  
  s <- silhouette( as.integer(TE001_subclustering$seurat_clusters) , euc_distance )
  
  
  pdf(file.path(figures_path, "step_2_silhouette.pdf") , width = 6 , height = 8)
  fviz_silhouette(s,print.summary = FALSE)			
  dev.off()
  
  TE001_subclustering <-  RunUMAP(TE001_subclustering, dims = 1:npcs)
  
  TE001_subclustering@reductions$umap@cell.embeddings <- as.matrix( TE001_subclustering@meta.data[c("UMAP_1_scanpy","UMAP_2_scanpy")] )
    colnames(TE001_subclustering@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  
  pdf(file.path(figures_path,"step_2_clustering_umap.pdf"), width = 12, height = 6)  # Adjust width and height as needed
  DimPlot(TE001_subclustering, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), 
          reduction="umap", pt.size=1.5 )
  dev.off()
}

# Save clustering result
{ 
  #saveRDS(object=TE001_subclustering, file=file.path(DataFolder, "TE001_iter_subclustering_gene_expr.rds"))
}



# Define heatmap function
{
  gExprHeatmap=function(genes,expression_to_plot,Clusters,colors){
    require(pheatmap)
    annotation_columns<-data.frame(Clusters)
    rownames(annotation_columns)<-colnames(expression_to_plot)
    colnames(annotation_columns)<-c("clusters")
    
    anno_colors <- list(clusters = colors)
    names(anno_colors[[1]]) <- levels(annotation_columns$cluster)
    
    cellType_Order<-unique(Clusters)
    idx<-c()
    for (i in 1:length(cellType_Order)){
      idx<-c(idx, which(annotation_columns[,1]==cellType_Order[i]))
    }
    
    expression_to_plot=expression_to_plot[,idx]
    annotation_columns=data.frame(annotation_columns[idx,])
    rownames(annotation_columns)<-colnames(expression_to_plot)
    colnames(annotation_columns)=c("clusters")
    
    quantile_breaks <- function(xs, n = 10) {
      breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
      breaks[!duplicated(breaks)]
    }
    mat_breaks <- quantile_breaks(as.matrix(apply(expression_to_plot,1,function(x){(x-mean(x))/sd(x)})), n = 30)
    
    library(pheatmap)
    file=c("Gene_Expression_Heatmap")
    p<-pheatmap(expression_to_plot, cluster_rows=FALSE, 
                show_rownames=TRUE, cluster_cols=FALSE, 
                annotation_col=annotation_columns, breaks=mat_breaks, 
                color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),
                fontsize_row = 8, show_colnames = FALSE, annotation_colors = anno_colors,
                scale="row")
    pdf(file= paste(file, "pdf", sep="."), height = 10, width=15)
    print(p)
    no_show<-dev.off()
    
  }
  
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
  
  
  #pdf(file.path(figures_path,"top_10_heatmap.pdf"), width = 10, height = 10)  # Adjust width and height as needed
  #gExprHeatmap(genes, mat_to_show, Clusters, colors)
  #dev.off()
 
  
  
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
  
  
 
  
  
  col_fun = colorRamp2(c(-2, 0, 2), c("deepskyblue3", "white", "brown3"))
  
  top10_hm <- Heatmap( mat_to_show,
                     col = col_fun ,
                     use_raster = FALSE ,
                     column_split = df_annot.df$cluster_id,
                     heatmap_legend_param = list(color_bar = "continuous") ,	
                     # show_heatmap_legend = FALSE , 
                     # clustering_method_columns = "ward.D2" ,			
                     cluster_columns = FALSE ,
                     cluster_rows = FALSE ,
                     show_column_names = FALSE ,
                     # cluster_columns = columns_dend ,
                     name = "Gene Expression" ,  
                     row_names_gp = gpar(fontsize = 10) , column_names_gp = gpar(fontsize = 8) ,
                     row_title = "Gene Expression", 
                     row_title_side = "left" , row_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  ht_list = df_annot_cols_cluster_score %v% df_annot_cols %v% top10_hm
  
  pdf(file.path(figures_path,"top_10_heatmap.pdf"), width = 10, height = 12)  # Adjust width and height as needed
  p <- draw(ht_list, column_km = 1 , 
            heatmap_legend_side = "right", 
            annotation_legend_side = "bottom" )
  dev.off()  
}


# Heatmap for relevant genes (cytotrace genes)
{
  
  # mat_cytotrace <- as.matrix(TE001_subclustering@assays$SCT@scale.data[cytotrace_markers,])
  # 
  # cyto_hm <- Heatmap( mat_cytotrace,
  #                      col = col_fun ,
  #                      use_raster = FALSE ,
  #                      column_split = df_annot.df$cluster_id,
  #                      heatmap_legend_param = list(color_bar = "continuous") ,	
  #                      # show_heatmap_legend = FALSE , 
  #                      # clustering_method_columns = "ward.D2" ,			
  #                      cluster_columns = FALSE ,
  #                      cluster_rows = FALSE ,
  #                      show_column_names = FALSE ,
  #                      # cluster_columns = columns_dend ,
  #                      name = "Gene Expression" ,  
  #                      row_names_gp = gpar(fontsize = 10) , column_names_gp = gpar(fontsize = 8) ,
  #                      row_title = "Gene Expression", 
  #                      row_title_side = "left" , row_title_gp = gpar(fontsize = 10, fontface = "bold")
  # )
  # 
  # cyto_list = df_annot_cols_cluster_score %v% df_annot_cols %v% cyto_hm
  # p <- draw(cyto_list, column_km = 1 , 
  #           heatmap_legend_side = "right", 
  #           annotation_legend_side = "bottom" )
  # 
  # pdf(file.path(figures_path,"cytotrace_heatmap.pdf"), width = 10, height = 10)  # Adjust width and height as needed
  # gExprHeatmap(cytotrace_markers, mat_cytotrace, Clusters, colors)
  # dev.off()
  
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

