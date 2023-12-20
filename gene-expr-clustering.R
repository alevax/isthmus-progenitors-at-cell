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
  
running_clustering <- FALSE # do not re-run clustering
cat("'runnig_clustering' is set to: ", running_clustering)

# optimal clustering parameters previously determined
# step 1
opt_param <- c(0.02213701, 7) # res NN 
# step 2 (subclustering)
optim_par_list <- matrix(c(0.02839557, 27, 0.05340340, 4, 0.04845794, 6),byrow=TRUE,nrow=3,ncol=2)
colnames(optim_par_list) <- c("res", "NN")
rownames(optim_par_list) <- c("0", "2", "1")


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
  
 
  PCs<-as.data.frame(TE001$pca@cell.embeddings[,1:10])
  euc_distance<-dist(PCs, method="euclidean")
  
  s <- silhouette( as.integer(TE001$seurat_clusters) , euc_distance )
  fviz_silhouette(s,print.summary = FALSE)			
  
  
  npcs <- dim( TE001@reductions$pca@cell.embeddings )[2]
  
  TE001 <-  RunUMAP(TE001, dims = 1:npcs)
  
  DimPlot(TE001, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), reduction="umap" )
  
  
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
  
  PCs<-as.data.frame(TE001_subclustering$pca@cell.embeddings[,1:10])
  euc_distance<-dist(PCs, method="euclidean")
  
  TE001_subclustering@meta.data$seurat_clusters <- as.factor(TE001_subclustering@meta.data$seurat_clusters)
  
  s <- silhouette( as.integer(TE001_subclustering$seurat_clusters) , euc_distance )
  fviz_silhouette(s,print.summary = FALSE)			
  
  
  npcs <- dim( TE001_subclustering@reductions$pca@cell.embeddings )[2]
  
  TE001_subclustering <-  RunUMAP(TE001_subclustering, dims = 1:npcs)
  
  TE001_subclustering@reductions$umap@cell.embeddings <- as.matrix( TE001_subclustering@meta.data[c("UMAP_1_scanpy","UMAP_2_scanpy")] )
    colnames(TE001_subclustering@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
  
  DimPlot(TE001_subclustering, group.by=c("seurat_clusters", "iter_cluster_id_with_paneth"), 
          reduction="umap", pt.size=1.5 )
  
}

# Save clustering result
{ 
  
  saveRDS(object=TE001_subclustering, file=file.path(DataFolder, "TE001_iter_subclustering_gene_expr.rds"))
  
  }