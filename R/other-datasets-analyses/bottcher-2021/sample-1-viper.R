## 
# Single cell Analysis of Ctrl1 from Bottcher, et al. 2021 (VIPER Analysis) 
# --------------------------------------------------------------------------------------------------
# system.time({ source("sources/bottcher-2021/sample-1-viper.R") })

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

set.seed(666)

retainOnlyTFs <- function(network) {
	
	print_msg_info(">>> Filtering using my lists")
	x <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat" , col_names = FALSE )
	y <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat" , col_names = FALSE )
	# z <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat" , col_names = FALSE )
	
	selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))	
	
	index <- names(network) %in% selected_regulators
	return( network[index] )
	
}

# my_sample_id <- "Ctrl1"
my_sample_id <- "Ctrl2"
# create_workspace("Ctrl1-viper-with-sct-TE001-ref-model")
create_workspace(paste0(my_sample_id,"-viper"))

library(logger)
log_threshold(DEBUG)
log_threshold(DEBUG,index=2)

log_filename <- file.path(reports.dir,"log.txt")
log_appender(appender = appender_file(log_filename),index=2)
# log_formatter(formatter_glue)
log_layout(layout_glue_colors)
log_warn('--- Start of Analysis ----')

require(tidyverse)
# filename <- "~/Clouds/Dropbox/Data/isc/bottcher-2021/ctrl1/ctrl1-cpm.rds"
# log_info('Readiimg CPM matrix: {filename}')
# fs::file_info(filename)
# cpm_matrix <- readRDS(filename)
# dim(cpm_matrix)

library(Seurat)
filename <- paste0("~/Clouds/Dropbox/Data/isc/bottcher-2021/",my_sample_id,"/",my_sample_id,"-seurat-analysis-data.rds")
log_info('Reading scale.data matrix: {filename}')
fs::file_info(filename)
sc_data <- readRDS(filename)
dim(sc_data)

# tmp <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")

sc_data <- sc_data %>% 
	SCTransform( 
		# variable.features.n = 2000 , 
							 # reference.SCT.model = tmp@assays$SCT@SCTModel.list$model1 , 
							 # residual.features = rownames(tmp@assays$SCT@SCTModel.list$model1@feature.attributes) ,
								residual.features = rownames(sc_data@assays$RNA@counts) ,
							 return.only.var.genes = TRUE ,
							 do.scale = TRUE , do.center = TRUE , 
							 seed.use = 666 )

my_ges <- as.matrix(sc_data@assays$SCT@scale.data)
dim(my_ges)

# # filename <- "~/Clouds/Dropbox/Data/isc/TE001/networks/4-lists/TE001-metacells_unPruned.rds"
# filename <- "~/Clouds/Dropbox/Data/isc/TE001/networks/"
# log_info('Reading network: {filename}')
# isc_network <- readRDS(filename)

print_msg_info(">>> >> Statically assigning networks to clusters")
{
	# clusters_networks.list <- vector("list",3)
	clusters_networks.list <- list()
	clusters_networks.list[["0"]] <- readRDS(file.path("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
	clusters_networks.list[["1"]] <- readRDS(file.path("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c2_unPruned.rds"))
	clusters_networks.list[["2"]] <- readRDS(file.path("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c3and4_unPruned.rds"))
	str(clusters_networks.list,1)
}

library(viper)
clusters_networks.list <- lapply(clusters_networks.list,pruneRegulon,50)
clusters_networks.list <- lapply(clusters_networks.list,retainOnlyTFs)

viper_method <- "none"
log_info("Signature for VIPER: {viper_method}")
log_info("Runnig VIPER analysis")
vpmat <- viper( my_ges ,
								regulon = clusters_networks.list , method = viper_method , mvws = 1)
log_info("Runnig VIPER analysis [completed]")

# x <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-stouffer-7-clusters.rds")
# colnames(x) <- c("stem-1","stem-2","absorptive-1","absorptive-2","secretory-1","secretory-2","secretory-3")
stouffer_table <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-stouffer-8-clusters.rds")
colnames(stouffer_table) <- c("stem-1","stem-2","absorptive-1","absorptive-2","secretory-1","secretory-2","secretory-3","z_paneth")

faked_reg <- generateRegulonObjectFromProteinActivityMatrix(stouffer_table,n_top = 25)
isc_enrichment_matrix <- aREA(vpmat,faked_reg,minsize = 5)$nes

my_func <- function(x) { -log10( pnorm( x , lower.tail = F )*2 ) }
x <- my_func(isc_enrichment_matrix)
x <- ifelse(x<5,0,x)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
x <- range01(x)
# rowMax(x)
# y <- apply(x,1,range01) %>% t()
# x <- apply(x,1,function(z)z/sum(z)) %>% t()
isc_enrichment_matrix <- apply(x,2,function(z)z/sum(z))

library(Seurat)
vp_seurat <- CreateSeuratObject( counts = vpmat ,
																 assay = "VIPER" )
rownames(vp_seurat@meta.data) <- colnames(vp_seurat)

vp_seurat@assays$VIPER@scale.data=as.matrix(vp_seurat@assays$VIPER@data) # TODO: Scale vpmat
vp_seurat@meta.data$nCount_VIPER <- 1000
vp_seurat@meta.data$nFeature_VIPER <- nrow(vpmat)

# vp_seurat <- RunPCA(vp_seurat,seed.use = 666,features = rownames(vp_seurat))
most_var_feats <- rownames(vp_seurat)[ order(rowVars(vp_seurat@assays$VIPER@scale.data),decreasing = T)[1:100] ]
vp_seurat <- RunPCA(vp_seurat,seed.use = 666,features = most_var_feats )
vp_seurat <- RunUMAP( vp_seurat ,
											assay = "VIPER" ,
											reduction.name = "densmap" ,
											densmap = TRUE ,
											dims = 1:30 ,
											# features = selected_features ,
											n.components = 3L ,
											reduction = "pca" ,
											metric = "euclidean" , 
											verbose = FALSE , seed.use = 666 )

umap_1 <- vp_seurat@reductions$densmap@cell.embeddings[,"UMAP_1"]
umap_2 <- vp_seurat@reductions$densmap@cell.embeddings[,"UMAP_2"]

dim(vpmat)

clustering_assignment <- apply(isc_enrichment_matrix,2,function(i) rownames(isc_enrichment_matrix)[which.max(i)],simplify = TRUE)
clustering_assignment <- unlist(clustering_assignment)
table(clustering_assignment)

my_tibble <- tibble( cell_id = names(clustering_assignment) , cluster_id = clustering_assignment )
table(my_tibble$cluster_id)

my_tibble$umap_1 = umap_1[ match( my_tibble$cell_id , names(umap_1) ) ]
my_tibble$umap_2 = umap_2[ match( my_tibble$cell_id , names(umap_2) ) ]

tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
														"#8c564b","#e377c2","#bcbd22","#7f7f7f","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
																"brown","pink","olive","gray","cyan")
# tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
# 									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
# names(tab10_palette) <- c("blue",'orange',"green","red","purple",
# 													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(my_tibble$cluster_id))
my_palette <- my_palette[!is.na(names(my_palette))]

p <- ggplot(my_tibble) +
	geom_point( data = my_tibble , 
							# aes(x,y), col="grey75" , size = 0.5 ) +
							aes(umap_1,umap_2,color=cluster_id) , size = 1 ) +
	theme_light() +
	scale_color_manual(values = my_palette)

pdf( file.path(reports.dir, paste0(my_sample_id,"-cluster-umap.pdf")) , width = 5 , height = 4 )
	print(p)
dev.off()

cyto_table <- read_delim( paste0("~/Clouds/Dropbox/Data/isc/bottcher-2021/",my_sample_id,"/",my_sample_id,"-cytotrace.csv"),"\t")
my_tibble <- left_join( my_tibble , cyto_table )

## Plotting Cytotrace as boxplot ----
p <- ggplot( my_tibble , aes(cluster_id, cytotrace_score , 
															 fill = cluster_id
) ) +
	geom_boxplot() +
	# scale_fill_brewer(palette = "Set1") +
	scale_fill_manual(values = my_palette) +
	scale_x_discrete( position = "top" ) +
	theme_light() +
	theme( axis.text.x = element_text(angle = 45, hjust = 0, face = "bold") )

pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-boxplot-cytotrace.pdf" )) , width = 4 , height = 3 )
print(p)
dev.off()		

p <- ggplot(my_tibble) +
	geom_point( data = my_tibble , 
							# aes(x,y), col="grey75" , size = 0.5 ) +
							aes(umap_1,umap_2,color=cytotrace_score) , size = 1 ) +
	theme_light() +
	scale_color_viridis()

pdf( file.path(reports.dir,paste0(my_sample_id,"-cytotrace-umap.pdf")) , width = 5 , height = 4 )
print(p)
dev.off()

log_warn('--- END of Analysis ----')
