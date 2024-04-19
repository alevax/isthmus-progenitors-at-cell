
# system.time({ source("sources/figures/differential-density-using-reference-as-null-TE001-vs-TE002.R") })

library(reticulate)
# use_python('/Applications/anaconda3/bin/python')
use_python('~/opt/anaconda3/envs/single-cell-env/bin/python')

library(Seurat)
library(scales)
library(viridis)

source("../vaxtools/R/utils.R")

my_seed <- 666
set.seed(my_seed)

n_tiles_per_axis <- 100
n_iter <- 100

create_workspace("diff-kde-irra-densumap-TE001-vs-TE002")
vp_seurat <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE002-viper.rds")

## Enrichment Analysis for Gathering Clustering Labels ----
vp_baseline <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-with-paneth.rds")
table(vp_baseline@meta.data$pas_cluster_id)
my_clusters <- factor( vp_baseline@meta.data$pas_cluster_id )
names(my_clusters) <- vp_baseline@meta.data$cell_id
stouffer_table <- doStoufferFromClusters( as.matrix(vp_baseline@assays$VIPER@scale.data) , my_clusters )
isc_markers.regulon <- generateRegulonObjectFromProteinActivityMatrix(stouffer_table,n_top = 25)

library(viper)
vpmat <- as.matrix(vp_seurat@assays$VIPER@scale.data)
isc_enrichment.mat <- aREA(vpmat,isc_markers.regulon)$nes

get_max_names <- function(mat) {
	row_index <- apply(mat,2,which.max)
	return(rownames(mat)[row_index])
}
clustering_solution_from_enrichment <- get_max_names(isc_enrichment.mat)
vp_seurat$cluster_id <- clustering_solution_from_enrichment
table(vp_seurat$cluster_id)

table(vp_seurat$sample_id,vp_seurat$cluster_id)
# create_workspace("diff-kde-haplo")
# vp_seurat <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE005-viper.rds")

# create_workspace("diff-kde-abla")
# vp_seurat <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE006-viper.rds")

# pca_1 <- vp_seurat@reductions$pca@cell.embeddings[,"PC_1"]
# pca_2 <- vp_seurat@reductions$pca@cell.embeddings[,"PC_2"]

# write_csv(vp_seurat@reductions$pca@cell.embeddings %>% as.data.frame(),"~/Downloads/embedding-for-heeju.csv")

vp_seurat <- RunUMAP( vp_seurat ,
											assay = "VIPER" ,
											reduction.name = "densmap" ,
											densmap = TRUE ,
											dims = 1:30 ,
											# features = selected_features ,
											n.components = 3L ,
											reduction = "pca" ,
											metric = "euclidean" , 
											verbose = FALSE , seed.use = my_seed )

pca_1 <- vp_seurat@reductions$densmap@cell.embeddings[,"densmap_1"]
pca_2 <- vp_seurat@reductions$densmap@cell.embeddings[,"densmap_2"]
pca_3 <- vp_seurat@reductions$densmap@cell.embeddings[,"densmap_3"]

# pca_1 <- vp_seurat@reductions$umap@cell.embeddings[,"UMAP_1"]
# pca_2 <- vp_seurat@reductions$umap@cell.embeddings[,"UMAP_2"]

pca_1 <- rescale(pca_1,to = c(-1,1))
pca_2 <- rescale(pca_2,to = c(-1,1))

limit_vector <- c( min(pca_1) , max(pca_1) ,
									 min(pca_2) , max(pca_2) )

index <- vp_seurat@meta.data$sample_id == "control"
pca_1_control <- pca_1[index]
pca_2_control <- pca_2[index]

index <- vp_seurat@meta.data$sample_id == "treatment"
pca_1_treatment <- pca_1[index]
pca_2_treatment <- pca_2[index]

library(MASS)

getX <- function(x,y,n_iter=100,n=25,limit_vector=NULL) {
	stopifnot(length(x)==length(y))
	if (is.null(limit_vector))
	{
		limit_vector <- c( min(x) , max(x) ,
											 min(y) , max(y))
	}

	res <- list()
	for ( i in 1:n_iter )
	{
		# index <- sample(1:length(x),size = 0.3*length(x) )
		index <- sample(1:length(x),size = 500 )
		res[[i]] <- kde2d( x = x[index] , y = y[index] , n = n , lims = limit_vector ) 
	}

	p <- array( unlist(lapply( res , `[[` , "z" )) , dim = c(n,n,n_iter) )
	mean <- apply( p , 1:2 , mean )
	std <- apply( p , 1:2 , sd )
	
	std <- ifelse( std < 1 , 1 , std )
	
	ret <- list(mean=mean,std=std)
	
	return(ret)
	
}

set.seed(666)

control_kde_total <- kde2d(pca_1_control,pca_2_control , n = n_tiles_per_axis )
str(control_kde_total,1)
treatment_kde_total <- kde2d(pca_1_treatment,pca_2_treatment , n = n_tiles_per_axis )
str(treatment_kde_total,1)

# limit_vector <- c( min(c(treatment_kde_total$x,control_kde_total$x)) , max(c(treatment_kde_total$x,control_kde_total$x)) ,
									 # min(c(treatment_kde_total$y,control_kde_total$y)) , max(c(treatment_kde_total$y,control_kde_total$y)) )

control <- getX( x = pca_1_control ,
								 y = pca_2_control , 
								 limit_vector = limit_vector ,
								 n_iter = n_iter , n = n_tiles_per_axis )

treatment <- getX( x = pca_1_treatment ,
								 y = pca_2_treatment,
								 limit_vector = limit_vector ,
								 n_iter = n_iter , n = n_tiles_per_axis )

m_ctrl_test <- ( control$mean - control$mean ) / control$std
m_ctrl_mean <- control$mean
m_ctrl_std <- control$std
m_diff_mean <- treatment$mean
m_diff <- ( treatment$mean - control$mean ) / control$std
# m_diff <- ( treatment$mean - control$mean )
# m_diff <- treatment$mean


# colnames(m_ctrl_test) <- 1:ncol(m_ctrl_test)
# my_tibble <- m_ctrl_test %>% as.data.frame() %>% rownames_to_column("y") %>% 
# 	pivot_longer( cols = -y , names_to = "x" ) %>% as_tibble()
# my_tibble$sample <- "ctrl_test"
# 
# colnames(m_diff) <- 1:ncol(m_diff)
# x <- m_diff %>% as.data.frame() %>% rownames_to_column("y") %>% 
# 	pivot_longer( cols = -y , names_to = "x" ) %>% as_tibble()
# x$sample <- "treatment_diff"
# 
# my_tibble <- bind_rows(my_tibble,x)
# 
# p <- ggplot( my_tibble , aes( x , y , z = value , fill = value ) ) +
# 	geom_tile() +
# 	# scale_fill_viridis() +
# 	scale_fill_gradient( c(lowest,0,highest) ,  c("blue","white","red") ) +
# 	facet_wrap(~sample,ncol = 2) +
# 	# geom_contour_filled() +
# 	geom_contour() +
# 	theme_void()
# 
# pdf( file.path(reports.dir,"diff-kde.pdf") , width = 8 , height = 4.5 )
# 	print(p)
# dev.off()

# tibble( control_test = m_ctrl_test ) ,
# 				control_mean = m_ctrl_mean ,
# 				control_std = m_ctrl_std ,
# 				treatment_mean = m_diff_mean , 
# 				treatemnt_diff = m_diff )

# m_diff <- ifelse( m_diff < 2 & m_diff > -2 , 0 , m_diff )

library(circlize)
library(ComplexHeatmap)

lowest <- min(c(m_ctrl_mean,m_diff))
highest <- max(c(m_ctrl_mean,m_diff))

# my_colors <- colorRamp2( seq(lowest,highest,length.out=10) ,  viridis(10) )
my_colors <- colorRamp2( c(lowest,0,highest) ,  c("blue","white","red") )


## Control mean and std ----
h_left <- ComplexHeatmap::Heatmap( m_ctrl_mean , col = my_colors , 
																	 name = "Control" ,
																		 cluster_columns = FALSE , cluster_rows = FALSE )
h_right <- ComplexHeatmap::Heatmap( m_ctrl_std , col = my_colors , 
																		name = "Treatment" ,
																			cluster_columns = FALSE , cluster_rows = FALSE )
p <- h_left + h_right

pdf( file.path(reports.dir,"control-mean-and-std-diff-kde.pdf") , width = 8 , height = 3.5 )
	print(p)
dev.off()

## Controlvs treatment ----
h_left <- ComplexHeatmap::Heatmap( m_ctrl_test , col = my_colors , 
																	 name = "Control" ,
																	 cluster_columns = FALSE , cluster_rows = FALSE )
h_right <- ComplexHeatmap::Heatmap( m_diff , col = my_colors , 
																		name = "Treatment" ,
																		cluster_columns = FALSE , cluster_rows = FALSE )
p <- h_left + h_right

pdf( file.path(reports.dir,"control-vs-treatment-diff-kde.pdf") , width = 8 , height = 3.5 )
print(p)
dev.off()

## Control mean vs treatment ----
h_left <- ComplexHeatmap::Heatmap( m_ctrl_mean , col = my_colors , 
																	 name = "Control" ,
																	 cluster_columns = FALSE , cluster_rows = FALSE )
h_right <- ComplexHeatmap::Heatmap( m_diff , col = my_colors , 
																		name = "Treatment" ,
																		cluster_columns = FALSE , cluster_rows = FALSE )
p <- h_left + h_right

pdf( file.path(reports.dir,"control-mean-vs-treatment-diff-kde.pdf") , width = 8 , height = 3.5 )
print(p)
dev.off()

## Control mean vs treatment VIRIDIS ----
h_left <- ComplexHeatmap::Heatmap( m_ctrl_mean , 
																	 name = "Control" ,
																	 col = viridis(10) , use_raster = TRUE ,
																	 cluster_columns = FALSE , cluster_rows = FALSE )
h_right <- ComplexHeatmap::Heatmap( m_diff , 
																		name = "Treatment" ,
																		col = viridis(10) , use_raster = TRUE ,
																		cluster_columns = FALSE , cluster_rows = FALSE )
p <- h_left + h_right

pdf( file.path(reports.dir,"control-mean-vs-treatment-diff-kde-viridis.pdf") , width = 8 , height = 3.5 )
print(p)
dev.off()

## Gathering Clustering Labels ----
x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
table(x$cluster_id)
stopifnot( identical( names(pca_1_control) , names(pca_2_control) ) )
index <- match( gsub("control_","",names(pca_1_control)) , x$sample_id )
df <- data.frame( x = pca_1_control , y = pca_2_control , cluster_id = x$cluster_id[index])
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(df$cluster_id))
my_palette <- my_palette[!is.na(names(my_palette))]

## With Countour ----
p_c <- ggplot(df) +
	geom_point( data = df , 
							# aes(x,y), col="grey75" , size = 0.5 ) +
							aes(x,y,color=cluster_id) , size = 1 ) +
	theme_bw() +
	scale_color_manual(values = my_palette) +
	geom_contour(data=data.frame( x=rep( control_kde_total$x , n_tiles_per_axis ) ,
																y=rep( control_kde_total$y , each=n_tiles_per_axis ),
																z = as.vector( m_ctrl_mean ) ),
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="black" ,size=0.5)
	# guides(color = guide_legend(nrow = 1)) +
	# theme(legend.position = "bottom")

p_t <- ggplot() +
	geom_point( data = data.frame( x = pca_1_treatment , y = pca_2_treatment ) , 
							aes(x,y), col="grey75" , size = 1 ) +
	theme_bw() +
	geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																z = as.vector( ifelse( m_diff < 0 , 0 , m_diff ) ) ) ,
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="red" ,size=0.5) +
	geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																z = -1*as.vector( ifelse( m_diff > 0 , 0 , m_diff ) ) ) ,
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="blue" ,size=0.5) 

p <-  p_c + p_t 

pdf( file.path(reports.dir,"diff-kde-contour.pdf") , width = 9 , height = 3.5 )
	print(p)
dev.off()

# ## FastClust Solution Gathering Clustering Labels ----
# x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-STEM-iter-clustering-table.csv")
# table(x$cluster_id)
# stopifnot( identical( names(pca_1_control) , names(pca_2_control) ) )
# index <- match( gsub("control_","",names(pca_1_control)) , x$sample_id )
# df <- data.frame( x = pca_1_control , y = pca_2_control , cluster_id = x$cluster_id[index])
# df$cluster_id[ is.na(df$cluster_id) ] <- "z_other"
# # tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
# # 									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
# # names(tab10_palette) <- c("blue",'orange',"green","red","purple",
# # 													"brown","pink","gray","olive","cyan")
# # my_palette <- tab10_palette
# # names(my_palette) <- names(table(df$cluster_id))
# # my_palette <- my_palette[!is.na(names(my_palette))]


## Ingest Solution Gathering Clustering Labels ----
x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
table(x$cluster_id)
stopifnot( identical( names(pca_1_control) , names(pca_2_control) ) )
index <- match( gsub("control_","",names(pca_1_control)) , x$sample_id )
df <- data.frame( x = pca_1_control , y = pca_2_control , cluster_id = x$cluster_id[index])
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(df$cluster_id))
my_palette <- my_palette[!is.na(names(my_palette))]

## With Countour ----
p_c <- ggplot(df) +
	geom_point( data = df , 
							# aes(x,y), col="grey75" , size = 0.5 ) +
							aes(x,y,color=cluster_id) , size = 1 ) +
	theme_bw() +
	scale_color_manual(values = my_palette) +
	# scale_color_brewer(palette = "Set1") +
	geom_contour(data=data.frame( x=rep( control_kde_total$x , n_tiles_per_axis ) ,
																y=rep( control_kde_total$y , each=n_tiles_per_axis ),
																z = as.vector( m_ctrl_mean ) ),
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="black" ,size=0.5)
# guides(color = guide_legend(nrow = 1)) +
# theme(legend.position = "bottom")

df <- data.frame( x = pca_1_treatment , y = pca_2_treatment )
x <- read_delim("~/Clouds/Dropbox/Data/isc/rad-metadata-ingest-cytoperm.csv",",")
df$cluster_id <- x$iter_cluster_id[ match( gsub("treatment_","",rownames(df)) , x$cell_id ) ]
df <- df[ !is.na(df$cluster_id) , ]

tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(df$cluster_id))
my_palette <- my_palette[!is.na(names(my_palette))]

p_t <- ggplot() +
	geom_point( data = df , 
							aes(x,y,color=cluster_id), size = 1 ) +
	scale_color_manual(values = my_palette) +	
	theme_bw() +
	geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																z = as.vector( ifelse( m_diff < 0 , 0 , m_diff ) ) ) ,
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="red" ,size=0.5) +
	geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																z = -1*as.vector( ifelse( m_diff > 0 , 0 , m_diff ) ) ) ,
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="blue" ,size=0.5) 

p <-  p_c + p_t 

pdf( file.path(reports.dir,"diff-kde-contour-ingest.pdf") , width = 9 , height = 3.5 )
print(p)
dev.off()

## Cytotrace plots ----
x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
table(x$cluster_id)
stopifnot( identical( names(pca_1_control) , names(pca_2_control) ) )
index <- match( gsub("control_","",names(pca_1_control)) , x$sample_id )
df <- data.frame( x = pca_1_control , y = pca_2_control , cluster_id = x$cluster_id[index])
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(df$cluster_id))
my_palette <- my_palette[!is.na(names(my_palette))]

x <- read_delim("~/Clouds/Dropbox/Data/isc/TE001/TE001-cytotrace.csv","\t")
df$cytotrace <- x$cytotrace_score[ match( gsub("control_","",rownames(df)) , x$cell_id ) ]

## With Countour ----
df <- df %>% arrange(cytotrace)
p_c <- ggplot(df) +
	geom_point( data = df , 
							# aes(x,y), col="grey75" , size = 0.5 ) +
							aes(x,y,color=cytotrace) , size = 1 ) +
	theme_bw() +
	# scale_color_manual(values = my_palette) +
	# scale_color_brewer(palette = "Set1") +
	scale_color_viridis() +
	geom_contour(data=data.frame( x=rep( control_kde_total$x , n_tiles_per_axis ) ,
																y=rep( control_kde_total$y , each=n_tiles_per_axis ),
																z = as.vector( m_ctrl_mean ) ),
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="black" ,size=0.5)
# guides(color = guide_legend(nrow = 1)) +
# theme(legend.position = "bottom")

# df <- data.frame( x = pca_1_treatment , y = pca_2_treatment )
# x <- read_delim("~/Clouds/Dropbox/Data/isc/TE002/TE002-cytotrace.csv","\t")
# df$cytotrace <- x$cytotrace_score[ match( gsub("treatment_","",rownames(df)) , x$cell_id ) ]
# df <- df %>% arrange(cytotrace)
df <- data.frame( x = pca_1_treatment , y = pca_2_treatment )
x <- read_delim("~/Clouds/Dropbox/Data/isc/rad-metadata-ingest-cytoperm.csv",",")
df$cytotrace <- x$cyto_perm[ match( gsub("treatment_","",rownames(df)) , x$cell_id ) ]
df <- df %>% arrange(cytotrace)

df <- df[ !is.na(df$cytotrace) , ]

p_t <- ggplot() +
	geom_point( data = df , 
							aes(x,y,color=cytotrace) , size = 1 ) +
	theme_bw() +
	scale_color_viridis() +
	geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																z = as.vector( ifelse( m_diff < 0 , 0 , m_diff ) ) ) ,
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="red" ,size=0.5) +
	geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																z = -1*as.vector( ifelse( m_diff > 0 , 0 , m_diff ) ) ) ,
							 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="blue" ,size=0.5) 

p <-  p_c + p_t 

pdf( file.path(reports.dir,"diff-kde-contour-cytotrace.pdf") , width = 9 , height = 3.5 )
print(p)
dev.off()

# Plot a list of genes ----
genes_to_plot <- c("Clu","Lgr5","Mki67","Alpi","Dll1","Atho1","Axin2","Nkx2.2","Lyz1","Mmp7","Sox9","Gfi1")

for (a_gene in genes_to_plot)
{
	print_msg_info("Plotting Gene " , a_gene)
	
	x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
	table(x$cluster_id)
	stopifnot( identical( names(pca_1_control) , names(pca_2_control) ) )
	index <- match( gsub("control_","",names(pca_1_control)) , x$sample_id )
	df <- data.frame( x = pca_1_control , y = pca_2_control , cluster_id = x$cluster_id[index])
	tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
										 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
	names(tab10_palette) <- c("blue",'orange',"green","red","purple",
														"brown","pink","gray","olive","cyan")
	my_palette <- tab10_palette
	names(my_palette) <- names(table(df$cluster_id))
	my_palette <- my_palette[!is.na(names(my_palette))]
	
	x <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-cpm.rds")
	if ( !(a_gene %in% rownames(x)) )
		next ;
	x <- log(x+1)
	x <- x[a_gene,] %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
	df[,a_gene] <- x$.[ match( gsub("control_","",rownames(df)) , x$cell_id ) ]
	# df <- df %>% arrange(a_gene)
	df <- df[order(df[,a_gene]),]
	
	## With Countour ----
	p_c <- ggplot(df) +
		geom_point( data = df , 
								# aes(x,y), col="grey75" , size = 0.5 ) +
								aes_string("x","y",color=a_gene) , size = 1 ) +
		theme_bw() +
		# scale_color_manual(values = my_palette) +
		# scale_color_brewer(palette = "Set1") +
		scale_color_viridis() +
		geom_contour(data=data.frame( x=rep( control_kde_total$x , n_tiles_per_axis ) ,
																	y=rep( control_kde_total$y , each=n_tiles_per_axis ),
																	z = as.vector( m_ctrl_mean ) ),
								 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="black" ,size=0.5)
	# guides(color = guide_legend(nrow = 1)) +
	# theme(legend.position = "bottom")
	
	df <- data.frame( x = pca_1_treatment , y = pca_2_treatment )
	x <- readRDS("~/Clouds/Dropbox/Data/isc/TE002/TE002-cpm.rds")
		if ( !(a_gene %in% rownames(x)) )
		next ;
	x <- log(x+1)
	x <- x[a_gene,] %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
	df[,a_gene] <- x$.[ match( gsub("treatment_","",rownames(df)) , x$cell_id ) ]
	# df <- df %>% arrange(a_gene)
	df <- df[order(df[,a_gene]),]
	
	p_t <- ggplot() +
		geom_point( data = df , 
								aes_string("x","y",color=a_gene) , size = 1 ) +
		theme_bw() +
		scale_color_viridis() +
		geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																	y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																	z = as.vector( ifelse( m_diff < 0 , 0 , m_diff ) ) ) ,
								 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="red" ,size=0.5) +
		geom_contour(data=data.frame( x=rep( treatment_kde_total$x , n_tiles_per_axis ) ,
																	y=rep( treatment_kde_total$y , each=n_tiles_per_axis ),
																	z = -1*as.vector( ifelse( m_diff > 0 , 0 , m_diff ) ) ) ,
								 aes(x=x,y=y,z=z, alpha=after_stat(level)), col="blue" ,size=0.5) 
	
	p <-  p_c + p_t 
	
	pdf( file.path(reports.dir, paste0("diff-kde-contour-",a_gene,".pdf") ), width = 9 , height = 3.5 )
	print(p)
	dev.off()
	
	
}

## Dot Plot Gene Expression ----
print_msg_info(">>> Printing dot plot for TE001")
{
	
	x <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-cpm.rds")
	x <- log(x+1)
	
	gex_and_metadata <- x %>%
		as.data.frame() %>%
		rownames_to_column("gene") %>%
		pivot_longer(cols = !c("gene") , names_to = "cell_id" , values_to = "expression" )
	gex_and_metadata <- gex_and_metadata %>% pivot_wider(names_from = "gene",values_from = "expression")
	
	x <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
	table(x$cluster_id)
	gex_and_metadata <- left_join( gex_and_metadata , x , by=c("cell_id"="sample_id") )
	
	tibble2plot <- gex_and_metadata
	table(tibble2plot$cluster_id)
	tibble2plot <- tibble2plot[ !is.na(tibble2plot$cluster_id) , ]
	
	my_threshold <- 0.1
	gene_list.available <- genes_to_plot[ genes_to_plot %in% colnames(tibble2plot) ]
	df_summary <- tibble2plot %>%
		dplyr::select( cluster_id,gene_list.available  ) %>%
		pivot_longer(cols = gene_list.available , names_to = "gene",values_to = "expr")
	# df_summary$score <- lgt(df_summary$score)
	df_summary <- df_summary %>%
		group_by(cluster_id, gene) %>%
		summarise(Avg = mean(expr),
							Pct = sum(expr >= my_threshold) / length(expr) * 100)
	df_summary$gene_f <- factor(x = df_summary$gene , levels = sort(gene_list.available) , ordered = TRUE)
	# df_summary$cluster_id <- factor(x = df_summary$cluster_id ,
	# levels = rev(c("stem-1","stem-2",
	# "absorbitive-1","absorbitive-2",
	# "secretory-1","secretory-2","secretory-3")) , ordered = TRUE)
	
	dot_plot_2 <- ggplot(df_summary, aes(x=gene_f, y=cluster_id)) +
		geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
		scale_size("% detected", range = c(0,6)) +
		scale_fill_gradient2(low="gray",high = "tomato1",mid="white",midpoint=1,
											 guide = guide_colorbar(ticks.colour = "black",
											 											 frame.colour = "black"),
											 name = "Average\nexpression") +		
		# scale_fill_gradientn(colours = viridis(5),
		# 										 guide = guide_colorbar(ticks.colour = "black",
		# 										 frame.colour = "black"),
		# 										 name = "Average\nExpression\n") +												 
		# scale_fill_gradient2(low="steelblue1",high = "tomato1",mid="white",midpoint=my_threshold,
												 # guide = guide_colorbar(ticks.colour = "black",
												 # frame.colour = "black"),
												 # name = "Average\nactivity\n-Log(p)") +
		ylab("Cluster") + xlab("") +
		# facet_grid(~pas_cluster_id_f) +
		theme_bw() +
		theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
					axis.text.y = element_text(size=12, color="black"),
					axis.title = element_text(size=14) ,
					panel.grid.major = element_blank(), panel.grid.minor = element_blank()
		)
	
	pdf( file = file.path(reports.dir,"TE001-expression-gene-list-dotplot.pdf") , width = 5 , height = 2.75 )
	print(dot_plot_2)
	dev.off()
	
}
print_msg_info(">>> Printing dot plot for TE002")
{
	
	x <- readRDS("~/Clouds/Dropbox/Data/isc/TE002/TE002-cpm.rds")
	x <- log(x+1)
	
	gex_and_metadata <- x %>%
		as.data.frame() %>%
		rownames_to_column("gene") %>%
		pivot_longer(cols = !c("gene") , names_to = "cell_id" , values_to = "expression" )
	gex_and_metadata <- gex_and_metadata %>% pivot_wider(names_from = "gene",values_from = "expression")
	
	x <- read_delim("~/Clouds/Dropbox/Data/isc/rad-metadata-ingest-cytoperm.csv",",")
	table(x$iter_cluster_id)
	gex_and_metadata <- left_join( gex_and_metadata , x )
	
	tibble2plot <- gex_and_metadata
	table(tibble2plot$iter_cluster_id)
	tibble2plot <- tibble2plot[ !is.na(tibble2plot$iter_cluster_id) , ]
	
	my_threshold <- 0.1
	gene_list.available <- genes_to_plot[ genes_to_plot %in% colnames(tibble2plot) ]
	df_summary <- tibble2plot %>%
		dplyr::select( iter_cluster_id,gene_list.available  ) %>%
		pivot_longer(cols = gene_list.available , names_to = "gene",values_to = "expr")
	# df_summary$score <- lgt(df_summary$score)
	df_summary <- df_summary %>%
		group_by(iter_cluster_id, gene) %>%
		summarise(Avg = mean(expr),
							Pct = sum(expr >= my_threshold) / length(expr) * 100)
	df_summary$gene_f <- factor(x = df_summary$gene , levels = sort(gene_list.available) , ordered = TRUE)
	# df_summary$cluster_id <- factor(x = df_summary$cluster_id ,
	# levels = rev(c("stem-1","stem-2",
	# "absorbitive-1","absorbitive-2",
	# "secretory-1","secretory-2","secretory-3")) , ordered = TRUE)
	
	dot_plot_2 <- ggplot(df_summary, aes(x=gene_f, y=iter_cluster_id)) +
		geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
		scale_size("% detected", range = c(0,6)) +
		scale_fill_gradient2(low="gray",high = "tomato1",mid="white",midpoint=1,
											 guide = guide_colorbar(ticks.colour = "black",
											 											 frame.colour = "black"),
											 name = "Average\nexpression") +		
		# scale_fill_gradientn(colours = viridis(5),
		# 										 guide = guide_colorbar(ticks.colour = "black",
		# 										 frame.colour = "black"),
		# 										 name = "Average\nExpression\n") +												 
		# scale_fill_gradient2(low="steelblue1",high = "tomato1",mid="white",midpoint=my_threshold,
												 # guide = guide_colorbar(ticks.colour = "black",
												 # frame.colour = "black"),
												 # name = "Average\nactivity\n-Log(p)") +
		ylab("Cluster") + xlab("") +
		# facet_grid(~pas_cluster_id_f) +
		theme_bw() +
		theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
					axis.text.y = element_text(size=12, color="black"),
					axis.title = element_text(size=14) ,
					panel.grid.major = element_blank(), panel.grid.minor = element_blank()
		)
	
	pdf( file = file.path(reports.dir,"TE002-expression-gene-list-dotplot.pdf") , width = 5 , height = 2.75 )
	print(dot_plot_2)
	dev.off()
	
}

## 3D Plots ----
tibble_to_plot <- tibble( cell_id = vp_seurat$cell_id , 
													cluster_id = vp_seurat$cluster_id , 
													condition = vp_seurat$sample_id ,
													UMAP_1 = vp_seurat@reductions$densmap@cell.embeddings[,"densmap_1"] ,
													UMAP_2 = vp_seurat@reductions$densmap@cell.embeddings[,"densmap_2"] ,
													UMAP_3 = vp_seurat@reductions$densmap@cell.embeddings[,"densmap_3"] )

tibble_to_plot$cluster_id <- factor(tibble_to_plot$cluster_id)
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#bcbd22","#7f7f7f","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","olive","gray","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(tibble_to_plot$cluster_id))
my_palette <- my_palette[!is.na(names(my_palette))]

library(rgl)
# my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(tibble_to_plot$cluster_id) )
my_colors <- my_palette
tibble_to_plot$cluster_color <- my_colors[ as.numeric(tibble_to_plot$cluster_id) ]			

plot3d( 
	x = tibble_to_plot$UMAP_1, y = tibble_to_plot$UMAP_2, z = tibble_to_plot$UMAP_3,
	col = tibble_to_plot$cluster_color , 
	# type = 's',
	radius = 10,
	xlab="UMAP 1", ylab="UMAP 2", zlab="UMAP 3" )

my_colors <- c("blue","red")
tibble_to_plot$cluster_color <- my_colors[ as.numeric(factor(tibble_to_plot$condition)) ]

plot3d( 
	x = tibble_to_plot$UMAP_1, y = tibble_to_plot$UMAP_2, z = tibble_to_plot$UMAP_3,
	col = tibble_to_plot$cluster_color , 
	# type = 's',
	radius = 10,
	xlab="UMAP 1", ylab="UMAP 2", zlab="UMAP 3" )

