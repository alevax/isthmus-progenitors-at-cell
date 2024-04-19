
# system.time({source("sources/figures/differential-density-using-permutation-tests.R")})

library(Seurat)

source("../vaxtools/R/utils.R")

n_tiles_per_axis <- 100
n_iter <- 100

create_workspace("TE001-vs-TE002-diff-kde-irra")
vp_seurat <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE002-viper.rds")

# create_workspace("diff-kde-haplo")
# vp_seurat <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE005-viper.rds")

# create_workspace("diff-kde-abla")
# vp_seurat <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE006-viper.rds")

# pca_1 <- vp_seurat@reductions$pca@cell.embeddings[,"PC_1"]
# pca_2 <- vp_seurat@reductions$pca@cell.embeddings[,"PC_2"]

# write_csv(vp_seurat@reductions$pca@cell.embeddings %>% as.data.frame() %>% rownames_to_column("sample_id"),"~/Downloads/embedding-for-aaron.csv")

pca_1 <- vp_seurat@reductions$umap@cell.embeddings[,"UMAP_1"]
pca_2 <- vp_seurat@reductions$umap@cell.embeddings[,"UMAP_2"]

library(scales)
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
		index <- sample( 1:length(x),size = length(pca_1_treatment) )
		m1 <- kde2d( x = x[index] , y = y[index] , n = n , lims = limit_vector ) 
		index <- sample( 1:length(x),size = length(pca_1_control) )
		m2 <- kde2d( x = x[index] , y = y[index] , n = n , lims = limit_vector ) 
		
		res[[i]] = abs(m1$z - m2$z)
	}

	ret <- array( unlist( res ) , dim = c(n,n,n_iter) )
	# mean <- apply( p , 1:2 , mean ))
	
	return(ret)
	
}

set.seed(666)

control_kde_total <- kde2d(pca_1_control,pca_2_control , n = n_tiles_per_axis , lims = limit_vector )
str(control_kde_total,1)
treatment_kde_total <- kde2d(pca_1_treatment,pca_2_treatment , n = n_tiles_per_axis , lims = limit_vector )
str(treatment_kde_total,1)
test_difference = abs(treatment_kde_total$z - control_kde_total$z)
str(test_difference,1)

n_iter <- 1e3
system.time({ control_difference <- getX( x = pca_1 ,
								 y = pca_2 , 
								 limit_vector = limit_vector ,
								 n_iter = n_iter , n = n_tiles_per_axis )
})
str(control_difference,1)

m <- matrix(data = 1,nrow = n_tiles_per_axis,ncol = n_tiles_per_axis)
for ( i in 1:dim(control_difference)[1] )
{
	for ( j in 1:dim(control_difference)[2] )
	{
		.diff <- sum( test_difference[i,j] > control_difference[i,j,] ) / length(control_difference[i,j,])
		if (.diff>0) { 
			m[i,j] = .diff 
		}
			
	}
}

m_adj <- apply( m , 1:2 , function(x) p.adjust(x,method = "BY") )
# m_adj <- p.adjust(as.vector(m),method = "BH")
sum( m_adj < 0.01 )

ComplexHeatmap::Heatmap(m_adj,col=viridis(10,direction = -1),cluster_columns = F,cluster_rows = F)
ComplexHeatmap::Heatmap(control_kde_total$z,cluster_columns = F,cluster_rows = F)
ComplexHeatmap::Heatmap(treatment_kde_total$z,cluster_columns = F,cluster_rows = F)


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

lowest <- min(c(m_ctrl,m_diff))
highest <- max(c(m_ctrl,m_diff))

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
							aes(x,y,color=cluster_id) , size = 0.5 ) +
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
							aes(x,y), col="grey75" , size = 0.5 ) +
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

pdf( file.path(reports.dir,"diff-kde-contour.pdf") , width = 9 , height = 4 )
	print(p)
dev.off()

