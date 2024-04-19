
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
library(Rmagic)
library(Seurat)

source("../vaxtools/R/utils.R")
# init_python_on_laptop()

# create_workspace(run_dir = "TE001-phate-analysis")

# ## GEX PHATE ----
# print_msg_info(">>> PHATE Analysis on Gene Expression")
# {
# 	my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")
# 	
# 	x <- as.matrix(my_seurat.obj@assays$RNA@counts)
# 	index <- grepl( "mt-" , rownames(x) )
# 	rownames(x)[index]
# 	x <- x[!index,]
# 	dim(x)
# 	
# 	x <- library.size.normalize(x)
# 	x <- sqrt(x)
# 	dim(x)
# 	
# 	# run PHATE
# 	x_PHATE <- phate( t(x) , knn = 51 , 
# 										t = 100 ,
# 										knn.dist.method = "cosine" ,
# 										# knn.dist.method = "euclidean" ,
# 										n.jobs = 5 , verbose = TRUE , seed = 666 )
# 	# str(x_PHATE,1)
# 	# x_PHATE$embedding
# 	df <- x_PHATE$embedding %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
# 	stopifnot( identical( df$cell_id , colnames(x) ) )
# 	
# 	df$cluster_id <- my_seurat.obj$seurat_clusters
# 	
# 	df$Lgr5 <- x["Lgr5",]
# 	df$Olfm4 <- x["Olfm4",]
# 	df$Top2a <- x["Top2a",]
# 	df$Krt19 <- x["Krt19",]
# 	
# 	# l <- levels(df$cluster_id)
# 	# cluster_colors <- c(brewer.pal(9,"Set1"),brewer.pal(6,"Set2"))[1:length(l)]
# 	# names(cluster_colors) <- l
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
# 		scale_fill_discrete() +
# 		# scale_fill_brewer(palette = "Set1") +
# 		# coord_fixed(ratio = 1) +
# 		theme_light()
# 	
# 	pdf( file.path(reports.dir,"TE001-ges-phate-clusters.pdf"))
# 		print(p)
# 	dev.off()
# 	
# 	p <- ggplot( df , aes(PHATE1, PHATE2 , 
# 												fill = Lgr5
# 												# fill = Top2a
# 												# fill = cluster_id
# 	) ) +
# 		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 		# labs(color="Mpo") +
# 		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
# 		scale_fill_viridis() +
# 		# scale_fill_discrete() +
# 		# scale_fill_brewer(palette = "Set1") +
# 		# coord_fixed(ratio = 1) +
# 		theme_light()
# 	
# 	pdf( file.path(reports.dir,"TE001-ges-phate-Lgr5.pdf"))
# 		print(p)
# 	dev.off()	
# 	
# 	p <- ggplot( df , aes(PHATE1, PHATE2 , 
# 												# fill = Lgr5
# 												fill = Top2a
# 												# fill = cluster_id
# 	) ) +
# 		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 		# labs(color="Mpo") +
# 		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
# 		scale_fill_viridis() +
# 		# scale_fill_discrete() +
# 		# scale_fill_brewer(palette = "Set1") +
# 		# coord_fixed(ratio = 1) +
# 		theme_light()
# 	
# 	pdf( file.path(reports.dir,"TE001-ges-phate-Top2a.pdf"))
# 		print(p)
# 	dev.off()		
# 	
# }

## PAS PHATE ----
print_msg_info(">>> PHATE Analysis on Protein Activity")
{
	# my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-viper-analysis-data.rds")
	my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-viper-analysis-with-metacell-data.rds")
	my_metadata.filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-analysis-with-metacell-metadata.rds"
	my_metadata.obj <- readRDS(my_metadata.filename)
	stopifnot( identical( colnames(my_seurat.obj) , my_metadata.obj$sample_id ) )
	
	x <- as.matrix(my_seurat.obj@assays$VIPER@data)
	
	# run PHATE
	x_PHATE <- phate( t(x) ,
										knn = 51 , 
										# ndim = 3 ,
										t = 100 , 
										npca = 10 ,
										# knn.dist.method = "cosine" ,
										knn.dist.method = "euclidean" ,
										n.jobs = 5 , verbose = TRUE , seed = 666 )
	# str(x_PHATE,1)
	# x_PHATE$embedding
	df <- x_PHATE$embedding %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
	stopifnot( identical( df$cell_id , colnames(x) ) )
	df$cluster_id <- my_seurat.obj$seurat_clusters
	
	df$sc_entropy <- my_metadata.obj$sc_entropy_pas
	df$stemness_index <- my_metadata.obj$stemness_index
	
	# df$Lgr5 <- x["Lgr5",]
	# df$Top2a <- x["Top2a",]
	# df$Krt19 <- x["Krt19",]
	
	# l <- levels(df$cluster_id)
	# cluster_colors <- c(brewer.pal(9,"Set1"),brewer.pal(6,"Set2"))[1:length(l)]
	# names(cluster_colors) <- l
	
	# df_3 <- df %>% pivot_longer( cols = contains("PHATE") , names_to = "phate_coordinates" , values_to = "phate_values" ) 
	
	# p <- ggplot( df , aes( x = PHATE2 , y = PHATE3 , 
	# 											# color = Lgr5
	# 											# fill = Top2a
	# 											fill = cluster_id
	# ) ) +
	# 	# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
	# 	# labs(color="Mpo") +
	# 	geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
	# 	# scale_fill_viridis() +
	# 	# scale_fill_discrete() +
	# 	scale_fill_brewer(palette = "Set1") +
	# 	# coord_fixed(ratio = 1) +
	# 	# facet_wrap(~phate_coordinates) +
	# 	theme_light()
	# 
	# pdf( file.path(reports.dir,"TE001-pas-phate-2-3.pdf"))
	# 	print(p)
	# dev.off()	
	
	p <- ggplot( df , aes(PHATE1, PHATE2 , 
												# color = Lgr5
												# fill = Top2a
												fill = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		theme_light()
	
	pdf( file.path(reports.dir,"TE001-pas-phate.pdf"))
		print(p)
	dev.off()
	
	
	# # https://www.rayshader.com/index.html
	# # devtools::install_github("tylermorganwall/rayshader")
	# require("rayshader")
	# plot_gg(p,multicore=TRUE,width=5,height=5,scale=250,windowsize=c(1400,866),
	# 				zoom = 0.65, phi = 30)
	# render_snapshot()
	# # Create the plot
	# p <- plot_ly(
	# 	df, x = ~wt, y = ~hp, z = ~qsec, 
	# 	color = ~am, colors = c('#BF382A', '#0C4B8E')
	# ) %>%
	# 	add_markers() %>%
	# 	layout(
	# 		scene = list(xaxis = list(title = 'Weight'),
	# 								 yaxis = list(title = 'Gross horsepower'),
	# 								 zaxis = list(title = '1/4 mile time'))
	# 	)
	
	
}	
	
## Box Plots of sc_entropy on PHATE ----
print_msg_info(">>> Plots of Entropy on PAS")
{
	
	p <- ggplot( df %>% arrange(sc_entropy) , aes(PHATE1, PHATE2 , 
																								# color = Lgr5
																								# fill = Top2a
																								fill = sc_entropy
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
		scale_fill_viridis(option = "inferno") +
		# scale_fill_discrete() +
		# scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		theme_light()
	
	pdf( file.path(reports.dir,"TE001-pas-phate-sc-entropy.pdf"))
		print(p)
	dev.off()	
	
	p <- ggplot( df %>% arrange(stemness_index) , aes(PHATE1, PHATE2 , 
																										# color = Lgr5
																										# fill = Top2a
																										fill = stemness_index
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
		scale_fill_viridis(option = "viridis") +
		# scale_fill_discrete() +
		# scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		theme_light()
	
	pdf( file.path(reports.dir,"TE001-pas-phate-stemness-index.pdf"))
		print(p)
	dev.off()		
	
	# levels(df$cluster_id)[1] <- "Krt19"
	# levels(df$cluster_id)[2] <- "Lgr5_Prolif"
	# levels(df$cluster_id)[3] <- "Proliferative"
	# levels(df$cluster_id)[4] <- "Atoh1"
	# levels(df$cluster_id)[5] <- "Lgr5"
	# levels(df$cluster_id)[6] <- "Olfm4"
	# levels(df$cluster_id)[7] <- "Dclk1"
	
	p <- ggplot( df , aes(PHATE1, PHATE2 , 
												# color = Lgr5
												# fill = Top2a
												fill = cluster_id 
												# groups = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_point(shape = 21, colour = "grey50",size = 1, stroke = 0.25) +
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		coord_fixed(ratio = 1) +
		# theme_void() +
		theme_light() +
		facet_wrap(~cluster_id)
	
	pdf( file.path(reports.dir,"TE001-pas-phate-clusters-facets.pdf"))
		print(p)
	dev.off()	
	
	df <- df %>% 
		group_by(cluster_id) %>% 
		add_tally(name="cluster_size") %>%
		ungroup()

	total_cells <- df$cluster_size %>% unique() %>% sum()
	df <- df %>% mutate( cluster_cell_perc = round(cluster_size / total_cells * 100 , 2) )
	
	p <- ggplot( df , aes(cluster_id, sc_entropy , 
												# color = sc_entropy ,
												# fill = Top2a
												fill = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_boxplot() +
		geom_text( aes(x = cluster_id, y = 2, 
									 label = paste0(cluster_cell_perc,"%")),
							 colour="black", size = 3) +		
		geom_text( aes(x = cluster_id, y = 2.05 , 
									 label = paste0(cluster_size,"\ncells")),
							 colour="black", size = 3) +				
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		# coord_fixed(ratio = 5) +
		# theme_void() +
		scale_x_discrete( position = "top" ) +
		theme_light() +
		theme( axis.text.x = element_text(angle = 45, hjust = 0, face = "bold") )
	
	pdf( file.path(reports.dir,"TE001-pas-phate-clusters-boxplots-entropy.pdf") , width = 6 )
		print(p)
	dev.off()		
	
	p <- ggplot( df , aes(cluster_id, stemness_index , 
												# color = sc_entropy ,
												# fill = Top2a
												fill = cluster_id
	) ) +
		# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
		# labs(color="Mpo") +
		geom_boxplot() +
		geom_text( aes(x = cluster_id, y = 1.01, 
									 label = paste0(cluster_cell_perc,"%")),
							 colour="black", size = 3) +		
		geom_text( aes(x = cluster_id, y = 1.05 , 
									 label = paste0(cluster_size,"\ncells")),
							 colour="black", size = 3) +				
		# scale_fill_viridis() +
		# scale_fill_discrete() +
		scale_fill_brewer(palette = "Set1") +
		# coord_fixed(ratio = 5) +
		# theme_void() +
		scale_x_discrete( position = "top" ) +
		theme_light() +
		theme( axis.text.x = element_text(angle = 45, hjust = 0, face = "bold") )
	
	pdf( file.path(reports.dir,"TE001-pas-phate-clusters-boxplots-stemness.pdf") , width = 6 )
		print(p)
	dev.off()	
	
	print_msg_warn(">>> Updating metadata file")
	{
		my_metadata.obj <- left_join( my_metadata.obj , df , by = c("sample_id"="cell_id") )
		write_csv( my_metadata.obj , "~/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-with-metacell-all-analysis.csv")
	}
}

# filename <- "/Users/afpvax/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-with-metacell-all-analysis-processed.csv"
# file.info(filename)
# x <- read_csv(filename)
# 
# x <- x %>% dplyr::rename(cell_id=X1)
# dim(x)
# 
# vpmat <- as.matrix(my_seurat.obj@assays$VIPER@data)
# 
# x <- x %>% arrange( latent_time )
# 
# index <- match( x$cell_id , colnames(vpmat) )
# vpmat <- vpmat[,index]
# stopifnot( identical( x$cell_id , colnames(vpmat) ) )
# 
# x
# y <- t(vpmat) %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
# df <- left_join(x,y)
# 
# # df_to_plot <- df %>% filter(cluster_id.x %in% c(6,5,2) )
# df_to_plot <- df 
# 
# ggplot( df_to_plot , 
# 				aes(x=latent_time,y=gene_Lgr5,fill=cluster_id.y,shape=cluster_id.y) ) +
# 	geom_point(shape = 21, colour = "grey50",size = 1, stroke = 0.25) +
# 	scale_fill_brewer(palette = "Set1") +
# 	# geom_smooth(method = "lm") +
# 	theme_light()
# ggplot( df_to_plot , 
# 				aes(x=latent_time,y=Nolc1,fill=cluster_id.y,shape=cluster_id.y) ) +
# 	geom_point(shape = 21, colour = "grey50",size = 1, stroke = 0.25) +
# 	scale_fill_brewer(palette = "Set1") +
# 	geom_smooth(method = "lm",inherit.aes = TRUE) +
# 	theme_light()
# ggplot( df_to_plot , 
# 				aes(x=stemness_index.x,y=Nolc1,fill=cluster_id.y,shape=cluster_id.y) ) +
# 	geom_point(shape = 21, colour = "grey50",size = 1, stroke = 0.25) +
# 	scale_fill_brewer(palette = "Set1") +
# 	geom_smooth(method = "lm",inherit.aes = TRUE) +
# 	theme_light()
# 





