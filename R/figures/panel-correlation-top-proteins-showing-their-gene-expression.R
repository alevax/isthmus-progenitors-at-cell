
# system.time({ source("sources/figures/correlation-analysis-cytotrace-mrs-and-genes.R") })

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

library(ComplexHeatmap)
# lt = lapply(1:20, function(x) cumprod(1 + runif(1000, -x/100, x/100)) - 1)
# ha = rowAnnotation(foo = anno_horizon(lt))	
# ha	
# draw(ha)

# x <- vpmat.tibble %>% filter(cytotrace > 0.8)
# # p <- x %>% group_by(protein) %>% group_map( ~ broom::tidy(lm(cytotrace~viper_score,data = .x) ) )
# p <- x %>% group_by(protein) %>% group_map( function(.x,.y) slope = lm(cytotrace~viper_score,data = .x)$coefficient[[2]] )
# head( sort(unlist(p),decreasing = F) )


# vpmat.tibble %>% filter(cytotrace > 0.95) %>% group_by(protein) %>% summarise(vp_mean=mean(viper_score)) %>% arrange(desc(vp_mean))


require(zoo)
# prots_corr_vector

## Protein Activity Correlation Plot ----
N <- 15
# selected_markers <- sort( doStouffer(mat[,1:25]) , decreasing = TRUE )[1:25]
selected_markers <- names(tail(sort(prots_corr_vector),N))

selected_markers <- c(rev(selected_markers),"Id3","Hnf4g","Atoh1","Spdef","Neurod1")

df_to_plot <- expmat.tibble %>% filter( gene %in% selected_markers )
# df_to_plot <- df_to_plot %>% dplyr::rename("cluster_id"="pas_cluster_id")
df_to_plot <- df_to_plot %>% arrange(desc(cytotrace))

df_to_plot <- left_join(df_to_plot,clustering_table,by=c("cell_id"="sample_id"))

df_to_plot <- df_to_plot %>% dplyr::rename("protein"="gene")

my_list <- vector(mode = "list",length = length(table(df_to_plot$protein)))
names(my_list) <- names(table(df_to_plot$protein))

# index <- match( prots_corr_table$protein[1:N] , names(my_list) )
index <- match( selected_markers , names(my_list) )
my_list <- my_list[index]

for ( a_protein in names(my_list) )
{
	x <- df_to_plot %>% dplyr::filter(protein == a_protein) %>% arrange(desc(cytotrace)) %>% pull(expression)
	my_list[[a_protein]] <- rollapply( x , width = 50, by = 10 , FUN = mean, align = "left")
	# my_list[[a_protein]] <- rescale(my_list[[a_protein]],to = c(0,1))
}

str(my_list,1)
mat <- do.call(rbind,my_list)
dim(mat)

x <- df_to_plot %>% dplyr::filter(protein == a_protein) %>% arrange(desc(cytotrace)) %>% pull(cytotrace)
cyto_rolled <- rollapply( x , width = 50, by = 10 , FUN = mean, align = "left")

colnames(mat) <- paste0( "foo_" , 1:ncol(mat) )
names(cyto_rolled) <- colnames(mat)

col_fun = colorRamp2(quantile(cyto_rolled, seq(0, 1, by = 0.1)), viridis(11))
ha_cyto = HeatmapAnnotation(cytotrace = cyto_rolled, col = list(cytotrace = col_fun), 
														height = unit(3, "mm") ,
														show_legend = FALSE , annotation_name_side = "left" )


require(circlize)

h_corr <- Heatmap( mat , 
									 width = unit(2, "in") , 
									 height = unit(2, "in") ,
									 use_raster = TRUE , raster_quality = 10 ,
									 top_annotation = ha_cyto ,
									 # column_title = "Protein Activity"
									 heatmap_legend_param = list(title = "Gene Expression") ,
									 col = circlize::colorRamp2(c(0, 3), c("white", "darkorange")),
									 border = FALSE , row_gap = unit(2, "mm") ,
									 row_names_side = "left" , show_column_names = FALSE ,
									 cluster_rows = FALSE ,
									 cluster_columns = FALSE, 
									 row_names_gp = gpar(fontsize = 6) )
ht_list <- h_corr

pdf( file.path(reports.dir,"correlation-plot-figure-one.pdf") , width = 6 , height = 4.75 )
	draw(ht_list, ht_gap = unit(2, "mm"))
dev.off()

# source("../vaxtools/R/utils.R")
# source("../vaxtools/R/cross-species-utils.R")
# 
# preppi_table <- read_csv( "~/Clouds/Dropbox/Data/tables/preppi_final600_modified.txt")
# stem_proteins_human <- stem_proteins %>% mouse_to_human()
# 
# preppi_table <- preppi_table %>% filter(symbol_prot1 %in% stem_proteins_human)
# 
# my_percentile <- 0.95
# dim(preppi_table)
# preppi_table <- preppi_table %>% filter( final_score > quantile(preppi_table$final_score,my_percentile) )
# dim(preppi_table)

# umap.tibble <- tibble( cell_id = vpmat_seurat@meta.data$cell_id ,
# 	umap_1 = vpmat_seurat@meta.data$UMAP_1.PAS ,
# 	umap_2 = vpmat_seurat@meta.data$UMAP_2.PAS )

## VIPER better than GEX ----
# umap.tibble <- read_delim("~/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-ingest-with-cluster-ids.csv")
# umap.tibble <- read_delim("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
umap.tibble <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-umap-and-clusters-for-paper.csv")
umap.tibble$cell_id <- umap.tibble$cell_id...1
umap.tibble <- umap.tibble %>% dplyr::select(cell_id,UMAP_1=UMAP_1_scanpy,UMAP_2=UMAP_2_scanpy)

# selected_features <- names(tail(sort(prots_corr_vector),16))
features_table <- readxl::read_excel("~/Clouds/Dropbox/Data/isc/TE001/Supplementary Table 4 (Protein Activity 7-cluster Solution).xlsx")
mat <- features_table %>% tibble2matrix()
selected_features <- Reduce( union , apply( mat , 2 , function(x) { x <- sort(x,decreasing = T) ; names(c(head(x,4))) } ) )
# selected_features <- apply( mat , 2 , function(x) { x <- sort(x,decreasing = T) ; names(c(head(x,5))) } )
# selected_features <- make.unique(selected_features)
length(selected_features)
# selected_features <- names(tail(sort(prots_corr_vector),16))

x <- expmat.tibble %>% filter( gene %in% selected_features)
df_expr <- left_join(x,umap.tibble)
df_expr$cluster <- vpmat.tibble$pas_cluster_id[ match( df_expr$cell_id , vpmat.tibble$cell_id ) ]
df_expr$expression_orig <- df_expr$expression
df_expr$expression <- ifelse( df_expr$expression > 5 , 5 , df_expr$expression )
df_expr$gene_f <- factor(x = df_expr$gene , levels = selected_features , ordered = TRUE)

library(ggrastr)
umap_plot_1 <- ggplot( data = df_expr , aes( UMAP_1 , UMAP_2 , colour = expression ) ) +
	# geom_point( alpha = 0.75 , stroke = 0.5 , size = 0.5 ) +
	geom_point_rast( shape = 16 ,  alpha = 0.5 , size = 0.15 , raster.dpi = 150 ) +
	scale_color_gradient2(aesthetics = "colour" , midpoint = 1 , 
												low="steelblue1",high = "tomato1", mid="white") +
	# scale_color_gradient2(low="white",high = "tomato1") +
	xlab( paste0( "UMAP 1") ) +
	ylab( paste0( "UMAP 2") ) +
	labs(colour='Log(cpm)') +
	ggtitle("Gene Expression") +
	# guides(color=guide_legend(title="log(cpm)")) +
	facet_wrap(~gene_f , nrow = 7 , ncol = 4 ) +
	theme_void() +
	theme( axis.line = element_line(colour = 'black', size = 0.5) ,
				 axis.text = element_text(face = "bold" , size = 8 ) ,
				 axis.title = element_text(face = "bold" , size = 8 ) ,
				 title = element_text(face = "bold" , size = 15 , hjust = 0.5) ,
				 plot.margin = margin(0.25, 0.25, 0.25, 0.25, "in") ,
				 plot.title = element_text(hjust = 0.5 , margin=margin(30,30,30,30) ) , 
				 panel.spacing.x=unit(0.1, "in") , panel.spacing.y=unit(0.1, "in") ,
				 strip.text = element_text(size = 12, face = "bold",hjust = 0.5) 
	)

pdf( file = file.path(reports.dir,"expression-top-markers-umap.pdf") , width = 7 , height = 9 )
	print(umap_plot_1)
dev.off()

x <- vpmat.tibble %>% filter( protein %in% selected_features )
df_prot <- left_join(x,umap.tibble)
df_prot$protein_f <- factor(x = df_prot$protein , levels = selected_features , ordered = TRUE)
df_prot$viper_score_orig <- df_prot$viper_score
df_prot$viper_score <- ifelse( df_prot$viper_score > 5 , 5 , df_prot$viper_score )
df_prot$viper_score <- ifelse( df_prot$viper_score < -5 , -5 , df_prot$viper_score )
# df_prot$viper_score <- ifelse( df_prot$viper_score < 1.964 , 0 , df_prot$viper_score )

my_threshold <- 5
df_prot$score = -log10( pnorm( df_prot$viper_score_orig , lower.tail = F )*2 )
df_prot$score <- ifelse( df_prot$score > 30 , 30 , df_prot$score )
df_prot$score <- ifelse( df_prot$score < my_threshold , 0 , df_prot$score )


library(ggrastr)
umap_plot_2 <- ggplot( data = df_prot , aes( UMAP_1 , UMAP_2 , color = score ) ) +
	# geom_point( alpha = 0.75 , stroke = 0.5 , size = 0.5 ) +
	geom_point_rast( alpha = 0.5 , stroke = 0.15 , size = 0.15 , raster.dpi = 150 ) +
	# scale_color_gradient2(low="steelblue",high = "darkred",mid="white",midpoint=0) +
	scale_color_gradient2(low="steelblue1",high = "tomato1",mid="white",midpoint=my_threshold) +
	xlab( paste0( "UMAP 1") ) +
	ylab( paste0( "UMAP 2") ) +
	labs(color='Activity') +
	ggtitle("Protein Activity") +
	# guides(color=guide_legend(title="log(cpm)")) +
	facet_wrap(~protein_f , nrow = 7 , ncol = 4 ) +
	theme_void() +
	theme( axis.line = element_line(colour = 'black', size = 0.5) ,
				 axis.text = element_text(face = "bold" , size = 8 ) ,
				 axis.title = element_text(face = "bold" , size = 8 ) ,
				 title = element_text(face = "bold" , size = 15 , hjust = 0.5) ,
				 plot.margin = margin(0.25, 0.25, 0.25, 0.25, "in") ,
				 plot.title = element_text(hjust = 0.5 , margin=margin(30,30,30,30) ) , 
				 panel.spacing.x=unit(0.1, "in") , panel.spacing.y=unit(0.1, "in") ,
				 strip.text = element_text(size = 12, face = "bold",hjust = 1,vjust = 1) 
	)

pdf( file = file.path(reports.dir,"activity-top-markers-umap.pdf") , width = 7 , height = 9 )
	print(umap_plot_2)
dev.off()

library(cowplot)
p <- plot_grid(umap_plot_1, umap_plot_2, nrow = 1, ncol = 2)
pdf( file = file.path(reports.dir,"activity-vs-expression-top-markers-umap.pdf") , width = 13 , height = 8 )
	print(p)
dev.off()


## Dot plot ----

markers_table <- tibble( cluster = rep( c("stem-1","stem-2",
															 "absorbitive-1","absorbitive-2",
															 "secretory-1","secretory-2","secretory-3") , each = 4 ) ,
												feature = selected_features )

df_summary <- df_expr %>%
	group_by(cluster, gene) %>%
	summarise(Avg = mean(expression),
						Pct = sum(expression > 0) / length(expression) * 100)
df_summary$gene_f <- factor(x = df_summary$gene , levels = selected_features , ordered = TRUE)
df_summary$cluster_f <- factor(x = df_summary$cluster , 
															 levels = rev(c("stem-1","stem-2",
															 					 "absorbitive-1","absorbitive-2",
															 					 "secretory-1","secretory-2","secretory-3")) , ordered = TRUE)
# df_summary$marker4cluster <- markers_table$cluster[ match( df_summary$gene , markers_table$feature ) ]

dot_plot_1 <- ggplot(df_summary, aes(x=gene_f, y=cluster_f)) +
	geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
	scale_size("% detected", range = c(0,6)) +
	# scale_fill_gradientn(colours = viridisLite::mako(100),
	# scale_fill_gradientn(colours = viridisLite::viridis(10),
	scale_fill_gradient2(low="steelblue1",high = "tomato1",mid="white",midpoint=1,
											 guide = guide_colorbar(ticks.colour = "black",
											 											 frame.colour = "black"),
											 name = "Average\nexpression") +
	ylab("Cluster") + xlab("") +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
				axis.text.y = element_text(size=12, color="black"),
				axis.title = element_text(size=14) ,
				panel.grid.major = element_blank(), panel.grid.minor = element_blank()
	)

pdf( file = file.path(reports.dir,"expression-top-markers-umap-dotplot.pdf") , width = 12 , height = 4 )
	print(dot_plot_1)
dev.off()

## pnorm conversion ----
df_prot$score = -log10( pnorm( df_prot$viper_score_orig , lower.tail = F )*2 )
df_prot$score <- ifelse( df_prot$score > 30 , 30 , df_prot$score )
my_threshold <- 5

df_summary <- df_prot %>%
	group_by(pas_cluster_id, protein) %>%
	summarise(Avg = mean(score),
						Pct = sum(score >= my_threshold) / length(score) * 100)	
	# summarise(Avg = mean(viper_score),
	# 					Pct = sum(viper_score > 2) / length(viper_score) * 100)
df_summary$protein_f <- factor(x = df_summary$protein , levels = selected_features , ordered = TRUE)
df_summary$pas_cluster_id_f <- factor(x = df_summary$pas_cluster_id , 
															 levels = rev(c("stem-1","stem-2",
															 					 "absorbitive-1","absorbitive-2",
															 					 "secretory-1","secretory-2","secretory-3")) , ordered = TRUE)

dot_plot_2 <- ggplot(df_summary, aes(x=protein_f, y=pas_cluster_id_f)) +
	geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
	scale_size("% detected", range = c(0,6)) +
	# scale_fill_gradientn(colours = viridisLite::mako(100),
	scale_fill_gradient2(low="steelblue1",high = "tomato1",mid="white",midpoint=my_threshold,
											 guide = guide_colorbar(ticks.colour = "black",
											 											 frame.colour = "black"),
											 name = "Average\nactivity\n-Log(p)") +
	ylab("Cluster") + xlab("") +
	# facet_grid(~pas_cluster_id_f) +
	theme_bw() +
	theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
				axis.text.y = element_text(size=12, color="black"),
				axis.title = element_text(size=14) ,
				panel.grid.major = element_blank(), panel.grid.minor = element_blank()
	)

pdf( file = file.path(reports.dir,"activity-top-markers-umap-dotplot.pdf") , width = 12 , height = 4 )
	print(dot_plot_2)
dev.off()

p <- plot_grid(dot_plot_1, dot_plot_2, nrow = 2, ncol = 1)
pdf( file = file.path(reports.dir,"activity-vs-expression-top-markers-dotplot.pdf") , width = 12 , height = 9 )
	print(p)
dev.off()

pdf( file = file.path(reports.dir,"activity-vs-expression-top-markers-dotplot-small.pdf") , width = 9 , height = 6 )
	print(p)
dev.off()

