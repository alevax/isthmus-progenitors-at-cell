

# source("sources/figures/lineage-mrs-three-dimension-other-samples.R")

library(tidyverse)
library(viper)
library(viridis)

source("../vaxtools/R/utils.R")

create_workspace("lineage-mrs-three-dimension-TE002")

filename <- "~/Clouds/Dropbox/Data/isc/TE002/TE002-seurat-metaviper-analysis-data.rds"
vpmat_seurat <- readRDS(filename)
vpmat <- vpmat_seurat@assays$VIPER@scale.data
dim(vpmat)
viper_clusters <- vpmat_seurat$pas_cluster_id
table(viper_clusters)

# Three cluster solution ----
stouffer_table <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-stouffer-with-metacell.csv")
stouffer_table <- stouffer_table %>% tibble2matrix()
head(stouffer_table)

isc_markers.regulon <- generateRegulonObjectFromProteinActivityMatrix(stouffer_table)
isc_enrichment.mat <- aREA(vpmat,isc_markers.regulon)$nes

my_tibble <- isc_enrichment.mat %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()

# Seven cluster solution ----
stouffer_table <- readxl::read_excel("~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-iter-clustering-viper-analysis-with-metacell.xlsx")
stouffer_table <- stouffer_table %>% tibble2matrix()
head(stouffer_table)

isc_markers.regulon <- generateRegulonObjectFromProteinActivityMatrix(stouffer_table)
isc_enrichment.mat <- aREA(vpmat,isc_markers.regulon)$nes

rownames(isc_enrichment.mat) <- paste0("iter_",rownames(isc_enrichment.mat))

my_tibble_2 <- isc_enrichment.mat %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
my_tibble <- left_join( my_tibble , my_tibble_2 )

cytotrace_table <- read_delim("~/Clouds/Dropbox/Data/isc/TE002/TE002-cytotrace.csv","\t")
cytotrace_table$cluster_id <- NULL
my_tibble <- left_join( my_tibble , cytotrace_table )

clustering_table <- readRDS("~/Clouds/Dropbox/Data/isc/TE002/TE002-seurat-metaviper-analysis-metadata.rds")
my_tibble <- left_join( my_tibble , clustering_table , by = c("cell_id"="sample_id"))

filename <- "~/Clouds/Dropbox/Data/isc/TE002/TE002-stem-mrs-enriched-corr.csv"
x <- read_csv(filename)
my_tibble <- left_join( my_tibble , x , by = c("cell_id"="cell_id"))

filename <- "~/Clouds/Dropbox/Data/isc/TE002/TE002-stem-mrs-gene-enriched-corr.csv"
x <- read_csv(filename)
my_tibble <- left_join( my_tibble , x , by = c("cell_id"="cell_id"))

df <- my_tibble

nas_cells <- is.na(df$cytotrace_score)
print_msg_warn(">>> Found CT SCORE NAs: " , sum(nas_cells) ) 
df <- df[!nas_cells,]

df <- df %>% dplyr::rename("ct_score"="cytotrace_score")


# df$cluster_id <- factor(df$cluster_id)
# 
# my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) )
# df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]

# library(car)
# scatter3d( x = df$Cluster_2,
# 					 y = df$Cluster_3,
# 					 z = df$ct_score,
# 					 groups = df$cluster_id ,
# 					 surface = FALSE ,
# 					 surface.col = RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) ) )


# library(rgl)
# plot3d(
# 	x = df$Cluster_2, y = df$Cluster_3, z = df$ct_score,
# 	col = df$cluster_color ,
# 	# type = 's',
# 	radius = 10,
# 	xlab="Absorptive", ylab="Secretory", zlab="Differentiation Status" )
# rgl.snapshot(filename = file.path(reports.dir,"3dscatter.png"))
# 
# writeWebGL( filename = file.path(reports.dir,"3dscatter.html") ,  width=1000, height=1000)

# filename <- file.path("~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data.rds")
# vpmat_seurat <- readRDS(filename)
# vpmat <- vpmat_seurat@assays$VIPER@scale.data
# dim(vpmat)

## Fixing palette colors ----
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(df$cluster_id))

# df$stem_mrs_corr <- ifelse( df$stem_mrs_corr > -2 & df$stem_mrs_corr < 2 , 0 , df$stem_mrs_corr )

## Plotting CytoTrace ----
print_msg_info("## Plotting CytoTrace ----")
{
	p1 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	p2 <- ggplot(df,aes(Cluster_2,ct_score,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Inferred Cell Potency (CT)")
	
	p3 <- ggplot(df,aes(Cluster_3,ct_score,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Secretory Lineage MR Enrichment") +
		ylab("Inferred Cell Potency (CT)")
	
	p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=ct_score)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		scale_fill_viridis(option = "B") +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	library("gridExtra")
	plot_one <- arrangeGrob(p1, p2, p3, 
													ncol=2, nrow=2, widths=c(5, 5), heights=c(2.5, 2.5) )
	
	library("gridExtra")
	plot_two <- arrangeGrob(p1, p4,
													ncol=2, nrow=1, widths=c(5,5), heights=c(5))
	
	pdf( file.path(reports.dir,"3D-flat-plot-cyto.pdf") , width = 10, height = 6)
	grid::grid.draw(plot_one)
	dev.off()
	
	pdf( file.path(reports.dir,"lineage-flat-plot-cyto.pdf") , width = 12.5, height = 5)
	grid::grid.draw(plot_two)
	dev.off()
	
}

## Plotting MRs of Stemness ----
print_msg_info("## Plotting MRs of Stemness ----")
{
	
	p1 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	# p2 <- ggplot(df,aes(Cluster_2,ct_score,fill=cluster_id)) +
	p2 <- ggplot(df,aes(Cluster_2,stem_mrs_corr,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Inferred Cell Potency (MRs Corr)")
	
	# p3 <- ggplot(df,aes(Cluster_3,ct_score,fill=cluster_id)) +
	p3 <- ggplot(df,aes(Cluster_3,stem_mrs_corr,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Secretory Lineage MR Enrichment") +
		ylab("Inferred Cell Potency (MRs Corr)")
	
	# p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=ct_score)) +
	p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=stem_mrs_corr)) +	
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_viridis(option = "B") +
		scale_fill_gradient2(low = "steelblue" , high = "red" , mid = "white", midpoint = 0) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	library("gridExtra")
	plot_one <- arrangeGrob(p1, p2, p3, 
													ncol=2, nrow=2, widths=c(5, 5), heights=c(2.5, 2.5) )
	
	library("gridExtra")
	plot_two <- arrangeGrob(p1, p4,
													ncol=2, nrow=1, widths=c(5,5), heights=c(5))
	
	pdf( file.path(reports.dir,"3D-flat-plot-mrs-stemness.pdf") , width = 10, height = 6)
	grid::grid.draw(plot_one)
	dev.off()
	
	pdf( file.path(reports.dir,"lineage-flat-plot-mrs-stemness.pdf") , width = 12.5, height = 5)
	grid::grid.draw(plot_two)
	dev.off()
	
}

## Plotting Gene Enriched of Stemness ----
print_msg_info("## Plotting Gene Enriched of Stemness ----")
{
	
	p1 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	# p2 <- ggplot(df,aes(Cluster_2,ct_score,fill=cluster_id)) +
	p2 <- ggplot(df,aes(Cluster_2,stem_gene_vp_corr,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Inferred Cell Potency (Gene Corr)")
	
	# p3 <- ggplot(df,aes(Cluster_3,ct_score,fill=cluster_id)) +
	p3 <- ggplot(df,aes(Cluster_3,stem_gene_vp_corr,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Secretory Lineage MR Enrichment") +
		ylab("Inferred Cell Potency (Gene Corr)")
	
	# p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=ct_score)) +
	p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=stem_gene_vp_corr)) +	
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_viridis(option = "B") +
		scale_fill_gradient2(low = "steelblue" , high = "red" , mid = "white", midpoint = 0) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	library("gridExtra")
	plot_one <- arrangeGrob(p1, p2, p3, 
													ncol=2, nrow=2, widths=c(5, 5), heights=c(2.5, 2.5) )
	
	library("gridExtra")
	plot_two <- arrangeGrob(p1, p4,
													ncol=2, nrow=1, widths=c(5,5), heights=c(5))
	
	pdf( file.path(reports.dir,"3D-flat-plot-gene-stemness.pdf") , width = 10, height = 6)
	grid::grid.draw(plot_one)
	dev.off()
	
	pdf( file.path(reports.dir,"lineage-flat-plot-gene-stemness.pdf") , width = 12.5, height = 5)
	grid::grid.draw(plot_two)
	dev.off()
	
}

# df$cluster_id <- factor(df$cluster_id)
# # my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) )
# my_colors <- my_palette
# df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]
# 
# library(rgl)
# plot3d(
# 	x = df$Cluster_2, y = df$Cluster_3, z = df$stem_mrs_corr,
# 	col = df$cluster_color ,
# 	type = 's',
# 	radius = 0.2,
# 	xlab="Absorptive", ylab="Secretory", zlab="Differentiation Status" )
# rgl.snapshot(filename = file.path(reports.dir,"3dscatter.png"))
# 
# writeWebGL( filename = file.path(reports.dir,"3dscatter.html") ,  width=1000, height=1000)

# df %>% filter( df$stem_mrs_corr > 10 )


# df$C1 <- rescale( ifelse(df$Cluster_1<1.96,0,df$Cluster_1) )*100
# df$C2 <- rescale( ifelse(df$Cluster_2<1.96,0,df$Cluster_2) )*100
# df$C3 <- rescale( ifelse(df$Cluster_3<1.96,0,df$Cluster_3) )*100
df$C1 <- rescale( df$Cluster_1) *100
df$C2 <- rescale( df$Cluster_2) *100
df$C3 <- rescale( df$Cluster_3) *100

library(ggtern)
ggtern(data=df,aes(C1,C2,C3)) + 
	geom_mask() +
	geom_point(fill="red",shape=21,size=0.5) + 
	theme_bw() +
	theme_showarrows() +
	theme_clockwise()

