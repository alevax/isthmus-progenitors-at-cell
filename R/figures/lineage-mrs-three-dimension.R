# system.time({source("sources/figures/lineage-mrs-three-dimension.R")})

library(tidyverse)
library(viper)
library(viridis)

source("../vaxtools/R/utils.R")

create_workspace("lineage-mrs-three-dimension-TE001")

# filename <- file.path("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-viper-analysis-with-metacell-data.rds")
# filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data.rds"
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-with-paneth.rds"
vpmat_seurat <- readRDS(filename)
vpmat <- vpmat_seurat@assays$VIPER@scale.data
dim(vpmat)
three_clusters <- vpmat_seurat$pas_cluster_id
table(three_clusters)

# stouffer_table <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-stouffer-with-metacell.csv")
stouffer_table <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-stouffer-with-metacell-3-clusters.csv")
stouffer_table <- stouffer_table %>% tibble2matrix()
head(stouffer_table)

isc_markers.regulon <- generateRegulonObjectFromProteinActivityMatrix(stouffer_table,n_top = 25)
isc_enrichment.mat <- aREA(vpmat,isc_markers.regulon)$nes

my_tibble <- isc_enrichment.mat %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
cytotrace_table <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-cytotrace-table.csv")
cytotrace_table$cluster_id <- NULL
my_tibble <- left_join( my_tibble , cytotrace_table )

# clustering_table <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv")
# my_tibble <- left_join( my_tibble , clustering_table , by = c("cell_id"="sample_id"))
table(vpmat_seurat@meta.data$pas_cluster_id)
clustering_table <- tibble( cell_id = vpmat_seurat@meta.data$cell_id , cluster_id = vpmat_seurat@meta.data$pas_cluster_id )
my_tibble <- left_join( my_tibble , clustering_table , by = c("cell_id"="cell_id"))

filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-stem-mrs-enriched-corr.csv"
x <- read_csv(filename)
my_tibble <- left_join( my_tibble , x , by = c("cell_id"="cell_id"))
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-stem-mrs-gene-enriched-corr.csv"
x <- read_csv(filename)
my_tibble <- left_join( my_tibble , x , by = c("cell_id"="cell_id"))

filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-cpm.rds"
x <- readRDS(filename)
x <- x["Lgr5",] %>% as.data.frame() %>% rownames_to_column("cell_id")
colnames(x) <- c("cell_id","Lgr5_expr")
my_tibble <- left_join( my_tibble , x )

df <- my_tibble

nas_cells <- is.na(df$ct_score)
print_msg_warn(">>> Found CT SCORE NAs: " , sum(nas_cells) ) 
df <- df[!nas_cells,]


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
# tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
# 									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
# names(tab10_palette) <- c("blue",'orange',"green","red","purple",
# 													"brown","pink","gray","olive","cyan")
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
														"#8c564b","#e377c2","#bcbd22","#7f7f7f","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
																"brown","pink","olive","gray","cyan")													
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
	
	df$Lgr5_expr_log <- log10(df$Lgr5_expr+1)
	p3_Lgr5 <- ggplot(df %>% arrange(Lgr5_expr_log),aes(Cluster_3,ct_score,fill=Lgr5_expr_log)) +
		# geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25,stroke=0.15) +
		# scale_fill_brewer(palette = "Set1")+
		# scale_fill_viridis() +
		geom_rug(data = df %>% filter(Lgr5_expr > 0), mapping = aes(alpha=0.15)) +
		scale_fill_gradient2(low = "steelblue" , high = "red" , mid = "white", midpoint = 0) +
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
	
	plot_one <- arrangeGrob(p1, p2, p3_Lgr5, 
													ncol=2, nrow=2, widths=c(5, 5), heights=c(2.5, 2.5) )
	
	pdf( file.path(reports.dir,"lineage-flat-plot-cyto-Lgr5-expr.pdf") , width = 10, height = 6)
	grid::grid.draw(plot_one)
	dev.off()	
	
	# library(ggExtra)
	# p3_Lgr5_marginal <- ggMarginal(data = p3_Lgr5, type="density",margins = 'x')
	# pdf( file.path(reports.dir,"lineage-flat-plot-cyto-Lgr5-expr-marginal.pdf") , width = 4, height = 3)
	# grid::grid.draw(p3_Lgr5_marginal)
	# dev.off()	
	
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
	
	pdf( file.path(reports.dir,"3D-flat-plot-mrs-corr-stemness.pdf") , width = 10, height = 6)
	grid::grid.draw(plot_one)
	dev.off()
	
	pdf( file.path(reports.dir,"lineage-flat-plot-mrs-corr-stemness.pdf") , width = 12.5, height = 5)
	grid::grid.draw(plot_two)
	dev.off()
	
}

## Plotting Three Lineages of MRs Enrichment ----
print_msg_info("## Plotting Three Lineages of MRs Enrichment ----")
{
	
	p1 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	# p2 <- ggplot(df,aes(Cluster_2,ct_score,fill=cluster_id)) +
	p2 <- ggplot(df,aes(Cluster_2,Cluster_1,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Stemness MRs")
	
	# p3 <- ggplot(df,aes(Cluster_3,ct_score,fill=cluster_id)) +
	p3 <- ggplot(df,aes(Cluster_3,Cluster_1,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Secretory Lineage MR Enrichment") +
		ylab("Stemness MRs")
	
	# p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=ct_score)) +
	p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=Cluster_1)) +	
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


# df <- my_tibble %>% dplyr::select(C1=Cluster_1,C2=Cluster_2,C3=Cluster_3)
df <- aREA(vpmat,isc_markers.regulon)$nes %>% t() %>% as.data.frame()
head(df)

write_csv(my_tibble %>% dplyr::select(cell_id,progenitors_lineage=Cluster_1,absorptive_lineage=Cluster_2,secretor_lineage=Cluster_3),"~/Clouds/Dropbox/Data/isc/TE001/TE001-mrs-enrichment-table-3D-plot.csv")

scale_to_1 <- function(x) { scales::rescale( x , to = c(0, 1) ) }
scale_to_gaussian <- function(x) { scales::rescale( x , to = c(-3, 3) ) }
# scale_to_100 <- function(x) { scales::rescale( x , to = c(0, 100) ) }
my_func <- function(x) { -log10( pnorm( x , lower.tail = F )*2 ) }


df_log10 <- apply(df , 1 , my_func ) %>% t()

my_tibble$Cluster_1_mlog10 <- my_tibble$Cluster_1 %>% my_func()
my_tibble$Cluster_2_mlog10 <- my_tibble$Cluster_2 %>% my_func()
my_tibble$Cluster_3_mlog10 <- my_tibble$Cluster_3 %>% my_func()

## Plotting Gene Enriched of Stemness -log10(p)----
print_msg_info("## Plotting Gene Enriched of Stemness -log10(p) ----")
{
	p1 <- ggplot(my_tibble,aes(Cluster_2_mlog10,Cluster_3_mlog10,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Secretory Lineage MR Enrichment")
	
	# p2 <- ggplot(df,aes(Cluster_2,ct_score,fill=cluster_id)) +
	p2 <- ggplot(my_tibble,aes(Cluster_2_mlog10,Cluster_1_mlog10,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Absorptive Lineage MR Enrichment") +
		ylab("Stem MR Enrichment")
	
	# p3 <- ggplot(df,aes(Cluster_3,ct_score,fill=cluster_id)) +
	p3 <- ggplot(my_tibble,aes(Cluster_3_mlog10,Cluster_1_mlog10,fill=cluster_id)) +
		geom_point(alpha=0.95,shape=21,color="gray50",size=1.25) +
		# scale_fill_brewer(palette = "Set1")+
		scale_fill_manual(values=my_palette) +
		theme_light() +
		xlab("Secretory Lineage MR Enrichment") +
		ylab("Stem MR Enrichment")
	
	# p4 <- ggplot(df,aes(Cluster_2,Cluster_3,fill=ct_score)) +
	p4 <- ggplot(my_tibble,aes(Cluster_2_mlog10,Cluster_3_mlog10,fill=Cluster_1_mlog10)) +	
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
	
	pdf( file.path(reports.dir,"3D-flat-plot-gene-stemness-mlog10.pdf") , width = 10, height = 6)
	grid::grid.draw(plot_one)
	dev.off()
	
	pdf( file.path(reports.dir,"lineage-flat-plot-gene-stemness-mlog10.pdf") , width = 12.5, height = 5)
	grid::grid.draw(plot_two)
	dev.off()
	
}

# # df[,c("C1","C2","C3")] <- apply(df[,c("C1","C2","C3")] , 1 , scale_to_100 ) %>% t()
# df$C1 <- ifelse(df$C1<5,0,df$C1)
# df$C2 <- ifelse(df$C2<5,0,df$C2)
# df$C3 <- ifelse(df$C3<5,0,df$C3)
# df <- apply(df,2,function(x)x/sum(x)) %>% as.data.frame()
# head(df)
# set.seed(666)
# # df$C1 <- ifelse(df$C1==0,runif(1,0,1),df$C1) 
# # df$C2 <- ifelse(df$C2==0,runif(1,0,1),df$C2) 
# # df$C3 <- ifelse(df$C3==0,runif(1,0,1),df$C3) 
# # df <- apply(df,1,function(x)x/sum(x)*100) %>% t() %>% as.data.frame()
# head(df)
# df$C1 <- scale_to_1( df$C1 )
# df$C2 <- scale_to_1( df$C2 )
# df$C3 <- scale_to_1( df$C3 )
# 
# # require(scales)
# # df$C1 <- rescale( ifelse( df$Cluster_1 < 0 , 0 , df$Cluster_1 ) )
# # df$C2 <- rescale( ifelse( df$Cluster_2 < 0 , 0 , df$Cluster_2 ) )
# # df$C3 <- rescale( ifelse( df$Cluster_3 < 0 , 0 , df$Cluster_3 ) )
# 
# # df[,c("C1","C2","C3")] <- apply( df[,c("C1","C2","C3")] , 1 , scale_to_1 ) %>% t()
# # df[,c("C1","C2","C3")] <- apply( df[,c("C1","C2","C3")] , 1 , norm_to_1 ) %>% t()
# 
# df <- apply(df,1,scale_to_gaussian) %>% t() %>% as.data.frame()
# df <- apply(df,2,scale_to_1) %>% as.data.frame()
# 
# library(ggtern)
# ggtern(data=df,aes(Cluster_1,Cluster_2,Cluster_3)) +
# # ggtern(data=data.frame(C1=1,C2=0.33,C3=0.33),aes(C1,C2,C3)) +
# 	geom_mask() +
# 	geom_point(fill="red",shape=21,size=0.5) + 
# 	# geom_density() +
# 	theme_bw() +
# 	theme_showarrows() +
# 	theme_clockwise()


# tmp <- data.frame(C1=c(100,20,30),C2=c(10,40,60) , C3=c(70,10,40))
# library(ggtern)
# ggtern(data=tmp,aes(C1,C2,C3)) + 
# 	geom_mask() +
# 	geom_point(fill=c("red","blue","green"),shape=21,size=5) + 
# 	theme_bw() +
# 	theme_showarrows() +
# 	theme_clockwise()

# install.packages('Ternary')
# https://cran.r-project.org/web/packages/Ternary/vignettes/Ternary.html


# # install.packages("archetypes")
# library(archetypes)
# # system.time({ res <- archetypes(vpmat,k = 3,verbose=TRUE) })
# 
# # system.time({ res_2 <- archetypes(vpmat,verbose=TRUE) })
# 
# df <- res$archetypes %>% t() %>% as.data.frame()
# colnames(df) <- c("C1","C2","C3")
# df$cell_id <- rownames(df)
# write_csv(df %>% dplyr::select(cell_id,C1,C2,C3),"~/Downloads/TE001-archetype-k3.csv")
# library(ggtern)
# ggtern(data=df,aes(C1,C2,C3)) +
# 	geom_mask() +
# 	geom_point(fill="red",shape=21,size=0.5) +
# 	# geom_density() +
# 	theme_bw() +
# 	theme_showarrows() +
# 	theme_clockwise()

# print_msg_info(">>> Make 3D Scatter Plot for UMAP")
# {
# 	library(car)
# 	# df <- vp_seurat@meta.data %>% dplyr::select(cluster_id=seurat_clusters,UMAP_1.PAS,UMAP_2.PAS,UMAP_3.PAS)
# 	df$cell_id <- rownames(df)
# 	df <- left_join( df , my_tibble %>% dplyr::select(cell_id,cluster_id) )
# 	df$cluster_id <- factor(df$cluster_id)
# 	
# 	library(rgl)
# 	my_colors <- RColorBrewer::brewer.pal(name = "Set1" , n = nlevels(df$cluster_id) )
# 	df$cluster_color <- my_colors[ as.numeric(df$cluster_id) ]
# 	
# 	plot3d(
# 		x = df$C1, y = df$C2, z = df$C3,
# 		col = df$cluster_color ,
# 		# type = 's',
# 		radius = 10,
# 		xlab="UMAP 1", ylab="UMAP 2", zlab="UMAP 3" )
# 	
# 	writeWebGL( filename = file.path(reports.dir,paste0(my_sample_id,"-3dscatter-UMAP.html") ),  width=1000, height=1000)
# }

