## 
# Correlation viper activity and Cytotrace
# -----------------------------------------
# system.time({ source("sources/figures/correlation-analysis-cytotrace-mrs-and-genes.R") })

library(tidyverse)
library(viper)
library(viridis)
library(Seurat)

source("../vaxtools/R/utils.R")

create_workspace("correlation-analysis")

set.seed(666)

## Cytotrace correlation at Protein Activity ----

filename <- file.path("~/Clouds/Dropbox/Data/isc/TE001/TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data.rds")
vpmat_seurat <- readRDS(filename)
vpmat <- vpmat_seurat@assays$VIPER@scale.data
dim(vpmat)
# three_clusters <- vpmat_seurat$pas_cluster_id
# table(three_clusters)
my_metadata <- vpmat_seurat@meta.data %>% dplyr::select(cytotrace=cytotrace_score.ges,pas_cluster_id) %>% 
	as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()

vpmat.tibble <- vpmat %>% as.data.frame() %>% rownames_to_column("protein") %>%
	pivot_longer(cols = -protein , names_to = "cell_id" , values_to = "viper_score" )

vpmat.tibble <- left_join( vpmat.tibble , my_metadata )
vpmat.tibble

cytotrace_threshold <- 0.95
vpmat.tibble %>% dplyr::filter(cytotrace > cytotrace_threshold) %>%
	group_by(cell_id) %>% summarise(n_distinct(cell_id)) %>% nrow()

library(broom)
my_cor_fun <- function(df) cor.test(df$viper_score,df$cytotrace,method="pea") %>% tidy()
# my_cor_fun <- function(df) print(str(df,1))
data_nest <- vpmat.tibble %>% 
	# dplyr::filter(cytotrace > cytotrace_threshold) %>%
	group_by(protein) %>% nest()

data_nest <- mutate(data_nest, model = map(data,my_cor_fun))
corr_pr <- dplyr::select(data_nest,-data) %>% unnest()
corr_pr

p_value_threshold <- 0.05
corr_pr <- mutate(corr_pr,sig=ifelse(p.value<p_value_threshold,"Sig.","Non Sig."))
corr_pr <- corr_pr %>% dplyr::filter(sig == "Sig.")
corr_pr <- corr_pr %>% dplyr::filter( p.value < p_value_threshold )
corr_pr <- corr_pr %>% arrange(desc(estimate))
corr_pr

library(xlsx)
filename <- "proteins-corr-table.xlsx"
write.xlsx(corr_pr %>% 
					 	dplyr::select(protein,estimate,p.value) %>%
					 	as.data.frame(), row.names = FALSE ,
					 file.path(reports.dir,filename))

## Cytotrace correlation at Gene Expression ----

filename <- file.path("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")
gene_seurat <- readRDS(filename)
# expmat <- gene_seurat@assays$SCT@data
expmat <- RelativeCounts(gene_seurat@assays$RNA@counts,scale.factor = 1e4) %>% as.matrix()
dim(expmat)

my_metadata <- gene_seurat@meta.data %>% dplyr::select(cytotrace=cytotrace_score.ges) %>% 
	as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()

expmat.tibble <- expmat %>% as.data.frame() %>% rownames_to_column("gene") %>%
	pivot_longer(cols = -gene , names_to = "cell_id" , values_to = "expression" )

expmat.tibble <- left_join( expmat.tibble , my_metadata )
expmat.tibble

cytotrace_threshold <- 0.95
expmat.tibble %>% filter(cytotrace > cytotrace_threshold) %>% 
	group_by(cell_id) %>% summarise(n_distinct(cell_id)) %>% nrow()

system.time({
	library(broom)
	my_cor_fun <- function(df) { cor.test(df$expression,df$cytotrace,method="pea") %>% tidy() }
	# my_cor_fun <- function(df) print(str(df,1))
	data_nest <- expmat.tibble %>% 
		# filter(gene == "Lgr4") %>%
		# filter(cytotrace > cytotrace_threshold) %>%
		group_by(gene) %>% nest()
	
	data_nest <- mutate(data_nest, model = map(data,my_cor_fun))
	corr_pr <- dplyr::select(data_nest,-data) %>% unnest()
	corr_pr
	
	p_value_threshold <- 0.05
	corr_pr <- mutate(corr_pr,sig=ifelse(p.value<p_value_threshold,"Sig.","Non Sig."))
	corr_pr <- corr_pr %>% filter(sig == "Sig.")
	corr_pr <- corr_pr %>% filter( p.value < p_value_threshold )
	corr_pr <- corr_pr %>% arrange(desc(estimate))
	corr_pr
	
	library(xlsx)
	filename <- "genes-corr-table.xlsx"
	write.xlsx(corr_pr %>% 
						 	dplyr::select(gene,estimate,p.value) %>%
						 	as.data.frame(), row.names = FALSE ,
						 file.path(reports.dir,filename))
	
})




## -----

genes_corr_table <- readxl::read_excel(file.path(reports.dir,"genes-corr-table.xlsx"))
genes_corr_vector <- genes_corr_table$estimate
names(genes_corr_vector) <- genes_corr_table$gene

prots_corr_table <- readxl::read_excel(file.path(reports.dir,"proteins-corr-table.xlsx"))
prots_corr_vector <- prots_corr_table$estimate
names(prots_corr_vector) <- prots_corr_table$protein

# C5_table <- readxl::read_excel("/Users/av2729/Workspace/isc-project/experiments/2022-02-17-stem-signature-C2-vs-C4-5-6/reports/stem-pas-signature-from-C5.xlsx")
# 
# C5_vector <- C5_table$nes
# names(C5_vector) <- C5_table$protein
# 
# tfs_in_common <- intersect( names(prots_corr_vector) , names(C5_vector) )
# x <- prots_corr_vector[tfs_in_common]
# y <- C5_vector[tfs_in_common]
# cor( x , y )
# 
# intersect( names( tail( sort(x) , 50 ) ) ,  names( tail( sort(y) , 50 ) ) )

## Plotting Top correlation markers ----
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")

clustering_table <- read_delim("~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv",",")

## Gene Expression Correlation Plot ----
N <- 10
selected_markers <- names(tail(sort(genes_corr_vector),10))
# selected_markers <- c("Ascl2","Lgr5","Tnfrsf19","Lgr4","Mki67","Notch1","Top2a","Ung","Brca1","Dll1")

df_to_plot <- expmat.tibble %>% filter( gene %in% selected_markers )
df_to_plot <- df_to_plot %>% arrange(desc(cytotrace))

df_to_plot <- left_join(df_to_plot,clustering_table,by=c("cell_id"="sample_id"))

df_to_plot$log2_expr <- log2(df_to_plot$expression+1)

my_palette <- tab10_palette
names(my_palette) <- names(table(df_to_plot$cluster_id))

p <- ggplot( df_to_plot , 
						 aes( x=cytotrace , y=log2_expr , color=cluster_id) ) +
						 # aes( x=cytotrace , y=log2_expr ) ) +
	# geom_point(shape = 21, colour = "grey50",size = 0.25, stroke = 0.25) +
	geom_point(shape = 21,size = 0.25, stroke = 0.25) +
	scale_color_manual(values = my_palette) +
	xlim(c(1,0)) +
	xlab("Cytotrace Score") + ylab("Log2(Expr)") +
	geom_smooth(method = "loess",inherit.aes = TRUE,size=0.5) +
	facet_wrap(~gene,ncol = 2,scales = "free_y") +	
	theme_void() +
	theme( axis.text=element_text(size=10) ,
				 axis.title.x = element_text(color="black", size=15, face="bold") ,
				 axis.title.y = element_text(color="black", size=15, face="bold",angle = 90) ,				 
				 strip.text = element_text(size = 20, face="bold") )

pdf( file.path(reports.dir,"top-genes-corr-cytotrace-ordered.pdf") , width = 8 , height = 8 )
	print(p)
dev.off()


## Protein Activity Correlation Plot ----
N <- 10
selected_markers <- names(tail(sort(prots_corr_vector),10))

df_to_plot <- vpmat.tibble %>% filter( protein %in% selected_markers )
# df_to_plot <- df_to_plot %>% dplyr::rename("cluster_id"="pas_cluster_id")
df_to_plot <- df_to_plot %>% arrange(desc(cytotrace))

df_to_plot <- left_join(df_to_plot,clustering_table,by=c("cell_id"="sample_id"))

my_palette <- tab10_palette
names(my_palette) <- names(table(df_to_plot$cluster_id))

p <- ggplot( df_to_plot , 
						 aes( x=cytotrace , y=viper_score , color=cluster_id) ) +
						 # aes( x=cytotrace , y=viper_score ) ) +
	geom_point(shape = 21, colour = "grey50",size = 0.25, stroke = 0.25) +
	# geom_point(shape = 21,size = 0.25, stroke = 0.25) +
	# scale_fill_brewer(palette = "Set1") +
	scale_color_manual(values = my_palette) +
	xlim(c(1,0)) +
	xlab("Cytotrace Score") + ylab("VIPER Score") +
	geom_smooth(method = "loess",inherit.aes = TRUE,size=0.5) +
	facet_wrap(~protein,ncol = 2,scales = "free_y") +	
	theme_void() +
	theme( axis.text=element_text(size=10) ,
				 axis.title.x = element_text(color="black", size=15, face="bold") ,
				 axis.title.y = element_text(color="black", size=15, face="bold",angle = 90) ,				 
				 strip.text = element_text(size = 20, face="bold") )

pdf( file.path(reports.dir,"top-proteins-corr-cytotrace-ordered.pdf") , width = 8 , height = 8 )
	print(p)
dev.off()

# library(scales)
# ggplot()+
# 	geom_tile(data=corr_pr,
# 						aes(protein,protein,fill=estimate),
# 						size=1,
# 						colour="white")+
# 	geom_tile(data=filter(corr_pr,sig=="Sig."),
# 						aes(protein,protein),
# 						size=1,
# 						colour="black",
# 						fill="transparent")+
# 	geom_text(data=corr_pr,
# 						aes(protein,protein,label=round(estimate,2),
# 								fontface=ifelse(sig=="Sig.","bold","plain")))+
# 	scale_fill_gradient2(breaks=seq(-1,1,0.2), low = "blue", mid = "white", high = "red")+
# 	labs(x="",y="",fill="",p.value="")+
# 	theme_minimal()+
# 	theme(panel.grid.major = element_blank(),
# 				panel.border = element_blank(),
# 				panel.background = element_blank(),
# 				axis.ticks = element_blank())

