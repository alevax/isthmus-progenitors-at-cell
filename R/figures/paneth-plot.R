
# system.time({source("sources/figures/paneth-plot.R")})

require(tidyverse)
library(Seurat)

source("../vaxtools/R/utils.R")
create_workspace("paneth-plots")

filename <- "~/Clouds/Dropbox/Data/isc/bottcher-2021/ctrl1/ctrl1-cpm.rds"
bottcher.cpm <- readRDS(filename)
dim(bottcher.cpm)
# bottcher.cpm <- sc_data@assays$RNA@counts %>% RelativeCounts(1e6) %>% as.matrix()

filename <- "~/Clouds/Dropbox/Data/isc/ayyaz-2019/C05/C05-cpm.rds"
ayyaz.cpm <- readRDS(filename)
dim(ayyaz.cpm)
# ayyaz.cpm <- sc_data@assays$RNA@counts %>% RelativeCounts(1e6) %>% as.matrix()

filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-cpm.rds"
te001.cpm <- readRDS(filename)
dim(te001.cpm)
# te001.cpm <- sc_data@assays$RNA@counts %>% RelativeCounts(1e6) %>% as.matrix()

# paneth_marker <- c("Dclk1")
paneth_marker <- c("Nupr1", "Dll4", "Mmp7", "Lyz1")
# paneth_marker <- c("Dll1")
# paneth_marker <- c("Defa22", "Defa24")
# paneth_marker <- c("Lyz1")
# paneth_marker <- c("Mptx1")
x <- bottcher.cpm %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::filter(gene %in% paneth_marker) %>% as.tibble() %>% 
	pivot_longer(cols=-c("gene") , names_to = "cell_id" , values_to = "expression" )
x$dataset <- "bottcher"
x$cell_id <- paste0(x$dataset,"_",x$cell_id)
x$log_expression <- log10(x$expression+1)

y <- ayyaz.cpm %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::filter(gene %in% paneth_marker) %>% as.tibble() %>%
	pivot_longer(cols=-c("gene") , names_to = "cell_id" , values_to = "expression" )
y$dataset <- "ayyaz"
y$cell_id <- paste0(y$dataset,"_",y$cell_id)
y$log_expression <- log10(y$expression+1)

z <- te001.cpm %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::filter(gene %in% paneth_marker) %>% as.tibble() %>%
	pivot_longer(cols=-c("gene") , names_to = "cell_id" , values_to = "expression" )
z$dataset <- "te001"
z$cell_id <- paste0(z$dataset,"_",z$cell_id)
z$log_expression <- log10(z$expression+1)

my_tibble <- do.call("rbind",list(x,y,z))

sum(x$log_expression > 2) / length(x$log_expression) * 100
sum(y$log_expression > 2) / length(y$log_expression) * 100
sum(z$log_expression > 2) / length(z$log_expression) * 100

res <- my_tibble %>% group_by(cell_id) %>% summarize( dataset , 
																											m_expression = mean(expression) ,
																											# g_up_counts = sum(log_expression > 2) ) %>%																											
																											g_up_counts = sum(expression != 0) ) %>%
	distinct()

res %>% group_by(dataset) %>% summarize(sum(g_up_counts==length(paneth_marker)))

df_to_plot <- res %>%
	group_by(dataset) %>%
	summarise(cnt = n() , paneth_counts = ifelse( g_up_counts == length(paneth_marker) , "Paneth" , "other" ) ) %>%
	mutate(paneth_freq = round( sum(paneth_counts=="Paneth") / cnt * 100, 3)) %>%
	mutate(other_freq = round( sum(paneth_counts=="other") / cnt * 100, 3)) %>%
	# mutate(freq = sum(paneth_counts=="Paneth")) %>%
	# mutate(cazzo = sum(cnt)) %>%
	dplyr::distinct() %>%
	arrange(desc(paneth_freq))

df_to_plot <- df_to_plot %>% pivot_longer(cols = c(paneth_freq,other_freq) , names_to = "is_paneth" , values_to = "percent" )

## Plotting Paneth as boxplot ----
p <- ggplot( df_to_plot , aes(dataset, percent , 
														 fill = is_paneth
) ) +
	geom_bar(position="fill", stat="identity") +
	scale_fill_brewer(palette = "Set2") +
	# geom_text(aes(y = percent, label = percent, group = is_paneth), color = "white") +
	# scale_fill_manual(values = my_palette) +
	# scale_x_discrete( position = "top" ) +
	theme_minimal() +
	theme( axis.text.x = element_text(angle = 90, hjust = 0, face = "bold") )

pdf( file = file.path( reports.dir , paste0("paneth-barplots.pdf")) , width = 2.5 , height = 3 )
print(p)
dev.off()	


sum( colSums( te001.cpm[paneth_marker,] > 0 ) == 4 )
paneth_cell_ids <- colnames(te001.cpm)[ colSums( te001.cpm[paneth_marker,] > 0 ) == 4 ]

my_metadata <- read.csv2("~/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-ingest-with-cluster-ids.csv",sep = ",") %>% as_tibble()
my_metadata$UMAP_1 <- as.numeric(my_metadata$UMAP_1)
my_metadata$UMAP_2 <- as.numeric(my_metadata$UMAP_2)
my_metadata$UMAP_3 <- as.numeric(my_metadata$UMAP_3)
my_metadata$iter_cluster_id_with_paneth <- my_metadata$iter_cluster_id
my_metadata$iter_cluster_id_with_paneth <- ifelse( my_metadata$cell_id %in% paneth_cell_ids , "z_paneth" , my_metadata$iter_cluster_id_with_paneth) 

write_csv(my_metadata,"~/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-ingest-with-cluster-ids-with-paneth-cluster.csv")

# my_metadata$iter_cluster_id_with_paneth <- factor(my_metadata$iter_cluster_id_with_paneth)

## Fixing palette colors ----
tab10_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
									 "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
names(tab10_palette) <- c("blue",'orange',"green","red","purple",
													"brown","pink","gray","olive","cyan")
my_palette <- tab10_palette
names(my_palette) <- names(table(my_metadata$iter_cluster_id_with_paneth))

my_palette <- my_palette[1:8]
my_palette[8] <- "black"

## Plotting Paneth as boxplot ----
p <- ggplot( my_metadata , aes(UMAP_1, UMAP_2 , 
														 color = iter_cluster_id_with_paneth
) ) +
	geom_point() +
	scale_color_manual(values = my_palette) +
	# geom_text(aes(y = percent, label = percent, group = is_paneth), color = "white") +
	# scale_fill_manual(values = my_palette) +
	# scale_x_discrete( position = "top" ) +
	theme_void() +
	theme( axis.text.x = element_text(angle = 90, hjust = 0, face = "bold") )

pdf( file = file.path( reports.dir , paste0("TE001-UMAP-inter-cluster-id-with-paneth.pdf")) , width = 6 , height = 6 )
	print(p)
dev.off()	

	isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
	my_sample_id <- "TE001"
	# seurat_viper_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-viper-analysis-with-metacell-data.rds") )
	# vp_seurat <- readRDS(seurat_viper_analysis_data_filename)
	# vp_seurat@assays$VIPER@scale.data %>% as.matrix()

	counts <- readRDS(file.path( isc_data_dir , paste0(my_sample_id,"-cpm.rds") ))
	print_msg_info(">>> >> Statically assigning networks to clusters")
	{
		# clusters_networks.list <- vector("list",3)
		clusters_networks.list <- list()
		clusters_networks.list[["0"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
		clusters_networks.list[["1"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c2_unPruned.rds"))
		clusters_networks.list[["2"]] <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c3and4_unPruned.rds"))
		str(clusters_networks.list,1)
	}
	library(viper)
	clusters_networks.list <- lapply(clusters_networks.list,pruneRegulon,50)
	
	x <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat") , col_names = FALSE )
	y <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat") , col_names = FALSE )
	z <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat") , col_names = FALSE )
	if( my_params.list$is_only_TFS )
	{
		selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))	
	} else {
		selected_regulators <- unique(c( x$X1 , y$X1 , z$X1 ))
	}
	print(length(selected_regulators))

	# TE001_network_full <- readRDS(file.path(isc_data_dir,"networks/TE001_Apr_29_2022_h15m03/TE001-metacells_unPruned.rds"))
	TE001_network_full <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/TE001_mc_unPruned.rds"))
	# TE001_network_full <- readRDS(file.path(isc_data_dir,"networks/TE001-full_unPruned.rds"))
		
	print_msg_info(">>> Total regulatory proteins in the actual network: " , length(TE001_network_full) )
	index <- names(TE001_network_full) %in% selected_regulators
	print_msg_info(">>> Total regulatory proteins to expect from selected ones: " , sum(index) )
 	TE001_network_full <- TE001_network_full[index]
	print_msg_info(">>> Total regulatory proteins recovered in the network: " , length(TE001_network_full) )
	TE001_network_full <- pruneRegulon(TE001_network_full,50)
	
	index <- colnames(counts) %in% paneth_cell_ids
	ges_m_samples <- rankNorm_with_ref( counts[,index] , ref_mat = counts[,!index] )
	
	# vp <- msviper( ges_m_samples  , TE001_network_full )
	# plot(vp,20)
	
	secretory_network <- clusters_networks.list[["2"]]
	print_msg_info(">>> Total regulatory proteins in the actual network: " , length(secretory_network) )
	index <- names(secretory_network) %in% selected_regulators
	print_msg_info(">>> Total regulatory proteins to expect from selected ones: " , sum(index) )
 	secretory_network <- secretory_network[index]
	print_msg_info(">>> Total regulatory proteins recovered in the network: " , length(secretory_network) )
	secretory_network <- pruneRegulon(secretory_network,50)
	
	vp <- msviper( ges_m_samples  , secretory_network )
	plot(vp,20)
	

pdf( file = file.path( reports.dir , paste0("TE001-paneth-mrs.pdf")) , width = 7 , height = 6 )
	plot(vp,20)
dev.off()	
