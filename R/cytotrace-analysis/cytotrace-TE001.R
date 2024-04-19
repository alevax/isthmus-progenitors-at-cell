
# source("sources/ermanno-data-analysis/cytotrace-analysis/cytotrace-TE001.R")

library(dplyr)
library(Seurat)

source("~/Downloads/CytoTRACE/R/zzz.R")
source("~/Downloads/CytoTRACE/R/CytoTRACE.R")
source("~/Downloads/CytoTRACE/R/plotCytoGenes.R")
source("~/Downloads/CytoTRACE/R/plotCytoTRACE.R")

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

my_sample_id <- "TE001"
create_workspace( "Cytotrace-TE001" )

getCountsLocally <- function(sample_id) {
	isc_data_dir <- file.path( "~/Clouds/Dropbox/Data/isc/" , sample_id )
	print_msg_info(">>> >> Opening connection with GEX Data File ...")
	seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(sample_id,"-seurat-analysis-data.rds") )
	my_isc.sdata <- readRDS(seurat_analysis_data_filename)	
	ret <- as.matrix(my_isc.sdata@assays$RNA@counts)
	print(dim(ret))
	return(ret)
}

## Comparison between TE001/TE005/TE006 ----

mat_TE001 <- getCountsLocally("TE001")

my_sample_list <- list("WT"=mat_TE001)
str(my_sample_list,1)

my_sample_list <- lapply( 1:length(my_sample_list) , 
													function(i) {
														ret <- my_sample_list[[i]]
														colnames(ret) <- paste0( names(my_sample_list)[i] , "_" , colnames(my_sample_list[[i]]) )
														return(ret)
													})

features_in_common <- Reduce( "intersect" , lapply(my_sample_list,rownames) )
print_msg_info( "Total number of features in common: # " , length(features_in_common) )
features_in_common <- sort(features_in_common)
my_big_mat <- do.call( cbind , lapply( my_sample_list , function(x) x[features_in_common,] ) )
print(dim(my_big_mat))

print_msg_info(">>> >> Running CytoTrace  ...")
cytotrace.data <- CytoTRACE(my_big_mat , enableFast = FALSE , ncores = 4 )

# str(cytotrace.data,1)
# cor( cytotrace.data$CytoTRACE , cytotrace.data$CytoTRACErank )
# plot( cytotrace.data$CytoTRACE , cytotrace.data$CytoTRACErank )
tail( sort( cytotrace.data$cytoGenes ) , 20 )

stemness_index <- getStemnessIndex(my_big_mat %>% mouse_to_human())	

my_tibble <- tibble( sample_id = gsub("(.*)_(.*)","\\1",names(cytotrace.data$CytoTRACE) ) ,
										 index = gsub("(.*)_(.*)","\\2",names(cytotrace.data$CytoTRACE) ),
										 cytotrace = cytotrace.data$CytoTRACE ,
										 stemness_index = stemness_index , 
)

write_tsv( my_tibble %>% filter(sample_id %in% "WT") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE001-only-cytotrace.tsv" )
# write_tsv( my_tibble %>% filter(sample_id %in% "Rad") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE002-cytotrace.tsv" )
# write_tsv( my_tibble %>% filter(sample_id %in% "Haplo") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE005-cytotrace.tsv" )
# write_tsv( my_tibble %>% filter(sample_id %in% "Ablation") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE006-cytotrace.tsv" )

TE001_metadata <- read_delim("~/Clouds/Dropbox/Data/isc/TE001-metadata-ingest.csv")

TE001_metadata$cell_id <- TE001_metadata$cell_id...1
TE001_metadata$cell_id...1 <- NULL
TE001_metadata$cell_id...2 <- NULL

plot( TE001_metadata$cytotrace , my_tibble$cytotrace[ match( TE001_metadata$cell_id , my_tibble$index  ) ] )

TE001_metadata <- left_join( TE001_metadata , my_tibble %>% dplyr::select(index,cytotrace,stemness_index) , 
					 by=c('cell_id'='index') , suffix = c("","_te001"))
# TE001_metadata$cytotrace_te001 <- my_tibble$cytotrace[ match( TE001_metadata$cell_id , my_tibble$index  ) ]

TE001_metadata$cytotrace_te001_rescaled <- ifelse(TE001_metadata$cytotrace_te001<0.95,0,1)

library(ggplot2)
library(viridis)
TE001_metadata <- TE001_metadata %>% arrange(cytotrace_te001)
p1 <- ggplot( data = TE001_metadata , aes( UMAP_1 , UMAP_2 , color = cytotrace_te001 ) ) +
	geom_point( alpha = 0.75 , size = 1 ,stroke = 0.25 ) +
	# geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
	scale_color_viridis(option = "B") +	
	# xlim(min(TE001_metadata$UMAP_1)-2,max(TE001_metadata$UMAP_1)+2) +
	# ylim(min(TE001_metadata$UMAP_2)-2,max(TE001_metadata$UMAP_2)+2) + 
	theme_void() +
	xlab("UMAP 1") +
	ylab("UMAP 2") +
	theme( 
		# axis.text.x = element_text(face = "plain",size = 9) , 
				 # axis.text.y = element_text(face = "plain",size = 9) ,
				 axis.title.x = element_text(face = "bold",size = 12) ,
				 axis.title.y = element_text(face = "bold",size = 12) 
				 # axis.line.x = element_line(size=0.75) ,
				 # axis.line.y = element_line(size=0.75) 
	)		

TE001_metadata <- TE001_metadata %>% arrange(cytotrace_te001)
p2 <- ggplot( data = TE001_metadata , aes( UMAP_1 , UMAP_2 , color = cytotrace_te001_rescaled ) ) +
	geom_point( alpha = 0.75 , size = 1 ,stroke = 0.25 ) +
	# geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
	scale_color_viridis(option = "B") +	
	# xlim(min(TE001_metadata$UMAP_1)-2,max(TE001_metadata$UMAP_1)+2) +
	# ylim(min(TE001_metadata$UMAP_2)-2,max(TE001_metadata$UMAP_2)+2) + 
	theme_void() +
	xlab("UMAP 1") +
	ylab("UMAP 2") +
	# coord_cartesian( xlim = c(min(TE001_metadata$UMAP_1)-2,max(TE001_metadata$UMAP_1)+2 ) , 
									 # ylim = c(min(TE001_metadata$UMAP_2)-2,max(TE001_metadata$UMAP_2)+2) )  +
	theme( 
		# axis.text.x = element_text(face = "plain",size = 9) , 
				 # axis.text.y = element_text(face = "plain",size = 9) ,
				 axis.title.x = element_text(face = "bold",size = 12) ,
				 axis.title.y = element_text(face = "bold",size = 12) 
				 # axis.line.x = element_line(size=0.75,lineend = "butt") ,
				 # axis.line.y = element_line(size=0.75) 
	)		

# library(cowplot)
# umap_plot <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
library(ggpubr)
umap_plot <- ggarrange(p1,p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
pdf( file = file.path( reports.dir , paste0(my_sample_id,"-pas-umap-cytotrace-on-TE001.pdf") ) , width = 8 , height = 4 )
	print(umap_plot)
dev.off()	

write_csv(TE001_metadata,"~/Clouds/Dropbox/Data/isc/TE001/TE001-metadata-ingest-with-cluster-ids.csv")

# p3 <- ggplot( data = TE001_metadata , aes( UMAP_1 , UMAP_2 , color = iter_cluster_id ) ) +
# 	geom_point( alpha = 0.75 , size = 1 ,stroke = 0.25 ) +
# 	# geom_density_2d(lineend = "butt", color = "black",size=0.25,alpha=0.75) +
# 	# scale_color_viridis(option = "B") +
# 	# xlim(min(TE001_metadata$UMAP_1)-2,max(TE001_metadata$UMAP_1)+2) +
# 	# ylim(min(TE001_metadata$UMAP_2)-2,max(TE001_metadata$UMAP_2)+2) +
# 	theme_void() +
# 	xlab("UMAP 1") +
# 	ylab("UMAP 2") +
# 	# coord_cartesian( xlim = c(min(TE001_metadata$UMAP_1)-2,max(TE001_metadata$UMAP_1)+2 ) ,
# 	# ylim = c(min(TE001_metadata$UMAP_2)-2,max(TE001_metadata$UMAP_2)+2) )  +
# 	theme(
# 		# axis.text.x = element_text(face = "plain",size = 9) ,
# 		# axis.text.y = element_text(face = "plain",size = 9) ,
# 		axis.title.x = element_text(face = "bold",size = 12) ,
# 		axis.title.y = element_text(face = "bold",size = 12)
# 		# axis.line.x = element_line(size=0.75,lineend = "butt") ,
# 		# axis.line.y = element_line(size=0.75)
# 	)
