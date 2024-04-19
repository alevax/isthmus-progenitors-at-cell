
# system.time({source("sources/figures/stemness-density-across-datasets.R")})

TE001.Seurat <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")
Ayyaz.Seurat <- readRDS("~/Clouds/Dropbox/Data/isc/ayyaz-2019/C05/C05-seurat-analysis-data.rds")
Bottcher.Seurat <- readRDS("~/Clouds/Dropbox/Data/isc/bottcher-2021/ctrl1/ctrl1-seurat-analysis-data.rds")

## CytoTRACE ----
set.seed(666)
get_cyto_distro <- function(seurat_obj,n_cells = 1000 , n_times = 1) {
	res <- rep( sample( seurat_obj@meta.data$cytotrace_score.ges , 
											size = n_cells , replace = TRUE ) ,
		 n_times )
	return(res)
}

my_tibble <- tibble( TE001 = get_cyto_distro(TE001.Seurat) ,
										 Ayyaz = get_cyto_distro(Ayyaz.Seurat) ,
										 Bottcher = get_cyto_distro(Bottcher.Seurat)
)

my_tibble <- my_tibble %>% pivot_longer(cols = colnames(my_tibble),names_to = "dataset",values_to = "CytoTRACE")

p <- ggplot( my_tibble , aes(x=CytoTRACE, group=dataset,fill=dataset) ) + 
	geom_density(alpha=.3) +
	scale_fill_brewer(palette="Accent",direction = -1) +
	xlab( paste0( "Dataset") ) +
	ylab( paste0( "CytoTRACE Density") ) +
	# coord_equal() +
	theme_minimal() +
	theme( legend.position = "right",
				 text = element_text(face="italic", colour="black" , size=10 ) ,
				 axis.title.x = element_text(face="bold", colour="black" , size=10 ) ,
				 axis.title.y = element_text(face="bold", colour="black" , size=10 ) ,
				 axis.text.x = element_text(face="italic", colour="black" , angle = 0 ,size = 10 ) ,
				 axis.text.y = element_text(face="italic", colour="black" , angle = 0 , size = 10 ) )	

pdf(file.path(reports.dir,"CytoTRACE-score-on-gene-expression-across-datasets.pdf") , width = 4 , height = 3 )
	print(p)
dev.off()

## CytoTRACE GCS ----
set.seed(666)
get_cyto_gcs_distro <- function(seurat_obj,n_cells = 1000 , n_times = 1) {
	res <- rep( sample( seurat_obj@meta.data$cytotrace_gcs.ges , 
											size = n_cells , replace = TRUE ) ,
							n_times )
	return(res)
}

my_tibble <- tibble( TE001 = get_cyto_gcs_distro(TE001.Seurat) ,
										 Ayyaz = get_cyto_gcs_distro(Ayyaz.Seurat) ,
										 Bottcher = get_cyto_gcs_distro(Bottcher.Seurat)
)

my_tibble <- my_tibble %>% pivot_longer(cols = colnames(my_tibble),names_to = "dataset",values_to = "CytoTRACE")

p <- ggplot( my_tibble , aes(x=CytoTRACE, group=dataset,fill=dataset) ) + 
	geom_density(alpha=.3) +
	scale_fill_brewer(palette="Accent",direction = -1) +
	xlab( paste0( "Dataset") ) +
	ylab( paste0( "CytoTRACE GCS Density") ) +
	# coord_equal() +
	theme_minimal() +
	theme( legend.position = "right",
				 text = element_text(face="italic", colour="black" , size=10 ) ,
				 axis.title.x = element_text(face="bold", colour="black" , size=10 ) ,
				 axis.title.y = element_text(face="bold", colour="black" , size=10 ) ,
				 axis.text.x = element_text(face="italic", colour="black" , angle = 0 ,size = 10 ) ,
				 axis.text.y = element_text(face="italic", colour="black" , angle = 0 , size = 10 ) )	

pdf(file.path(reports.dir,"CytoTRACE-gcs-score-on-gene-expression-across-datasets.pdf") , width = 4 , height = 3 )
print(p)
dev.off()

## Stemness Index ----
set.seed(666)
get_stem_distro <- function(seurat_obj,n_cells = 1000 , n_times = 1) {
	res <- rep( sample( seurat_obj@meta.data$stemness_index.ges ,
	# res <- rep( sample( getStemnessIndex( as.matrix( RelativeCounts(seurat_obj@assays$RNA@counts,1e6) ) %>% mouse_to_human() ) , 
											size = n_cells , replace = FALSE ) ,
							n_times )
	return(res)
}

my_tibble <- tibble( TE001 = get_stem_distro(TE001.Seurat) ,
										 Ayyaz = get_stem_distro(Ayyaz.Seurat) ,
										 Bottcher = get_stem_distro(Bottcher.Seurat)
)

my_tibble <- my_tibble %>% pivot_longer(cols = colnames(my_tibble),names_to = "dataset",values_to = "stemness_index")

p <- ggplot( my_tibble , aes(x=stemness_index, group=dataset,fill=dataset) ) + 
	geom_density(alpha=.3) +
	scale_fill_brewer(palette="Accent",direction = -1) +
	xlab( paste0( "Dataset") ) +
	ylab( paste0( "Stemness Index Density") ) +
	# coord_equal() +
	theme_minimal() +
	theme( legend.position = "right",
				 text = element_text(face="italic", colour="black" , size=10 ) ,
				 axis.title.x = element_text(face="bold", colour="black" , size=10 ) ,
				 axis.title.y = element_text(face="bold", colour="black" , size=10 ) ,
				 axis.text.x = element_text(face="italic", colour="black" , angle = 0 ,size = 10 ) ,
				 axis.text.y = element_text(face="italic", colour="black" , angle = 0 , size = 10 ) )	

pdf(file.path(reports.dir,"Stemness-index-on-gene-expression-across-datasets.pdf") , width = 4 , height = 2.5)
print(p)
dev.off()
