
# source("sources/figures/stem-signature.R")

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

create_workspace("stem-signature-C2-vs-C4-5-6")

require(tidyverse)
require(Seurat)

## Loading ATAC Clusters and CT scores ----

# atac_clusters <- read_csv("~/Downloads/ataclabels.csv") # Ermanno's analysis
# atac_clusters <- atac_clusters %>% dplyr::rename("cell_id"="...1","cluster"="Clusters")
atac_clusters <- read_csv("~/Clouds/Dropbox/Data/isc/TE013/ATAC-Seq/processed/TE013-archr-clusters.csv")
atac_clusters$cell_id <- gsub("TE013#","",atac_clusters$cell_id)

te013_table <- read_csv( "~/Clouds/Dropbox/Data/isc/TE013/TE013-cytotrace.csv" )

te013_table <- left_join(te013_table,atac_clusters)
table(te013_table$cluster_id)

sum(is.na(te013_table$cluster_id))
sum(is.na(te013_table$cytotrace_score))

table(te013_table$cluster_id)
te013_table <- te013_table[ !is.na(te013_table$cluster_id) , ]
te013_table <- te013_table[ !is.na(te013_table$cytotrace_score) , ]
table(te013_table$cluster_id)

p <- ggplot( te013_table , aes(cluster_id, cytotrace_score , 
															 # color = sc_entropy ,
															 # fill = Top2a
															 fill = cluster_id
) ) +
	geom_boxplot() +
	scale_fill_brewer(palette = "Set1") +
	scale_x_discrete( position = "top" ) +
	theme_light() +
	theme( axis.text.x = element_text(angle = 45, hjust = 0, face = "bold") )

pdf( file = file.path( reports.dir , "TE013-atac-boxplot-cytotrace.pdf" ) , width = 4 , height = 2.5 )
	print(p)
dev.off()	

## 43 cell stem signature ----
filename <- "~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-original-data.rds"
gene_seurat <- readRDS(filename)

dim(gene_seurat)
dim(te013_table)

gene_seurat$clusters_atac <- te013_table$cluster_id[ match( colnames(gene_seurat) , te013_table$cell_id ) ]
gene_seurat$cytotrace_score <- te013_table$cytotrace_score[ match( colnames(gene_seurat) , te013_table$cell_id ) ]
table(gene_seurat$clusters_atac)
table(te013_table$cluster_id)
sum(is.na(gene_seurat$clusters_atac))

# gene_seurat <- SetIdent(gene_seurat, value = "unk")
# gene_seurat <- SetIdent(gene_seurat, 
# 												cells = stem_cell_id, 
# 												value = "stem")
# gene_seurat <- SetIdent(gene_seurat, 
# 												cells = terminal_cell_id, 
# 												value = "terminal")
# table(Idents(gene_seurat))

gene_seurat <- SCTransform(gene_seurat, vst.flavor = "v2",seed.use = 666,
													 return.only.var.genes = FALSE , 
													 variable.features.n = nrow(gene_seurat@assays$RNA@data) )
gene_seurat <- PrepSCTFindMarkers(gene_seurat)

table(gene_seurat$clusters_atac)
summary(gene_seurat$mt_percent[gene_seurat$clusters_atac == "C6"])
summary(gene_seurat$mt_percent[gene_seurat$clusters_atac == "C5"])

summary(gene_seurat$cytotrace_score[gene_seurat$clusters_atac == "C6"])
summary(gene_seurat$cytotrace_score[gene_seurat$clusters_atac == "C5"])

summary(gene_seurat$nCount_RNA[gene_seurat$clusters_atac == "C6"])
summary(gene_seurat$nCount_RNA[gene_seurat$clusters_atac == "C5"])
summary(gene_seurat$nFeature_RNA[gene_seurat$clusters_atac == "C6"])
summary(gene_seurat$nFeature_RNA[gene_seurat$clusters_atac == "C5"])

# FindAllMarkers()
dge_table <- FindMarkers(gene_seurat, assay = "SCT", 
												 group.by = "clusters_atac",
												 test.use = "wilcox",
												 ident.1 = "C6",
												 ident.2 = c("C5","C7"),
												 verbose = TRUE, 
												 features = rownames(gene_seurat@assays$SCT@scale.data) ,
												 recorrect_umi = FALSE)
dge_table <- dge_table %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble()
# View(dge_table)
# dge_table$diff <- dge_table$pct.1 - dge_table$pct.2
# sum( dge_table$p_val_adj < 0.1 )
# nrow(dge_table)



## Gene Markers of the 2 stem ----
# filename <- "~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-original-data.rds"
filename <- "~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-analysis-data.rds"
fs::file_info(filename)
gene_seurat <- readRDS(filename)

nrow(te013_table)
te013_table <- te013_table[ complete.cases(te013_table) , ]
nrow(te013_table)

stem_cell_id <- te013_table %>% 
	dplyr::filter(cluster_id == "C5") %>%
		# dplyr::filter(cluster_id == "C2") %>%
	dplyr::filter(cytotrace_score > 0.9) %>%
	dplyr::select(cell_id) %>% pull()
length(stem_cell_id)

terminal_cell_id <- te013_table %>% 
	# dplyr::filter(cluster %in% c("C4","C5","C6") ) %>%
	# dplyr::filter(cluster_id %in% c("C4","C7") ) %>%
	# dplyr::filter(cytotrace_score < 0.2) %>%
	dplyr::filter(cluster_id %in% c("C1","C2","C3","C8") ) %>%
	dplyr::filter(cytotrace_score < 0.1) %>%
	# dplyr::filter(cluster %in% c("C1","C3") ) %>%
	# dplyr::filter(cluster %in% c("C1") ) %>%
	# dplyr::filter(cluster %in% c("C3") ) %>%
	# dplyr::filter(cytotrace_score < 0.5) %>%
	dplyr::select(cell_id) %>% pull()
length(terminal_cell_id)

gene_seurat <- SetIdent(gene_seurat, value = "unk")
gene_seurat <- SetIdent(gene_seurat, 
												cells = stem_cell_id, 
												value = "stem")
gene_seurat <- SetIdent(gene_seurat, 
												cells = terminal_cell_id, 
												value = "terminal")
table(Idents(gene_seurat))

gene_seurat <- SCTransform(gene_seurat, vst.flavor = "v2",seed.use = 666,
													 return.only.var.genes = FALSE , 
													 variable.features.n = nrow(gene_seurat@assays$RNA@data) )
gene_seurat <- PrepSCTFindMarkers(gene_seurat)

dge_table <- FindMarkers(gene_seurat, assay = "SCT", 
												 test.use = "wilcox",
												 ident.1 = "stem", 
												 ident.2 = "terminal",
												 verbose = FALSE, 
												 features = rownames(gene_seurat@assays$SCT@scale.data) , logfc.threshold = 0.01 ,
												 recorrect_umi = FALSE)
dge_table <- dge_table %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble()
dge_table$diff <- dge_table$pct.1 - dge_table$pct.2
sum( dge_table$p_val_adj < 0.1 )
nrow(dge_table)
# View(dge_table)

ges <- qnorm(dge_table$p_val/2 , lower.tail = F ) * sign(dge_table$avg_log2FC)
names(ges) <- dge_table$gene

## VIPER Analysis ----
isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
TE001_network_stem <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
x <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat") , col_names = FALSE )
y <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat") , col_names = FALSE )
z <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat") , col_names = FALSE )
if( TRUE )
{
	selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))	
} else {
	selected_regulators <- unique(c( x$X1 , y$X1 , z$X1 ))
}
print(length(selected_regulators))

require(viper)
TE001_network_stem <- TE001_network_stem[ names(TE001_network_stem) %in% selected_regulators]
class(TE001_network_stem) <- "regulon"
TE001_network_stem <- pruneRegulon(TE001_network_stem,50)

x <- as.matrix(cbind(ges,ges))
rownames(x) <- names(ges)
colnames(x) <- c("stem_signature","shit")
vp <- viper(x,regulon = TE001_network_stem,method="none")
dim(vp)

pwe <- do_pathways_with_aREA(vp %>% mouse_to_human())

dim(pwe)

df <- data.frame( protein = rownames(vp) , nes = vp[,1] )
df <- df[ order(df$nes,decreasing = T) , ]
writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-from-C5.xlsx"))
# writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-from-C2.xlsx"))

df <- data.frame( pathway = rownames(pwe) , nes = pwe[,1] )
df <- df[ order(df$nes,decreasing = T) , ]
writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-from-C5-pathways.xlsx"))
# writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-from-C2-pathways.xlsx"))




## Stem Signature of Stem clusters vs terminally differentiated ----
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-original-data.rds"
fs::file_info(filename)
gene_seurat <- readRDS(filename)

filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-iter-clustering-viper-analysis-with-metacell.tsv"
cluster_id.tibble <- read_csv(filename)
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-cytotrace.csv"
cytotrace.tibble <- read_delim(filename,"\t")

my_tibble <- left_join(cluster_id.tibble,cytotrace.tibble,by=c("sample_id"="cell_id"))

my_tibble <- my_tibble %>% dplyr::rename(cell_id=sample_id)

nrow(my_tibble)
my_tibble <- my_tibble[ complete.cases(my_tibble) , ]
nrow(my_tibble)

gene_seurat <- SCTransform(gene_seurat, vst.flavor = "v2",seed.use = 666,
													 return.only.var.genes = FALSE , 
													 variable.features.n = nrow(gene_seurat@assays$RNA@data) )
gene_seurat <- PrepSCTFindMarkers(gene_seurat)

stem_cell_id <- my_tibble %>% 
	dplyr::filter(cluster_id %in% c("stem-1","stem-2")) %>%
	# dplyr::filter(cluster_id == "C2") %>%
	dplyr::filter(cytotrace_score > 0.9975) %>%
	dplyr::select(cell_id) %>% pull()
length(stem_cell_id)

terminal_cell_id <- my_tibble %>%
	# dplyr::filter( cluster_id %in% c("stem-1","stem-2") ) %>%
	dplyr::filter( !( cluster_id %in% c("stem-1","stem-2")) ) %>%
	# dplyr::filter( !( cluster_id %in% c("secretory-1","secretory-2","secretory-3")) ) %>%
	# dplyr::filter( !( cluster_id %in% c("absorbitive-1","absorbitive-2")) ) %>%
	# dplyr::filter(cytotrace_score < 0.4) %>%
	dplyr::select(cell_id) %>% pull()
length(terminal_cell_id)

gene_seurat <- SetIdent(gene_seurat, value = "unk")
gene_seurat <- SetIdent(gene_seurat, 
												cells = stem_cell_id, 
												value = "stem")
gene_seurat <- SetIdent(gene_seurat, 
												cells = terminal_cell_id, 
												value = "terminal")
table(Idents(gene_seurat))

dge_table <- FindMarkers(gene_seurat, assay = "SCT", 
												 test.use = "wilcox",
												 ident.1 = "stem", 
												 ident.2 = "terminal",
												 verbose = FALSE, 
												 features = rownames(gene_seurat@assays$SCT@scale.data) , logfc.threshold = 0.01 ,
												 recorrect_umi = FALSE)
dge_table <- dge_table %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble()
dge_table$diff <- dge_table$pct.1 - dge_table$pct.2
sum( dge_table$p_val_adj < 0.1 )
nrow(dge_table)
# View(dge_table)

ges <- qnorm(dge_table$p_val/2 , lower.tail = F ) * sign(dge_table$avg_log2FC)
names(ges) <- dge_table$gene

## VIPER Analysis ----
isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
TE001_network_stem <- readRDS(file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds"))
x <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat") , col_names = FALSE )
y <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat") , col_names = FALSE )
z <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat") , col_names = FALSE )
if( TRUE )
{
	selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))	
} else {
	selected_regulators <- unique(c( x$X1 , y$X1 , z$X1 ))
}
print(length(selected_regulators))

require(viper)
TE001_network_stem <- TE001_network_stem[ names(TE001_network_stem) %in% selected_regulators]
class(TE001_network_stem) <- "regulon"
TE001_network_stem <- pruneRegulon(TE001_network_stem,50)

x <- as.matrix(cbind(ges,ges))
rownames(x) <- names(ges)
colnames(x) <- c("stem_signature","shit")
vp <- viper(x,regulon = TE001_network_stem,method="none")
dim(vp)

pwe <- do_pathways_with_aREA(vp %>% mouse_to_human())

dim(pwe)

df <- data.frame( protein = rownames(vp) , nes = vp[,1] )
df <- df[ order(df$nes,decreasing = T) , ]
writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-stem-vs-diff.xlsx"))
# writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-from-C2.xlsx"))

df <- data.frame( pathway = rownames(pwe) , nes = pwe[,1] )
df <- df[ order(df$nes,decreasing = T) , ]
writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-stem-vs-diff-pathways.xlsx"))
# writexl::write_xlsx( df , file.path(reports.dir,"stem-pas-signature-from-C2-pathways.xlsx"))

