
# system.time({ source("sources/figures/correlation-analysis-cytotrace-mrs-and-genes.R") })

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")

library(ComplexHeatmap)
# lt = lapply(1:20, function(x) cumprod(1 + runif(1000, -x/100, x/100)) - 1)
# ha = rowAnnotation(foo = anno_horizon(lt))	
# ha	
# draw(ha)

require(zoo)
# prots_corr_vector

## Protein Activity Correlation Plot ----
N <- 100
# selected_markers <- sort( doStouffer(mat[,1:25]) , decreasing = TRUE )[1:25]
selected_markers <- names(tail(sort(prots_corr_vector),N))

df_to_plot <- vpmat.tibble %>% filter( protein %in% selected_markers )
# df_to_plot <- df_to_plot %>% dplyr::rename("cluster_id"="pas_cluster_id")
df_to_plot <- df_to_plot %>% arrange(desc(cytotrace))

df_to_plot <- left_join(df_to_plot,clustering_table,by=c("cell_id"="sample_id"))

my_list <- vector(mode = "list",length = length(table(df_to_plot$protein)))
names(my_list) <- names(table(df_to_plot$protein))

index <- match( prots_corr_table$protein[1:N] , names(my_list) )
my_list <- my_list[index]

for ( a_protein in names(my_list) )
{
	x <- df_to_plot %>% dplyr::filter(protein == a_protein) %>% arrange(desc(cytotrace)) %>% pull(viper_score)
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
									 height = unit(10, "in") , 
									 use_raster = TRUE , raster_quality = 10 ,
									 top_annotation = ha_cyto ,
									 # column_title = "Protein Activity"
									 heatmap_legend_param = list(title = "Protein Activity") ,
									 col = circlize::colorRamp2(c(-5, 0, 5), c("deepskyblue3", "white", "brown3")),
									 border = FALSE , row_gap = unit(2, "mm") ,
									 row_names_side = "left" , show_column_names = FALSE ,
									 cluster_rows = FALSE ,
									 cluster_columns = FALSE, 
									 row_names_gp = gpar(fontsize = 6) )

# ha = rowAnnotation(foo = anno_horizon(my_list,normalize = TRUE , 
# 																			gap = unit(1, "mm") ,negative_from_top = FALSE) )	
# pdf( file.path(reports.dir,"panel-correlation.pdf") , height = 2 , width = 5)
# draw(ha)
# dev.off()

# ha = rowAnnotation(foo = anno_density(my_list, type = "heatmap", width = unit(6, "cm")))	
# pdf( file.path(reports.dir,"panel-correlation.pdf") , height = 2 , width = 5)
# draw(ha)
# dev.off()

# BiocManager::install("topGO")
# library(topGO)

# install.packages("enrichR")
library(enrichR)
# listEnrichrSites()
setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()
websiteLive <- TRUE
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# stem_proteins <- names( head( sort(prots_corr_vector,decreasing = TRUE) , N ) ) %>% sort()
stem_proteins <- selected_markers
length(stem_proteins)

dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021","Transcription_Factor_PPIs")
if (websiteLive) {
	enriched <- enrichr(stem_proteins, dbs)
}

# View( if (websiteLive) enriched[["GO_Biological_Process_2021"]] ) 
# View( if (websiteLive) enriched[["GO_Molecular_Function_2021"]] )
# View( if (websiteLive) enriched[["GO_Cellular_Component_2021"]] )
# View( if (websiteLive) enriched[["Transcription_Factor_PPIs"]] )
if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

# dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021","Transcription_Factor_PPIs")
# if (websiteLive) {
# 	enriched <- enrichr(enriched[["Transcription_Factor_PPIs"]]$Term , dbs)
# }

# View(preppi_table %>% filter(symbol_prot1 == "ETV5") )

x <- enriched[["GO_Biological_Process_2021"]]
# x <- enriched[["GO_Molecular_Function_2021"]]
x <- x %>% arrange(Adjusted.P.value)
# x <- x[1:50,]
x <- x %>% dplyr::select(Term,protein=Genes)
x <- x %>% separate_rows(protein, convert = TRUE)
x$presence <- 1
x <- x %>% pivot_wider(names_from = protein,values_from = presence,values_fill = 0)
mat_go <- tibble2matrix(x) %>% t()
length(rownames(mat_go))
# mat_go <- mat_go[,1:round(N/4+5)]
mat_go <- mat_go[,1:100]

mat_go <- mat_go[ , names( sort( colSums(mat_go), decreasing = T ) ) ]

# tmp <- rownames(mat_go) %>% human_to_mouse(na.rm = F)
# na_index <- is.na(tmp)
# rownames(mat_go)[na_index]

# print_msg_warn("FACKING STUPID FIX converting ZRSR1 to ZRSR2")
# rownames(mat_go) <- ifelse(rownames(mat_go) == "Zrsr1","Zrsr2",rownames(mat_go))
# missing_prots <- stem_proteins[ !(stem_proteins %in% rownames(mat_go)) ]

# m <- matrix(0,ncol = ncol(mat_go) , nrow = length(missing_prots) , dimnames = list(missing_prots,colnames(mat_go)))

# mat_go <- rbind(mat_go,m)

print_msg_warn("**** VERY BAD HUMAN/MOUSE GENE CONVERSION: WE MISS ATAD2 in the canonical conversion approach ****")
library(stringr)
rownames(mat_go) <- str_to_title( tolower( rownames(mat_go) ))

index <- match( rownames(mat) , rownames(mat_go) )
mat_go <- mat_go[index,]

df_annot_go_count <- columnAnnotation( n_go = anno_barplot(
	colSums(mat_go,na.rm = T) ,
	bar_width = 1,
	gp = gpar(col = "white", fill = "#FFE200"),
	border = FALSE, 
	baseline = 0 ,
	axis_param = list(direction = "reverse") ,
	# axis_param = list(at = c(0, 5e5, 1e6, 1.5e6),
	#     labels = c("0", "500k", "1m", "1.5m")),
	# width = unit(2, "cm") 
	height = unit(2, "in")
	)
)    

h_go <- Heatmap( mat_go , width = unit(10, "in") ,
								 height = unit(10, "in") , 
								 row_gap = unit(5, "mm") , show_heatmap_legend = FALSE ,
									 col = circlize::colorRamp2(c(1,0), c("gray20", "gray98")), na_col = "gray98",
									 border = FALSE , rect_gp = gpar(col = "white", lwd = 2) ,
									 cluster_rows = FALSE ,
								 cluster_columns = FALSE , 
								 column_names_side = "top" , show_row_names = FALSE , 
								 bottom_annotation = df_annot_go_count ,
								 row_names_gp = gpar(fontsize = 8) , column_names_gp = gpar(fontsize = 5) , column_names_rot = 60 )


ht_list <- h_corr + h_go

# pdf( file.path(reports.dir,"correlation-plot-figure-one.pdf") , width = 6 , height = 4.75 )
pdf( file.path(reports.dir,"correlation-plot-figure-one.pdf") , width = 20 , height = 16 )
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


chromatin_index <- grepl("histon|chromatin|chromosome organization" , colnames(mat_go) , ignore.case = T )
sum(rowSums(mat_go[,chromatin_index]) > 0 ,na.rm = TRUE)
dna_damage_index <- grepl("DNA damage response|chromatin" , colnames(mat_go) , ignore.case = T )
sum(rowSums(mat_go[,dna_damage_index]) > 0,na.rm = TRUE)
cell_cycle_index <- grepl("DNA replication|DNA repair" , colnames(mat_go) , ignore.case = T )
sum(rowSums(mat_go[,cell_cycle_index]) > 0,na.rm = TRUE)



