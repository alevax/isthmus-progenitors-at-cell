## 
# Plotting Analysis for Cell Paper
# --------------------------------
# system.time({source("sources/isc-multiome/multiome-plots-for-paper.R")})

source("../vaxtools/R/utils.R")
my_sample_id <- "TE013"
create_workspace( isWorkspaceToClean = FALSE , 
									experiments_dir = "experiments" , 
									# run_dir = "multiome-plots-for-paper-with-ingest" )
									run_dir = "multiome-plots-for-paper" )

library(logger)
# Using logger library
log_threshold(DEBUG)
log_threshold(DEBUG,index=2)

log_filename <- file.path(reports.dir,"log.txt")
log_appender(appender = appender_file(log_filename),index=2)
# log_formatter(formatter_glue)
log_layout(layout_glue_colors)
log_warn('--- Start of Analysis ----')
# demo(colors, package = 'logger', echo = FALSE)	

log_info("Working on Sample: " , my_sample_id)

library(tidyverse)
library(hdf5r)
library(Signac)
library(Seurat)

library(future)
options(future.globals.maxSize = 100 * 1024 ^ 3) # for 50 Gb RAM
plan()
library(parallel)
plan("multicore", workers = 8)

log_info(">>> Loading Signac file")
signac_filename <- "~/Clouds/Dropbox/Data/isc/TE013/TE013-signac-subset.rds"
signac_obj <- readRDS(signac_filename)
# signac_obj@meta.data

log_info(">>> Loading CytoTRACE Data")
cyto_te013 <- read_delim("~/Clouds/Dropbox/Data/isc/TE013/TE013-cytotrace.csv",delim = ",")

peaks_counts.mat <- as.matrix(signac_obj@assays$peaks@counts)

my_peaks <- colSums(peaks_counts.mat > 0)
peaks_tbl <- tibble( cell_id = names(my_peaks) ,
										 open_peaks = my_peaks )

my_tibble <- right_join( cyto_te013 , peaks_tbl )
my_tibble

# Calculate correlation and p-value
correlation_test <- cor.test(my_tibble$cytotrace_score, my_tibble$open_peaks)
cor_value <- correlation_test$estimate
p_value <- correlation_test$p.value

# Create the plot
p <- ggplot(my_tibble, aes(x=cytotrace_score, y=open_peaks)) +
    geom_point(size=0.005, alpha=0.5) +
    coord_trans(y="log10") +
    geom_smooth(method = "glm", formula = y~x, color="red",
                method.args = list(family = gaussian(link = 'log'))) +
    xlim(c(0,1)) +  # Adjust this line to set the x-axis limits as needed
    xlab("CytoTRACE Score") +
    ylab("log10(total peaks)") +
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = sprintf("Correlation: %.2f, p-value: %.2e", cor_value, p_value))

# Save the plot to a PDF
pdf(file.path(reports.dir,"log-corr-number-of-peaks-vs-cytotrace.pdf"), width = 5, height = 3)
print(p)
dev.off()

# cor.test(my_tibble$cytotrace_score,my_tibble$open_peaks,method = "pea")

# peaks_counts.mat <- as.matrix(signac_obj@assays$peaks@counts)
peaks_annotation <- ClosestFeature(signac_obj,regions = rownames(peaks_counts.mat))
dim(peaks_counts.mat)
dim(peaks_annotation)

regions_in_common <- intersect( peaks_annotation$query_region , rownames(peaks_counts.mat) )
peaks_annotation <- peaks_annotation[ match( regions_in_common , peaks_annotation$query_region ) , ]
peaks_counts.mat <- peaks_counts.mat[ match( regions_in_common , rownames(peaks_counts.mat) ) , ]
stopifnot(identical( peaks_annotation$query_region , rownames(peaks_counts.mat) ))

# We thank the reviewer for the suggestion. We now show lineage-associated gene accessibility scores and compare 
# their variability to Lrg4/5 expression to hopefully show that low to high Lgr5 expression correlate 
# with low-to-high lineage-specific accessibility scores, while our new markers do the opposite.


## Open VIPER clustering object ----
filename <- "~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-viper-analysis.rds"
vp_seurat <- readRDS(filename)

# ## With Ingest
# filename <- "~/Clouds/Dropbox/Data/isc/multiome-metadata-ingest.csv"
# md_ingest <- read_csv(filename)
# 
# cells_in_common <- intersect( colnames(vp_seurat) , md_ingest$cell_id )
# index <- match( colnames(vp_seurat) , md_ingest$cell_id )
# md_ingest <- md_ingest[index,]
# table(md_ingest$cluster_id)
# 
# # setdiff( colnames(vp_seurat) , md_ingest$cell_id )
# # setdiff( md_ingest$cell_id , colnames(vp_seurat))
# vp_seurat$cluster_id_aREA_assigned_main_lineages <- md_ingest$cluster_id

## with OncoMatch-like enrichment
table(vp_seurat$cluster_id_aREA_assigned)
vp_seurat$cluster_id_aREA_assigned_main_lineages <- "nothing"
vp_seurat$cluster_id_aREA_assigned_main_lineages <- ifelse( grepl( "secret|pane" , vp_seurat$cluster_id_aREA_assigned ) , "secretory" , vp_seurat$cluster_id_aREA_assigned_main_lineages )
vp_seurat$cluster_id_aREA_assigned_main_lineages <- ifelse( grepl( "abso" , vp_seurat$cluster_id_aREA_assigned ) , "absorptive" , vp_seurat$cluster_id_aREA_assigned_main_lineages )
vp_seurat$cluster_id_aREA_assigned_main_lineages <- ifelse( grepl( "stem" , vp_seurat$cluster_id_aREA_assigned ) , "stem" , vp_seurat$cluster_id_aREA_assigned_main_lineages )
table(vp_seurat$cluster_id_aREA_assigned_main_lineages)

vp_tibble <- tibble( cell_id = colnames(vp_seurat) ,
										 cluster_id = vp_seurat$cluster_id_aREA_assigned_main_lineages )

## Gene names involved in differentiations (lLineage Correlation Plot) ----
secretory_gene_list <- c("Atoh1","Dll1","Neurod1","Dclk1","Lyz1")
absorptive_gene_list <- c("Alpi","Fabp2","Krt20")
progenitors_gene_list <- c("Atad2","Stmn1","Smarca5","Dek")

lineages_to_plot <- list(secretory=secretory_gene_list,absorptive=absorptive_gene_list,progenitor=progenitors_gene_list)

for (i in 1:length(lineages_to_plot)) {
	
	gene_list <- lineages_to_plot[[i]]
	lineage_name <- names(lineages_to_plot[i])
	
	index <- apply(sapply(peaks_annotation$gene_name, 
												function(word) grepl(word,gene_list)), 2, any)
	
	# peaks_annotation[index,]
	sum(index)
	# index <- index[ peaks_annotation$distance == 0 ]
	# sum(index)
	my_peaks <- colSums(peaks_counts.mat[index,] >0)
	peaks_tbl <- tibble( cell_id = names(my_peaks) ,
											 open_peaks = my_peaks )
	
	my_tibble <- right_join( cyto_te013 , peaks_tbl )
	my_tibble
	
	my_tibble <- right_join( vp_tibble , my_tibble )
	
	# for (a_cluster in names(table(vp_seurat$cluster_id_aREA_assigned))) {
	for (a_cluster in names(table(vp_seurat$cluster_id_aREA_assigned_main_lineages))) {
		tibble2plot <- my_tibble %>% dplyr::filter( cluster_id %in% a_cluster )
		print(nrow(tibble2plot))
		
		# # To prepare for log transformation in plotting
		# tibble2plot$open_peaks <- ifelse( tibble2plot$open_peaks == 0 , 0.999 , tibble2plot$open_peaks )
		
		# Removing 0s because we don't know whether it's dropout or actual signal
		tibble2plot <- tibble2plot[ tibble2plot$open_peaks > 0 , ]
		
		# Calculate correlation and p-value
		cor_test <- cor.test(tibble2plot$cytotrace_score, tibble2plot$open_peaks , 
												 method = "pea")
		cor_value <- cor_test$estimate
		p_value <- cor_test$p.value
		
		p <- ggplot(tibble2plot, aes(x=cytotrace_score,y=open_peaks)) +
			geom_point(size=0.005,alpha=0.5) +
			# scale_y_continuous(trans='log10') +
			# coord_trans(y="log10") +
			geom_smooth(method = "glm", formula = y~x, color="red" ) +
			annotate("text", x = Inf, y = Inf, label = sprintf("R: %.2f, p: %.5f", cor_value, p_value), 
							 hjust = -2.1, vjust = 1.1)  +
			xlim(c(1,0)) +
			xlab("CytoTRACE Score") +
			ylab("ATAC Fragments") +
			theme_minimal()
		
		pdf( file.path(reports.dir,paste0(lineage_name,"-lineage-genes-total-peaks-vs-cytotrace-on-cluster-",a_cluster,".pdf") ) , width = 5 , height = 3)
		print(p)
		dev.off()
	}
	
}

## Correlation between peaks and expression ----
# gene_list_for_expression <- c("Lgr4","Lgr5","Mki67",secretory_gene_list,absorptive_gene_list,progenitors_gene_list,"Fgfbp1")
gene_list_for_expression <- c("Lgr4","Lgr5","Mki67",secretory_gene_list,absorptive_gene_list,progenitors_gene_list)
seurat_obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE013/TE013-seurat-analysis-data.rds")
seurat_obj
seurat_obj <- seurat_obj %>% subset( cells = colnames(peaks_counts.mat) )

gex.cpm <- Seurat::RelativeCounts(seurat_obj@assays$RNA@counts,scale.factor = 1e6)
gex.cpm.log10 <- log10(gex.cpm+1)

gex_tibble <- gex.cpm.log10[gene_list_for_expression,] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(cols = !c("gene") , names_to = "cell_id" , values_to = "expression" )
gex_tibble <- gex_tibble %>% pivot_wider(names_from = "gene",values_from = "expression") 

tmp <- tibble( cell_id = colnames(seurat_obj) , stemness = seurat_obj@meta.data$stemness_index.ges )
gex_tibble <- left_join(gex_tibble,tmp)

## ATAC Tibble remodeling
index <- apply(sapply(peaks_annotation$gene_name, 
											function(word) grepl(word,secretory_gene_list)), 2, any)
# index <- index[ peaks_annotation$distance == 0 ]
my_peaks_secretory <- colSums(peaks_counts.mat[index,] > 0)

index <- apply(sapply(peaks_annotation$gene_name, 
											function(word) grepl(word,absorptive_gene_list)), 2, any)
# index <- index[ peaks_annotation$distance == 0 ]
my_peaks_absorptive <- colSums(peaks_counts.mat[index,] > 0)

index <- apply(sapply(peaks_annotation$gene_name, 
											function(word) grepl(word,progenitors_gene_list)), 2, any)
# index <- index[ peaks_annotation$distance == 0 ]
my_peaks_progenitor <- colSums(peaks_counts.mat[index,] > 0)

my_peaks_total <- colSums(peaks_counts.mat > 0)

peaks_tbl <- tibble( cell_id = names(my_peaks) ,
										 peaks_secretory = my_peaks_secretory ,
										 peaks_absorptive = my_peaks_absorptive ,
										 peaks_progenitor = my_peaks_progenitor ,
										 peaks_total = my_peaks_total )
my_tibble <- right_join( cyto_te013 , peaks_tbl )
my_tibble <- right_join( vp_tibble , my_tibble )

multimodal_tibble <- left_join( gex_tibble , my_tibble )

library(corrplot)
library(gridGraphics)

for (a_cluster in names(table(vp_seurat$cluster_id_aREA_assigned_main_lineages))) {
	
	my_mat <- multimodal_tibble %>% dplyr::filter(cluster_id == a_cluster) %>% dplyr::select(-cell_id,-cluster_id) %>% as.matrix()
	M <- cor(my_mat,method = "pea")
	
	# mat : is a matrix of data
	# ... : further arguments to pass to the native R cor.test function
	cor.mtest <- function(mat, ...) {
		mat <- as.matrix(mat)
		n <- ncol(mat)
		p.mat<- matrix(NA, n, n)
		diag(p.mat) <- 0
		for (i in 1:(n - 1)) {
			for (j in (i + 1):n) {
				tmp <- cor.test(mat[, i], mat[, j], ...)
				p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
			}
		}
		colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
		p.mat
	}
	# matrix of the p-value of the correlation
	p.mat <- cor.mtest(my_mat,method="pea")
	head(p.mat[, 1:5])
	
	col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA") %>% rev() )
	corrplot(M, method="color", col=col(200),  
					 type="upper", 
					 # order="hclust",	
					 addCoef.col = "black", # Add coefficient of correlation ,
					 tl.col="black", tl.srt=45, #Text label color and rotation
					 # Combine with significance
					 p.mat = p.mat, sig.level = 0.1, insig = "blank", number.cex = 0.6 ,
					 order = 'hclust' ,
					 # hide correlation coefficient on the principal diagonal
					 diag=TRUE 
	)
	
	grid.echo()
	p <- grid.grab()
	
	pdf( file.path(reports.dir,paste0("corr-mat-atac+rna-",a_cluster,".pdf")) , height = 10 , width = 10 )
	grid.draw(p)
	dev.off()
	
}

## Correlation Plot across all the dataset ----
my_mat <- multimodal_tibble %>% dplyr::select(-cell_id,-cluster_id) %>% as.matrix()
M <- cor(my_mat,method = "pea")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(my_mat,method = "pea")
head(p.mat[, 1:5])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA") %>% rev() )
corrplot(M, method="color", col=col(200),  
				 type="upper", 
				 # order="hclust",	
				 addCoef.col = "black", # Add coefficient of correlation ,
				 tl.col="black", tl.srt=45, #Text label color and rotation
				 # Combine with significance
				 p.mat = p.mat, sig.level = 0.01, insig = "blank", number.cex = 0.6 ,
				 order = 'FPC' ,
				 # hide correlation coefficient on the principal diagonal
				 diag=TRUE 
)
library(gridGraphics)
grid.echo()
p <- grid.grab()

pdf( file.path(reports.dir,paste0("corr-mat-atac+rna-all-data.pdf")) , height = 10 , width = 10 )
grid.draw(p)
dev.off()







# for (i in 1:length(lineages_to_plot)) {
# 	
# 	gene_list <- lineages_to_plot[[i]]
# 	lineage_name <- names(lineages_to_plot[i])
# 	
# 	# ATAC Tibble remodeling
# 	index <- apply(sapply(peaks_annotation$gene_name, 
# 												function(word) grepl(word,gene_list)), 2, any)
# 	
# 	my_peaks <- colSums(peaks_counts.mat[index,])
# 	peaks_tbl <- tibble( cell_id = names(my_peaks) ,
# 											 open_peaks = my_peaks )
# 	my_tibble <- right_join( cyto_te013 , peaks_tbl )
# 	my_tibble <- right_join( vp_tibble , my_tibble )
# 	
# 	multimodal_tibble <- left_join( gex_tibble , my_tibble )
# 	
# 	for (a_cluster in names(table(multimodal_tibble$cluster_id))) {
# 		tibble2plot <- multimodal_tibble %>% dplyr::filter( cluster_id %in% a_cluster )
# 		print(nrow(tibble2plot))
# 		
# 		# # Removing 0s because we don't know whether it's dropout or actual signal
# 		# tibble2plot <- tibble2plot[ tibble2plot$open_peaks > 0 , ]
# 		
# 		# Calculate correlation and p-value
# 		cor_test <- cor.test(tibble2plot$cytotrace_score, tibble2plot$open_peaks , method = "spe")
# 		cor_value <- cor_test$estimate
# 		p_value <- cor_test$p.value
# 		p <- ggplot(tibble2plot, aes(x=cytotrace_score,y=open_peaks)) +
# 			geom_point(size=0.005,alpha=0.5) +
# 			# scale_y_continuous(trans='log10') +
# 			# coord_trans(y="log10") +
# 			geom_smooth(method = "glm", formula = y~x, color="red" ) +
# 			annotate("text", x = Inf, y = Inf, label = sprintf("R: %.2f, p: %.5f", cor_value, p_value), 
# 							 hjust = -2.1, vjust = 1.1)  +
# 			xlim(c(1,0)) +
# 			xlab("CytoTRACE Score") +
# 			ylab("ATAC Fragments") +
# 			theme_minimal()
# 		
# 		pdf( file.path(reports.dir,paste0(lineage_name,"-lineage-genes-total-peaks-vs-cytotrace-on-cluster-",a_cluster,".pdf") ) , width = 5 , height = 3)
# 		print(p)
# 		dev.off()
# 	}
# }
