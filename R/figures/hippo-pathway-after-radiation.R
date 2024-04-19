
# system.time({ source("sources/figures/hippo-pathway-after-radiation.R") })

source("../vaxtools/R/utils.R")
source("../vaxtools/R/cross-species-utils.R")
source("../vaxtools/R/interactome_handler.R")

create_workspace("hippo")

require(tidyverse)
library(Seurat)
seurat_obj.expr <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE002-gex.rds")
seurat_obj.expr
table(seurat_obj.expr@meta.data$sample_id)

seurat_obj.expr_control <- seurat_obj.expr %>% subset( sample_id == "control" )
seurat_obj.expr_treatment <- seurat_obj.expr %>% subset( sample_id == "treatment" )


seurat_obj.expr_control <- seurat_obj.expr_control %>% subset( cytotrace_score > 0.9 )
seurat_obj.expr_control
seurat_obj.expr_treatment <- seurat_obj.expr_treatment %>% subset( cytotrace_score > 0.9 )
seurat_obj.expr_treatment

	selected_cells <- colnames(seurat_obj.expr_control)
	Idents(seurat_obj.expr_control) <- "unk"
	seurat_obj.expr_control <- SetIdent( seurat_obj.expr_control , cells = selected_cells , "stem" )
	table(Idents(seurat_obj.expr_control))
	
	selected_cells <- colnames(seurat_obj.expr_treatment)
	Idents(seurat_obj.expr_treatment) <- "unk"
	seurat_obj.expr_treatment <- SetIdent( seurat_obj.expr_treatment , cells = selected_cells , "stem-response" )
	table(Idents(seurat_obj.expr_treatment))
	
	selected_cells_t <- colnames(seurat_obj.expr_treatment)
	selected_cells_c <- colnames(seurat_obj.expr_control)
	Idents(seurat_obj.expr) <- "unk"
	seurat_obj.expr <- SetIdent( seurat_obj.expr , cells = selected_cells_c , "stem" )
	seurat_obj.expr <- SetIdent( seurat_obj.expr , cells = selected_cells_t , "stem-response" )
	table(Idents(seurat_obj.expr))
	
	seurat_obj.expr <- SCTransform(seurat_obj.expr, vst.flavor = "v2",seed.use = 666,
																 return.only.var.genes = FALSE , 
																 variable.features.n = nrow(seurat_obj.expr@assays$RNA@data) )	
	seurat_obj.expr <- PrepSCTFindMarkers(seurat_obj.expr)
	markers_TvsC <- FindMarkers(seurat_obj.expr,ident.1 = "stem-response",ident.2 = "stem" ,
															features = rownames(seurat_obj.expr@assays$SCT@data) ,logfc.threshold = 0 )
	markers_TvsC <- as_tibble(markers_TvsC %>% as.data.frame() %>% rownames_to_column("gene") )
	markers_TvsC
	
	correctMarkersTable <- function(table) {
		
		min_p <- min( table$p_val[ table$p_val != 0 ] )
		table$p_val_corr <- ifelse( table$p_val == 0 , min_p , table$p_val )
		
		min_p_adj <- min( table$p_val_adj[ table$p_val_adj != 0 ] )
		table$p_val_adj_corr <- ifelse( table$p_val_adj == 0 , min_p_adj , table$p_val_adj )
		
		table$p_val_adj_mlog10 <- -log10(table$p_val_adj_corr)
		
		return(table)
	}
	
	markers_TvsC <- correctMarkersTable(markers_TvsC)
	markers_TvsC$is_dge <- ifelse(markers_TvsC$p_val_adj < 0.1,"Yes","No")
	fc_threshold <- 1
	
	markers_TvsC$is_a_validated_gene <- ifelse(markers_TvsC$gene %in% c("Clu","Areg","Ly6a"),"Yes","No")
	
	plotVolcano <- function(my_df) {
		require(ggplot2)
		require(ggrepel)
		my_plot <- ggplot( data = my_df , 
											 aes( x = avg_log2FC , y = p_val_adj_mlog10 , color = "grey50" ) ) +
			
			geom_point( alpha = 0.5 , stroke = 1 , size = 1 , color = "gray" , shape = 20 ) +
					geom_point( data = my_df %>% filter( is_dge == "No" ) , alpha = 0.75 , stroke = 1 , size = 1 , color = "gray" , fill = "white" , shape = 16 ) +
					geom_point( data = my_df %>% filter( is_dge == "Yes" & avg_log2FC > fc_threshold ) , alpha = 0.75 , stroke = 1 , size = 1 , color = "brown3", fill ="brown3" , shape = 16 ) +
					geom_point( data = my_df %>% filter( is_dge == "Yes" & avg_log2FC < -fc_threshold ) , alpha = 0.75 , stroke = 1 , size = 1 , color = "deepskyblue3", fill ="deepskyblue3" , shape = 16 ) +
					geom_point( data = my_df %>% filter( is_a_validated_gene == "Yes") , 
											color = "darkorange" , alpha = 0.75 , stroke = 2 , size = 3 , fill = "white" , shape = 1 )  +
			# geom_text_repel( data = my_df %>% filter( p_val_adj_mlog10 > 20 || is_a_validated_gene == "Yes" ) ,
			geom_text_repel( data = my_df %>% filter( is_a_validated_gene == "Yes" ) ,			
											 aes( label = gene ) ,
											 max.iter = 3e3 , max.overlaps = 20 ,
											 box.padding = unit(0.5, 'lines') , point.padding = unit(0.5, 'lines') ,
											 size = 3 , color = "gray25" , segment.color = 'gray75',
											 # nudge_x = ifelse( my_df %>%
											 # 										# filter( is_me2Dn_me3Up == "Yes" ) %>%
											 # 										dplyr::select(avg_log2FC) %>% pull() > 0, 1.5, -1.5)
			) +
			xlim(c( -1*max(abs(my_df$avg_log2FC))-2 , + max(abs(my_df$avg_log2FC))+2 )) +
			ylim(c( 0 , + max(my_df$p_val_adj_mlog10)+2 )) +
			xlab( paste0( "LogFC") ) +
			ylab( paste0( "-log10(Adj P Value)") ) +
			# ggtitle( paste0( file_tag ) ) +
			theme_minimal() +
			theme( 
				plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
				axis.title.x = element_text(size=15,face="bold") , # family="Open Sans") ,
				axis.title.y = element_text(size=15,face="bold") ) # family="Open Sans") )
		
		return(my_plot)
	}
	
		my_plot <- plotVolcano(markers_TvsC)
		pdf( file.path(reports.dir,paste0("TvsC-dp-treatment-vs-dp-control-volcano-plot.pdf")) )
		print(my_plot)
		dev.off()
		
		markers_TvsC$pct_diff <- markers_TvsC$pct.1 - markers_TvsC$pct.2
		markers_TvsC <- markers_TvsC %>% arrange(desc(pct_diff))
		
		filename <- file.path(reports.dir,paste0("gene-expression-irradiation-response-markers.xlsx"))
		# markers_TvsC$sample_id <- file_tag
		writexl::write_xlsx(markers_TvsC,filename)
		
		
		# log_info("Running Pathway Enrichment Analysis")
		source("../vaxtools/R/pathway-analysis.R")
		
		my_signature <- qnorm( markers_TvsC$p_val/2, lower.tail = FALSE ) * sign( markers_TvsC$avg_log2FC )
		names(my_signature) <- markers_TvsC$gene
		
		x <- as.matrix(my_signature %>% mouse_to_human())
		x <- cbind(x,x)
		colnames(x) <- c("signature","two")
		my_pathways <- do_pathways_with_aREA( x , threshold = 0 )
		my_pathways["REACTOME_SIGNALING_BY_HIPPO",]
		
		index <- !grepl("HALLMARK" , rownames(my_pathways) )
		my_pathways <- my_pathways[index,]
		
		message(">>> Histogram Plot for Pathway Analysis")
		{
			require(scales)
			colnames(my_pathways)
			tmp <- my_pathways[,1]
			my_tibble <- tibble( pathway = names(tmp) , score = tmp )
			a <- my_tibble %>% slice_max(score,n=20) %>% arrange(desc(score))
			b <- my_tibble %>% slice_min(score,n=20) %>% arrange(desc(score))
			my_tibble <- bind_rows(a,b)
			
			my_tibble$pathway <- strtrim(my_tibble$pathway,40)
			
			my_tibble$pathway_factor <- factor(my_tibble$pathway,levels = rev(my_tibble$pathway) , ordered = TRUE )
			
			p <- 	ggplot( data = my_tibble , aes(x = pathway_factor , y = score , fill = score ) ) +
				geom_bar( stat = "identity" , position = position_stack(reverse = TRUE) ) +
				scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0 ) +
				coord_flip() +
				# geom_hline( yintercept = -1.96 , color = "red", size = 1 ) +
				# geom_hline( yintercept = +1.96 , color = "red", size = 1 ) +
				theme_light() +
				theme( legend.position = "right" , 
							 axis.title.x = element_text(face="bold", colour="#990000" , size=10) ,
							 axis.text.y  = element_text(angle=0, size = 5 ,hjust = 1)
				) + 
				ylab( "Normalized Enrichment Score") + 
				xlab( "Pathways") +
				ggtitle( "Pathways" )
			
			pdf( file.path( reports.dir , paste0("pathway-analysis-at-ges.pdf" ) ) , width = 7 , height = 4.5 )
				print(p)
			dev.off()			
		}		
		
	networks.list <- list()
	networks.list[["c1"]] <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds")
	networks.list[["c2"]] <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c2_unPruned.rds")
	networks.list[["c3and4"]] <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c3and4_unPruned.rds")
	str(networks.list,1)
	
		retainOnlyTFs <- function(network) {
		
		print_msg_info(">>> Filtering using my lists")
		x <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat" , col_names = FALSE )
		y <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat" , col_names = FALSE )
		# z <- read_csv( "~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat" , col_names = FALSE )
		
		selected_regulators <- unique(c( x$X1 , y$X1 , c("Rspo1","Olfm4","Krt19","Lgr5","Tnfrsf19","Chga")))	
		
		index <- names(network) %in% selected_regulators
		return( network[index] )
		
	}

	networks.list <- lapply( networks.list , retainOnlyTFs )
	networks.list <- lapply( networks.list , pruneRegulon , 50 )
		
	x <- as.matrix(my_signature)
	x <- cbind(x,x)
	colnames(x) <- c("signature","two")	
	vp <- viper( x , networks.list[["c1"]] , method = "none" )[,1]
	x <- as.matrix(vp %>% mouse_to_human())
	x <- cbind(x,x)
	colnames(x) <- c("signature","two")
	my_pathways <- do_pathways_with_aREA( x , threshold = 0 )
	# my_pathways["REACTOME_SIGNALING_BY_HIPPO",]	
	
	index <- !grepl("HALLMARK" , rownames(my_pathways) )
	my_pathways <- my_pathways[index,]
	
	x <- x[,1,drop=FALSE] %>% as.data.frame() %>% rownames_to_column("protein") %>% as_tibble()
	
	filename <- file.path(reports.dir,paste0("mrs-rad-vs-wt-top-cyto-rank-list.xlsx"))
	# markers_TvsC$sample_id <- file_tag
	writexl::write_xlsx(x,filename)
		
	message(">>> Histogram Plot for Pathway Analysis")
	{
		require(scales)
		colnames(my_pathways)
		tmp <- my_pathways[,1]
		my_tibble <- tibble( pathway = names(tmp) , score = tmp )
		a <- my_tibble %>% slice_max(score,n=20) %>% arrange(desc(score))
		b <- my_tibble %>% slice_min(score,n=20) %>% arrange(desc(score))
		my_tibble <- bind_rows(a,b)
		
		my_tibble$pathway <- strtrim(my_tibble$pathway,40)
		
		my_tibble$pathway_factor <- factor(my_tibble$pathway,levels = rev(my_tibble$pathway) , ordered = TRUE )
		
		p <- 	ggplot( data = my_tibble , aes(x = pathway_factor , y = score , fill = score ) ) +
			geom_bar( stat = "identity" , position = position_stack(reverse = TRUE) ) +
			scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0 ) +
			coord_flip() +
			# geom_hline( yintercept = -1.96 , color = "red", size = 1 ) +
			# geom_hline( yintercept = +1.96 , color = "red", size = 1 ) +
			theme_light() +
			theme( legend.position = "right" , 
						 axis.title.x = element_text(face="bold", colour="#990000" , size=10) ,
						 axis.text.y  = element_text(angle=0, size = 5 ,hjust = 1)
			) + 
			ylab( "Normalized Enrichment Score") + 
			xlab( "Pathways") +
			ggtitle( "Pathways" )
		
		pdf( file.path( reports.dir , paste0("pathway-analysis-at-pas.pdf" ) ) , width = 7 , height = 4.5 )
		print(p)
		dev.off()			
	}		
		
	vp_table <- vp %>% as.data.frame() %>% rownames_to_column("protein") %>% as_tibble()
	vp_table <- vp_table %>% dplyr::rename("score"=".") %>% arrange(desc(score))
	vp_table$protein <- factor( vp_table$protein , levels = vp_table$protein , ordered = TRUE )
	vp_table$protein_order <- 1:nrow(vp_table)
	
	my_threshold <- 2
	
	my_plot <- ggplot( data = vp_table , aes( x = protein_order , y = score ) ) +
		
		geom_point( data = vp_table %>% filter( score < my_threshold ) , alpha = 0.5 , stroke = 1 , size = 1 , color = "gray" , shape = 20 ) +
		geom_point( data = vp_table %>% filter( score >= my_threshold ) , alpha = 0.5 , stroke = 1 , size = 3 , color = "orange" , shape = 20 ) +
		
      geom_text_repel( data = vp_table %>% filter( score >=  my_threshold) , mapping = aes( label = protein ) , 
                       # nudge_x = -3.5, 
                       box.padding = 0.5,
                       # nudge_y = 1,
      								 nudge_x = 300,
                       segment.curvature = -0.1,
                       segment.ncp = 3,
                       segment.angle = 20 
      ) +					
		ylim(-10,max(vp_table$score+1)) +
		xlab( paste0( "Scored Motifs") ) +
		ylab( paste0( "log10 Variability") ) +
		ggtitle( paste0( "Variability Plot") ) +
		theme_minimal() +
		theme( 
			plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
			axis.title.x = element_text(size=15,face="bold") , # family="Open Sans") ,
			axis.title.y = element_text(size=15,face="bold") ) # family="Open Sans") )
	
	# print(my_plot)
	
	# report_dir <- file.path("experiments/2020-01-21-analysis-on-atac-seq-chromVar/reports/")
	pdf(file.path(reports.dir,"variability-plot.pdf") , width = 4 , height = 24 )
		print(my_plot)
	dev.off()
	
	
	
	
	
	markers_TvsUnk <- FindMarkers(seurat_obj.expr,ident.1 = "stem-response",ident.2 = "unk" ,
															features = rownames(seurat_obj.expr@assays$SCT@data) ,logfc.threshold = 0 )
	markers_TvsUnk <- as_tibble(markers_TvsUnk %>% as.data.frame() %>% rownames_to_column("gene") )
	markers_TvsUnk
	markers_TvsUnk.signature <- qnorm( markers_TvsUnk$p_val/2, lower.tail = FALSE ) * sign( markers_TvsUnk$avg_log2FC )
	names(markers_TvsUnk.signature) <- markers_TvsUnk$gene
	
	markers_CvsUnk <- FindMarkers(seurat_obj.expr,ident.1 = "stem",ident.2 = "unk" ,
																features = rownames(seurat_obj.expr@assays$SCT@data) ,logfc.threshold = 0 )
	markers_CvsUnk <- as_tibble(markers_CvsUnk %>% as.data.frame() %>% rownames_to_column("gene") )
	markers_CvsUnk
	markers_CvsUnk.signature <- qnorm( markers_CvsUnk$p_val/2, lower.tail = FALSE ) * sign( markers_CvsUnk$avg_log2FC )
	names(markers_CvsUnk.signature) <- markers_CvsUnk$gene
	
	genes_in_common <- intersect( names(markers_TvsUnk.signature) , names(markers_CvsUnk.signature) )
	
	markers_TvsUnk.signature <- markers_TvsUnk.signature[genes_in_common]
	markers_CvsUnk.signature <- markers_CvsUnk.signature[genes_in_common]
	
	stopifnot(identical(names(markers_TvsUnk.signature),names(markers_CvsUnk.signature)))
	
	x <- cbind(TvsUnk=markers_TvsUnk.signature,CvsUnk=markers_CvsUnk.signature)
	vp <- viper( x , networks.list[["c1"]] , method = "none" )
	df <- as.matrix(vp %>% mouse_to_human())

	x <- df %>% as.data.frame() %>% rownames_to_column("protein") %>% as_tibble()
	x$diff <- x$TvsUnk - x$CvsUnk
	x <- x %>% arrange(desc(diff))
	x$diff_2 <- ifelse( x$TvsUnk > 1.5 & x$CvsUnk < -1.5 | x$CvsUnk > 1.5 & x$TvsUnk < -1.5 , "marked" , "no" )
	table(x$diff_2)
	
	filename <- file.path(reports.dir,paste0("mrs-homeo-injury-sorted.xlsx"))
	# markers_TvsC$sample_id <- file_tag
	writexl::write_xlsx(x,filename)
	
	my_pathways <- do_pathways_with_aREA( df , threshold = 0 )
	# my_pathways <- do_pathways_with_aREA( df )
	# my_pathways["REACTOME_SIGNALING_BY_HIPPO",]	
	
	index <- !grepl("HALLMARK" , rownames(my_pathways) )
	my_pathways <- my_pathways[index,]
	
	x <- my_pathways %>% as.data.frame() %>% rownames_to_column("pathway") %>% as_tibble()
	x$diff <- x$TvsUnk - x$CvsUnk
	x <- x %>% arrange(desc(diff))
	x$diff_2 <- ifelse( x$TvsUnk > 1.5 & x$CvsUnk < -1.5 | x$CvsUnk > 1.5 & x$TvsUnk < -1.5 , "marked" , "no" )
	table(x$diff_2)
	
	filename <- file.path(reports.dir,paste0("pathways-homeo-injury-sorted.xlsx"))
	# markers_TvsC$sample_id <- file_tag
	writexl::write_xlsx(x,filename)
	
	
	