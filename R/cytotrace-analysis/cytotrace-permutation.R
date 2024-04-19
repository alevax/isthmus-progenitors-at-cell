

		## CytoTrace Analysis ----
		print_msg_info(">>> CytoTrace Analysis")
		{
			source("~/Downloads/CytoTRACE/R/zzz.R")
			source("~/Downloads/CytoTRACE/R/CytoTRACE.R")
			source("~/Downloads/CytoTRACE/R/plotCytoGenes.R")
			source("~/Downloads/CytoTRACE/R/plotCytoTRACE.R")
			
			set.seed(666)
			
			# sc_data <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE001-TE002-gex.rds")
			sc_data <- readRDS("~/Clouds/Dropbox/Data/isc/preprocessed/anchoring-TE005-TE006-gex.rds")
			table(sc_data$sample_id)
			
			sc_control <- sc_data %>% subset( sample_id == "control" )
			sc_treatment <- sc_data %>% subset( sample_id == "treatment" )
			sc_control <- sc_control@assays$RNA@counts %>% as.matrix()
			sc_treatment <- sc_treatment@assays$RNA@counts %>% as.matrix()
			dim(sc_control)
			dim(sc_treatment)
			
			n_cells <- 5000
			sc_control <- sc_control[ , sample(colnames(sc_control) , size = n_cells,replace = FALSE) ]
			sc_treatment <- sc_treatment[ , sample(colnames(sc_treatment) , size = n_cells,replace = FALSE) ]
			
			dim(sc_control)
			dim(sc_treatment)
			
			n_cells <- ncol(sc_treatment)
			
			n_iter <- 10
			cyto_mat <- matrix(0,nrow = n_iter,ncol = ncol(sc_treatment) , dimnames = list(1:n_iter,colnames(sc_treatment)) )
			
			for (i in 1:n_iter)
			{
				sc_control_subset <- sc_control[ , sample(colnames(sc_control) , size = n_cells,replace = FALSE) ]
				dim(sc_control_subset)
				dim(sc_treatment)
				
				mat <- cbind(sc_control_subset,sc_treatment)
				dim(mat)
				
				print_msg_info(">>> >> Running CytoTrace iteration: " , i)
				cytotrace.data <- CytoTRACE(mat , enableFast = FALSE , ncores = 10 )
				
				cyto_mat[i,] <- cytotrace.data$CytoTRACE[ match( colnames(cyto_mat) , names(cytotrace.data$CytoTRACE) ) ]
			}
			
			cyto_final <- apply( cyto_mat , 2 , mean )

			names(cyto_final) <- gsub("treatment_","",names(cyto_final))
			
			# my_ingest <- read_csv("~/Clouds/Dropbox/Data/isc/rad-metadata-ingest.csv")
			my_ingest <- read_csv("~/Clouds/Dropbox/Data/isc/ablation-metadata-ingest.csv")
			
			my_ingest$cyto_perm <- cyto_final[ match( my_ingest$cell_id , names(cyto_final) ) ]
			
			# write_csv(my_ingest,"~/Clouds/Dropbox/Data/isc/rad-metadata-ingest-cytoperm.csv")
			write_csv(my_ingest,"~/Clouds/Dropbox/Data/isc/ablation-metadata-ingest-cytoperm.csv")
			
			ggplot( my_ingest , aes( x = iter_cluster_id , y = cyto_perm , fill = iter_cluster_id) ) +
				geom_point() + 
				geom_boxplot() +
				scale_fill_brewer(palette = "Set3") + 
				theme_light()
			
		}	
		
		