##
# Analysis of Ayyaz 2019 samples with SingleR
# ---------------------------------------------
# system.time({source("sources/ayyaz-data-analysis/ayyaz-SingleR.R")})

source("../vaxtools/R/utils.R")

create_workspace( isWorkspaceToClean = FALSE , 
									experiments_dir = "experiments" , 
									run_dir = "ayyaz-SingleR" )

## ---- C05 ----

sc_data <- readRDS("/Users/afpvax/Library/CloudStorage/Dropbox/Data/isc/ayyaz-2019/C05/C05-seurat-original-data.rds")
# sc_data <- readRDS("/Users/afpvax/Library/CloudStorage/Dropbox/Data/isc/ayyaz-2019/C07/C07-seurat-original-data.rds")

library(SingleR)
library(celldex)
bped.se <- BlueprintEncodeData()
bped.se

my_mat <- as.matrix(sc_data@assays$RNA@counts)
system.time({
	sR_BP_predicted_labels <- SingleR( test = my_mat %>% mouse_to_human(), 
																		 ref = bped.se, 
																		 # fine.tune = FALSE,prune=T,tune.thresh=0.1,
																		 # genes = "sd" ,
																		 # assay.type.test=1,
																		 labels = bped.se$label.main)
})
table(sR_BP_predicted_labels$labels)

library(ggplot2)
hist_plot <- ggplot(data.frame(x = sR_BP_predicted_labels$labels), aes(x)) +
	geom_histogram(stat = "count" , fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
	theme_minimal() +
	labs(
		title = "C05",
		x = "Cell types",
		y = "Frequency"
	) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.line = element_line(color = "black"),
		axis.text = element_text(color = "black"),
		axis.text.x = element_text(color = "black" , angle = 90 , vjust = 0.5, hjust=1 ),
		axis.title = element_text(color = "black", size = 12),
		plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
		legend.position = "none"
	)

pdf( file.path(reports.dir , "SingleR-Ayyaz-C05.pdf") , width = 3 , height = 4 )
	print(hist_plot)
dev.off()

ncol(sc_data@assays$RNA@counts)
sum(sc_data@assays$RNA@counts["Ptprc",] > 0 )
sum( sc_data@assays$RNA@counts["Ptprc",] > 0 | sc_data@assays$RNA@counts["Cd8a",] > 0 )
sum( ( sc_data@assays$RNA@counts["Ptprc",] > 0 | sc_data@assays$RNA@counts["Cd8a",] > 0 ) & sc_data@assays$RNA@counts["Epcam",] == 0 )

ncol(sc_data@assays$RNA@counts) - sum( sc_data@assays$RNA@counts["Ptprc",] > 0 | sc_data@assays$RNA@counts["Cd8a",] > 0 )

sum( ( sc_data@assays$RNA@counts["Ptprc",] > 0 & sc_data@assays$RNA@counts["Epcam",] == 0 ) )

## ---- C07 ----

sc_data <- readRDS("/Users/afpvax/Library/CloudStorage/Dropbox/Data/isc/ayyaz-2019/C07/C07-seurat-original-data.rds")
# sc_data <- readRDS("/Users/afpvax/Library/CloudStorage/Dropbox/Data/isc/ayyaz-2019/C07/C07-seurat-original-data.rds")

library(SingleR)
library(celldex)
bped.se <- BlueprintEncodeData()
bped.se

my_mat <- as.matrix(sc_data@assays$RNA@counts)
system.time({
	sR_BP_predicted_labels <- SingleR( test = my_mat %>% mouse_to_human(), 
																		 ref = bped.se, 
																		 # fine.tune = FALSE,prune=T,tune.thresh=0.1,
																		 # genes = "sd" ,
																		 # assay.type.test=1,
																		 labels = bped.se$label.main)
})
table(sR_BP_predicted_labels$labels)

library(ggplot2)
hist_plot <- ggplot(data.frame(x = sR_BP_predicted_labels$labels), aes(x)) +
	geom_histogram(stat = "count" , fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
	theme_minimal() +
	labs(
		title = "C07",
		x = "Cell types",
		y = "Frequency"
	) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		axis.line = element_line(color = "black"),
		axis.text = element_text(color = "black"),
		axis.text.x = element_text(color = "black" , angle = 90 , vjust = 0.5, hjust=1 ),
		axis.title = element_text(color = "black", size = 12),
		plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
		legend.position = "none"
	)

pdf( file.path(reports.dir , "SingleR-Ayyaz-C07.pdf") , width = 3 , height = 4 )
print(hist_plot)
dev.off()

ncol(sc_data@assays$RNA@counts)
sum(sc_data@assays$RNA@counts["Epcam",] > 0 )
sum(sc_data@assays$RNA@counts["Ptprc",] > 0 )
sum( sc_data@assays$RNA@counts["Ptprc",] > 0 | sc_data@assays$RNA@counts["Cd8a",] > 0 )
sum( ( sc_data@assays$RNA@counts["Ptprc",] > 0 | sc_data@assays$RNA@counts["Cd8a",] > 0 ) & sc_data@assays$RNA@counts["Epcam",] == 0 )

ncol(sc_data@assays$RNA@counts) - sum( sc_data@assays$RNA@counts["Ptprc",] > 0 | sc_data@assays$RNA@counts["Cd8a",] > 0 )

sum( ( sc_data@assays$RNA@counts["Ptprc",] > 0 & sc_data@assays$RNA@counts["Epcam",] == 0 ) )

