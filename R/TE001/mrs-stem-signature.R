##
# Generating MRs signature of Instestinal Stem Cells (ISC)
# -------------------------------------------------------
# source("sources/ermanno-data-analysis/TE001/mrs-stem-signature.R")

source("../vaxtools/R/utils.R")

create_workspace("stem-signatures")

library(logger)
log_threshold(DEBUG)
log_threshold(DEBUG,index=2)

log_filename <- file.path(reports.dir,"log.txt")
log_appender(appender = appender_file(log_filename),index=2)
# log_formatter(formatter_glue)
log_layout(layout_glue_colors)
log_warn('--- Start of Analysis ----')
# log_info("is_viper_with_subset_of_cells: {is_viper_with_subset_of_cells}")

require(tidyverse)
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-cpm.rds"
fs::file_info(filename)
cpm_matrix <- readRDS(filename)
dim(cpm_matrix)
filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-cytotrace.csv"
fs::file_info(filename)
cytotrace_table <- read_delim(filename,"\t")

stem_supertop_cell_ids <- cytotrace_table %>% 
	filter(cytotrace_score >= 0.995) %>%
	pull(cell_id)
length(stem_supertop_cell_ids)

stem_cell_ids <- cytotrace_table %>% 
	filter(cytotrace_score >= 0.95) %>%
	pull(cell_id)
length(stem_cell_ids)

progenitor_cell_ids <- cytotrace_table %>% 
	filter(cytotrace_score < 0.85) %>%
	filter(cytotrace_score > 0.8) %>%
	pull(cell_id)
length(progenitor_cell_ids)

terminal_cell_ids <- cytotrace_table %>% 
	filter(cytotrace_score < 0.05) %>%
	filter(cytotrace_score > 0) %>%
	pull(cell_id)
length(terminal_cell_ids)

require(viper)
stem_network <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression/TE001_c1_unPruned.rds")
# stem_network <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/TE001_mc_unPruned.rds")

tfs <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/tf-mus-current-symbol.dat" , col_names = FALSE )$X1
cotfs <- read_csv("~/Clouds/Dropbox/Data/isc/TE001/networks/run-with-metacell-from-gene-expression-one-network-only/lists/cotf-mus-current-symbol.dat" , col_names = FALSE )$X1
# z <- read_csv( file.path(isc_data_dir,"networks/run-with-metacell-from-gene-expression-one-network-only/lists/surface-mus-current-symbol.dat") , col_names = FALSE )

index <- names(stem_network) %in% c(tfs,cotfs)
sum(index)
stem_network <- stem_network[index]
stem_network <- pruneRegulon(stem_network,50)

log_info("Total TF/coTFs: {length(stem_network)}")

# --

log_info("Signature for VIPER: stem vs progenitors")
# log_cpm_matrix <- log(cpm_matrix)
vp_stem <- viper( cpm_matrix[,stem_cell_ids] ,
									cpm_matrix[,progenitor_cell_ids] ,
									regulon = stem_network , method = "mad")
log_info("Running VIPER using MAD as signature generation")

# vs <- viperSignature(cpm_matrix[,stem_cell_ids],
# 										 cpm_matrix[,progenitor_cell_ids],
# 										 method = "zscore")
# vp_stem <- viper(vs,regulon = stem_network , method = "none")

mrs_list <- doStouffer(vp_stem)
mrs_list <- sort(mrs_list,decreasing = T)
head(mrs_list,50)
mrs_table <- tibble(MR=names(mrs_list),score=mrs_list)
xlsx::write.xlsx( mrs_table , file.path(reports.dir,"mrs-stem-of-stem-vs-progenitor.xlsx" ) )

# --

log_info("Signature for VIPER: stem vs terminal")
vp_stem <- viper( cpm_matrix[,stem_cell_ids] ,
									cpm_matrix[,terminal_cell_ids] ,
									regulon = stem_network , method = "mad")
log_info("Running VIPER using MAD as signature generation")

# vs <- viperSignature(cpm_matrix[,stem_cell_ids],
# 										 cpm_matrix[,progenitor_cell_ids],
# 										 method = "zscore")
# vp_stem <- viper(vs,regulon = stem_network , method = "none")

mrs_list <- doStouffer(vp_stem)
mrs_list <- sort(mrs_list,decreasing = T)
head(mrs_list,50)
mrs_table <- tibble(MR=names(mrs_list),score=mrs_list)
xlsx::write.xlsx( mrs_table , file.path(reports.dir,"mrs-stem-of-stem-vs-terminal.xlsx" ) )

# --

log_info("Signature for VIPER: super stem vs progenitor")
vp_stem <- viper( cpm_matrix[,stem_supertop_cell_ids] ,
									cpm_matrix[,progenitor_cell_ids] ,
									regulon = stem_network , method = "mad")
log_info("Running VIPER using MAD as signature generation")

# vs <- viperSignature(cpm_matrix[,stem_cell_ids],
# 										 cpm_matrix[,progenitor_cell_ids],
# 										 method = "zscore")
# vp_stem <- viper(vs,regulon = stem_network , method = "none")

mrs_list <- doStouffer(vp_stem)
mrs_list <- sort(mrs_list,decreasing = T)
head(mrs_list,50)
mrs_table <- tibble(MR=names(mrs_list),score=mrs_list)
xlsx::write.xlsx( mrs_table , file.path(reports.dir,"mrs-stem-of-superstem-vs-progenitor.xlsx" ) )

# seurat_obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")
# sct_matrix <- seurat_obj@assays$SCT@scale.data %>% as.matrix()
# 
# vp_stem <- viper( sct_matrix[,stem_cell_ids] , 
# 									sct_matrix[,progenitor_cell_ids] , 
# 									regulon = stem_network , method = "mad")
# 
# mrs_list <- doStouffer(vp_stem)
# mrs_list <- sort(mrs_list,decreasing = T)
# head(mrs_list,50)

log_warn('--- END of Analysis ----')

