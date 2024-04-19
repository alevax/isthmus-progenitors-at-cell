
# source("sources/ermanno-data-analysis/cytotrace-analysis/cytotrace-across-samples.R")

library(dplyr)
library(Seurat)

source("~/Downloads/CytoTRACE/R/zzz.R")
source("~/Downloads/CytoTRACE/R/CytoTRACE.R")
source("~/Downloads/CytoTRACE/R/plotCytoGenes.R")
source("~/Downloads/CytoTRACE/R/plotCytoTRACE.R")

source("../vaxtools/R/utils.R")


create_workspace( "cytotrace-and-stemness-across-samples" )

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
mat_TE002 <- getCountsLocally("TE002")
mat_TE005 <- getCountsLocally("TE005")
mat_TE006 <- getCountsLocally("TE006")

my_sample_list <- list("WT"=mat_TE001,"Rad"=mat_TE002,"Haplo"=mat_TE005,"Ablation"=mat_TE006)
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

write_tsv( my_tibble %>% filter(sample_id %in% "WT") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE001-cytotrace.tsv" )
write_tsv( my_tibble %>% filter(sample_id %in% "Rad") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE002-cytotrace.tsv" )
write_tsv( my_tibble %>% filter(sample_id %in% "Haplo") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE005-cytotrace.tsv" )
write_tsv( my_tibble %>% filter(sample_id %in% "Ablation") , "~/Clouds/Dropbox/Data/isc/cytotrace-and-stemmness/TE006-cytotrace.tsv" )

# ## Comparison between TE001/TE002 ----	
# 
# mat_TE001 <- getCountsLocally("TE001")
# mat_TE002 <- getCountsLocally("TE002")
# 
# my_sample_list <- list("WT"=mat_TE001,"Rad"=mat_TE002)
# str(my_sample_list,1)
# 
# my_sample_list <- lapply( 1:length(my_sample_list) , 
# 													function(i) {
# 														ret <- my_sample_list[[i]]
# 														colnames(ret) <- paste0( names(my_sample_list)[i] , "_" , colnames(my_sample_list[[i]]) )
# 														return(ret)
# 													})
# 
# features_in_common <- Reduce( "intersect" , lapply(my_sample_list,rownames) )
# print_msg_info( "Total number of features in common: # " , length(features_in_common) )
# features_in_common <- sort(features_in_common)
# my_big_mat <- do.call( cbind , lapply( my_sample_list , function(x) x[features_in_common,] ) )
# print(dim(my_big_mat))
# 
# print_msg_info(">>> >> Running CytoTrace  ...")
# cytotrace.data <- CytoTRACE(my_big_mat , enableFast = FALSE , ncores = 4 )
# 
# # str(cytotrace.data,1)
# cor( cytotrace.data$CytoTRACE , cytotrace.data$CytoTRACErank )
# # plot( cytotrace.data$CytoTRACE , cytotrace.data$CytoTRACErank )
# tail( sort( cytotrace.data$cytoGenes ) , 20 )
# 
# my_tibble <- tibble( sample_id = gsub("(.*)_(.*)","\\1",names(cytotrace.data$CytoTRACE) ) ,
# 										 index = gsub("(.*)_(.*)","\\2",names(cytotrace.data$CytoTRACE) ),
# 										 cytotrace = cytotrace.data$CytoTRACE
# )
# 
# write_tsv( my_tibble %>% filter(sample_id %in% "WT") , "~/Clouds/Dropbox/Data/isc/cytotrace/irradiation/TE001-cytotrace.tsv" )
# write_tsv( my_tibble %>% filter(sample_id %in% "Rad") , "~/Clouds/Dropbox/Data/isc/cytotrace/irradiation/TE002-cytotrace.tsv" )

