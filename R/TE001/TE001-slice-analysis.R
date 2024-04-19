# devtools::install_github("xu-lab/SLICE")

# loading SLICE functions
library(SLICE)
require(Seurat)

message(">>> SLICE Analysis at Gene Expression")
{
	my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")
	my_metadata.filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-analysis-with-metacell-metadata.rds"
	my_metadata.obj <- readRDS(my_metadata.filename)
	stopifnot( identical( colnames(my_seurat.obj) , my_metadata.obj$sample_id ) )
	
	# require(os)
	file.info(my_metadata.filename)
	
	require(pryr)
	object_size(my_seurat.obj)

	x <- as.matrix( my_seurat.obj@assays$RNA@counts %>% RelativeCounts(scale.factor = 1e6) )
	# x <- as.matrix( my_seurat.obj@assays$SCT@scale.data )
	dim(x)
	stopifnot( identical( colnames(x) , names(my_metadata.obj$cluster_id) ))

	sc.obj = construct( exprmatrix = as.data.frame(x),
											cellidentity = my_metadata.obj$cluster_id )
	# cellidentity = my_metadata.obj$gene_expression_clusters )

	km = cor(t(x),method="pearson")
	dim(km)
	### In entropy calculation, Gene expression = positive viper scores
	sc.obj <- getEntropy( sc.obj,
												km = km , # use the pre-computed kappa similarity matrix of mouse genes
												calculation = "bootstrap",# choose the bootstrap calculation
												B.num = 20 , # default: 100 iterations
												exp.cutoff = 0 , # the threshold for expressed genes
												B.size = 1000 ,  # the size of bootstrap sample
												clustering.k = floor(sqrt(1000/2)),# the number of functional clusters
												# clustering.k = 50 ,
												# clustering.k = 7,# the number of functional clusters
												random.seed = 42 )

	my_metadata.obj$sc_entropy_ges <- sc.obj@entropies$scEntropy.bootstrap
	saveRDS(my_metadata.obj,my_metadata.filename)
	plotEntropies(sc.obj)
}

message(">>> SLICE Analysis at Protein Activity")
{
	# my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-viper-analysis-data.rds")
	my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-viper-analysis-with-metacell-data.rds")
	my_metadata.filename <- "~/Clouds/Dropbox/Data/isc/TE001/TE001-viper-analysis-with-metacell-metadata.rds"
	my_metadata.obj <- readRDS(my_metadata.filename)
	stopifnot( identical( colnames(my_seurat.obj) , my_metadata.obj$sample_id ) )
	
	require(pryr)
	object_size(my_seurat.obj)
	
	x <- as.matrix(my_seurat.obj@assays$VIPER@data)
	
	table(my_seurat.obj$pas_cluster_id)
	sc.obj = construct( exprmatrix = as.data.frame(x), 
											cellidentity = my_seurat.obj$pas_cluster_id)
	km = cor(t(x),method="pearson")
	dim(km)
	### In entropy calculation, Gene expression = positive viper scores
	sc.obj <- getEntropy( sc.obj, 
												km = km , # use the pre-computed kappa similarity matrix of mouse genes
												calculation = "bootstrap",# choose the bootstrap calculation
												B.num = 20 , # default: 100 iterations
												exp.cutoff = 0 , # the threshold for expressed genes
												B.size = 1000 ,  # the size of bootstrap sample
												clustering.k = floor(sqrt(1000/2)),# the number of functional clusters
												# clustering.k = 50 ,
												# clustering.k = 7,# the number of functional clusters
												random.seed = 42 )
	
	my_metadata.obj$sc_entropy_pas <- sc.obj@entropies$scEntropy.bootstrap
	saveRDS(my_metadata.obj,my_metadata.filename)
	plotEntropies(sc.obj)
	
}
# no_metacell <- x.entropy$entropy1
# names(no_metacell) <- rownames(x.entropy)
# head(no_metacell)
# 
# with_metacell <- x.entropy$entropy1
# names(with_metacell) <- rownames(x.entropy)
# head(with_metacell)
# 
# plot(no_metacell,with_metacell)

# plotEntropies(sc.obj)

# df <- my_metadata.obj$sc_entropy %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
# stopifnot( identical( df$cell_id , colnames(x) ) )
# df$cluster_id <- my_seurat.obj$seurat_clusters
# df$Lgr5 <- x["Lgr5",]
# df$Top2a <- x["Top2a",]
# df$Krt19 <- x["Krt19",]

# l <- levels(df$cluster_id)
# cluster_colors <- c(brewer.pal(9,"Set1"),brewer.pal(6,"Set2"))[1:length(l)]
# names(cluster_colors) <- l
# p <- ggplot( df , aes(entropy1, scEntropy.bootstrap , 
# 											# color = Lgr5
# 											# fill = Top2a
# 											fill = cluster_id
# ) ) +
# 	# geom_point(aes(PHATE1, PHATE2, color=x$Mpo)) +
# 	# labs(color="Mpo") +
# 	geom_point(shape = 21, colour = "white",size = 2, stroke = 0.25) +
# 	# scale_fill_viridis() +
# 	# scale_fill_discrete() +
# 	scale_fill_brewer(palette = "Set1") +
# 	coord_fixed(ratio = 1) +
# 	theme_light()
# 
# pdf( file.path(reports.dir,"TE001-pas-slice.pdf"))
# 	print(p)
# dev.off()

# # use a variable to store the path to the data directory; getwd() functions returns the full path of current working directory
# #data.dir <- paste(getwd(),"/data/", sep="")
# 
# # use a variable to store the name of the dataset and analysis
# data.name <- "TE001-gene-expression"
# 
# cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Start SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))
# 
# # a string describing the analysis and will be attached to the names of all the files generated during the execution
# context_str = paste("SLICE-", data.name, "-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="") 
# 
# 
# ###########################################################################
# #####                                                               #######
# #####               Loading data and preprocessing                  #######
# #####                                                               #######
# ###########################################################################
# 
# 
# # loading the FB data set; expression profiles with same gene symbol annotations have been averaged
# data(FB)
# 
# library(pryr)
# library("BisqueRNA") # SeuratToExpressionSet
# library(Seurat)
# 
# # isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/TE001/")
# my_seurat.obj <- readRDS("~/Clouds/Dropbox/Data/isc/TE001/TE001-seurat-analysis-data.rds")
# object_size(my_seurat.obj)
# my_expr_set <- SeuratToExpressionSet(my_seurat.obj,version = "v3",delimiter = "",position = 1)
# # the expressions, cell info, and gene info are stored as an Biobase::ExpressionSet object
# es <- FB
# 
# # remove ERCC genes and Ribosomal genes
# ercc.genes <- grep("^ERCC-", rownames(fData(es)), value = TRUE)
# rb.genes <- grep("^Rpl|^Rps|^Mrpl|^Mrps", rownames(fData(es)), value = TRUE)
# es <- es[which(!(rownames(fData(es)) %in% c(rb.genes, ercc.genes))), ]
# 
# 
# ###########################################################################
# #####                                                               #######
# #####                   Initlializing SLICE object                  #######
# #####                                                               #######
# ###########################################################################
# 
# 
# # use the expression matrix of FB cells (n=79) to create a SLICE object
# # use the predicted fibroblast subtype information from the original analysis (encoded in the Type attribute of the pData) as the original identity
# # set projname to describe this analysis
# sc <- construct(exprmatrix=as.data.frame(exprs(es)), 
# 								cellidentity=factor(pData(es)$Type, levels=c("PMP","IF1","IF2","MyoF","MFB")),
# 								projname=context_str)
# 
# ###########################################################################
# #####                                                               #######
# #####     Measuring cell differentiation state using scEntropy      #######
# #####                                                               #######
# ###########################################################################
# 
# 
# # loading the pre-computed mouse gene-gene Kappa similarity matrix
# # The matrix will be loaded into a variable named "km"
# data(mm_kappasim,verbose = TRUE)
# 
# # bootstrap calculation of scEntropy
# sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of mouse genes
# 								 calculation="bootstrap",               # choose the bootstrap calculation
# 								 B.num=100,                             # 100 iterations
# 								 exp.cutoff=1,                          # the threshold for expressed genes
# 								 B.size=1000,                           # the size of bootstrap sample
# 								 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
# 								 random.seed=201602)                    # set the random seed to reproduce the results in the paper
# 
# 
# # plot the entropy; will be save to a pdf named with "entropies" suffix in the working directory
# plotEntropies(sc)
# 
# 
# ###########################################################################
# #####                                                               #######
# #####                      Lineage reconstruction                   #######
# #####                                                               #######
# ###########################################################################
# 
# 
# # perform PCA using the expression of predicted signature genes and with a high variance (greater than 8 variance in log2 FPKM)
# # The predicted signature genes are stored in FB.sig in FB.rda in the data folder
# sc <- getRDS(sc, genes.use=FB.sig, method="pca", num_dim=2, log.base=2, do.center=TRUE, do.scale=FALSE, use.cor=TRUE, min.var=8, min.cells=0)
# 
# # use the graph-based method in SLICE to infer lineage model
# # results will be plotted on screen and saved in the PDF files with file names containing "lineage-g"
# sc <- getLineageModel(sc, lm.method="graph",                          # select graph-based method
# 											model.type="tree",                              # infer mst based lineage model
# 											reverse=F,                                      # infer differentiation lineage
# 											ss.method="pcst", ss.threshold=0.25,            # use linear prize collecting tree algorithm to find core cell sets
# 											wiring.threshold=function(mst) max(mst)*1.2,    # the threshold for local wiring
# 											community.method="louvain")                     # use the louvain algorithm to detect cell clusters
# 
# # marker genes for lineage validation
# markers <- c("Fn1", "Tcf21", "Vcam1","Nr3c1", "Pdgfra", "Myocd", "Actg2", "Myh11", "Foxm1", "Top2a")
# 
# 
# # use the shortest-path method in SLICE to reconstruct cell transitional path following C1->C4
# # results will be plotted on screen and saved in the PDF files with file names containing "path-sp"
# sc <- getTrajectories(sc, method="sp", start=1, end=4, network="mst", NN.threshold=0.8, NN.type="all", do.plot=T)
# 
# # extract and visualize expression profiles of marker genes in the reconstructed shrotest-path transitional path
# # results will be plotted on screen and saved in the PDF files with file names containing "path-sp" and "profiles"
# sc <- getProfiles(sc, trajectory.type="sp", genes.use=markers)
# 
# 
# if (TRUE) { # This is to reproduce the differentiation expression and temporal pattern analysis of the PMP->SM branch in Figure 6. Might take ~15 minutes
# 	
# 	# extract expression profiles of all genes in the reconstructed shrotest-path transitional path
# 	# no plot will be generated due to the large number of genes
# 	sc <- getProfiles(sc, trajectory.type="sp")
# 	
# 	# detect lineage depenent differentially expressed genes
# 	# criteria: expressed (>thresh.exp) in at least 2 (thresh.cell) cells, log2 variance >=0.5 (thresh.var), and FDR <0.1 (diff.thresh)
# 	df.sm <- getDiffGenes(sc@profiles, diff.thresh=0.1, log.base=2, thresh.exp=1, thresh.cell=2, thresh.var=0.5)
# 	
# 	# identy temporal patterns of differentially expressed genes
# 	# results will be plotted on screen and saved in a PDF file in the working directory
# 	pn.sm <- getPatterns(df.sm[[1]]$exprs.fit.diff, c.method="pam", k=3, plot.filename = paste(context_str, "path-diffgenes-patterns.sm.pdf", sep=""))
# 	
# 	# save differential analysis results
# 	pn.sm.df <- data.frame(gene=names(pn.sm$clustering), pattern=as.numeric(pn.sm$clustering), check.names=FALSE)
# 	diff.sm <- df.sm[[1]]$diff
# 	diff.sm <- diff.sm[which(diff.sm$significant==1), ]
# 	write.table(cbind(pn.sm.df, diff.sm), file=paste(context_str, "path-diffgenes-patterns.sm.txt", sep=""), sep="\t", col.names=T, row.names=F)
# }
# 
# 
# # use the shortest-path method in SLICE to reconstruct cell transitional path following C1->C3->C2
# # results will be plotted on screen and saved in the PDF files with file names containing "path-sp"
# sc <- getTrajectories(sc, method="sp", start=1, end=2, network="mst", NN.threshold=0.8, NN.type="all", do.plot=T)
# 
# # extract and visualize expression profiles of marker genes in the reconstructed shrotest-path transitional path
# # results will be plotted on screen and saved in the PDF files with file names containing "path-sp" and "profiles"
# sc <- getProfiles(sc, trajectory.type="sp", genes.use=markers)
# 
# 
# if (TRUE) { # This is to reproduce the differentiation expression and temporal pattern analysis of the PMP->MFB branch in Figure 6. Might take ~15 minutes
# 	
# 	# extract expression profiles of all genes in the reconstructed shrotest-path transitional path
# 	# no plot will be generated due to the large number of genes
# 	sc <- getProfiles(sc, trajectory.type="sp")
# 	
# 	# detect lineage depenent differentially expressed genes
# 	# criteria: expressed (>thresh.exp) in at least 2 (thresh.cell) cells, log2 variance >=0.5 (thresh.var), and FDR <0.1 (diff.thresh)
# 	df.mfb <- getDiffGenes(sc@profiles, diff.thresh=0.1, log.base=2, thresh.exp=1, thresh.cell=2, thresh.var=0.5)
# 	
# 	# identy temporal patterns of differentially expressed genes
# 	# results will be plotted on screen and saved in a PDF file in the working directory
# 	pn.mfb <- getPatterns(df.mfb[[1]]$exprs.fit.diff, c.method="pam", k=3, plot.filename = paste(context_str, "path-diffgenes-patterns.mfb.pdf", sep=""))
# 	
# 	# save differential analysis results
# 	pn.mfb.df <- data.frame(gene=names(pn.mfb$clustering), pattern=as.numeric(pn.mfb$clustering), check.names=FALSE)
# 	diff.mfb <- df.mfb[[1]]$diff
# 	diff.mfb <- diff.mfb[which(diff.mfb$significant==1), ]
# 	write.table(cbind(pn.mfb.df, diff.mfb), file=paste(context_str, "path-diffgenes-patterns.mfb.txt", sep=""), sep="\t", col.names=T, row.names=F)
# 	
# }
# 
# 
# 
# # use the clustering-based method in SLICE to infer lineage model
# # results will be plotted on screen and saved in the PDF files with file names containing "lineage-c"
# sc <- getLineageModel(sc.obj, lm.method="clustering", model.type="tree", reverse=F, ss.threshold=0.25, cluster.method="pam", k=4)
# 
# # use the principal-curve method in SLICE to reconstruct cell transitional path following C2->C3
# # results will be plotted on screen and saved in the PDF files with file names containing "path-pc"
# sc <- getTrajectories(sc, method="pc", start=2, end=3, do.trim=T,  do.plot=T, stepwise=T)
# 
# # extract and visualize expression profiles of marker genes in the reconstructed principal-curve transitional path
# # results will be plotted on screen and saved in the PDF files with file names containing "path-pc" and "profiles"
# sc <- getProfiles(sc, trajectory.type="pc", genes.use=markers)
# 
# # use the principal-curve method in SLICE to reconstruct cell transitional path following C2->C4->C1
# # results will be plotted on screen and saved in the PDF files with file names containing "path-pc"
# sc <- getTrajectories(sc, method="pc", start=2, end=1, do.trim=T,  do.plot=T)
# 
# # extract and visualize expression profiles of marker genes in the reconstructed principal-curve transitional path
# # results will be plotted on screen and saved in the PDF files with file names containing "path-pc" and "profiles"
# sc <- getProfiles(sc, trajectory.type="pc", genes.use=markers)
# 
# 
# cat(paste("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n   Finish SLICE analysis of ", data.name, "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n", sep=""))
