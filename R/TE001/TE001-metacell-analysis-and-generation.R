
# source("sources/ermanno-data-analysis/TE001/TE001-metacell-analysis-and-generation.R")

source("../vaxtools/R/utils.R")

library(metacell)
library(Seurat)

set.seed(42)

knn_nn <- 51
metacell_min_size <- 5

create_workspace("metacell-with-metacell-package")
scdb_init(processed_data.dir, force_reinit=T)
scfigs_init(reports.dir)

my_sample_id <- "TE001"
isc_data_dir <- file.path("~/Clouds/Dropbox/Data/isc/",my_sample_id)
my_path <- file.path(isc_data_dir,"/filtered_feature_bc_matrix/unzipped/")
mat <- mcell_import_scmat_10x(mat_nm = "TE001", 
															matrix_fn = "~/Clouds/Dropbox/Data/isc//TE001/filtered_feature_bc_matrix/unzipped/matrix.mtx" ,
															genes_fn = "~/Clouds/Dropbox/Data/isc//TE001/filtered_feature_bc_matrix/unzipped/features.tsv" ,
															cells_fn = "~/Clouds/Dropbox/Data/isc//TE001/filtered_feature_bc_matrix/unzipped/barcodes.tsv" )

# seurat_analysis_data_filename <- file.path( isc_data_dir , paste0(my_sample_id,"-seurat-analysis-data.rds") )
# sc_data <- readRDS(seurat_analysis_data_filename)

mat = scdb_mat("TE001")
print(dim(mat@mat))
mcell_plot_umis_per_cell("TE001")

nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^Igj", nms, v=T),
                grep("^Igh",nms,v=T),
                grep("^Igk", nms, v=T),
                grep("^Igl", nms, v=T))
immuno_genes = c(grep("^Ptprc$", nms, v=T))

bad_genes <- unique(c(grep("^mt-", nms, v=T),grep("^Mt-", nms, v=T), grep("^Mtmr", nms, v=T), grep("^Mtnd", nms, v=T),"Neat1","Tmsb4x", "Tmsb10", ig_genes , immuno_genes ))
# bad_genes
mcell_mat_ignore_genes(new_mat_id="TE001", mat_id="TE001", bad_genes, reverse=F)
mcell_mat_ignore_small_cells("TE001", "TE001", 1000)

mat = scdb_mat("TE001")
print(dim(mat@mat))

mcell_add_gene_stat(gstat_id="TE001", mat_id="TE001", force=T)

mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="TE001", T_vm=0.05, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="TE001", T_tot=200, T_top3=2)
mcell_plot_gstats(gstat_id="TE001", gset_id="test_feats")

mcell_add_cgraph_from_mat_bknn(mat_id="TE001",
                gset_id = "test_feats",
                graph_id="test_graph",
                K=knn_nn,
                dsamp=F)
mcell_coclust_from_graph_resamp(
                coc_id="test_coc500",
                graph_id="test_graph",
                min_mc_size=metacell_min_size,
                p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
                coc_id="test_coc500",
                mat_id= "TE001",
                mc_id= "test_mc",
                K=knn_nn, min_mc_size=metacell_min_size, alpha=2)

mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "TE001", T_lfc=3)

mcell_mc_split_filt(new_mc_id="test_mc",
            mc_id="test_mc",
            mat_id="TE001",
            T_lfc=3, plot_mats=F)

mcell_gset_from_mc_markers(mc_id="test_mc",gset_id="test_markers")
marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
mc_colorize("test_mc", marker_colors=marks_colors)
mcell_mc_plot_marks(mc_id="test_mc", gset_id="test_markers", mat_id="TE001")
mc = scdb_mc("test_mc")
table(mc@colors)

mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc", graph_id="test_graph",ignore_mismatch = TRUE)
#> comp mc graph using the graph test_graph and K 20
#> Missing coordinates in some cells that are not ourliers or ignored - check this out! (total 1 cells are missing, maybe you used the wrong graph object? first nodes CCCTCCTAGTCGCCGT
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="test_2dproj")

mc_hc = mcell_mc_hclust_confu(mc_id="test_mc",graph_id="test_graph")
mc_sup = mcell_mc_hierarchy(mc_id="test_mc",mc_hc=mc_hc, T_gap=0.04)
mcell_mc_plot_hierarchy(mc_id="test_mc",
                   graph_id="test_graph",
                    mc_order=mc_hc$order,
                    sup_mc = mc_sup,
                    width=2800, heigh=2000, min_nmc=2)

x <- mcell_mc_cell_homogeneity(mc_id="test_mc",graph_id="test_graph")

mc = scdb_mc("test_mc")


