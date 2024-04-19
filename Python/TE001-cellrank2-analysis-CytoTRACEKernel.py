python3

### a) Import packages and data
# a.1) setup path to data-containing folder and savings
h5ad_path = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001-h5ad/"
figures_dir = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/"
# a.2) packages
import sys
import cellrank as cr
import scanpy as sc
import numpy as np
import pandas as pd
import os
from matplotlib import rc_context
import matplotlib.pyplot as plt
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

import warnings
warnings.simplefilter("ignore", category=UserWarning)

cytotrace_markers = ['Smarca5','Rbbp7','Tcerg1','Hnrnpd','Hmg20b','Nelfe','Ube2i','Etv5','Ubn1','Mbd3','Dek','Maz',
                     'Itgb3bp','Ilf2','Pa2g4'] # Id3','Hnf4g','Atoh1','Spdef','Neurod1' markers upregulated in cytotrace (Fig 1e) 

# a.3) load counts data (exported from Seurat@RNA assay)
counts_h5ad = h5ad_path + "TE001-counts.h5ad"
adata = sc.read_h5ad(counts_h5ad) # load in object


# a.4) load metadata for TE001 
metadata_csv = h5ad_path + "TE001-metadata-umap-and-clusters-for-paper.csv"
metadata = pd.read_csv(metadata_csv)

# a.5) process metadata in adata
specified_columns = ["cell_id", "nCount_RNA", "nFeature_RNA", "mt_percent", "cytotrace_score.ges",
                     "cytotrace_gcs.ges", "S.Score", "G2M.Score", "Phase", "seurat_clusters", "singleR_labels", "stemness_index.ges"]

adata.obs = adata.obs[specified_columns]
cells_to_analyze = metadata['cell_id'] # cells to analyze
adata = adata[adata.obs_names.isin(cells_to_analyze)] # subset cells to analyze in adata

adata.obs = pd.merge(adata.obs, metadata, on='cell_id', how='left') # merge metadata and include into counts object
 
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category') # clusters as categorical variable

# a.6) set UMAP coordinates to those obtained at protein activty
umap_coordinates = np.array(adata.obs.loc[:, ['UMAP_1_scanpy','UMAP_2_scanpy']]) 
adata.obsm['X_umap'] = umap_coordinates


# a.7) Include metadata of terminal states for CellRank analysis
adata.obs['terminal_states'] = adata.obs['iter_cluster_id_with_paneth']
adata.obs['terminal_states'].iloc[adata.obs['terminal_states'].isin(["stem-1","stem-2"])] = np.nan

print("adata contains the counts for the TE001 dataset")
#counts = adata.raw.to_adata() # the adata already contains the counts matrix

# a.8) display specific marker genes
log_expression = adata.copy()
sc.pp.normalize_total(log_expression, target_sum=1e4)
sc.pp.log1p(log_expression)

sc.pl.umap(log_expression,color=["Lgr4","Lgr5"], use_raw=False, cmap='viridis',add_outline=True)


#####################


### b) Preprocess the data and UMAP visualization
sc.tl.pca(adata, random_state=0)
sc.pp.neighbors(adata, random_state=0)
sc.pl.umap(adata, color=["seurat_clusters","iter_cluster_id_with_paneth","cytotrace","stemness_index"], ncols=2, add_outline=True)

#####################

### c) CytoTRACE kernel
# c.1) Setup kernel
print("Working with CytoTRACE kernel")
from cellrank.kernels import CytoTRACEKernel
import scvelo as scv
# CytoTRACE by default uses imputed data - a simple way to compute
# k-NN imputed data is to use scVelo's moments function.
# However, note that this function expects `spliced` counts because
# it's designed for RNA velocity, so we're using a simple hack here:
if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
    adata.layers['spliced'] = adata.X
    adata.layers['unspliced'] = adata.X
scv.pp.moments(adata) # hack for CytoTRACEkernel

ctk = CytoTRACEKernel(adata) # initialize the CellRank2 kernel

# c.2) compute transition matrix
ctk = ctk.compute_cytotrace().compute_transition_matrix(threshold_scheme="soft",nu=0.5) # compute transition matrix

figures_dir_CytoTRACE = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/CR2_CytoTRACEKernel/"
os.mkdir(figures_dir_CytoTRACE)

ctk_pseudotime_figure = figures_dir_CytoTRACE + "CytoTRACE_pseudotime.pdf"
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['ct_score', 'ct_pseudotime'], show=False, add_outline=True)
    plt.savefig(ctk_pseudotime_figure)

ctk_pseudotime_vln = figures_dir_CytoTRACE + "CytoTRACE_pseudotime_vln.pdf"
with rc_context({'figure.figsize': (8.5, 8.5)}):
    sc.pl.violin(adata, keys=["ct_pseudotime"], groupby="iter_cluster_id_with_paneth", rotation=90, show=False)
    plt.savefig(ctk_pseudotime_vln)

# c.3) Simulate a random walk on the Markov chain implied by the transition matrix 
ctk_rw_figure = figures_dir_CytoTRACE + "CytoTRACE_random_walk.pdf"
ctk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs=None,
    legend_loc="right",
    dpi=100,
    save=ctk_rw_figure,
    figsize=(3,3)
)

# c.4) visualize the transition matrix
differentiation_figure = figures_dir_CytoTRACE + "CytoTRACE_differentiation_ges_clusters.png"
ctk.plot_projection(basis="umap", color="seurat_clusters", 
                    legend_loc="right", save=differentiation_figure, show=False)


differentiation_figure = figures_dir_CytoTRACE + "CytoTRACE_differentiation_pa_clusters.png"
ctk.plot_projection(basis="umap", color="iter_cluster_id_with_paneth", 
                    legend_loc="right", save=differentiation_figure, show=False)


# d.5) Check terminal states from annotations
sc.pl.embedding(adata, basis="umap", color="terminal_states", add_outline=True)
#####################


### d) Connectivity kernel
# d.1) Setup kernel 
print("Working with Connectivity kernel")
from cellrank.kernels import ConnectivityKernel 

# d.2) compute transition matrix
ck = ConnectivityKernel(adata).compute_transition_matrix()

# d.3) Simulate a random walk on the Markov chain implied by the transition matrix 
ck_rw_figure = figures_dir + "Connectivity_random_walk.pdf"
ck.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs=None,
    legend_loc="right",
    dpi=100,
    save=ck_rw_figure
)

# d.4) Check terminal states from annotations
sc.pl.embedding(adata, basis="umap", color="terminal_states", add_outline=True)



###################


### e) Advanced estimator-based analysis of the CytoTRACE kernel with estimators 
from cellrank.estimators import GPCCA
from cellrank.estimators import CFLARE

# e.1) compute estimator
g_ctk = cr.estimators.GPCCA(ctk)
g_ctk.compute_schur(n_components=100) # compute Schur decomposition
eig_figure = figures_dir_CytoTRACE + "GPCCA_eig.pdf"
g_ctk.plot_spectrum(real_only=True, n=None, show_all_xticks=False, save=eig_figure, figsize=(20,5))

g_ctk.compute_macrostates(n_states=7, cluster_key="iter_cluster_id_with_paneth")
macrostates_figure = figures_dir_CytoTRACE + "GPCCA_macrostates.pdf"
g_ctk.plot_macrostates(which="all", legend_loc="right", s=100, save=macrostates_figure, show=False)

macrostates_figure_composition = figures_dir_CytoTRACE + "GPCCA_macrostates_composition.pdf"
g_ctk.plot_macrostate_composition(key="iter_cluster_id_with_paneth", figsize=(8.5,4), show=False, save=macrostates_figure_composition) # composition of each macrostate

coarse_T_figure = figures_dir_CytoTRACE + "GPCCA_coarse_T.pdf"
g_ctk.plot_coarse_T(annotate=True, save=coarse_T_figure) # plot transition matrix

# e.2) compute terminal states
g_ctk.predict_terminal_states()
terminal_states_figure = figures_dir_CytoTRACE + "GPCCA_terminal_states.pdf"
g_ctk.plot_macrostates(which="terminal", legend_loc="right", s=100, save=terminal_states_figure, show=False)

# e.3) predict initial states
g_ctk.predict_initial_states(allow_overlap=True)
g_ctk.plot_macrostates(which="initial", s=100)

# e.4) verify the initial score by aggregating across cytotrace markers

sc.tl.score_genes(
    adata, gene_list=cytotrace_markers, score_name="initial_score", use_raw=False
) # compute a score in scanpy by aggregating across a few cytotrace markers
# collect macrostates to AnnData
adata.obs["macrostates"] = g_ctk.macrostates 
adata.uns["macrostates_colors"] = g_ctk.macrostates_memberships.colors
# visualize macrostates via initial_scores
sc.pl.violin(adata, keys="initial_score", groupby="macrostates", rotation=90)


#########################################################
# f) Estimating Fate Probabilities 

# f.1) compute fate probabilities towards identified terminal states
g_ctk.compute_fate_probabilities()

# f.2) visualize fate probabilities on UMAP
fate_probability_figure = figures_dir_CytoTRACE + "Fate_probability.pdf"
g_ctk.plot_fate_probabilities(same_plot=False, vmax=1,show=False, save=fate_probability_figure)

# f.3) visualize fate probabilities on a circular projection
fate_circular_figure = figures_dir_CytoTRACE + "Fate_circular.pdf"
cr.pl.circular_projection(adata, keys="iter_cluster_id_with_paneth", legend_loc="right", save=fate_circular_figure)

# f.4) aggregate fate probabilities and visualize how they are committed towards selected cell types   
progenitor_states = ["stem-1", "stem-2"] # select stem-1 and stem-2 as the progenitor states to aggregate their probabilities
fate_committed_vln_figure = figures_dir_CytoTRACE + "Fate_committed_vln.pdf"
cr.pl.aggregate_fate_probabilities(
    adata,
    mode="violin",
    lineages=["absorbitive-2", "stem-1_2", "secretory-1", "secretory-2", "z_paneth"], # "secretory-3"
    cluster_key="iter_cluster_id_with_paneth",
    clusters=progenitor_states,
    save=fate_committed_vln_figure
)


#########################################################
# g) Uncover driver genes

sc.pl.violin(adata, keys=["ct_pseudotime"], groupby="iter_cluster_id_with_paneth", 
             order=["stem-2","stem-1","absorbitive-1","absorbitive-2","secretory-1","secretory-2","secretory-3","z_paneth"], rotation=90)

# g.1) Fit GAM model
GAM_model = cr.models.GAMR(adata, n_knots=6, smoothing_penalty=10.0)

gene_trends_figure = figures_dir_CytoTRACE + "gene_trends.pdf"
# g.2) Plot gene trends 
cr.pl.gene_trends(
    adata,
    model=GAM_model,
    #data_key="magic_imputed_data",
    #genes=["Lgr4","Lgr5","Smarca5","Rbbp7","Tcerg1","Hnrnpd","Hmg20b","Nelfe","Ube2i","Etv5","Ubn1","Mbd3"],
    genes=cytotrace_markers + ["Id3","Hnf4g","Atoh1","Spdef","Neurod1","Lgr4", "Lgr5", "Atad2"],
    same_plot=True,
    ncols=4,
    time_key="ct_pseudotime",
    figsize=(20,12),
    hide_cells=True,
    legend_loc=None,
    #return_figure=True,
    weight_threshold=(1e-3, 1e-3),
    save=gene_trends_figure
)


# g.3) compute putative driver genes for the stem-1, stem-2 trajectories
drivers = g_ctk.compute_lineage_drivers(lineages=["secretory-1", "secretory-2", "secretory-3"])

# g.3.1) plot heatmap for secretory-1
gene_cascade_figure = figures_dir_CytoTRACE + "gene_cascade_secretory-1.pdf"
cr.pl.heatmap(
    adata,
    model=GAM_model,  # use the GAM model from before
    lineages="secretory-1",
    cluster_key="iter_cluster_id_with_paneth",
    show_fate_probabilities=True,
    genes=drivers.head(40).index,
    time_key="ct_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save=gene_cascade_figure
)

# g.3.2) plot heatmap for secretory-2
gene_cascade_figure = figures_dir_CytoTRACE + "gene_cascade_secretory-2.pdf"
cr.pl.heatmap(
    adata,
    model=GAM_model,  # use the GAM model from before
    lineages="secretory-2",
    cluster_key="iter_cluster_id_with_paneth",
    show_fate_probabilities=True,
    genes=drivers.head(40).index,
    time_key="ct_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save=gene_cascade_figure
)

# g.3.2) plot heatmap for secretory-3
gene_cascade_figure = figures_dir_CytoTRACE + "gene_cascade_secretory-3.pdf"
cr.pl.heatmap(
    adata,
    model=GAM_model,  # use the GAM model from before
    lineages="secretory-3",
    cluster_key="iter_cluster_id_with_paneth",
    show_fate_probabilities=True,
    genes=drivers.head(40).index,
    time_key="ct_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save=gene_cascade_figure
)

#driver_clusters = ["stem-1", "stem-2", "absorbitive-1", "absorbitive-2", "secretory-2", "secretory-3"] # "secretory-1"

## g.1) compute driver genes for all trajectories
#driver_df = g_ctk.compute_lineage_drivers(cluster_key="iter_cluster_id_with_paneth", clusters=driver_clusters, use_raw=False) # compute driver genes
#annotated_genes = ['Lgr4','Lgr5'] # define set of genes to annotate

#genes_oi = {
 #   "annotated_genes": annotated_genes
#}

## make sure all of these exist in AnnData
#assert [
#    gene in adata.var_names for genes in genes_oi.values() for gene in genes
#], "Did not find all genes"

#adata.var["mean expression"] = adata.X.A.mean(axis=0) # compute mean gene expression across all cells

## visualize in a scatter plot
#g_ctk.plot_lineage_drivers_correlation(
 #   lineage_x="stem-2",
 #   lineage_y="absorbitive-2",
 #   adjust_text=True,
 #   gene_sets=genes_oi,
 #   color="mean expression",
 #   legend_loc="none",
 #   figsize=(5, 5),
 #   dpi=150,
 #   fontsize=9,
 #   size=50,
 #   use_raw=False
#)


#sc.pl.embedding(
#    adata,
#    basis="umap",
#    color=["iter_cluster_id_with_paneth", "ct_pseudotime"],
#    color_map="gnuplot",
#    legend_loc="on data",
#)







#############################################

### f) Analyze Connectivity kernel with estimators
# # f.1) Setup estimator 
print("Analysing Connectivity kernel with GPCCA")
gpcca_ck = GPCCA(ck)

# f.2) Fit the estimator to compute macrostates of cellular dynamics
gpcca_ck.fit(n_states=8, cluster_key="seurat_clusters")
gpcca_ck.plot_macrostates(which="all")

# f.3) Predict initial states
gpcca_ck.predict_initial_states() 
gpcca_ck.plot_macrostates(which="initial")

# f.4) Predict terminal states
gpcca_ck.predict_terminal_states(method="top_n", n_states=3)
gpcca_ck.plot_macrostates(which="terminal")

# f.5) Compute and plot fate probabilities 
gpcca_ck.compute_fate_probabilities()
gpcca_ck.plot_fate_probabilities(legend_loc="right")

# f.6) Compute 
cr.pl.circular_projection(adata, keys="seurat_clusters", legend_loc="right")

# f.6) Inferring putative driver genes for any of these trajectories
#mono_drivers = g.compute_lineage_drivers(lineages="Mono_1_1")
#mono_drivers.head(10)

# f.7) Visualizing expression trends
model_ck = cr.models.GAMR(adata)


##############################################################################








cr.logging.print_versions()




##### Standard way to setup and run the estimator
# e.1) Setup estimator
#print("Analysing CytoTRACE kernel with GPCCA")
#gpcca_ctk = GPCCA(ctk)

# e.2) Fit the estimator to compute macrostates of cellular dynamics
#gpcca_ctk.fit(n_states=[2, 20], cluster_key="iter_cluster_id_with_paneth") # find macrostates in the interval [2, 20]
#gpcca_ctk.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100) # shown n_cells most associated with each macrostate


# e.3) Predict initial states
#gpcca_ctk.predict_initial_states() 
#gpcca_ctk.plot_macrostates(which="initial", legend_loc="right", s=100)

# e.4) Predict terminal states
#gpcca_ctk.predict_terminal_states()
#gpcca_ctk.plot_macrostates(which="terminal", legend_loc="right", s=10q0)
#gpcca_ctk.plot_macrostates(which="terminal", discrete=False)

# e.5) Plot coarse grain transition matrix 
#ctk_tm_figure = figures_dir + "ctk_transition_matrix.pdf"
#gpcca_ctk.plot_coarse_T(save=ctk_tm_figure)


# e.6) Compute and plot fate probabilities
#gpcca_ctk.compute_fate_probabilities()
#gpcca_ctk.plot_fate_probabilities(legend_loc="right")

#cr.pl.circular_projection(adata, keys="seurat_clusters", legend_loc="right")


# e.6) Inferring putative driver genes for any of these trajectories
#mono_drivers = g.compute_lineage_drivers(lineages="Mono_1_1")
#mono_drivers.head(10)

# e.7) Visualizing expression trends
#model_ctk = cr.models.GAMR(adata)

#########################