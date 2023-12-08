python3

### a) Import packages and data
# a.1) setup path to data-containing folder and savings and parameters
h5ad_path = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE001-h5ad/"
figures_dir = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/"

n_macro_CytoTRACE = 7 # number of macrostates 
n_macro_Connectivity= 7 # number of macrostates
# a.2) packages
import sys
import cellrank as cr
import scanpy as sc
import numpy as np
import pandas as pd
import os
from matplotlib import rc_context
import matplotlib.pyplot as plt
import scvelo as scv

sc.settings.set_figure_params(frameon=False, dpi=100)
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

import warnings
warnings.simplefilter("ignore", category=UserWarning)

cytotrace_markers = ['Smarca5','Rbbp7','Tcerg1','Hnrnpd','Hmg20b','Nelfe','Ube2i','Etv5','Ubn1','Mbd3','Dek','Maz',
                     'Itgb3bp','Ilf2','Pa2g4'] # Id3','Hnf4g','Atoh1','Spdef','Neurod1' markers upregulated in cytotrace (Fig 1e) 

# a.3) load counts data (exported from Seurat@RNA assay)
counts_h5ad = h5ad_path + "TE001-counts.h5ad"
adata = sc.read_h5ad(counts_h5ad) # load in object

# a.4) load .loom file with RNA velocity annalysis in it
counts_loom = h5ad_path + "TE001.loom"
ldata = scv.read(counts_loom, cache=True)

# a.5) merge loom file with gene expression adata object
adata = scv.utils.merge(adata,ldata)

# a.6) load metadata for TE001 
metadata_csv = h5ad_path + "TE001-metadata-umap-and-clusters-for-paper.csv"
metadata = pd.read_csv(metadata_csv)

# a.7) process metadata in adata
specified_columns = ["cell_id", "nCount_RNA", "nFeature_RNA", "mt_percent", "cytotrace_score.ges",
                     "cytotrace_gcs.ges", "S.Score", "G2M.Score", "Phase", "seurat_clusters", "singleR_labels", "stemness_index.ges"]

adata.obs = adata.obs[specified_columns]
cells_to_analyze = metadata['cell_id'] # cells to analyze
adata = adata[adata.obs['cell_id'].isin(cells_to_analyze)] # subset cells to analyze in adata

adata.obs = pd.merge(adata.obs, metadata, on='cell_id', how='left') # merge metadata and include into counts object

adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category') # clusters as categorical variable

# a.8) set UMAP coordinates to those obtained at protein activty
umap_coordinates = np.array(adata.obs.loc[:, ['UMAP_1_scanpy','UMAP_2_scanpy']]) 
adata.obsm['X_umap'] = umap_coordinates


# a.9) Include metadata of terminal states for CellRank analysis
adata.obs['terminal_states'] = adata.obs['iter_cluster_id_with_paneth']
adata.obs['terminal_states'].iloc[adata.obs['terminal_states'].isin(["stem-1","stem-2"])] = np.nan

print("adata contains the counts and RNA velocity for the TE001 dataset")
#counts = adata.raw.to_adata() # the adata already contains the counts matrix

# a.8) display specific marker genes
log_expression = adata.copy()
sc.pp.normalize_total(log_expression, target_sum=1e4)
sc.pp.log1p(log_expression)

sc.pl.umap(log_expression,color=["Lgr4","Lgr5"], use_raw=False, cmap='viridis',add_outline=True)


#####################


### b) Preprocess the data and UMAP visualization
scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
)
sc.tl.pca(adata, random_state=0)
sc.pp.neighbors(adata,  n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

### c) Run scVelo's dynamical model to estimate model parameters
# c.1) get the parameters
from comm import create_comm
scv.tl.recover_dynamics(adata, n_jobs=8)

# c.2) compute the actual velocities
scv.tl.velocity(adata, mode="dynamical")


##########################################################################################################
##########################################################################################################
### CellRank2 analysis with VelocityKernel
##########################################################################################################
##########################################################################################################

### c) Velocity kernel
# c.1) Setup kernel
print("Working with Velocity kernel")
# The VelocityKernel computes transition probabilities Tij to each cell j in the neighborhood of i,
# by quantifying how much the velocity vector vi of cell i points towards each of its nearest 
# neighbors.

vk = cr.kernels.VelocityKernel(adata) # initialize the CellRank2 kernel

# c.2) compute transition probability
vk.compute_transition_matrix(model="stochastic") # compute transition matrix

figures_dir_Velocity = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/CR2_VelocityKernel/"

if os.path.exists(figures_dir_Velocity):
    print("'figures_dir_Velocity' directory already exists")
else:
    os.mkdir(figures_dir_Velocity)

# c.3) Simulate a random walk on the Markov chain implied by the transition matrix 
vk_rw_figure = figures_dir_Velocity + "Velocity_random_walk.pdf"
vk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs=None,
    legend_loc="right",
    dpi=100,
    save=vk_rw_figure,
    figsize=(3,3)
)


vk_rw_figure = figures_dir_Velocity + "Velocity_random_walk_stem-1.pdf"
vk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={'iter_cluster_id_with_paneth': ["stem-1"]},
    legend_loc="right",
    dpi=100,
    save=vk_rw_figure
)

# sampling cells randomly from stem-2 population
vk_rw_figure = figures_dir_Velocity + "Velocity_random_walk_stem-2.pdf"
vk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={'iter_cluster_id_with_paneth': ["stem-2"]},
    legend_loc="right",
    dpi=100,
    save=vk_rw_figure
)


# c.4) visualize the transition matrix
differentiation_figure = figures_dir_Velocity + "Velocity_differentiation_ges_clusters.png"
vk.plot_projection(basis="umap", color="seurat_clusters", 
                    legend_loc="right", save=differentiation_figure, show=False)


differentiation_figure = figures_dir_Velocity + "Velocity_differentiation_pa_clusters.png"
vk.plot_projection(basis="umap", color="iter_cluster_id_with_paneth", 
                    legend_loc="right", save=differentiation_figure, show=False)


# c.5) Check terminal states from annotations
sc.pl.embedding(adata, basis="umap", color="terminal_states", add_outline=True)
#####################


###################


### d) Advanced estimator-based analysis of the CytoTRACE kernel with estimators 
from cellrank.estimators import GPCCA
from cellrank.estimators import CFLARE

# d.1) compute estimator
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

# d.2) compute terminal states
g_ctk.predict_terminal_states()
terminal_states_figure = figures_dir_CytoTRACE + "GPCCA_terminal_states.pdf"
g_ctk.plot_macrostates(which="terminal", legend_loc="right", s=100, save=terminal_states_figure, show=False)

# d.3) predict initial states
g_ctk.predict_initial_states(allow_overlap=True)
g_ctk.plot_macrostates(which="initial", s=100)

# d.4) verify the initial score by aggregating across cytotrace markers

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



##############################################################################



##########################################################################################################
##########################################################################################################
### CellRank2 analysis with Connectivity Kernel
##########################################################################################################
##########################################################################################################



### h) Connectivity kernel
# h.1) Setup kernel 
print("Working with Connectivity kernel")
from cellrank.kernels import ConnectivityKernel 


figures_dir_Connectivity = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/CR2_ConnectivityKernel/" # where to save figures for connectivity kernel
if os.path.exists(figures_dir_Connectivity):
    print("'figures_dir_Connectivity' directory already exists")
else:
    os.mkdir(figures_dir_Connectivity)

# h.2) compute transition matrix
ck = ConnectivityKernel(adata).compute_transition_matrix()


# h.3) Simulate a random walk on the Markov chain implied by the transition matrix 
# sampling cells randomly among all clusters
ck_rw_figure = figures_dir_Connectivity + "Connectivity_random_walk.pdf"
ck.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs=None,
    legend_loc="right",
    dpi=100,
    save=ck_rw_figure
)

# h.4) visualize the transition matrix
# sampling cells randomly from stem-1 population
ck_rw_figure = figures_dir_Connectivity + "Connectivity_random_walk_stem-1.pdf"
ck.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={'iter_cluster_id_with_paneth': ["stem-1"]},
    legend_loc="right",
    dpi=100,
    save=ck_rw_figure
)

# sampling cells randomly from stem-2 population
ck_rw_figure = figures_dir_Connectivity + "Connectivity_random_walk_stem-2.pdf"
ck.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={'iter_cluster_id_with_paneth': ["stem-2"]},
    legend_loc="right",
    dpi=100,
    save=ck_rw_figure
)

# h.4) visualize the transition matrix
differentiation_figure = figures_dir_Connectivity + "Connectivity_differentiation_ges_clusters.png"
ck.plot_projection(basis="umap", color="seurat_clusters", 
                    legend_loc="right", save=differentiation_figure, show=False)


differentiation_figure = figures_dir_Connectivity + "Connectivity_differentiation_pa_clusters.png"
ck.plot_projection(basis="umap", color="iter_cluster_id_with_paneth", 
                    legend_loc="right", save=differentiation_figure, show=False)



# h.5) Check terminal states from annotations
sc.pl.embedding(adata, basis="umap", color="terminal_states", add_outline=True)






#############################################

### i) Analyze Connectivity kernel with estimators 
# # i.1) Setup estimator 
from cellrank.estimators import GPCCA
from cellrank.estimators import CFLARE
print("Analysing Connectivity kernel with GPCCA")

# i.1) compute estimator
g_ck = GPCCA(ck)
g_ck.compute_schur(n_components=100) # compute Schur decomposition
eig_figure = figures_dir_Connectivity + "GPCCA_eig.pdf"
g_ck.plot_spectrum(real_only=True, n=None, show_all_xticks=False, save=eig_figure, figsize=(20,5))

g_ck.compute_macrostates(n_states=n_macro_Connectivity, cluster_key="iter_cluster_id_with_paneth")
macrostates_figure = figures_dir_Connectivity + "GPCCA_macrostates.pdf"
g_ck.plot_macrostates(which="all", legend_loc="right", s=100, save=macrostates_figure, show=False)

macrostates_figure_composition = figures_dir_Connectivity + "GPCCA_macrostates_composition.pdf"
g_ck.plot_macrostate_composition(key="iter_cluster_id_with_paneth", figsize=(8.5,4), show=False, save=macrostates_figure_composition) # composition of each macrostate


coarse_T_figure = figures_dir_Connectivity + "GPCCA_coarse_T.pdf"
g_ck.plot_coarse_T(annotate=True, save=coarse_T_figure) # plot transition matrix


# i.2) compute terminal states
g_ck.predict_terminal_states()
terminal_states_figure = figures_dir_Connectivity + "GPCCA_terminal_states.pdf"
g_ck.plot_macrostates(which="terminal", legend_loc="right", s=100, save=terminal_states_figure, show=False)


# i.3) predict initial states
g_ck.predict_initial_states(allow_overlap=True)
g_ck.plot_macrostates(which="initial", s=100)



#########################################################
# j) Estimating Fate Probabilities for Connectivity kernel

# j.1) compute fate probabilities towards identified terminal states
g_ck.compute_fate_probabilities()

# j.2) visualize fate probabilities on UMAP
fate_probability_figure = figures_dir_Connectivity + "Fate_probability.pdf"
g_ck.plot_fate_probabilities(same_plot=False, vmax=1,show=False, save=fate_probability_figure)

# j.3) visualize fate probabilities on a circular projection
fate_circular_figure = figures_dir_Connectivity + "Fate_circular.pdf"
cr.pl.circular_projection(adata, keys="iter_cluster_id_with_paneth", legend_loc="right", save=fate_circular_figure)

# j.4) aggregate fate probabilities and visualize how they are committed towards selected cell types   
progenitor_states = ["stem-1", "stem-2"] # select stem-1 and stem-2 as the progenitor states to aggregate their probabilities
fate_committed_vln_figure = figures_dir_Connectivity + "Fate_committed_vln.pdf"
ck_terminal_states = g_ck.terminal_states.unique().dropna().to_list()
cr.pl.aggregate_fate_probabilities(
    adata,
    mode="violin",
    lineages=ck_terminal_states, # "secretory-3"
    cluster_key="iter_cluster_id_with_paneth",
    clusters=progenitor_states,
    save=fate_committed_vln_figure
)



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