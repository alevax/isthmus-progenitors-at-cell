### a) Import packages and data
# Setup path to data-containing folder and savings and parameters

# a.1) setup path to data-containing folder and savings and parameters
h5ad_path = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/TE005/"
figures_dir = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/"
figures_dir_CytoTRACE = "/Volumes/ac_lab_scratch/lz2841/ics-rebuttal/figures/CR2_CytoTRACEKernel/TE005/"

n_macro_CytoTRACE = 8 # number of macrostates 


#Import packages and set markers of interest
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


#Load counts data and metadata for TE006. Setup UMAP coordinates for UMAP visualization 
# a.3) load counts data 
counts_h5ad = h5ad_path + "TE005-cells-original-data.h5ad"
adata = sc.read_h5ad(counts_h5ad) # load in object

# a.4) load metadata for TE005 
metadata_csv = h5ad_path + "TE005-metadata-ingest-with-cluster-ids-with-paneth-cluster.csv"
metadata = pd.read_csv(metadata_csv)
metadata = metadata.rename(columns={'index': 'cell_id'})

metadata_csv_umap = h5ad_path + "haplo-metadata-ingest.csv" 
metadata_umap = pd.read_csv(metadata_csv_umap)
metadata_umap = metadata_umap[['cell_id','UMAP_1','UMAP_2']]

# a.5) process metadata in adata
cells_to_analyze = metadata['cell_id'] # cells to analyze
adata = adata[adata.obs_names.isin(cells_to_analyze)] # subset cells to analyze in adata
adata.obs['cell_id'] = adata.obs_names

adata.obs = pd.merge(adata.obs, metadata, on='cell_id', how='left') # merge metadata and include into counts object
adata.obs = pd.merge(adata.obs, metadata_umap, on='cell_id', how='left').set_index('cell_id') # merge metadata and include into counts object

#adata.obs_names = adata.obs['cell_id']

# a.6) set UMAP coordinates to those obtained at protein activty
umap_coordinates = np.array(adata.obs.loc[:, ['UMAP_1','UMAP_2']]) 
adata.obsm['X_umap'] = umap_coordinates


# a.7) Include metadata of terminal states for CellRank analysis
adata.obs['terminal_states'] = adata.obs['iter_cluster_id_with_paneth']
adata.obs['terminal_states'].iloc[adata.obs['terminal_states'].isin(["stem-1","stem-2"])] = np.nan

print("adata contains the counts for the TE005 dataset")
#counts = adata.raw.to_adata() # the adata already contains the counts matrix

#Show Lgr4 and Lgr5 marker genes
# a.8) display specific marker genes
log_expression = adata.copy()
sc.pp.normalize_total(log_expression, target_sum=1e4)
sc.pp.log1p(log_expression)

print("UMAP showing Lgr4 and Lgr5 expression")
sc.pl.umap(log_expression,color=["Lgr4","Lgr5"], use_raw=False, cmap='viridis',add_outline=True, show=False)


# add log1p-transformed expression to selected genes to the adata  
log_selected_df = log_expression[:,['Lgr4','Lgr5','Atad2','Mki67']].to_df()
log_selected_df = log_selected_df.rename(columns={'Lgr4':'log(Lgr4)', 'Lgr5':'log(Lgr5)', 'Atad2':'log(Atad2)', 'Mki67':'log(Mki67)'})
adata.obs = pd.merge(adata.obs, log_selected_df, left_index=True, right_index=True, how='left') # merge metadata and include into log_selected_df



#Show clusters at gene expression, clusters at protein activity, precomputed cytotrace score (with Cytotrace, not CellRank) and stemness index
sc.pl.umap(adata, color=["iter_cluster_id_with_paneth","cytotrace","stemness_index"], add_outline=True, show=False)

#Preprocess data for CellRank2 analysis
### b) Preprocess the data 
print("Preprocessing counts matrix for CellRank 2 analysis")
sc.pp.filter_genes(adata, min_cells=10)
sc.tl.pca(adata, random_state=0)
sc.pp.neighbors(adata, random_state=0)

### CellRank2 analysis with CytoTRACE Kernel
# Import CellRank2 kernel, compute CytoTRACE score with CellRank (sparse), and compute transition matrix for CellRank2 analysis


##########################################################################################################
##########################################################################################################
### CellRank2 analysis with CytoTRACE Kernel
##########################################################################################################
##########################################################################################################

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

# Initialize folder collecting CytoTRACE figures
if os.path.exists(figures_dir_CytoTRACE):
    print("'figures_dir_CytoTRACE' directory already exists")
else:
    os.mkdir(figures_dir_CytoTRACE)

ctk_pseudotime_figure = figures_dir_CytoTRACE + "CytoTRACE_pseudotime.pdf"
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['ct_score', 'ct_pseudotime'], show=False, add_outline=True)
    plt.savefig(ctk_pseudotime_figure)

ctk_pseudotime_vln = figures_dir_CytoTRACE + "CytoTRACE_pseudotime_vln.pdf"
with rc_context({'figure.figsize': (8.5, 8.5)}):
    sc.pl.violin(adata, keys=["ct_pseudotime"], groupby="iter_cluster_id_with_paneth", rotation=90, show=False)
    plt.savefig(ctk_pseudotime_vln)

#Simulate random walk on the Markov Chain implied by the transition matrix. Starting cells are selected at random. 
# In the first figure, starting cells are selected randomly from all clusters. In figures 2 and 3 they are randomly sampled from the 'stem-1' and 'stem-2' clusters. 
# 100 trajectories are simulated in each Random Walk. Black dots = cells of departure; yellow dots = cells of arrival.

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

# c.4) Simulate a random walk on the Markov chain implied by the transition matrix 
# sampling cells randomly among all clusters
ctk_rw_figure = figures_dir_CytoTRACE + "CytoTRACE_random_walk_stem-1.pdf"
ctk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={'iter_cluster_id_with_paneth':'stem-1'},
    legend_loc="right",
    dpi=100,
    save=ctk_rw_figure
)

# c.5) visualize the transition matrix
# sampling cells randomly from stem-1 population
ctk_rw_figure = figures_dir_CytoTRACE + "CytoTRACE_random_walk_stem-2.pdf"
ctk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={'iter_cluster_id_with_paneth': "stem-2"},
    legend_loc="right",
    dpi=100,
    save=ctk_rw_figure
)


#Visualize the projected Transition Probability matrix on the UMAP, with clusters colored at gene expression (panel 1) and protein activity (panel 2).
# Also, show what we consider as the most differentiated states in the dataset (panel 3).
# c.6) visualize the transition matrix
differentiation_figure = figures_dir_CytoTRACE + "CytoTRACE_differentiation_pa_clusters.png"
ctk.plot_projection(basis="umap", color="iter_cluster_id_with_paneth", 
                    legend_loc="right", save=differentiation_figure, show=True)


# c.7) Check terminal states from annotations
annotated_terminal_states_figure = figures_dir_CytoTRACE + "annotated_terminal_states.pdf"
sc.pl.embedding(adata, basis="umap", color="terminal_states", add_outline=True, title="Selected (putatively most differentiated) cell clusters") 
#####################


# CellRank 2 (advanced) estimator-baed analysis of the Transition Matrix and calculation of the terminal states. 
# Show the real value of the eigenvalues from Schur's decomposition (Figure 1), cell-type distrbution of the terminal states (Figure 2) and 
# visualization of the coarse-grain transition matrix (Figure 3). 
# To change the number of selected macrostate, change `n_macro_CytoTRACE` in the first cell of this notebook.

### d) Advanced estimator-based analysis of the CytoTRACE kernel with estimators 
from cellrank.estimators import GPCCA
from cellrank.estimators import CFLARE

# d.1) compute estimator
g_ctk = cr.estimators.GPCCA(ctk)
g_ctk.compute_schur(n_components=100) # compute Schur decomposition
eig_figure = figures_dir_CytoTRACE + "GPCCA_eig.pdf"
g_ctk.plot_spectrum(real_only=True, n=None, show_all_xticks=False, save=eig_figure, figsize=(20,5))

g_ctk.compute_macrostates(n_states=n_macro_CytoTRACE, cluster_key="iter_cluster_id_with_paneth")
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

#Collect the inferred terminal states in a list for ease of manipulation.

inferred_terminal_states = g_ctk.terminal_states.unique().dropna().to_list()

# Estimate Fate probabilities towards the identified terminal states


#########################################################
# f) Estimating Fate Probabilities 

# f.1) compute fate probabilities towards identified terminal states
g_ctk.compute_fate_probabilities()

# f.2) visualize fate probabilities on UMAP
fate_probability_figure = figures_dir_CytoTRACE + "Fate_probability.pdf"
g_ctk.plot_fate_probabilities(same_plot=False, vmax=1,show=True, save=fate_probability_figure)

# f.3) visualize fate probabilities on a circular projection
fate_circular_figure = figures_dir_CytoTRACE + "Fate_circular.pdf"
cr.pl.circular_projection(adata, keys="iter_cluster_id_with_paneth", legend_loc="right", save=fate_circular_figure)

fate_circular_figure = figures_dir_CytoTRACE + "Fate_circular_genes.pdf"
cr.pl.circular_projection(adata, keys=["log(Lgr4)", "log(Lgr5)", "log(Atad2)", "log(Mki67)"], cmap = "Oranges", ncols=2, legend_loc="right", save=fate_circular_figure)


# f.4) aggregate fate probabilities and visualize how they are committed towards selected cell types   
progenitor_states = ["stem-1", "stem-2"] # select stem-1 and stem-2 as the progenitor states to aggregate their probabilities
fate_committed_vln_figure = figures_dir_CytoTRACE + "Fate_committed_vln.pdf"
cr.pl.aggregate_fate_probabilities(
    adata,
    mode="violin",
    lineages=inferred_terminal_states,
    cluster_key="iter_cluster_id_with_paneth",
    clusters=progenitor_states,
    save=fate_committed_vln_figure
)


# Show Cytotrace-based pseudotime and compute lineage drivers towards all the terminal states

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
drivers = g_ctk.compute_lineage_drivers(lineages=inferred_terminal_states)


for state in inferred_terminal_states:
    gene_cascade_figure = figures_dir_CytoTRACE + "gene_cascade_" + state + ".pdf"
    # g.3.1) plot heatmap for the given state
    cr.pl.heatmap(
        adata,
        model=GAM_model,  # use the GAM model from before
        lineages=state,
        cluster_key="iter_cluster_id_with_paneth",
        show_fate_probabilities=True,
        genes=drivers.head(40).index,
        time_key="ct_pseudotime",
        figsize=(12, 10),
        show_all_genes=True,
        weight_threshold=(1e-3, 1e-3),
        save=gene_cascade_figure
    )

cr.logging.print_versions()
