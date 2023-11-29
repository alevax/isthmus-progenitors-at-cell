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
from matplotlib import rc_context
import matplotlib.pyplot as plt
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

import warnings
warnings.simplefilter("ignore", category=UserWarning)

# a.3) load counts data (exported from Seurat@RNA assay)
counts_h5ad = h5ad_path + "TE001-counts.h5ad"
adata = sc.read_h5ad(counts_h5ad) # load in object

adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category') # clusters as categorical variable


print("adata contains the counts for the TE001 dataset")
#counts = adata.raw.to_adata() # the adata already contains the counts matrix

#####################


### b) Preprocess the data and UMAP visualization
sc.tl.pca(adata, random_state=0)
sc.pp.neighbors(adata, random_state=0)
sc.pl.umap(adata, color=["seurat_clusters","singleR_labels","cytotrace_score.ges","stemness_index.ges"])

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
ctk = ctk.compute_cytotrace().compute_transition_matrix() # compute transition matrix

ctk_pseudotime_figure = figures_dir + "CytoTRACE_pseudotime.pdf"
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['ct_score', 'ct_pseudotime'], show=False)
    plt.savefig(ctk_pseudotime_figure)

ctk_pseudotime_vln = figures_dir + "CytoTRACE_pseudotime_vln.pdf"
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.violin(adata, keys=["ct_pseudotime"], groupby="seurat_clusters", rotation=90, show=False)
    plt.savefig(ctk_pseudotime_vln)

# c.3) Simulate a random walk on the Markov chain implied by the transition matrix 
ctk_rw_figure = figures_dir + "CytoTRACE_random_walk.pdf"
ctk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs=None,
    legend_loc="right",
    dpi=100,
    save=ctk_rw_figure,
    show=False
)

# c.4) visualize the transition matrix
differentiation_figure = figures_dir + "CytoTRACE_differentiation.pdf"
ctk.plot_projection(basis="umap", color="seurat_clusters", 
                    legend_loc="right", save=differentiation_figure, show=False)

# c.5) 
sc.pl.umap(adata, color="ct_terminal_states")
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


###################


### e) Analyze CytoTRACE kernel with estimators and compute fate probabilities 
from cellrank.estimators import GPCCA
from cellrank.estimators import CFLARE

# e.1) Setup estimator
print("Analysing CytoTRACE kernel with GPCCA")
gpcca_ctk = GPCCA(ctk)

# e.2) Fit the estimator to compute macrostates of cellular dynamics
gpcca_ctk.fit(n_states=[2, 20], cluster_key="seurat_clusters") # find macrostates in the interval [2, 20]
gpcca_ctk.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100) # shown n_cells most associated with each macrostate


# e.3) Predict initial states
gpcca_ctk.predict_initial_states() 
gpcca_ctk.plot_macrostates(which="initial", legend_loc="right", s=100)

# e.4) Predict terminal states
gpcca_ctk.predict_terminal_states()
gpcca_ctk.plot_macrostates(which="terminal", legend_loc="right", s=100)
gpcca_ctk.plot_macrostates(which="terminal", discrete=False)

# e.5) Plot coarse grain transition matrix 
ctk_tm_figure = figures_dir + "ctk_transition_matrix.pdf"
gpcca_ctk.plot_coarse_T(save=ctk_tm_figure)


# e.6) Compute and plot fate probabilities
gpcca_ctk.compute_fate_probabilities()
gpcca_ctk.plot_fate_probabilities(legend_loc="right")

cr.pl.circular_projection(adata, keys="seurat_clusters", legend_loc="right")


# e.6) Inferring putative driver genes for any of these trajectories
#mono_drivers = g.compute_lineage_drivers(lineages="Mono_1_1")
#mono_drivers.head(10)

# e.7) Visualizing expression trends
model_ctk = cr.models.GAMR(adata)

#########################

### f) Analyze Connectivity kernel with estimators
# # f.1) Setup estimator 
print("Analysing Connectivity kernel with GPCCA")
gpcca_ck = GPCCA(ck)

# f.2) Fit the estimator to compute macrostates of cellular dynamics
gpcca_ck.fit(n_states=20, cluster_key="seurat_clusters")
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

cr.pl.circular_projection(adata, keys="seurat_clusters", legend_loc="right")

# f.6) Inferring putative driver genes for any of these trajectories
#mono_drivers = g.compute_lineage_drivers(lineages="Mono_1_1")
#mono_drivers.head(10)

# f.7) Visualizing expression trends
model_ck = cr.models.GAMR(adata)


#################################
# g) Advanced estimator-based analysis of the CytoTRACE kernel
ctk2 = cr.estimators.GPCCA(ctk)
ctk2.compute_schur()
ctk2.plot_spectrum(real_only=True)

ctk2.compute_macrostates(n_states=12, cluster_key="seurat_clusters")
ctk2.plot_macrostates(which="all", legend_loc="right", s=100)
ctk2.plot_macrostate_composition(key="seurat_clusters", figsize=(7,4)) # composition of each macrostate
ctk2.plot_coarse_T(annotate=False) # plot transition matrix







cr.logging.print_versions()
