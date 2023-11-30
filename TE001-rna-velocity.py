python3

#########################################
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
import scvelo as scv
import warnings
warnings.simplefilter("ignore", category=UserWarning)

# scvelo settings
scv.set_figure_params() # to visualize RNA velocity plots
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# a.3) load .loom file with RNA velocity annalysis in it
counts_loom = h5ad_path + "TE001.loom"

ldata = scv.read(counts_loom, cache=True)

# a.4) load counts data (exported from Seurat@RNA assay)
counts_h5ad = h5ad_path + "TE001-counts.h5ad"
counts = scv.read(counts_h5ad) # load in object


print("counts contains the counts for the TE001 dataset")
#counts = adata.raw.to_adata() # the adata already contains the counts matrix

# a.5) load protein activity data (exported from Seurat@VIPER assay)
pa_h5ad = h5ad_path + "TE001-subnetworks-one-signature-seurat-viper-analysis-with-metacell-data-with-paneth.h5ad"
pa = scv.read(pa_h5ad) # load in object


# a.6) merge loom file with protein activity adata object
adata = scv.utils.merge(counts,ldata)


# a.7) load metadata for TE001 
metadata_csv = h5ad_path + "TE001-metadata-ingest-with-cluster-ids-with-paneth-cluster.csv"
metadata = pd.read_csv(metadata_csv)

# a.8) process metadata in adata
specified_columns = ["cell_id", "nCount_RNA", "nFeature_RNA", "mt_percent", "cytotrace_score.ges",
                     "cytotrace_gcs.ges", "S.Score", "G2M.Score", "Phase", "seurat_clusters", "singleR_labels", "stemness_index.ges",
                     "sample_batch", "TotalUMIs","initial_size_unspliced", "initial_size_spliced", "initial_size"]

adata.obs = adata.obs[specified_columns]
adata.obs = pd.merge(adata.obs, metadata, on='cell_id', how='left') # merge metadata and include into counts object
 
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category') # clusters as categorical variable

# a.7) set UMAP coordinates to those obtained at protein activty
umap_coordinates = np.array(adata.obs.loc[:, ['UMAP_1','UMAP_2','UMAP_3']]) 

adata.obsm['X_umap'] = umap_coordinates



# a.7) show proportions of spliced/unspliced counts
scv.pl.proportions(adata)


###########################################

# b.1) Data preprocessing
# preprocessing consists of gene selection by detection (with a minimum number of counts)
# and highly variability (dispersion), normalizing every cell by its total size and 
# logarithmizing X. Filtering and normalization is applied to spliced and unspliced counts
# and X. Logarithmizing is only applied to X. If X is already preprocessed from former analysis,
# it will not be touched.

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
sc.pp.log1p(adata)

# b.2) Compute first (means) and second order moments (uncentered variances)
# means: needed for deterministic velocity estimation
# variances: needed for stochastic velocity estimation 

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


###########################

# c) Estimating RNA velocity 
# c.1) Computing RNA velocities 
scv.tl.velocity(adata) 

# c.2) Compute transition probabilities 
from comm import create_comm
scv.tl.velocity_graph(adata) # generate velocity_graph and add it to layer 
#scv.utils.get_transition_matrix(adata) # to get the transition matrix



# c.3) Project the velocities onto the embedding
# --as streamlines
velocity_streamline_fig = figures_dir + "velocity_streamline.pdf"
scv.pl.velocity_embedding_stream(adata,basis='umap', color='seurat_clusters', 
                                 save=velocity_streamline_fig, show=False)

# --as fine grain velocity embedding
velocity_fine_grain_fig = figures_dir + "velocity_fine_grain.pdf"
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=1.5, 
                          dpi=120, color='seurat_clusters',
                          save=velocity_fine_grain_fig, show=False)





# GET FULL DYNAMICAL MODEL 