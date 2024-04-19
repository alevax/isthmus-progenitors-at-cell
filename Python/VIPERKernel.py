import numpy as np
from cellrank.kernels import PrecomputedKernel
from cellrank.kernels.mixins import ConnectivityMixin


from tqdm import tqdm

def _get_MRs_ranked_array(vpmat):
    # For each sample, we calculate index arragement that would sort the vector
    sorted_order_array = np.argsort(-vpmat.values)
    # We then get the MRs ranked for each sample by indexing with this sorted order
    mrs_ranked_array = np.array(vpmat.columns)[sorted_order_array]
    return mrs_ranked_array

def _get_VIPER_transition_mat_with_mrs_knn(vpmat, mrs_ranked_array, knn_array, n_mrs = 100, symmetric_mrs = True):
    transition_matrix = np.zeros(knn_array.shape)
    n_samples = knn_array.shape[0]

    # If using top and bottom MRs
    if symmetric_mrs:
        n_mrs_half = int(n_mrs/2)
        mr_indices = np.concatenate([np.arange(n_mrs_half), np.arange(-n_mrs_half, 0)])
    # If using only top MRs
    else:
        mr_indices = np.arange(n_mrs)

    # Compute the probabilities for each sample
    for i in tqdm(range(n_samples)):
        # Get the KNNs and top MRs of sample i
        knn_i = np.nonzero(knn_array[i,:])[0]
        top_mrs_i = mrs_ranked_array[i, mr_indices]

        # Compute candidate MR scores for each KNN
        # Take the difference in activity between each KNN and sample i
        vpmat_i = vpmat.loc[vpmat.index[i], top_mrs_i]
        vpmat_knn = vpmat.loc[vpmat.index[knn_i], top_mrs_i]
        diffs = vpmat_knn.values - vpmat_i.values

        # Correct the difference of bottom MRs by multiplying by -1 (their sign)
        diffs = diffs * np.sign(vpmat_i.values)

        # Sum diff values across all MRs to get a score measuring a
        # strengthening in the MRs defining sample i's cellular identity
        scores = np.sum(diffs, axis=1)

        # Normalize the scores into probability values
        probs = scores + abs(np.min(scores))
        probs = probs/np.sum(probs)

        # Each sample gives us a row of the matrix
        transition_matrix[i, knn_i] = probs
    return transition_matrix

def _get_VIPER_transition_mat(adata, n_mrs = 100, symmetric_mrs = True, conn_key = 'connectivities'):
    vpmat = adata.to_df().copy()
    # Get an array of the top MRs for each sample
    mrs_ranked_array = _get_MRs_ranked_array(vpmat)

    try:
        knn_array = adata.obsp[conn_key].toarray()
    except KeyError:
        raise KeyError('connectivities are missing from adata.obsp["' + str(conn_key) + '"]. Run sc.pp.neighbors on adata to generate them.')

    # Compute the transition matrix
    transition_matrix = _get_VIPER_transition_mat_with_mrs_knn(vpmat, mrs_ranked_array, knn_array, n_mrs, symmetric_mrs)
    return transition_matrix

class VIPERKernel(PrecomputedKernel, ConnectivityMixin):

    def __init__(self, adata, n_mrs = 100, symmetric_mrs = True, conn_key = 'connectivities'):
        """\
        Kernel which computes directed transition probabilities using the
        activity of the top candidate master regulators (MRs) of each cell.
        Just like any other CellRank kernel, this kernel can be used for
        downstream analysis, including initializing an estimator and computing
        initial and terminal states, fate probabilities, and driver genes.

        Parameters
        ----------
        adata
            An anndata.AnnData or pd.DataFrame containing protein activity (NES),
            where rows are observations/samples (e.g. cells or groups) and
            columns are features (e.g. proteins or pathways). A connectivities
            matrix (e.g. generated using sc.pp.neighbors) must be located in
            adata.obsp[conn_key].
        n_mrs (default: 100)
            Number of top MRs from each sample to use to compute probabilities
            of it transitioning to each of its nearest-neighbors.
        symmetric_mrs (default: True)
            Whether to use n_mrs/2 top MRs and n_mrs/2 bottom MRs instead of
            using only n_mrs top MRs.
        conn_key (default: 'connectivities')
            Where the connectivitiy matrix is located within adata.obsp.

        References
        ----------
        Lange, M., Bergen, V., Klein, M., Setty, M., Reuter, B., Bakhti, M., ... &
        Theis, F. J. (2022). CellRank for directed single-cell fate mapping. Nature
        methods, 19(2), 159-170.
        """

        # Compute transition matrix
        transition_matrix = _get_VIPER_transition_mat(adata, n_mrs, symmetric_mrs, conn_key)

        # Assign connectivity characteristics to the kernel
        self._conn_key = conn_key
        self._conn = adata.obsp[self._conn_key]

        # Setup kernel using CellRank's PrecomputedKernel
        super().__init__(transition_matrix, adata = adata.copy())
