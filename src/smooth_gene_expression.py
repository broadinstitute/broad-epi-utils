import numpy as np
from scipy.sparse import csr_matrix

# 1. Get precomputed distances (sparse matrix)
D = adata.obsp["distances"]  # shape: (n_cells, n_cells)

# 2. Turn distances into similarity weights using a Gaussian kernel
#    (small distance = higher weight)
sigma = np.mean(D.data)
W = D.copy()
W.data = np.exp(-W.data**2 / (2 * sigma**2))  # apply Gaussian kernel

# 3. Row-normalize the weights so they sum to 1
row_sums = np.array(W.sum(axis=1)).flatten()
W_norm = W.multiply(1 / row_sums[:, np.newaxis])  # shape: (n_cells, n_cells)


# If working with just the marker genes

# 1. Get the gene expression (as a flat array)
gene = "" # Put marker gene name here
x = adata[:, gene].X
x = x.toarray().flatten() if hasattr(x, "toarray") else x.flatten()

# 2. Compute smoothed gene expression
x_smooth = W_norm.dot(x)  # shape: (n_cells,)

# 3. Store result
adata.obs[f"{gene}_smoothed"] = x_smooth


# If working with all genes

# 1. Get the gene expression (as a sparse matrix)
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

# 2. Compute smoothed gene expression
X_smoothed = W_norm @ X

# 3. Store result
adata.layers["X_smoothed"] = X_smoothed

# 4. Access the smoothed expression for a specific gene
adata[:, gene].layers["X_smoothed"]



# parameterized function to use with scanpy/snapatac2

def smooth_layer_from_knn(adata, genes=None, layer_name="smooth", sigma=None, n_neighbors=None):
    """
    Smooth gene expression using KNN with Gaussian kernel from precomputed distances.

    Parameters
    ----------
    adata : AnnData
        Must have `obsp["distances"]` populated.
    genes : list or None
        List of gene names to smooth. If None, smooth all genes.
    layer_name : str
        Name of the new layer to store smoothed matrix in AnnData.
    sigma : float or None
        Bandwidth for Gaussian kernel. If None, uses mean of non-zero distances.
    n_neighbors : int or None
        Number of neighbors to keep. If None, use all neighbors. This takes a while to compute for all neighbors - suggestion is take to closest 50.
    """

    D = adata.obsp["distances"].tocsr() 

    if sigma is None:
        sigma = np.mean(D.data)

    W = D.copy()
    W.data = np.exp(-W.data**2 / (2 * sigma**2))

    # keep only top n_neighbors if specified
    if n_neighbors is not None:
        for i in range(W.shape[0]): #for every row, take the first and last non-zero value, sort, set sorted neighbors distance past n_neighbors to 0, and then remove them
            row_start = W.indptr[i]
            row_end = W.indptr[i+1]
            row_data = W.data[row_start:row_end]
            if len(row_data) > n_neighbors:
                idx = np.argsort(row_data)[-n_neighbors:]
                mask = np.zeros_like(row_data, dtype=bool)
                mask[idx] = True
                W.data[row_start:row_end][~mask] = 0
        W.eliminate_zeros()

    row_sums = np.array(W.sum(axis=1)).flatten()
    W = W.multiply(1 / row_sums[:, np.newaxis])

    X = adata[:, genes].X if genes is not None else adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    #Compute smoothed gene expression
    X_smooth = W @ X

    # Save as new layer in AnnData object
    adata.layers[layer_name] = X_smooth