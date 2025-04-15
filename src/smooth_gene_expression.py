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
