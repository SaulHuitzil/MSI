"""
Hierarchical clustering utilities
"""

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform


def hierarchical_clustering(cov_matrix, method="average"):
    """
    Perform hierarchical clustering on covariance matrix
    
    Args:
        cov_matrix: Covariance matrix as numpy array
        method: Linkage method ('average', 'complete', 'ward', etc.)
    
    Returns:
        Z: Linkage matrix
        dist: Distance matrix
        corr_matrix: Correlation matrix
    """
    # Convert to correlation matrix
    std_devs = np.sqrt(np.diag(cov_matrix))
    corr_matrix = cov_matrix / np.outer(std_devs, std_devs)
    
    # Define distance: d_ij = 1 - |corr_ij|
    dist = 1 - np.abs(corr_matrix)
    
    # Condensed distance matrix for linkage
    dist_condensed = squareform(dist, checks=False)
    
    # Clustering
    Z = linkage(dist_condensed, method=method)
    
    return Z, dist, corr_matrix


def get_clustered_order(Z, variable_names):
    """
    Get the order of variables from hierarchical clustering
    
    Args:
        Z: Linkage matrix from hierarchical clustering
        variable_names: List of variable names
    
    Returns:
        ordered_labels: List of variable names in clustered order
        ordered_indices: List of indices in clustered order
    """
    dendro = dendrogram(Z, labels=variable_names, no_plot=True)
    ordered_labels = dendro["ivl"]
    ordered_indices = [variable_names.index(label) for label in ordered_labels]
    
    return ordered_labels, ordered_indices
