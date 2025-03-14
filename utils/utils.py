"""
utils.py

This file contains utility functions adapted from the pymodulon repository:
    https://github.com/SBRG/pymodulon

Please refer to the original repository for additional context and licensing details.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import warnings
from tqdm import tqdm
from sklearn.decomposition import FastICA, PCA
from scipy import sparse, stats
from sklearn.cluster import DBSCAN
from sklearn.exceptions import EfficiencyWarning

def DrawOntologyPlot(p_value, generatio, figsize=(7.1, 8.5), save_path=None):
    """
    Draw an ontology plot with circle sizes determined by 'generatio'
    and colors mapped from 'p_value'. The plot is saved if save_path is provided.
    """
    N, M = p_value.shape
    ylabels = p_value.index
    xlabels = p_value.columns

    minimum_circlesize_base = 0.3

    x, y = np.meshgrid(np.arange(M), np.arange(N))
    s = generatio.values + minimum_circlesize_base
    c = p_value

    fig, ax = plt.subplots(figsize=figsize)

    R = s / (((1.0 + minimum_circlesize_base) * 1.4) / 2)
    circles = [plt.Circle((j, i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]

    col = PatchCollection(circles, array=c.to_numpy().flatten(), cmap="plasma")
    ax.add_collection(col)

    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.grid(which='minor')

    fig.colorbar(col, label='p value', shrink=0.4, aspect=10, pad=0.02, anchor=(0.5, 0.1))
    plt.xticks(rotation=90)

    for area in [0.1, 0.5, 1.0]:
        ax.scatter([], [], c='k', alpha=0.7, s=area / (((1.0 + 0.1) * 1.4) / 2), label=str(area))
    plt.legend(scatterpoints=1, labelspacing=1, title='Generatio', loc='upper left', bbox_to_anchor=(1.05, 0.9))
    plt.xlim(-0.5, M - 0.5)
    plt.ylim(N - 0.5, -0.5)

    for i in range(4, len(ax.get_xticks()), 4):
        ax.axvline(x=i - 0.5, color='black', linestyle='--')

    ax.set_aspect('equal')

    if save_path is not None:
        plt.savefig(save_path, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.show()


def perform_ICA(rnaseq_df):
    """
    Perform Independent Component Analysis on the provided RNA-seq dataframe.
    """
    # Determine the number of components with PCA
    pca = PCA().fit(rnaseq_df.transpose())
    pca_var = np.cumsum(pca.explained_variance_ratio_)
    k_comp = np.where(pca_var > 0.99)[0][0] + 1

    n_genes, m_samples = rnaseq_df.shape
    print("Data: {} genes x {} samples".format(n_genes, m_samples))
    print("Found {} dimensions from PCA".format(k_comp))

    # Perform ICA runs
    n_runs = 100
    S_list, A_list = [], []
    for _ in tqdm(range(n_runs)):
        ica = FastICA(whiten="arbitrary-variance", max_iter=int(1e10), tol=1e-7, n_components=k_comp)
        S = pd.DataFrame(ica.fit_transform(rnaseq_df), index=rnaseq_df.index)
        A = pd.DataFrame(ica.mixing_, index=rnaseq_df.columns)
        S_list.append(S)
        A_list.append(A)

    # Cluster the ICA runs to find stable clusters
    labels = cluster_ica_components(S_list, n_runs)
    n_clusters = max(labels) + 1
    print(f"Found {n_clusters} clusters")

    # Aggregate ICA components by cluster
    S_final, A_final, df_stats = aggregate_components(S_list, A_list, labels, n_clusters)

    # Sort components by explained variance
    return sort_components_by_variance(S_final, A_final, df_stats, n_clusters)


def cluster_ica_components(S_list, n_runs):
    """
    Cluster ICA components from multiple runs to find stable clusters.
    """
    # Create a sparse distance matrix for the components
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=EfficiencyWarning)
        D = create_sparse_distance_matrix(S_list)
        D = sparse.csr_matrix((1 - D.data, D.indices, D.indptr))
        
        dbscan = DBSCAN(eps=0.05, min_samples=int(0.5 * n_runs), metric="precomputed")
        labels = dbscan.fit_predict(D)
    return labels


def create_sparse_distance_matrix(S_list):
    """
    Create a sparse distance matrix from a list of ICA component DataFrames.
    """
    n_runs = len(S_list)
    sparse_dist_list = {}
    for i in range(n_runs):
        for j in range(i, n_runs):
            S1, S2 = S_list[i], S_list[j]
            dist = abs(np.dot(S1.T, S2))
            dist[dist < 0.5] = 0
            sparse_dist = sparse.coo_matrix(np.clip(dist, 0, 1))
            sparse_dist_list[(i, j)] = sparse_dist

    # Assemble the block matrix
    block = []
    for i in range(n_runs):
        col = []
        for j in range(n_runs):
            if i <= j:
                mat = sparse_dist_list[(i, j)]
                col.append(mat)
            else:
                mat = sparse_dist_list[(j, i)]
                col.append(mat.T)
        block.append(col)
    D = sparse.bmat(block, "csr")
    return D


def aggregate_components(S_list, A_list, labels, n_clusters):
    """
    Aggregate ICA components by clustering them into n_clusters.
    """
    n_components = S_list[0].shape[1]

    S_bins = {i: [] for i in range(n_clusters)}
    A_bins = {i: [] for i in range(n_clusters)}

    for i, label in enumerate(labels):
        if label != -1:
            S_bins[label].append(S_list[i // n_components][i % n_components])
            A_bins[label].append(A_list[i // n_components][i % n_components])

    S_final = pd.DataFrame(columns=range(n_clusters), index=S_list[0].index)
    A_final = pd.DataFrame(columns=range(n_clusters), index=A_list[0].index)
    df_stats = pd.DataFrame(columns=["S_mean_std", "A_mean_std", "count"], index=range(n_clusters))

    for label in range(n_clusters):
        S_clust, A_clust = S_bins[label], A_bins[label]
        Svec0, Avec0 = S_clust[0], A_clust[0]

        if abs(min(Svec0)) > max(Svec0):
            Svec0, Avec0 = -Svec0, -Avec0

        S_single, A_single = [Svec0], [Avec0]

        for j in range(1, len(S_clust)):
            Svec, Avec = S_clust[j], A_clust[j]
            if stats.pearsonr(Svec, Svec0)[0] > 0:
                S_single.append(Svec)
                A_single.append(Avec)
            else:
                S_single.append(-Svec)
                A_single.append(-Avec)

        # Add centroid of cluster to final S and A matrices
        S_final[label] = np.mean(S_single, axis=0)
        A_final[label] = np.mean(A_single, axis=0)

        # Calculate statistics for components in each cluster
        df_stats.loc[label, "S_mean_std"] = np.std(S_single, axis=0).mean()
        df_stats.loc[label, "A_mean_std"] = np.std(A_single, axis=0).mean()
        df_stats.loc[label, "count"] = len(S_single)

    return S_final, A_final, df_stats


def sort_components_by_variance(S_final, A_final, df_stats, n_clusters):
    """
    Sort the aggregated components by explained variance.
    """
    explained_variances = []
    for label in range(n_clusters):
        k_ica_comp_matrix = np.outer(S_final.iloc[:, label], A_final.iloc[:, label])
        explained_variance = np.var(k_ica_comp_matrix, axis=0).sum()
        explained_variances.append(explained_variance)

    indices_sorted = np.argsort(explained_variances)[::-1]
    S_final_sorted = pd.DataFrame(S_final.iloc[:, indices_sorted].values,
                                  columns=[i for i in range(n_clusters)],
                                  index=S_final.index)
    A_final_sorted = pd.DataFrame(A_final.iloc[:, indices_sorted].values,
                                  columns=[i for i in range(n_clusters)],
                                  index=A_final.index)
    df_stats_sorted = df_stats.iloc[indices_sorted].reset_index(drop=True)

    return S_final_sorted, A_final_sorted, df_stats_sorted
