import numpy as np
import pandas as pd
import networkx as nx


def _rankdata_2d(values):
    temp = np.argsort(values, axis=1)
    ranks = np.empty_like(temp, dtype=float)
    row_index = np.arange(values.shape[0])[:, None]
    ranks[row_index, temp] = np.arange(1, values.shape[1] + 1).astype(float)
    return ranks


def prepare_data_driven(
    expression_data: pd.DataFrame,
    cor_threshold: float = 0.9,
):
    expression_filtered = expression_data.loc[expression_data.std(axis=1) > 0]

    ranked = _rankdata_2d(expression_filtered.values).astype(float)
    ranked -= ranked.mean(axis=1, keepdims=True)
    norms = np.sqrt((ranked ** 2).sum(axis=1, keepdims=True)) + 1e-10
    ranked /= norms
    corr_matrix = ranked @ ranked.T
    genes = expression_filtered.index.tolist()

    rows_idx, cols_idx = np.where(np.triu(np.abs(corr_matrix) > cor_threshold, k=1))

    graph = nx.Graph()
    base_weights = {}
    for i, j in zip(rows_idx, cols_idx):
        u = genes[i]
        v = genes[j]
        weight = float(abs(corr_matrix[i, j]))
        graph.add_edge(u, v, weight=weight)
        base_weights[(u, v)] = weight

    used_genes = list(graph.nodes)
    expression_filtered = expression_data.loc[expression_data.index.isin(used_genes)]
    return graph, base_weights, expression_filtered
