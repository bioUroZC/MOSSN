import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix


def normalize_sample_id(sample_id: str) -> str:
    return sample_id.replace("-", "_")[:12]


def iqr_normalize(series: pd.Series) -> pd.Series:
    median = series.median()
    iqr = series.quantile(0.75) - series.quantile(0.25)
    if iqr == 0:
        iqr = 1e-6
    return (series - median) / iqr


def filter_links_and_expression(links: pd.DataFrame, expression_data: pd.DataFrame):
    genes_in_links = pd.unique(links[["protein1", "protein2"]].values.ravel())
    common_genes = expression_data.index.intersection(genes_in_links)
    link_filtered = links[
        links["protein1"].isin(common_genes) & links["protein2"].isin(common_genes)
    ].drop_duplicates(subset=["protein1", "protein2"])

    used_genes = pd.unique(link_filtered[["protein1", "protein2"]].values.ravel())
    expression_filtered = expression_data.loc[expression_data.index.isin(used_genes)]
    return link_filtered, expression_filtered


def build_graph_from_links(links: pd.DataFrame, uniform_weight: float = 1.0):
    graph = nx.Graph()
    base_weights = {}
    for _, row in links.iterrows():
        weight = float(uniform_weight)
        u = row["protein1"]
        v = row["protein2"]
        graph.add_edge(u, v, weight=weight)
        base_weights[(u, v)] = weight
    return graph, base_weights


def filter_graph_and_expression(graph: nx.Graph, expression_data: pd.DataFrame):
    common_genes = expression_data.index.intersection(pd.Index(graph.nodes))
    filtered_graph = graph.subgraph(common_genes).copy()
    expression_filtered = expression_data.loc[expression_data.index.isin(filtered_graph.nodes)]
    return filtered_graph, expression_filtered
def reset_graph_edge_attributes(graph: nx.Graph, base_weights: dict):
    for (u, v), weight in base_weights.items():
        if graph.has_edge(u, v):
            graph[u][v]["weight"] = weight
            graph[u][v]["adjusted_weight"] = weight
            graph[u][v]["correction_score"] = 1.0


def build_transition_from_graph(graph: nx.Graph):
    nodes = list(graph.nodes)
    node_to_index = {node: i for i, node in enumerate(nodes)}
    row, col, edge_weights = [], [], []

    for u, v, data in graph.edges(data=True):
        weight = data["adjusted_weight"]
        i = node_to_index[u]
        j = node_to_index[v]
        row.extend([i, j])
        col.extend([j, i])
        edge_weights.extend([weight, weight])

    transition = csr_matrix((edge_weights, (row, col)), shape=(len(nodes), len(nodes)))
    row_sums = np.array(transition.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    transition = transition.multiply(1 / row_sums[:, None])
    return transition, nodes, node_to_index


def apply_expression_correction(graph: nx.Graph, base_weights: dict, expr_norm: pd.Series, gamma: float = 2.0):
    for u, v, data in graph.edges(data=True):
        zu = expr_norm.get(u, 0)
        zv = expr_norm.get(v, 0)
        activation = zu + zv
        sigmoid = 1 / (1 + np.exp(-activation))
        correction_score = max(1 + gamma * (sigmoid - 0.5), 0)
        base_weight = base_weights.get((u, v), base_weights.get((v, u), 0))
        data["adjusted_weight"] = base_weight * correction_score
        data["correction_score"] = correction_score


def seed_vector_from_expression(sample_expression: pd.Series, nodes: list, quantile: float = 0.9):
    threshold = sample_expression.quantile(quantile)
    seed_nodes = set(sample_expression[sample_expression > threshold].index).intersection(nodes)
    p0 = np.array([1.0 if node in seed_nodes else 0.0 for node in nodes])
    if p0.sum() == 0:
        p0[:] = 1.0
    return p0 / p0.sum()


def run_rwr(transition, p0, rwr_alpha: float = 0.3, max_iter: int = 50, tol: float = 1e-4):
    p = p0.copy()
    for _ in range(max_iter):
        p_new = rwr_alpha * p0 + (1 - rwr_alpha) * transition.T @ p
        if np.linalg.norm(p - p_new) < tol:
            return p_new
        p = p_new
    return p


def rank_normalize(values):
    ranks = np.argsort(np.argsort(values))
    return ranks / max(len(values) - 1, 1)


def edge_table_from_graph(sample_id: str, graph: nx.Graph, base_weights: dict, importance, node_to_index: dict):
    rows = []
    for u, v, data in graph.edges(data=True):
        rows.append(
            {
                "Sample": sample_id,
                "Node1": u,
                "Node2": v,
                "BaseWeight": base_weights.get((u, v), base_weights.get((v, u), 0)),
                "FinalWeight": round(
                    data["adjusted_weight"] * (importance[node_to_index[u]] + importance[node_to_index[v]]),
                    5,
                ),
            }
        )
    return pd.DataFrame(rows)
