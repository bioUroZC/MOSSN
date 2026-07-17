import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix

from ._utils import iqr_normalize, normalize_sample_id, rank_normalize, run_rwr


def prepare_data_direct(links: pd.DataFrame, omic_data: dict, direct_omics, uniform_weight: float = 1.0):
    all_omics = {"EXP": omic_data["EXP"]}
    for omic in direct_omics:
        if omic in omic_data:
            all_omics[omic] = omic_data[omic]

    all_omics = {k: v.copy().rename(columns=normalize_sample_id) for k, v in all_omics.items()}

    genes_in_links = pd.unique(links[["protein1", "protein2"]].values.ravel())
    common_genes = pd.Index(genes_in_links)
    for table in all_omics.values():
        common_genes = common_genes.intersection(table.index)

    filtered_links = links[
        links["protein1"].isin(common_genes) & links["protein2"].isin(common_genes)
    ].drop_duplicates(subset=["protein1", "protein2"])

    used_genes = pd.unique(filtered_links[["protein1", "protein2"]].values.ravel())
    common_genes = common_genes.intersection(used_genes)

    common_samples = None
    for table in all_omics.values():
        columns = table.columns
        common_samples = columns if common_samples is None else common_samples.intersection(columns)

    all_omics = {
        k: v.loc[list(common_genes), list(common_samples)].sort_index(axis=0).sort_index(axis=1)
        for k, v in all_omics.items()
    }

    graph = nx.Graph()
    weights = {}
    for _, row in filtered_links.iterrows():
        u = row["protein1"]
        v = row["protein2"]
        graph.add_edge(u, v, weight=float(uniform_weight))
        weights[(u, v)] = float(uniform_weight)

    for omic in direct_omics:
        if omic not in all_omics:
            continue
        for gene in common_genes:
            omic_node = f"{omic}:{gene}"
            graph.add_node(omic_node)
            graph.add_edge(gene, omic_node, weight=1.0)
            weights[(gene, omic_node)] = 1.0

    return graph, weights, all_omics, list(common_genes)


def prepare_data_direct_coupled(
    links: pd.DataFrame,
    omic_data: dict,
    coupled_omics,
    uniform_weight: float = 1.0,
):
    return prepare_data_direct(
        links=links,
        omic_data=omic_data,
        direct_omics=coupled_omics,
        uniform_weight=uniform_weight,
    )


def run_direct_single_sample(
    sample_id: str,
    graph,
    base_weights: dict,
    omic_data: dict,
    exp_genes,
    direct_omics,
    dynamic_cross: bool = True,
    rwr_alpha: float = 0.3,
):
    alpha_mod = 1.0
    gamma = 2.0

    z_norms = {}
    for omic, data in omic_data.items():
        if sample_id in data.columns:
            z_norms[omic] = iqr_normalize(data[sample_id]).to_dict()

    z_exp = z_norms.get("EXP", {})

    for u, v, data in graph.edges(data=True):
        base = base_weights.get((u, v), base_weights.get((v, u), 0))
        u_is_exp = ":" not in u
        v_is_exp = ":" not in v

        if u_is_exp and v_is_exp:
            zu = z_exp.get(u, 0)
            zv = z_exp.get(v, 0)
            sigmoid = 1 / (1 + np.exp(-alpha_mod * (zu + zv)))
            data["adjusted_weight"] = max(base * (1 + gamma * (sigmoid - 0.5)), 0)
        else:
            gene = u if u_is_exp else v
            omic_node = v if u_is_exp else u
            omic = omic_node.split(":")[0]
            z_other = z_norms.get(omic, {})

            if dynamic_cross:
                zo = z_other.get(gene, 0)
                ze = z_exp.get(gene, 0)
                sigmoid = 1 / (1 + np.exp(-(zo * ze)))
                data["adjusted_weight"] = max(base * (1 + gamma * (sigmoid - 0.5)), 0)
            else:
                data["adjusted_weight"] = base

    nodes = list(graph.nodes)
    idx = {node: i for i, node in enumerate(nodes)}
    rows, cols, vals = [], [], []
    for u, v, data in graph.edges(data=True):
        weight = data["adjusted_weight"]
        i = idx[u]
        j = idx[v]
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([weight, weight])

    transition = csr_matrix((vals, (rows, cols)), shape=(len(nodes), len(nodes)))
    row_sums = np.array(transition.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    transition = transition.multiply(1 / row_sums[:, None])

    seed_nodes = set()
    exp_data = omic_data["EXP"]
    if sample_id in exp_data.columns:
        q90 = exp_data[sample_id].quantile(0.9)
        seed_nodes.update(gene for gene in exp_data[sample_id][exp_data[sample_id] > q90].index if gene in graph.nodes)
    for omic in direct_omics:
        if omic not in omic_data or sample_id not in omic_data[omic].columns:
            continue
        series = omic_data[omic][sample_id]
        q90 = series.quantile(0.9)
        seed_nodes.update(
            f"{omic}:{gene}" for gene in series[series > q90].index if f"{omic}:{gene}" in graph.nodes
        )

    p0 = np.array([1.0 if node in seed_nodes else 0.0 for node in nodes])
    if p0.sum() == 0:
        p0[:] = 1.0
    p0 /= p0.sum()

    p = run_rwr(transition, p0, rwr_alpha=rwr_alpha)
    exp_idx = [idx[gene] for gene in exp_genes if gene in idx]
    p_exp_sub = p[exp_idx]
    norm_sub = rank_normalize(p_exp_sub)

    importance = np.zeros(len(nodes))
    for rank_i, node_i in enumerate(exp_idx):
        importance[node_i] = norm_sub[rank_i]

    rows_out = []
    for u, v, data in graph.edges(data=True):
        if ":" in u or ":" in v:
            continue
        base = base_weights.get((u, v), base_weights.get((v, u), 0))
        final_weight = data["adjusted_weight"] * (importance[idx[u]] + importance[idx[v]])
        rows_out.append(
            {
                "Sample": sample_id,
                "Node1": u,
                "Node2": v,
                "BaseWeight": base,
                "FinalWeight": round(final_weight, 5),
            }
        )
    return pd.DataFrame(rows_out)


def run_direct_coupled_single_sample(
    sample_id: str,
    graph,
    base_weights: dict,
    omic_data: dict,
    exp_genes,
    coupled_omics,
    dynamic_cross: bool = True,
    rwr_alpha: float = 0.3,
):
    return run_direct_single_sample(
        sample_id=sample_id,
        graph=graph,
        base_weights=base_weights,
        omic_data=omic_data,
        exp_genes=exp_genes,
        direct_omics=coupled_omics,
        dynamic_cross=dynamic_cross,
        rwr_alpha=rwr_alpha,
    )
