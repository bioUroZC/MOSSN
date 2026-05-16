from typing import Optional

import numpy as np
import pandas as pd
import networkx as nx

from ._utils import (
    apply_expression_correction,
    build_graph_from_links,
    build_transition_from_graph,
    edge_table_from_graph,
    extract_base_weights,
    filter_links_and_expression,
    filter_graph_and_expression,
    iqr_normalize,
    rank_normalize,
    reset_graph_edge_attributes,
    run_rwr,
    seed_vector_from_expression,
)


def prepare_data(
    expression_data: pd.DataFrame,
    links: Optional[pd.DataFrame] = None,
    graph: Optional[nx.Graph] = None,
    base_weights: Optional[dict] = None,
    weight_key: str = "score",
    use_prior: bool = True,
    uniform_weight: float = 1.0,
):
    if links is None and graph is None:
        raise ValueError("Either `links` or `graph` must be provided.")
    if links is not None and graph is not None:
        raise ValueError("Provide only one of `links` or `graph`, not both.")

    if links is not None:
        filtered_links, expression_filtered = filter_links_and_expression(links, expression_data)
        if use_prior:
            graph_obj, base_weights_obj = build_graph_from_links(filtered_links, weight_key=weight_key)
        else:
            graph_obj, base_weights_obj = build_graph_from_links(filtered_links, default_weight=uniform_weight)
        return graph_obj, base_weights_obj, expression_filtered

    filtered_graph, expression_filtered = filter_graph_and_expression(graph, expression_data, weight_key=weight_key)
    base_weights_obj = base_weights if base_weights is not None else extract_base_weights(filtered_graph, weight_key=weight_key)
    return filtered_graph, base_weights_obj, expression_filtered


def run_single_sample(
    sample_id: str,
    graph,
    base_weights: dict,
    expression_data: pd.DataFrame,
    gamma: float = 2.0,
    rwr_alpha: float = 0.3,
    seed_quantile: float = 0.9,
    use_seed: bool = True,
    use_rwr: bool = True,
    use_correction: bool = True,
):
    reset_graph_edge_attributes(graph, base_weights)
    sample_expression = expression_data[sample_id]
    if use_correction:
        expr_norm = iqr_normalize(sample_expression)
        apply_expression_correction(graph, base_weights, expr_norm, gamma=gamma)

    if not use_rwr:
        rows = []
        for u, v, data in graph.edges(data=True):
            rows.append(
                {
                    "Sample": sample_id,
                    "Node1": u,
                    "Node2": v,
                    "BaseWeight": base_weights.get((u, v), base_weights.get((v, u), 0)),
                    "FinalWeight": round(data["adjusted_weight"], 5),
                }
            )
        return pd.DataFrame(rows)

    transition, nodes, node_to_index = build_transition_from_graph(graph)
    if use_seed:
        p0 = seed_vector_from_expression(sample_expression, nodes, quantile=seed_quantile)
    else:
        p0 = np.ones(len(nodes)) / len(nodes)
    importance = rank_normalize(run_rwr(transition, p0, rwr_alpha=rwr_alpha))

    if use_correction:
        return edge_table_from_graph(sample_id, graph, base_weights, importance, node_to_index)

    rows = []
    for u, v in graph.edges():
        rows.append(
            {
                "Sample": sample_id,
                "Node1": u,
                "Node2": v,
                "BaseWeight": base_weights.get((u, v), base_weights.get((v, u), 0)),
                "FinalWeight": round(
                    base_weights.get((u, v), base_weights.get((v, u), 0))
                    * (importance[node_to_index[u]] + importance[node_to_index[v]]),
                    5,
                ),
            }
    )
    return pd.DataFrame(rows)


def run_samples(
    graph,
    base_weights: dict,
    expression_data: pd.DataFrame,
    sample_ids=None,
    gamma: float = 2.0,
    rwr_alpha: float = 0.3,
    seed_quantile: float = 0.9,
    use_seed: bool = True,
    use_rwr: bool = True,
    use_correction: bool = True,
):
    if sample_ids is None:
        sample_ids = list(expression_data.columns)

    results = []
    for sample_id in sample_ids:
        results.append(
            run_single_sample(
                sample_id=sample_id,
                graph=graph,
                base_weights=base_weights,
                expression_data=expression_data,
                gamma=gamma,
                rwr_alpha=rwr_alpha,
                seed_quantile=seed_quantile,
                use_seed=use_seed,
                use_rwr=use_rwr,
                use_correction=use_correction,
            )
        )

    if not results:
        return pd.DataFrame(columns=["Sample", "Node1", "Node2", "BaseWeight", "FinalWeight"])
    return pd.concat(results, ignore_index=True)
