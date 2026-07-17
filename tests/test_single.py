import pandas as pd
import networkx as nx

from mossn import (
    prepare_data,
    run_single_sample,
)


def make_links():
    return pd.DataFrame(
        {
            "protein1": ["A", "A", "B", "D"],
            "protein2": ["B", "C", "C", "E"],
        }
    )


def make_expression():
    return pd.DataFrame(
        {
            "s1": [4.2, 7.1, 3.5, 2.0, 1.0],
            "s2": [5.3, 6.4, 2.8, 2.5, 0.8],
        },
        index=["A", "B", "C", "D", "E"],
    )


def assert_edge_table(table, sample_id):
    assert not table.empty
    assert set(table.columns) == {"Sample", "Node1", "Node2", "BaseWeight", "FinalWeight"}
    assert (table["Sample"] == sample_id).all()


def test_clean_main_api_runs():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample(
        "s1",
        graph,
        base_weights,
        expression,
        gamma=1.5,
        rwr_alpha=0.2,
        seed_quantile=0.8,
    )
    assert_edge_table(table, "s1")
    assert set(base_weights.values()) == {1.0}


def test_main_api_supports_uniform_weights():
    graph, base_weights, expression = prepare_data(
        expression_data=make_expression(),
        links=make_links(),
        uniform_weight=2.0,
    )
    table = run_single_sample("s1", graph, base_weights, expression)
    assert_edge_table(table, "s1")
    assert set(table["BaseWeight"]) == {2.0}


def test_link_scores_are_ignored():
    links = make_links().assign(score=[0.8, 0.5, 0.9, 0.4])
    graph, base_weights, expression = prepare_data(
        expression_data=make_expression(),
        links=links,
    )
    table = run_single_sample("s1", graph, base_weights, expression)
    assert_edge_table(table, "s1")
    assert set(table["BaseWeight"]) == {1.0}


def test_main_api_supports_no_seed():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample("s1", graph, base_weights, expression, use_seed=False)
    assert_edge_table(table, "s1")


def test_main_api_supports_no_rwr():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample("s1", graph, base_weights, expression, use_rwr=False)
    assert_edge_table(table, "s1")


def test_main_api_supports_no_correction():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample("s1", graph, base_weights, expression, use_correction=False)
    assert_edge_table(table, "s1")


def test_clean_main_api_accepts_custom_graph():
    custom_graph = nx.Graph()
    custom_graph.add_edge("A", "B", weight=0.7)
    custom_graph.add_edge("B", "C", weight=0.6)
    custom_graph.add_edge("D", "E", weight=0.4)

    graph, base_weights, expression = prepare_data(
        expression_data=make_expression(),
        graph=custom_graph,
    )
    table = run_single_sample("s1", graph, base_weights, expression)
    assert_edge_table(table, "s1")
    assert graph.number_of_nodes() == 5
    assert len(base_weights) == graph.number_of_edges()
    assert set(base_weights.values()) == {1.0}


def test_uniform_pipeline_uses_constant_base_weights():
    graph, base_weights, expression = prepare_data(
        expression_data=make_expression(),
        links=make_links(),
        uniform_weight=1.0,
    )
    table = run_single_sample("s1", graph, base_weights, expression)
    assert_edge_table(table, "s1")
    assert table["BaseWeight"].nunique() == 1
    assert table["BaseWeight"].iloc[0] == 1.0


def test_no_seed_pipeline_runs():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample("s1", graph, base_weights, expression, use_seed=False)
    assert_edge_table(table, "s1")


def test_no_rwr_pipeline_returns_adjusted_weights():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample("s1", graph, base_weights, expression, use_rwr=False)
    assert_edge_table(table, "s1")
    assert (table["FinalWeight"] >= 0).all()


def test_no_corr_pipeline_runs():
    graph, base_weights, expression = prepare_data(expression_data=make_expression(), links=make_links())
    table = run_single_sample("s1", graph, base_weights, expression, use_correction=False)
    assert_edge_table(table, "s1")
    assert (table["FinalWeight"] >= 0).all()
