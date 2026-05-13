import pandas as pd

from mossn import prepare_data_driven, run_samples, run_single_sample


def make_expression():
    return pd.DataFrame(
        {
            "s1": [1.0, 1.1, 5.0, 0.5],
            "s2": [2.0, 2.1, 4.0, 0.7],
            "s3": [3.0, 3.2, 3.0, 0.9],
            "s4": [4.0, 4.1, 2.0, 1.1],
        },
        index=["G1", "G2", "G3", "G4"],
    )


def test_prepare_data_driven_filters_zero_variance_rows():
    expression = make_expression().copy()
    expression.loc["CONST"] = [2.0, 2.0, 2.0, 2.0]
    graph, base_weights, filtered = prepare_data_driven(expression, cor_threshold=0.5)
    assert "CONST" not in filtered.index
    assert filtered.shape[0] == 4
    assert graph.number_of_nodes() > 0
    assert len(base_weights) == graph.number_of_edges()


def test_prepare_data_driven_runs():
    expression = make_expression()
    graph, base_weights, filtered = prepare_data_driven(expression, cor_threshold=0.5)
    assert graph.number_of_nodes() > 0
    assert graph.number_of_edges() > 0
    assert len(base_weights) == graph.number_of_edges()


def test_prepare_data_driven_and_single_sample_run():
    graph, base_weights, expression = prepare_data_driven(make_expression(), cor_threshold=0.5)
    table = run_single_sample("s1", graph, base_weights, expression)
    assert not table.empty
    assert set(table.columns) == {"Sample", "Node1", "Node2", "BaseWeight", "FinalWeight"}


def test_prepare_data_driven_and_batch_run():
    graph, base_weights, expression = prepare_data_driven(make_expression(), cor_threshold=0.5)
    table = run_samples(graph, base_weights, expression, sample_ids=["s1", "s2"])
    assert not table.empty
    assert set(table["Sample"]) == {"s1", "s2"}
