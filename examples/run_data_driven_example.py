from mossn import prepare_data_driven, run_samples
from mossn.example_data import load_example_expression


def main():
    expression_data = load_example_expression()

    graph, base_weights, expression_data = prepare_data_driven(
        expression_data=expression_data,
        cor_threshold=0.9,
    )
    sample_id = expression_data.columns[0]
    edge_table = run_samples(
        graph=graph,
        base_weights=base_weights,
        expression_data=expression_data,
        sample_ids=[sample_id],
    )

    print(f"sample_id: {sample_id}")
    print(f"nodes: {graph.number_of_nodes()}")
    print(f"edges: {graph.number_of_edges()}")
    print(edge_table.head())


if __name__ == "__main__":
    main()
