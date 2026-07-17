from mossn import prepare_data, run_single_sample
from mossn.example_data import load_example_expression, load_example_links


def main():
    links = load_example_links()
    expression_data = load_example_expression()

    graph, base_weights, expression_data = prepare_data(
        expression_data=expression_data,
        links=links,
    )
    sample_id = expression_data.columns[0]

    edge_table = run_single_sample(
        sample_id=sample_id,
        graph=graph,
        base_weights=base_weights,
        expression_data=expression_data,
    )

    print(f"sample_id: {sample_id}")
    print(f"nodes: {graph.number_of_nodes()}")
    print(f"edges: {graph.number_of_edges()}")
    print(edge_table.head())


if __name__ == "__main__":
    main()
