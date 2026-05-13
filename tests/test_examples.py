from mossn import (
    get_example_expression_path,
    get_example_links_path,
    load_example_expression,
    load_example_links,
)


def test_example_paths_exist():
    assert get_example_expression_path().endswith("tcga_blca_expression.csv")
    assert get_example_links_path().endswith("string_links.csv")


def test_example_data_loads():
    expression = load_example_expression()
    links = load_example_links()
    assert expression.shape[0] > 0
    assert expression.shape[1] > 0
    assert {"protein1", "protein2", "score"}.issubset(links.columns)
