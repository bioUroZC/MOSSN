from importlib.resources import files

import pandas as pd


def get_example_expression_path() -> str:
    return str(files("mossn.data").joinpath("tcga_blca_expression.csv"))


def get_example_links_path() -> str:
    return str(files("mossn.data").joinpath("string_links.csv"))


def load_example_expression() -> pd.DataFrame:
    return pd.read_csv(get_example_expression_path(), index_col=0)


def load_example_links() -> pd.DataFrame:
    return pd.read_csv(get_example_links_path(), index_col=0)
