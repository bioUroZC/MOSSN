import pandas as pd

from mossn import prepare_data_direct_coupled, run_direct_coupled_single_sample


def make_links():
    return pd.DataFrame(
        {
            "protein1": ["A", "A", "B"],
            "protein2": ["B", "C", "C"],
            "score": [0.8, 0.5, 0.9],
        }
    )


def make_omics():
    exp = pd.DataFrame(
        {
            "S-1": [4.2, 7.1, 3.5],
            "S-2": [5.3, 6.4, 2.8],
        },
        index=["A", "B", "C"],
    )
    cnv = pd.DataFrame(
        {
            "S-1": [0.1, 0.8, 0.2],
            "S-2": [0.3, 0.6, 0.4],
        },
        index=["A", "B", "C"],
    )
    return {"EXP": exp, "CNV": cnv}


def test_direct_coupled_pipeline_runs():
    graph, base_weights, omics, exp_genes = prepare_data_direct_coupled(
        make_links(),
        make_omics(),
        coupled_omics=["CNV"],
    )
    sample_id = sorted(omics["EXP"].columns)[0]
    table = run_direct_coupled_single_sample(
        sample_id=sample_id,
        graph=graph,
        base_weights=base_weights,
        omic_data=omics,
        exp_genes=exp_genes,
        coupled_omics=["CNV"],
    )
    assert not table.empty
    assert set(table.columns) == {"Sample", "Node1", "Node2", "BaseWeight", "FinalWeight"}
