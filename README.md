# MOSSN

[![PyPI version](https://img.shields.io/pypi/v/mossn.svg)](https://pypi.org/project/mossn/)
[![Python versions](https://img.shields.io/pypi/pyversions/mossn.svg)](https://pypi.org/project/mossn/)
[![CI](https://github.com/bioUroZC/MOSSN/actions/workflows/ci.yml/badge.svg)](https://github.com/bioUroZC/MOSSN/actions/workflows/ci.yml)

`mossn` packages the MOSSN algorithm for constructing sample-specific protein
interaction networks from gene expression data, together with a direct-coupled
multi-omics extension.

## Highlights

- Single-sample and multi-sample workflows for sample-specific network scoring
- Support for either a PPI links table or a custom `networkx.Graph`
- Direct-coupled multi-omics workflow for integrating matched omics layers
- Bundled example datasets for quick experiments and reproducible demos

## Method Overview

MOSSN builds a sample-specific interaction network by combining a background
protein-protein interaction graph with sample-level molecular measurements.
Every background edge starts with the same weight, following the uniform MOSSN
workflow. Input edge scores are not used.
In the single-omics workflow, the package follows four main steps:

1. Filter the background network and expression matrix to their shared genes.
2. Reweight each edge using IQR-normalized sample expression and the `gamma`
   correction strength.
3. Run random walk with restart from high-expression seed genes, or from a
   uniform prior when seed selection is disabled.
4. Combine edge-level correction and node-level importance into final
   sample-specific edge weights.

The direct-coupled multi-omics workflow adds cross-layer edges linking
expression nodes to other omics measurements for the same gene, then applies a
similar propagation scheme on the coupled graph.

## When To Use Which Workflow

- Use `prepare_data(...)` and `run_single_sample(...)` when you already have a
  curated PPI network and matched expression data.
- Use `prepare_data_direct_coupled(...)` when you have matched expression and
  additional omics layers such as CNV.

## Installation

`mossn` supports Python 3.9, 3.10, and 3.11. The GitHub Actions CI workflow
tests these versions on every push and pull request.

```bash
pip install mossn
```

For local development:

```bash
pip install -e .[test]
```

## Quick Start

```python
import pandas as pd
from mossn import prepare_data, run_single_sample

links = pd.DataFrame(
    {
        "protein1": ["A", "A", "B"],
        "protein2": ["B", "C", "C"],
    }
)

expression_data = pd.DataFrame(
    {
        "sample_1": [4.2, 7.1, 3.5],
        "sample_2": [5.3, 6.4, 2.8],
    },
    index=["A", "B", "C"],
)

graph, base_weights, expression_data = prepare_data(
    expression_data=expression_data,
    links=links,
)

edge_table = run_single_sample(
    sample_id="sample_1",
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
    gamma=2.0,
    rwr_alpha=0.3,
    seed_quantile=0.9,
)

print(edge_table.head())
```

To process every sample in the expression matrix:

```python
from mossn import prepare_data, run_samples

graph, base_weights, expression_data = prepare_data(
    expression_data=expression_data,
    links=links,
)

edge_table = run_samples(
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
)
```

## Core API

The recommended public API consists of these five functions:

- `prepare_data(...)`: build a single-omics background network from a links
  table or `networkx.Graph`
- `run_single_sample(...)`: infer a sample-specific network for one sample
- `run_samples(...)`: run the single-omics workflow across multiple samples
- `prepare_data_direct_coupled(...)`: construct a direct-coupled multi-omics
  graph
- `run_direct_coupled_single_sample(...)`: score one sample in the direct-coupled
  multi-omics setting

## Common Inputs

### PPI links table

The `links` table should contain:

- `protein1`
- `protein2`

Any additional columns, including edge scores, are ignored. Set
`uniform_weight` in `prepare_data(...)` to choose the common starting weight for
all background edges.

### Expression matrix

The expression matrix should use:

- rows as genes or proteins
- columns as sample IDs

## Main Parameters

The single-omics workflow exposes three main numeric parameters:

- `gamma`: strength of expression-based edge reweighting
- `rwr_alpha`: restart probability in random walk with restart
- `seed_quantile`: expression quantile used to define seed genes

It also exposes the following switches for sample-specific scoring:

- `use_seed=True`: use high-expression seed genes
- `use_rwr=True`: use random walk with restart
- `use_correction=True`: use expression-based edge correction

For example, to disable seed genes:

```python
edge_table = run_single_sample(
    sample_id="sample_1",
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
    use_seed=False,
)
```

## Output Format

The main scoring functions return a tidy edge table with these columns:

- `Sample`: sample identifier
- `Node1`: first endpoint of the interaction
- `Node2`: second endpoint of the interaction
- `BaseWeight`: background-network edge weight before sample-specific scoring
- `FinalWeight`: sample-specific edge score returned by MOSSN

This output format is designed to be easy to sort, filter, merge, and export
for downstream biological analysis.

## More Examples

### Use your own network

```python
import networkx as nx
from mossn import prepare_data, run_single_sample
from mossn.example_data import load_example_expression

graph = nx.Graph()
graph.add_edge("A", "B")
graph.add_edge("B", "C")

expression_data = load_example_expression()
graph, base_weights, expression_data = prepare_data(
    expression_data=expression_data,
    graph=graph,
)

edge_table = run_single_sample(
    sample_id=expression_data.columns[0],
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
    gamma=1.5,
    rwr_alpha=0.2,
    seed_quantile=0.8,
)
```

### Use the bundled example data

```python
from mossn import prepare_data, run_single_sample
from mossn.example_data import load_example_expression, load_example_links

links = load_example_links()
expression_data = load_example_expression()

graph, base_weights, expression_data = prepare_data(
    expression_data=expression_data,
    links=links,
)

edge_table = run_single_sample(
    sample_id=expression_data.columns[0],
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
)
```

### Direct-coupled multi-omics

```python
from mossn import prepare_data_direct_coupled, run_direct_coupled_single_sample

graph, base_weights, omic_data, exp_genes = prepare_data_direct_coupled(
    links=links,
    omic_data=omic_data,
    coupled_omics=["CNV"],
)

edge_table = run_direct_coupled_single_sample(
    sample_id=sample_id,
    graph=graph,
    base_weights=base_weights,
    omic_data=omic_data,
    exp_genes=exp_genes,
    coupled_omics=["CNV"],
)
```

## Example Data Helpers

The package includes:

- TCGA BLCA expression matrix
- STRING-derived PPI links

You can access them with:

```python
from mossn.example_data import (
    get_example_expression_path,
    get_example_links_path,
    load_example_expression,
    load_example_links,
)
```

## Optional Data-Driven Extension

When no background interaction network is available, the optional
`prepare_data_driven(...)` helper can infer one from expression correlations.
This is an auxiliary workflow rather than the recommended uniform PPI workflow.

```python
from mossn import prepare_data_driven, run_samples

graph, base_weights, expression_data = prepare_data_driven(
    expression_data=expression_data,
    cor_threshold=0.9,
)
edge_table = run_samples(
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
)
```

## Repository Guide

- [CONTRIBUTING.md](CONTRIBUTING.md): local development and pull request notes
- `.github/workflows/`: CI and release automation
- `.github/ISSUE_TEMPLATE/`: issue forms for bug reports, feature requests, and usage questions
- `examples/`: runnable scripts
- `tests/`: automated test suite
- `src/mossn/`: package source code

## License

`mossn` is distributed under a research-only license. Non-commercial research,
teaching, and evaluation use are allowed. Commercial use requires prior written
permission from the copyright holder.

## Development Status

The repository currently maintains CI coverage for Python 3.9, 3.10, and 3.11.
If you encounter a version-specific issue, please open a GitHub issue and
include your Python version and a minimal reproduction example.
