# mossn

`mossn` packages the MOSSN algorithm for constructing sample-specific protein
interaction networks from gene expression data, together with ablation variants
and a direct-coupled multi-omics extension.

The main single-omics API is:

- `prepare_data(...)`
- `run_single_sample(...)`
- `run_samples(...)`

## Features

- Sample-specific edge reweighting using gene-expression-derived correction
  scores.
- Random walk with restart (RWR) to estimate node importance per sample.
- User-tunable parameters in the main API:
  - `gamma`
  - `rwr_alpha`
  - `seed_quantile`
  - `use_seed`
  - `use_rwr`
  - `use_correction`
  - `use_prior`
- Support for either:
  - a user-provided PPI links table
  - a user-provided `networkx.Graph`
- Multi-omics extension:
  - direct-coupled cross-layer graph

## Installation

```bash
pip install mossn
```

For local development:

```bash
pip install -e .
```

## Quick Start

```python
import pandas as pd
from mossn import prepare_data, run_single_sample

links = pd.DataFrame(
    {
        "protein1": ["A", "A", "B"],
        "protein2": ["B", "C", "C"],
        "score": [0.8, 0.6, 0.9],
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
    use_prior=True,
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

If you want to run the full expression matrix sample by sample:

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

## Example Data

If you want to try the packaged BLCA and STRING example files:

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

## Main Parameters

The main single-omics workflow exposes three tunable parameters:

- `gamma`: strength of expression-based edge reweighting
- `rwr_alpha`: restart probability in random walk with restart
- `seed_quantile`: expression quantile used to define seed genes

It also exposes four logical switches:

- `use_seed=True`: use high-expression seed genes
- `use_rwr=True`: use random walk with restart
- `use_correction=True`: use expression-based edge correction
- `use_prior=True`: use input edge weights instead of uniform weights

For example, if you do not want to use seed genes:

```python
edge_table = run_single_sample(
    sample_id="sample_1",
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
    use_seed=False,
)
```

If you want to ignore prior edge weights and use a uniform network:

```python
graph, base_weights, expression_data = prepare_data(
    expression_data=expression_data,
    links=links,
    use_prior=False,
    uniform_weight=1.0,
)
```

## Use Your Own Network

You can also provide your own `networkx.Graph` instead of a links table. Edge
weights are read from the `weight` attribute by default.

```python
import networkx as nx
from mossn import prepare_data, run_single_sample
from mossn.example_data import load_example_expression

graph = nx.Graph()
graph.add_edge("A", "B", weight=0.8)
graph.add_edge("B", "C", weight=0.6)

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

## Bundled Example Data

The package includes the following example datasets:

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

## Input format

### PPI links

The `links` table must contain:

- `protein1`
- `protein2`
- `score`

### Expression matrix

The expression matrix must use:

- rows as genes or proteins
- columns as sample IDs

## Main API

The recommended public API consists of six core functions.

### 1. `prepare_data(...)`

Prepare a single-omics background network from either:

- a PPI links table
- a user-provided `networkx.Graph`

Returns:

- `graph`
- `base_weights`
- filtered `expression_data`

### 2. `run_single_sample(...)`

Run MOSSN for one sample.

Returns:

- one sample-specific edge-weight table

### 3. `run_samples(...)`

Run MOSSN across multiple samples in an expression matrix.

Returns:

- a combined edge-weight table for all requested samples

### 4. `prepare_data_driven(...)`

Construct a data-driven background network directly from the expression matrix
when no external prior network is available.

Returns:

- `graph`
- `base_weights`
- filtered `expression_data`

### 5. `prepare_data_direct_coupled(...)`

Prepare the direct-coupled multi-omics graph.

Returns:

- `graph`
- `base_weights`
- filtered `omic_data`
- `exp_genes`

### 6. `run_direct_coupled_single_sample(...)`

Run the direct-coupled multi-omics version for one sample.

Returns:

- one sample-specific edge-weight table

The current package release supports only the direct-coupled multi-omics
extension.

## Direct-Coupled Multi-Omics Example

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

## License

`mossn` is distributed under a research-only license. Non-commercial research,
teaching, and evaluation use are allowed. Commercial use requires prior written
permission from the copyright holder.

## Notes

- The package expects matched identifiers between the network and omics tables.
- The main single-sample API lets you set `gamma`, `rwr_alpha`,
  `seed_quantile`, `use_seed`, `use_rwr`, and `use_correction` directly.
- Sample-specific normalization uses median and interquartile range (IQR).
- Node importance is rank-normalized before computing final edge weights.
- The data-driven mode first infers a background graph from expression
  correlations when an external reference network is unavailable.

## Data-Driven Example

```python
from mossn import prepare_data_driven, run_samples
from mossn.example_data import load_example_expression

expression_data = load_example_expression()
graph, base_weights, expression_data = prepare_data_driven(
    expression_data=expression_data,
    cor_threshold=0.9,
)

edge_table = run_samples(
    graph=graph,
    base_weights=base_weights,
    expression_data=expression_data,
    sample_ids=[expression_data.columns[0]],
)
```
