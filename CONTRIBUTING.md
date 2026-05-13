# Contributing to MOSSN

Thanks for your interest in improving `mossn`.

## Development setup

```bash
git clone https://github.com/bioUroZC/MOSSN.git
cd MOSSN
pip install -e .[test]
```

## Running tests

```bash
pytest
```

## Project layout

- `src/mossn/`: package source code
- `src/mossn/data/`: bundled example datasets
- `tests/`: automated test suite
- `examples/`: runnable usage examples

## Pull requests

Please keep pull requests focused and include:

- a short explanation of the change
- tests for behavior changes when practical
- documentation updates if the public API or workflow changed

## Reporting issues

When opening an issue, it helps to include:

- the package version
- your Python version
- a minimal reproducible example
- the full traceback or error message
