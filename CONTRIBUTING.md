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

## Release workflow

The repository includes a tag-based GitHub Actions release workflow.

1. Update the version in `pyproject.toml`.
2. Update `src/mossn/__init__.py` to the same version.
3. Run `pytest`.
4. Commit the release changes.
5. Create and push a tag such as `v0.1.2`.

```bash
git tag v0.1.2
git push origin main --tags
```

When the tag reaches GitHub, the release workflow will:

- verify that the tag matches the package version
- run the test suite
- build the source distribution and wheel
- publish the package to PyPI
- create a GitHub Release and attach the built artifacts

## PyPI publishing setup

The release workflow is designed for PyPI Trusted Publishing.

Before the first automated release, configure `bioUroZC/MOSSN` as a trusted
publisher for the `mossn` project in PyPI. Once that is configured, no API
token needs to be stored in GitHub Actions secrets.

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
