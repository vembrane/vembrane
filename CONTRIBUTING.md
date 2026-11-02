# Contributing

## Setup

vembrane is managed with [uv](https://docs.astral.sh/uv/). From the repository root:
```sh
uv sync --locked
source .venv/bin/activate   # or prefix commands with `uv run`
vembrane --help
```

Dependencies are listed in `pyproject.toml` and locked in `uv.lock`.
`uv export` prints the locked environment in `requirements.txt` format.

## Checks

CI runs the following tox environments (Python >= 3.12, tested on 3.12):
```sh
uv run tox -e format-check   # ruff format --check
uv run tox -e lints          # ruff check
uv run tox -e mypy           # mypy
uv run tox -e 3.12           # pytest tests
```
Arguments after `--` are passed to pytest. To run pytest directly:
```sh
uv run pytest tests
uv run pytest tests -k 'filter-test_ann_impact'
```
Test ids look like `test_command[filter-test_ann_impact-cyvcf2-Context.NONE]`, so `-k 'filter-test_ann_impact'` selects all variants of one case.

Run all commands from the repository root: test cases reference files such as `tests/resources/context.py` by relative path.

## pre-commit

The hooks run `ruff --fix` and `ruff-format`, so formatting needs no separate commit:
```sh
uv run pre-commit install
```

## Commit messages

Use [Conventional Commits](https://www.conventionalcommits.org/) (`fix: ...`, `feat: ...`, `docs: ...`).
release-please derives version bumps and `CHANGELOG.md` from them, so do not edit `CHANGELOG.md` or the version in `vembrane/__init__.py` by hand.

## Test cases

Each case is a directory `tests/testcases/<subcommand>/<case>/` with
- `config.yaml`: `function` (the subcommand) plus the subcommand's arguments. Other keys map to CLI flags, with `_` replaced by `-` (e.g. `output_fmt` becomes `--output-fmt`).
- `test.vcf`: the input.
- the expected output: `expected.vcf`, `expected.tsv`, or `expected.<output_fmt>` for `structured` and `fhir`.

Every case runs against both backends, pysam and cyvcf2.
Look at an existing case of the same subcommand for the required keys.

## VCF specification

Where vembrane does not rely directly on a backend, it follows the current [VCF specification](https://github.com/samtools/hts-specs).
The backends do not define what is valid: they can accept more or less than the specification allows, so check the specification when in doubt or when they disagree.
Test data should be valid under it too.

## Documentation

`docs/<subcommand>.md` and `README.md` contain copies of the `--help` output. Update them when you change options.
