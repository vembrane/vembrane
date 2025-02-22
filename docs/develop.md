# Development

This project is managed using [uv](https://docs.astral.sh/uv/).
You can build and install vembrane using
```sh
uv sync --locked
source .venv/bin/activate
vembrane --help
```

## Dependencies

All packages that vembrane relies on can be found in the `dependencies` section of `pyproject.toml`.
The exact environment is usually locked (see `uv.lock` for details).

To see the dependencies in the `requirements.txt` format, type in the root directory
```sh
uv export
```

uv uses `virtualenvs` internally.
Either activate the virtualenv with
```sh
source .venv/bin/activate
```
or run a command `cmd` within the `venv` directly with
```sh
uv run cmd
```

## Testing
To run the tests, use
```sh
uv run tox -e $TARGET
```
where `$TARGET` is any of the following:
- `format-check`: check code formatting (with `ruff format --check`)
- `lints`: check code style (with `ruff check`)
- `mypy`: check type annotations (with `mypy`)
- `3.10`, `3.11`, `3.12`, `3.13`, … : run tests with the specified Python version (with `pytest`)


## `pre-commit` hooks

Since we enforce code formatting with `black` by checking for that in CI, we can avoid "fmt" commits by ensuring formatting is done upon comitting changes:
1. make sure `pre-commit` is installed on your machine / in your env (should be available in pip, conda, archlinux repos, ...)
2. run `pre-commit install`. This will activate pre-commit hooks to your _local_ .git

Now when calling `git commit`, your changed code will be formatted with `black` and `isort`, checked with`flake8`, get trailing whitespace removed and trailing newlines added (if needed).
