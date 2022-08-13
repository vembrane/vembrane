# Development

This project is managed using [poetry](https://python-poetry.org/).
You can build and install vembrane using
```sh
poetry install
```

## Dependencies

All packages that vembrane relies on can be found in the `[tool.poetry.dependencies]` section of `pyproject.toml`.
The exact environment is usually locked (see `poetry.lock` for details).


## `pre-commit` hooks

Since we enforce code formatting with `black` by checking for that in CI, we can avoid "fmt" commits by ensuring formatting is done upon comitting changes:
1. make sure `pre-commit` is installed on your machine / in your env (should be available in pip, conda, archlinux repos, ...)
2. run `pre-commit install`. This will activate pre-commit hooks to your _local_ .git

Now when calling `git commit`, your changed code will be formatted with `black` and `isort`, checked with`flake8`, get trailing whitespace removed and trailing newlines added (if needed).
