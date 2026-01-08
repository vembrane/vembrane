import argparse
import builtins
import csv
import json
import os
import tempfile
from enum import Enum, auto
from itertools import product, zip_longest
from pathlib import Path
from typing import Any, Callable

import pyarrow.parquet as pq
import pytest
import yaml

from vembrane import __version__, errors
from vembrane.backend.base import Backend
from vembrane.common import create_reader
from vembrane.modules import fhir, filter, sort, structured, table, tag

TEST_DIR = Path(__file__).parent
FILTER_CASES = TEST_DIR / "testcases/filter"
TABLE_CASES = TEST_DIR / "testcases/table"
TAG_CASES = TEST_DIR / "testcases/tag"
ANNOTATE_CASES = TEST_DIR / "testcases/annotate"
STRUCTURED_CASES = TEST_DIR / "testcases/structured"
FHIR_CASES = TEST_DIR / "testcases/fhir"
SORT_CASES = TEST_DIR / "testcases/sort"


def test_version():
    assert __version__ != "unknown"


def idfn(val):
    if isinstance(val, os.PathLike):
        return "-".join(os.path.normpath(val).split(os.sep)[-2:])


class Context(Enum):
    FILE = auto()
    STATEMENT = auto()
    NONE = auto()
    ERROR = auto()


def run_test_filter_tag_sort(args, config, path, tmp_out, backend):
    """Handles verification for `filter`, `tag`, and `sort` commands."""
    match args.command:
        case "filter":
            module = filter
        case "tag":
            module = tag
        case "sort":
            module = sort
        case _:
            raise ValueError(f"Unknown VCF command: {args.command}")

    expected = str(path / "expected.vcf")
    module.execute(args)

    with create_reader(tmp_out.name, backend=backend) as vcf_actual:
        with create_reader(expected, backend=backend) as vcf_expected:
            for r1, r2 in zip_longest(vcf_actual, vcf_expected):
                assert r1 == r2
            assert vcf_actual.header.get_generic("vembraneVersion") == __version__


def run_test_structured_fhir(args, config, path, tmp_out, backend):
    """Handles verification for `structured` and `fhir` commands."""
    expected_path = (path / "expected").with_suffix(f".{config['output_fmt']}")

    if args.command == "structured":
        structured.execute(args)
    else:
        fhir.execute(args)

    with open(expected_path) as e, open(tmp_out.name) as t:
        if config["output_fmt"] == "jsonl":
            t_out = [json.loads(line) for line in t]
            e_out = [json.loads(line) for line in e]
        elif config["output_fmt"] == "json":
            t_out = json.load(t)
            e_out = json.load(e)
        elif config["output_fmt"] == "yaml":
            t_out = yaml.safe_load(t)
            e_out = yaml.safe_load(e)
        else:
            raise ValueError(f"Unknown output format: {config['output_fmt']}")

    assert t_out == e_out


def run_test_table(args, config, path, tmp_out, backend):
    expected_path = str(path / "expected.tsv")

    args.output_fmt = "csv"
    table.execute(args)

    def get_stripped_content(filename):
        with open(filename) as f:
            return "".join(line for line in f if not line.startswith("##vembrane"))

    actual_csv = get_stripped_content(tmp_out.name)
    with open(expected_path) as f:
        expected_csv = f.read()

    assert actual_csv == expected_csv

    args.output_fmt = "parquet"
    with open(tmp_out.name, "wb") as f:
        f.truncate(0)

    table.execute(args)
    check_parquet_output(args, expected_path, tmp_out.name)


def check_parquet_output(args, expected_path, tmp_out_path):
    has_header = args.header != "none"

    with open(expected_path, newline="") as f:
        reader = csv.reader(f, delimiter=args.separator)
        if has_header:
            expected_header = next(reader)
        expected_rows = list(reader)

    t = pq.read_table(tmp_out_path)
    actual_header = t.column_names
    actual_rows_raw = t.to_pylist()

    if has_header:
        expected_header = [c.replace('""', '"') for c in expected_header]
        assert actual_header == expected_header
    else:
        # If no header, ensure column counts match
        assert len(actual_header) == len(expected_rows[0]) if expected_rows else True

    assert len(actual_rows_raw) == len(expected_rows)

    for i, (row_raw, row_expected) in enumerate(
        zip(actual_rows_raw, expected_rows, strict=False)
    ):
        if has_header:
            row_vals = [row_raw[col] for col in actual_header]
        else:
            row_vals = list(row_raw.values())

        for val_arrow, val_expected in zip(row_vals, row_expected, strict=False):
            norm_val = normalize_arrow_value(val_arrow)

            candidates = {str(norm_val)}
            if isinstance(norm_val, list):
                candidates.add(str(tuple(norm_val)))
                str_items = [str(v) for v in norm_val]
                candidates.add("&".join(str_items))

            assert val_expected in candidates, (
                f"Row {i} mismatch: Expected '{val_expected}', "
                f"but got Parquet value formatting to {candidates}"
            )


TEST_RUNNERS: dict[str, Callable] = {
    "filter": run_test_filter_tag_sort,
    "tag": run_test_filter_tag_sort,
    "sort": run_test_filter_tag_sort,
    "structured": run_test_structured_fhir,
    "fhir": run_test_structured_fhir,
    "table": run_test_table,
}


@pytest.mark.parametrize(
    "testcase,backend,context",
    product(
        (
            case_path.joinpath(d)
            for case_path in [
                FILTER_CASES,
                TABLE_CASES,
                TAG_CASES,
                ANNOTATE_CASES,
                STRUCTURED_CASES,
                FHIR_CASES,
                SORT_CASES,
            ]
            for d in os.listdir(case_path)
            if not d.startswith(".")
        ),
        (Backend.pysam, Backend.cyvcf2),
        (Context.NONE, Context.FILE, Context.STATEMENT, Context.ERROR),
    ),
    ids=idfn,
)
def test_command(testcase: Path, backend: Backend, context: Context | None):
    path = testcase
    vcf_path = path.joinpath("test.vcf")

    with open(path.joinpath("config.yaml")) as config_fp:
        config = yaml.load(config_fp, Loader=yaml.FullLoader)

    # emulate command-line command setup to use argparse afterwards
    cmd = config["function"]
    command = parse_command_config(cmd, config, vcf_path)

    match context:
        case Context.STATEMENT:
            command.extend(["--context", "import random; dummy_func = lambda x: x"])
        case Context.FILE:
            command.extend(["--context-file", "tests/resources/context.py"])
        case Context.ERROR:
            if "raises" not in config:
                command.extend(["--context-file", "tests/resources/context_error.py"])
                config["raises"] = "VembraneError"
            else:
                pytest.skip(
                    "Test expects other error, skipping context error injection."
                )
        case Context.NONE:
            if config.get("needs_context"):
                pytest.skip("Test needs context, skipping.")

    parser = construct_parser()

    try:
        args, unknown = parser.parse_known_args(command)
    # we have to catch and compare such exceptions here, because they
    # do not cause a SystemExit like the .execute() instances below
    except Exception as e:
        assert type(e).__name__ == config["raises"]
        return

    if "raises" in config:
        exception = config["raises"]
        try:
            exception = getattr(errors, exception)
        except AttributeError:
            exception = getattr(builtins, exception)

        with pytest.raises(SystemExit), pytest.raises(exception):
            if args.command == "filter":
                filter.execute(args)
            elif args.command == "table":
                table.execute(args)
            elif args.command == "tag":
                tag.execute(args)
            elif args.command == "structured":
                structured.execute(args)
            elif args.command == "fhir":
                fhir.execute(args)
            elif args.command == "sort":
                sort.execute(args)
            else:
                raise AssertionError("Unknown subcommand")
        return

    with tempfile.NamedTemporaryFile(mode="w+t") as tmp_out:
        args.output = tmp_out.name
        args.backend = backend

        runner = TEST_RUNNERS.get(args.command)
        if runner:
            runner(args, config, path, tmp_out, backend)
        else:
            raise ValueError(f"No verifier defined for command: {args.command}")


def construct_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest="command",
        description="valid subcommands",
        required=True,
    )
    filter.add_subcommmand(subparsers)
    table.add_subcommmand(subparsers)
    tag.add_subcommand(subparsers)
    structured.add_subcommand(subparsers)
    fhir.add_subcommand(subparsers)
    sort.add_subcommand(subparsers)
    return parser


def parse_command_config(cmd, config, vcf_path):
    if cmd in ("filter", "table"):
        command = [cmd, config["expression"], str(vcf_path)]
    elif cmd == "tag":
        command = [cmd]
        tags = config["tags"]
        for name, expr in tags.items():
            command += ["--tag", f"{name}={expr}"]
        command.append(str(vcf_path))
    elif cmd == "structured":
        template_path = vcf_path.parent / "template.yte.yaml"
        command = [cmd, str(template_path), str(vcf_path)]
    elif cmd == "fhir":
        command = [cmd, str(vcf_path), config["sample"], config["assembly"]]
        del config["sample"]
        del config["assembly"]
    elif cmd == "sort":
        command = [cmd, config["expression"], str(vcf_path)]
        del config["expression"]
    else:
        raise ValueError(f"Unknown subcommand {config['function']}")
    for key in config:
        if key == "function":
            continue
        if isinstance(config[key], str):
            command.append(config_key_to_arg(key))
            command.append(config[key])
        elif isinstance(config[key], bool):
            if config[key]:
                command.append(config_key_to_arg(key))
        elif isinstance(config[key], (int, float)):
            command.append(config_key_to_arg(key))
            command.append(str(config[key]))
        else:
            if isinstance(config[key], dict):
                for k, v in config[key].items():
                    command.append(config_key_to_arg(key))
                    command.append(f"{k}={v}")
            else:
                command.append(config_key_to_arg(key))
                for argument in config[key]:
                    command.append(argument)
    return command


def config_key_to_arg(key: str) -> str:
    """Convert a config key to an argument name."""
    return f"--{key.replace('_', '-')}"


def normalize_arrow_value(value: Any) -> Any:
    if value is None:
        return ""

    if isinstance(value, dict):
        keys = set(value.keys())

        # RangeTotal
        if keys == {"start", "stop", "total"}:
            start, stop, total = value["start"], value["stop"], value["total"]
            length = stop - start
            if length == 1:
                return f"number / total: {start} / {total}"
            else:
                return f"range / total: {start} - {stop - 1} / {total}"

        # NumberTotal
        if keys == {"number", "total"}:
            return f"number / total: {value['number']} / {value['total']}"

        # PosRange
        if keys == {"start", "end"}:
            start, end = value["start"], value["end"]
            if start is None or end is None:
                s_str = start if start is not None else ""
                e_str = end if end is not None else ""
                return f"(start: {s_str}, end: {e_str}, length: )"
            return f"(start: {start}, end: {end}, length: {end - start})"

        return str(value)

    if isinstance(value, list):
        return [normalize_arrow_value(v) for v in value]

    return value
