import argparse
import builtins
import json
import os
import tempfile
from enum import Enum, auto
from itertools import product, zip_longest
from pathlib import Path

import pytest
import yaml

from vembrane import __version__, errors
from vembrane.backend.base import Backend
from vembrane.common import create_reader
from vembrane.modules import fhir, filter, sort, structured, table, tag

FILTER_CASES = Path(__file__).parent.joinpath("testcases/filter")
TABLE_CASES = Path(__file__).parent.joinpath("testcases/table")
TAG_CASES = Path(__file__).parent.joinpath("testcases/tag")
ANNOTATE_CASES = Path(__file__).parent.joinpath("testcases/annotate")
STRUCTURED_CASES = Path(__file__).parent.joinpath("testcases/structured")
FHIR_CASES = Path(__file__).parent.joinpath("testcases/fhir")
SORT_CASES = Path(__file__).parent.joinpath("testcases/sort")


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
def test_command(testcase: os.PathLike, backend: Backend, context: Context | None):
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
                raise AssertionError from None
    else:
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_out:
            args.output = tmp_out.name
            args.backend = backend
            if (
                args.command == "filter"
                or args.command == "tag"
                or args.command == "sort"
            ):
                match args.command:
                    case "filter":
                        module = filter
                    case "tag":
                        module = tag
                    case "sort":
                        module = sort
                expected = str(path.joinpath("expected.vcf"))
                module.execute(args)
                with create_reader(tmp_out.name, backend=backend) as vcf_actual:
                    with create_reader(expected, backend=backend) as vcf_expected:
                        for r1, r2 in zip_longest(vcf_actual, vcf_expected):
                            assert r1 == r2

                        assert (
                            vcf_actual.header.get_generic("vembraneVersion")
                            == __version__
                        )
            elif args.command == "table":
                expected = str(path.joinpath("expected.tsv"))
                table.execute(args)
                t_out = "".join(
                    line for line in tmp_out if not line.startswith("##vembrane")
                )
                with open(expected) as e:
                    e_out = e.read()
                assert t_out == e_out
            elif args.command == "structured" or args.command == "fhir":
                expected = (path / "expected").with_suffix(f".{config['output_fmt']}")

                if args.command == "structured":
                    structured.execute(args)
                else:
                    fhir.execute(args)

                with open(expected) as e:
                    if config["output_fmt"] == "jsonl":
                        t_out = [json.loads(line) for line in tmp_out]
                        e_out = [json.loads(line) for line in e]
                    elif config["output_fmt"] == "json":
                        t_out = json.load(tmp_out)
                        e_out = json.load(e)
                    elif config["output_fmt"] == "yaml":
                        t_out = yaml.safe_load(tmp_out)
                        e_out = yaml.safe_load(e)
                assert t_out == e_out
            else:
                assert args.command in {
                    "filter",
                    "table",
                    "tag",
                    "structured",
                    "fhir",
                    "sort",
                }, "Unknown subcommand"


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
