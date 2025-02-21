import argparse
import builtins
import json
import os
import tempfile
from itertools import product, zip_longest
from pathlib import Path

import pytest
import yaml

from vembrane import __version__, errors
from vembrane.backend.base import Backend
from vembrane.common import create_reader
from vembrane.modules import filter, structured, table, tag

FILTER_CASES = Path(__file__).parent.joinpath("testcases/filter")
TABLE_CASES = Path(__file__).parent.joinpath("testcases/table")
TAG_CASES = Path(__file__).parent.joinpath("testcases/tag")
ANNOTATE_CASES = Path(__file__).parent.joinpath("testcases/annotate")
STRUCTURED_CASES = Path(__file__).parent.joinpath("testcases/structured")


def test_version():
    assert __version__ != "unknown"


def idfn(val):
    if isinstance(val, os.PathLike):
        return "-".join(os.path.normpath(val).split(os.sep)[-2:])


@pytest.mark.parametrize(
    "testcase,backend",
    product(
        (
            case_path.joinpath(d)
            for case_path in [
                FILTER_CASES,
                TABLE_CASES,
                TAG_CASES,
                ANNOTATE_CASES,
                STRUCTURED_CASES,
            ]
            for d in os.listdir(case_path)
            if not d.startswith(".")
        ),
        (Backend.pysam, Backend.cyvcf2),  # Backend.cyvcf2
    ),
    ids=idfn,
)
def test_command(testcase: os.PathLike, backend: Backend):
    path = testcase
    vcf_path = path.joinpath("test.vcf")

    with open(path.joinpath("config.yaml")) as config_fp:
        config = yaml.load(config_fp, Loader=yaml.FullLoader)

    # emulate command-line command setup to use argparse afterwards
    cmd = config["function"]
    command = parse_command_config(cmd, config, vcf_path)
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
            with pytest.raises(SystemExit), pytest.raises(exception):
                if args.command == "filter":
                    filter.execute(args)
                elif args.command == "table":
                    table.execute(args)
                elif args.command == "tag":
                    tag.execute(args)
                else:
                    raise AssertionError from None
        except AttributeError:
            exception = getattr(builtins, exception)
            with pytest.raises(exception):
                if args.command == "filter":
                    filter.execute(args)
                elif args.command == "table":
                    table.execute(args)
                elif args.command == "tag":
                    tag.execute(args)
                else:
                    raise AssertionError from None
    else:
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_out:
            args.output = tmp_out.name
            args.backend = backend
            if args.command == "filter" or args.command == "tag":
                module = filter if args.command == "filter" else tag
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
            elif args.command == "structured":
                expected = (path / "expected").with_suffix(f".{config['output_fmt']}")
                structured.execute(args)

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
    else:
        raise ValueError(f"Unknown subcommand {config['function']}")
    for key in config:
        if isinstance(config[key], str):
            command.append(f"--{key.replace('_', '-')}")
            command.append(config[key])
        else:
            if isinstance(config[key], dict):
                for k, v in config[key].items():
                    command.append(f"--{key.replace('_', '-')}")
                    command.append(f"{k}={v}")
            else:
                command.append(f"--{key.replace('_', '-')}")
                for argument in config[key]:
                    command.append(argument)
    return command
