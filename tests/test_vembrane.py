import argparse
import builtins
import os
import tempfile
from itertools import product, zip_longest
from pathlib import Path

import pytest
import yaml

from vembrane import __version__, errors
from vembrane.backend.base import Backend
from vembrane.common import create_reader
from vembrane.modules import filter, table, tag

CASES = Path(__file__).parent.joinpath("testcases")


def test_version():
    assert __version__ != "unknown"


@pytest.mark.parametrize(
    "testcase,backend",
    product(
        (d for d in os.listdir(CASES) if not d.startswith(".")),
        (Backend.pysam, ),  # TODO: include Backend.cyvcf2
    ),
)
def test_filter(testcase: os.PathLike, backend: Backend):
    path = CASES.joinpath(testcase)

    vcf_path = path.joinpath("test.vcf")

    with open(path.joinpath("config.yaml")) as config_fp:
        config = yaml.load(config_fp, Loader=yaml.FullLoader)

    # emulate command-line command setup to use argparse afterwards
    cmd = config["function"]
    if cmd in ("filter", "table"):
        command = [cmd, config["expression"], str(vcf_path)]
    elif cmd == "tag":
        command = [cmd]
        tags = config["tags"]
        for name, expr in tags.items():
            command += ["--tag", f"{name}={expr}"]
        command.append(str(vcf_path))
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

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest="command", description="valid subcommands", required=True
    )
    filter.add_subcommmand(subparsers)
    table.add_subcommmand(subparsers)
    tag.add_subcommand(subparsers)
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
            with pytest.raises(SystemExit):
                with pytest.raises(exception):
                    if args.command == "filter":
                        filter.execute(args)
                    elif args.command == "table":
                        table.execute(args)
                    elif args.command == "tag":
                        tag.execute(args)
                    else:
                        assert False
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
                    assert False
    else:
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_out:
            args.output = tmp_out.name
            if args.command == "filter" or args.command == "tag":
                module = filter if args.command == "filter" else tag
                expected = str(path.joinpath("expected.vcf"))
                module.execute(args)
                with create_reader(tmp_out.name, backend=backend) as vcf_actual:
                    with create_reader(expected, backend=backend) as vcf_expected:
                        for r1, r2 in zip_longest(vcf_actual, vcf_expected):
                            assert r1 == r2

                        vcf_actual.header.records == vcf_expected.header.records
                        # TODO: the records dont hold the generic entries anymore
                        # implement the vembraneVersion check
            elif args.command == "table":
                expected = str(path.joinpath("expected.tsv"))
                table.execute(args)
                t_out = "".join(
                    line for line in tmp_out if not line.startswith("##vembrane")
                )
                with open(expected, mode="r") as e:
                    e_out = e.read()
                assert t_out == e_out
            else:
                assert args.command in {"filter", "table"}, "Unknown subcommand"
