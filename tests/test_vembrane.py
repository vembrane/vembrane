import argparse
import builtins
from itertools import zip_longest
import os
from pathlib import Path
import tempfile

import pytest
import yaml
from cyvcf2 import VCF

from vembrane import errors, __version__
from vembrane.modules import filter, table


CASES = Path(__file__).parent.joinpath("testcases")


def test_version():
    assert __version__ != "unknown"


@pytest.mark.parametrize(
    "testcase", [d for d in os.listdir(CASES) if not d.startswith(".")]
)
def test_filter(testcase):
    path = CASES.joinpath(testcase)

    vcf_path = path.joinpath("test.vcf")

    with open(path.joinpath("config.yaml")) as config_fp:
        config = yaml.load(config_fp, Loader=yaml.FullLoader)

    # emulate command-line command setup to use argparse afterwards
    command = [config["function"], config["expression"], str(vcf_path)]
    for key in config:
        command.append(f"--{key.replace('_','-')}")
        if isinstance(config[key], str):
            command.append(config[key])
        else:
            for argument in config[key]:
                command.append(argument)

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest="command", description="valid subcommands", required=True
    )
    filter.add_subcommmand(subparsers)
    table.add_subcommmand(subparsers)
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
                    else:
                        assert False
        except AttributeError:
            exception = getattr(builtins, exception)
            with pytest.raises(exception):
                if args.command == "filter":
                    filter.execute(args)
                elif args.command == "table":
                    table.execute(args)
                else:
                    assert False
    else:
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_out:
            args.output = tmp_out.name
            if args.command == "filter":
                expected = str(path.joinpath("expected.vcf"))
                filter.execute(args)
                vcf_actual = VCF(tmp_out.name)
                vcf_expected = VCF(expected)
                for r1, r2 in zip_longest(vcf_actual, vcf_expected):
                    # FIXME cyvcf2 doesn't properly implement __eq__ for all attributes
                    # which is why we're using the string representation here instead
                    s1 = str(r1)
                    s2 = str(r2)
                    assert s1 == s2

                for r1, r2 in zip_longest(
                    vcf_actual.raw_header.split("\n"),
                    vcf_expected.raw_header.split("\n"),
                ):
                    if r1.startswith("##vembraneVersion="):
                        continue
                        # assert r1[18:] == __version__
                    elif r1.startswith("##vembraneCmd="):
                        assert r1[14:].startswith("vembrane ")
                    else:
                        continue
                        # assert r1 == r2

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
