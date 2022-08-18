import argparse
import sys
from itertools import chain, islice
from sys import stderr
from typing import Dict, Iterator, Set

from pysam.libcbcf import VariantFile, VariantHeader, VariantRecord

from .. import __version__
from ..common import check_expression
from ..errors import FilterAlreadyDefined, VembraneError
from ..representations import Environment
from .filter import DeprecatedAction, read_auxiliary


class TagExpressionList(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        values = values[0].strip().split(",")
        for value in values:
            key, value = value.strip().split("=")
            expr = check_expression(value)
            getattr(namespace, self.dest)[key] = expr


def add_subcommand(subparsers):
    parser = subparsers.add_parser("tag")
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "expression",
        help="Tag records using the FILTER field.",
        nargs=1,
        metavar="[TAG=EXPRESSION,]+",
        action=TagExpressionList,
    )
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )
    parser.add_argument(
        "--output-fmt",
        "-O",
        default="vcf",
        choices=["vcf", "bcf", "uncompressed-bcf"],
        help="Output format.",
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--aux",
        "-a",
        nargs=2,
        action="append",
        metavar=("NAME", "PATH"),
        default=[],
        help="Path to an auxiliary file containing a set of symbols",
    )
    parser.add_argument(
        "--overwrite-number",
        help="Deprecated. "
        "Use --overwrite-number-info or --overwrite-number-format instead.",
        action="deprecated",
    )
    parser.add_argument(
        "--overwrite-number-info",
        nargs=2,
        action="append",
        metavar=("FIELD", "NUMBER"),
        default=[],
        help="Overwrite the number specification for INFO fields "
        "given in the VCF header. "
        "Example: `--overwrite-number cosmic_CNT .`",
    )
    parser.add_argument(
        "--overwrite-number-format",
        nargs=2,
        action="append",
        metavar=("FIELD", "NUMBER"),
        default=[],
        help="Overwrite the number specification for FORMAT fields "
        "given in the VCF header. "
        "Example: `--overwrite-number-format DP 2`",
    )


def test_record(
    env: Environment,
    idx: int,
    record: VariantRecord,
    ann_key: str,
) -> (VariantRecord, bool):
    env.update_from_record(idx, record)
    if env.expression_annotations():
        # if the expression contains a reference to the ANN field
        # get all annotations from the record.info field
        # (or supply an empty ANN value if the record has no ANN field)
        try:
            annotations = record.info[ann_key]
        except KeyError:
            annotations = [""]

        #  … and check if the expression evaluates to true for any  of the annotations
        filtered = any(map(env.evaluate, annotations))
        return record, filtered
    else:
        # otherwise, the annotations are irrelevant w.r.t. the expression,
        # so we can omit them
        return record, env.evaluate()


def tag_vcf(
    vcf: VariantFile,
    expressions: Dict[str, str],
    ann_key: str,
    auxiliary: Dict[str, Set[str]] = {},
    overwrite_number: Dict[str, Dict[str, str]] = {},
) -> Iterator[VariantRecord]:

    # For each tag-expression pair, a different Environment must be used.
    envs = {
        tag: Environment(expression, ann_key, vcf.header, auxiliary, overwrite_number)
        for tag, expression in expressions.items()
    }

    record: VariantRecord
    for idx, record in enumerate(vcf):
        for tag, env in envs.items():
            record, keep = test_record(env, idx, record, ann_key)
            if keep:
                record.filter.add(tag)
        yield record


def execute(args):
    aux = read_auxiliary(args.aux)
    with VariantFile(args.vcf) as vcf:
        header: VariantHeader = vcf.header
        header.add_meta("vembraneVersion", __version__)
        # NOTE: If .modules.filter.execute might be used as a library function
        #       in the future, we should not record sys.argv directly below.
        header.add_meta(
            "vembraneCmd",
            "vembrane "
            + " ".join(
                "'" + arg.replace("'", '"') + '"' if " " in arg else arg
                for arg in sys.argv[1:]
            ),
        )

        overwrite_number = {
            "INFO": dict(args.overwrite_number_info),
            "FORMAT": dict(args.overwrite_number_format),
        }

        expressions = dict(args.expression)
        for tag, expr in expressions.items():
            for t, rec in vcf.header.filters.items():
                if t == tag:
                    e = FilterAlreadyDefined(tag)
                    print(e, file=stderr)
                    exit(1)
            vcf.header.add_meta(
                key="FILTER", items=[("ID", tag), ("Description", expr)]
            )

        records = tag_vcf(
            vcf,
            expressions,
            args.annotation_key,
            auxiliary=aux,
            overwrite_number=overwrite_number,
        )

        try:
            first_record = list(islice(records, 1))
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)

        records = chain(first_record, records)
        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
        with VariantFile(
            args.output,
            f"w{fmt}",
            header=header,
        ) as out:
            try:
                for record in records:
                    out.write(record)

            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
