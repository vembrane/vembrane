import re
import sys
from itertools import chain, islice
from sys import stderr
from typing import Dict, Iterator, Set, Tuple

from pysam.libcbcf import VariantFile, VariantHeader, VariantRecord

from .. import __version__
from ..common import (
    AppendKeyValuePair,
    AppendTagExpression,
    check_expression,
    normalize,
    read_auxiliary,
    single_outer,
    swap_quotes,
)
from ..errors import FilterAlreadyDefined, FilterTagNameInvalid, VembraneError
from ..representations import Environment
from .filter import DeprecatedAction


def add_subcommand(subparsers):
    parser = subparsers.add_parser("tag")
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "--tag",
        "-t",
        help="Tag records using the FILTER field.\n"
        "Note: tag names cannot contain `;` or whitespace and must not be '0'.\n"
        'Example: `--tag q_above_30="not (QUAL<=30)"`',
        nargs=1,
        metavar="TAG=EXPRESSION",
        action=AppendTagExpression,
        required=True,
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
        nargs=1,
        action=AppendKeyValuePair,
        metavar="NAME=PATH",
        default={},
        help="Path to an auxiliary file containing a set of symbols",
    )
    parser.add_argument(
        "--overwrite-number-info",
        nargs=2,
        action="append",
        metavar=("FIELD", "NUMBER"),
        default=[],
        help="Overwrite the number specification for INFO fields "
        "given in the VCF header. "
        "Example: `--overwrite-number cosmic_CNT=.`",
    )
    parser.add_argument(
        "--overwrite-number-format",
        nargs=2,
        action="append",
        metavar=("FIELD", "NUMBER"),
        default=[],
        help="Overwrite the number specification for FORMAT fields "
        "given in the VCF header. "
        "Example: `--overwrite-number-format DP=2`",
    )
    parser.add_argument(
        "--tag-mode",
        "-m",
        default="pass",
        choices=["pass", "fail"],
        help="Set, whether to tag records that pass the tag expression(s), "
        "or records that fail them."
        "By default, vembrane tags records for which the tag expression(s) pass. "
        "This allows for descriptive tag names such as `q_at_least_30`, "
        "which would correspond to an expression `QUAL >= 30`. "
        "However, the VCF specification (`v4.4`) defines tags to be set when a "
        "filter expression is failed, so vembrane also offers the `fail` mode.",
    )


def test_record(
    env: Environment,
    idx: int,
    record: VariantRecord,
    ann_key: str,
) -> Tuple[VariantRecord, bool]:
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
    invert: bool = False,
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
            if invert:
                keep = not keep
            if keep:
                record.filter.add(tag)
        yield record


def check_tag(tag: str):
    if re.search(r"^0$|[\s;]", tag):
        raise FilterTagNameInvalid(tag)


def execute(args) -> None:
    aux = read_auxiliary(args.aux)
    with VariantFile(args.vcf) as vcf:
        header: VariantHeader = vcf.header

        overwrite_number = {
            "INFO": dict(args.overwrite_number_info),
            "FORMAT": dict(args.overwrite_number_format),
        }

        expressions = dict(args.tag)
        for tag, expr in expressions.items():
            for t, rec in vcf.header.filters.items():
                if t == tag:
                    e = FilterAlreadyDefined(tag)
                    print(e, file=stderr)
                    exit(1)
            try:
                check_tag(tag)
            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
            expr = swap_quotes(expr) if single_outer(expr) else expr
            check_expression(expr)
            vcf.header.add_meta(
                key="FILTER", items=[("ID", tag), ("Description", expr)]
            )

        header.add_meta("vembraneVersion", __version__)
        header.add_meta(
            "vembraneCmd",
            "vembrane "
            + " ".join(normalize(arg) if " " in arg else arg for arg in sys.argv[1:]),
        )

        records = tag_vcf(
            vcf,
            expressions,
            args.annotation_key,
            auxiliary=aux,
            overwrite_number=overwrite_number,
            invert=(args.tag_mode == "fail"),
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
