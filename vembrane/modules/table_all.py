import csv
import sys
from collections.abc import Iterator
from enum import Enum
from sys import stderr
from typing import Any

from ..backend.base import VCFHeader, VCFReader, VCFRecord
from ..common import (
    add_common_arguments,
    create_reader,
    get_annotation_keys,
    smart_open,
)
from ..errors import VembraneError
from ..representations import FuncWrappedExpressionEnvironment
from .filter import DeprecatedAction


class NamingConvention(Enum):
    DICTIONARY = "dictionary"
    UNDERSCORE = "underscore"
    SLASH = "slash"


def add_subcommand(subparsers):
    parser = subparsers.add_parser("table-all")
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "vcf",
        help="The file containing the variants.",
        nargs="?",
        default="-",
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--separator",
        "-s",
        default="\t",
        metavar="CHAR",
        help="Define the field separator (default: \\t).",
    )
    parser.add_argument(
        "--naming-convention",
        metavar="CONVENTION",
        default="dictionary",
        type=NamingConvention,
        choices=list(NamingConvention),
        help="The naming convention to use for column names.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )
    add_common_arguments(parser)


def tableize_vcf(
    vcf: VCFReader,
    ann_key: str,
    naming_convention: NamingConvention = NamingConvention.DICTIONARY,
    overwrite_number: dict[str, dict[str, str]] | None = None,
) -> Iterator[tuple]:
    if overwrite_number is None:
        overwrite_number = {}

    kwargs: dict[str, Any] = dict()
    default_expression, _header = construct_expression_and_header(
        vcf.header, ann_key, naming_convention
    )
    expression = f"(({default_expression}) for SAMPLE in SAMPLES)"
    env = FuncWrappedExpressionEnvironment(expression, ann_key, vcf.header, **kwargs)

    record: VCFRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.expression_annotations():
            annotations = env.get_record_annotations(idx, record)
            for annotation in annotations:
                yield from env.table_row(annotation)
        else:
            yield from env.table_row()


def construct_expression_and_header(
    header: VCFHeader, ann_key: str, naming_convention: NamingConvention
) -> tuple[str, list[str]]:
    default_header = ["SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "ID"]
    default_expression = ", ".join(default_header)

    def format_key(field: str, key: str) -> str:
        match naming_convention:
            case NamingConvention.DICTIONARY:
                return f"{field}['{key}']"
            case NamingConvention.UNDERSCORE:
                return f"{field}_{key}"
            case NamingConvention.SLASH:
                return f"{field}/{key}"

    info_keys = list(header.infos.keys())
    if info_keys:
        info_keys_without_ann = [key for key in info_keys if key != ann_key]
        info_expr = ", ".join(f'INFO["{key}"]' for key in info_keys_without_ann)
        default_header += [format_key("INFO", key) for key in info_keys_without_ann]
        default_expression += ", " + info_expr

    format_keys = list(header.formats.keys())
    if format_keys:
        default_header += [format_key("FORMAT", key) for key in format_keys]
        format_expr = ", ".join(
            f'(FORMAT.get("{key}") or {{}}).get(SAMPLE, NA)' for key in format_keys
        )
        default_expression += ", " + format_expr

    annotation_keys = get_annotation_keys(header, ann_key)
    if annotation_keys:
        annotation_expr = ", ".join(f'ANN["{key}"]' for key in annotation_keys)
        default_header += [format_key(ann_key, key) for key in annotation_keys]
        default_expression += ", " + annotation_expr

    return default_expression, default_header


def get_row(row):
    if not isinstance(row, tuple):
        row = (row,)
    return row


def execute(args):
    overwrite_number = {
        "INFO": dict(args.overwrite_number_info),
        "FORMAT": dict(args.overwrite_number_format),
    }
    with create_reader(
        args.vcf,
        backend=args.backend,
        overwrite_number=overwrite_number,
    ) as vcf:
        _expr, header = construct_expression_and_header(
            vcf.header, args.annotation_key, args.naming_convention
        )
        rows = tableize_vcf(
            vcf,
            args.annotation_key,
            naming_convention=args.naming_convention,
            overwrite_number=overwrite_number,
        )

        try:
            with smart_open(args.output, "wt", newline="") as csvfile:
                writer = csv.writer(
                    csvfile,
                    delimiter=args.separator,
                    quoting=csv.QUOTE_MINIMAL,
                )
                writer.writerow(header)
                writer.writerows(get_row(row) for row in rows)
        except VembraneError as ve:
            print(ve, file=stderr)
            sys.exit(1)
