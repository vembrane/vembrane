import contextlib
import csv
import sys

from pysam.libcbcf import VariantFile, VariantRecord
from typing import Iterator, List, Optional
from sys import stderr

from ..common import check_expression
from ..errors import VembraneError
from ..representations import Environment

import ast


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("table")
    parser.add_argument(
        "expression",
        type=check_expression,
        help="The expression for the output.",
    )
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
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
        "--all",
        "-a",
        default=False,
        action="store_true",
        help="Do not filter duplicate entries.",
    )
    parser.add_argument(
        "--header",
        default="auto",
        metavar="TEXT",
        help='Override the automatically generated header. Provide "auto" (default) \
              to automatically generate the header from the expression. Provide a \
              comma separated string to manually set the header. Provide "none" to \
              disable any header output.',
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )


def tableize_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str,
) -> Iterator[tuple]:
    expression = f"({expression})"
    env = Environment(expression, ann_key, vcf.header)

    record: VariantRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.expression_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            try:
                annotations = record.info[ann_key]
            except KeyError:
                annotations = [""]
            for annotation in annotations:
                yield env.table(annotation)
        else:
            yield env.table()


def get_header(args) -> Optional[List[str]]:
    if args.header == "none":
        return
    elif args.header == "auto":
        header = args.expression
    else:
        header = args.header
    if "," not in header:
        return [header]
    else:
        # use the nodes of the first layer of the ast tree as header names
        elts = ast.parse(header).body[0].value.elts
        return [header[n.col_offset : n.end_col_offset] for n in elts]


def get_row(row):
    if not isinstance(row, tuple):
        row = (row,)
    return row


@contextlib.contextmanager
def smart_open(filename=None, *args, **kwargs):
    if filename and filename != "-":
        fh = open(filename, *args, **kwargs)
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def execute(args):
    with VariantFile(args.vcf) as vcf:
        rows = tableize_vcf(
            vcf,
            args.expression,
            args.annotation_key,
        )

        try:
            with smart_open(args.output, "wt", newline="") as csvfile:
                writer = csv.writer(csvfile, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
                if args.header:
                    writer.writerow(get_header(args))
                writer.writerows(get_row(row) for row in rows)
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)
