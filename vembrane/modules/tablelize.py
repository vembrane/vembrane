from pysam.libcbcf import VariantFile, VariantRecord
from typing import Iterator
from sys import stderr

from ..common import check_filter_expression
from ..errors import VembraneError
from ..representations import Environment

import ast


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("tablelize")
    parser.add_argument(
        "expression",
        type=check_filter_expression,
        help="Filter variants and annotations. If this removes all annotations, "
        "the variant is removed as well.",
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
        help='Override the automatically generated header. Provide "auto" (default) to automatically generate the header from the expression. Provide a comma separated string to manually set the header. Provide "none" to disable any header output.',
    )


def tablelize_vcf(
    vcf: VariantFile, expression: str, ann_key: str, keep_unmatched: bool = False,
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
                yield env.tablelize(annotation)
        else:
            yield env.tablelize()


def print_header(args):
    if args.header == "none":
        return
    elif args.header == "auto":
        header = args.expression
    else:
        header = args.header

    if not "," in header:
        print(f"#{header}")
    else:
        # print the nodes of the first layer of the ast tree as header names
        elts = ast.parse(header).body[0].value.elts
        header_fields = [header[n.col_offset : n.end_col_offset] for n in elts]
        header_string = "#" + args.separator.join(map(str.strip, header_fields))
        print(header_string)


def print_row(row, args):
    if not isinstance(row, tuple):
        row = (row,)
    out_string = args.separator.join(map(str, row))
    if args.all or out_string != print_row.last_string:
        print(*row, sep=args.separator)
        print_row.last_string = out_string


print_row.last_string = None


def execute(args):
    with VariantFile(args.vcf) as vcf:
        rows = tablelize_vcf(vcf, args.expression, args.annotation_key,)

        try:
            if args.header:
                print_header(args)
            for row in rows:
                print_row(row, args)
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)
