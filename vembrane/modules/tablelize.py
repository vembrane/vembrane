from pysam.libcbcf import VariantFile
from typing import Iterator

from .common import check_filter_expression
from ..errors import VembraneError
from ..representations import Environment


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


def tablelize_vcf(
    vcf: VariantFile, expression: str, ann_key: str, keep_unmatched: bool = False,
) -> Iterator[tuple]:

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


def execute(args):
    with VariantFile(args.vcf) as vcf:
        rows = tablelize_vcf(vcf, args.expression, args.annotation_key,)

        try:
            for row in rows:
                print(row)
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)
