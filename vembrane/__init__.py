__version__ = "0.1.0"

# import stuff we want to be available in eval by default:
import argparse
import math
import re
from inspect import getmembers, isbuiltin
from typing import Iterator

from pysam import VariantFile, VariantRecord

globals_whitelist = {
    **{
        "__builtins__": None,
        "__builtin__": None,
        "__file__": None,
        "__name__": None,
        "__doc__": None,
        "__package__": None,
    },
    **{mod.__name__: mod for mod in [any, all, min, max, re, list, dict, zip]},
    **dict(getmembers(math, predicate=isbuiltin)),
}


def filter_vcf(
    vcf: VariantFile, expression: str, vcf_info_field: str = "ANN"
) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    for name in header.info:
        env[name] = None

    annotation_keys = []
    info_field_found = False
    for rec in header.records:
        if rec.get("ID") == vcf_info_field:
            annotation_keys = list(
                map(str.strip, rec.get("Description").split("'")[1].split("|"))
            )
            info_field_found = True
            break

    if not info_field_found:
        raise ValueError(
            f'VCF info field "{vcf_info_field}" not found in header'
        )

    for record in vcf:
        for key in record.info:
            env[key] = record.info[key]
        ann = env.get(vcf_info_field, [])
        env["ANNO"] = dict(
            zip(
                annotation_keys,
                zip(*(list(map(str.strip, a.split("|"))) for a in ann)),
            )
        )
        if eval(expression, globals_whitelist, env):
            yield record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="The file containing the variants.")
    parser.add_argument(
        "expression", help="An expression to filter the variants."
    )
    parser.add_argument(
        "--vcf_info_field",
        metavar="key",
        default="ANN",
        help="The INFO key VEMBRANE expects to decode. (Default=ANN)",
    )

    args = parser.parse_args()
    expression = args.expression
    if ".__" in expression:
        raise ValueError("basic sanity check failed")

    with VariantFile(args.vcf) as vcf:
        with VariantFile("-", "w", header=vcf.header) as out:
            for record in filter_vcf(vcf, expression, args.vcf_info_field):
                out.write(record)
