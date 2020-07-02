__version__ = "0.1.0"

import argparse
import math
import re
from itertools import chain
from typing import Iterator, List, Dict

from pysam import VariantFile, VariantRecord, VariantHeader

globals_whitelist = {
    **{
        "__builtins__": {},
        "__builtin__": {},
        "__file__": None,
        "__name__": None,
        "__doc__": None,
        "__package__": None,
    },
    **{mod.__name__: mod for mod in [any, all, min, max, re, list, dict, zip]},
    **{name: mod for name, mod in vars(math).items() if not name.startswith("__")},
}


def get_annotation_keys(header: VariantHeader) -> List[str]:
    for rec in header.records:
        if rec.get("ID") == "ANN":
            return list(map(str.strip, rec.get("Description").split("'")[1].split("|")))
    return []


def parse_annotation_entry(entry: str) -> List[str]:
    return list(map(str.strip, entry.split("|")))


def filter_annotation_entries(
    entries: List[str], ann_filter_expression: str, annotation_keys: List[str], env: dict
) -> Iterator[str]:
    for entry in entries:
        env["ANN"] = dict(zip(annotation_keys, parse_annotation_entry(entry)))
        if eval(ann_filter_expression, globals_whitelist, env):
            yield entry


def filter_vcf(
    vcf: VariantFile, filter_expression: str, ann_filter_expression: str
) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    for name in header.info:
        env[name] = None

    annotation_keys = get_annotation_keys(header)

    for record in vcf:
        # obtain annotation entries
        ann = dict(
            zip(
                annotation_keys,
                zip(
                    *(
                        parse_annotation_entry(entry)
                        for entry in record.info.get("ANN", [])
                    )
                ),
            )
        )

        # setup filter expression env
        for key in record.info:
            if key != "ANN":
                env[key] = record.info[key]
        env["ANN"] = ann
        env["CHROM"] = record.chrom
        env["POS"] = record.pos
        env["REF"], env["ALT"] = chain(record.alleles)

        if not filter_expression or eval(filter_expression, globals_whitelist, env):
            if ann_filter_expression:
                # filter annotation entries
                ann = record.info.get("ANN")
                if ann:
                    filtered_ann = list(
                        filter_annotation_entries(
                            ann, ann_filter_expression, annotation_keys, env
                        )
                    )
                    if not filtered_ann:
                        # skip this record if filter removed all annotations
                        continue
                    record.info["ANN"] = filtered_ann

            yield record


def check_filter_expression(expression: str):
    if ".__" in expression:
        raise ValueError("basic sanity check failed")  # TODO: better error message


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="The file containing the variants.")
    parser.add_argument(
        "--filter-expression", help="An expression to filter the variants."
    )
    parser.add_argument(
        "--ann-filter-expression",
        help="Filter annotation entries. If this removes all annotations, "
        "the variant is removed as well.",
    )
    args = parser.parse_args()

    with VariantFile(args.vcf) as vcf:
        with VariantFile("-", "w", header=vcf.header) as out:
            for record in filter_vcf(
                vcf, args.filter_expression, args.ann_filter_expression
            ):
                out.write(record)
