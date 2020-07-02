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


def get_annotation_keys(header):
    for rec in header.records:
        if rec.get("ID") == "ANN":
            return list(map(str.strip, rec.get("Description").split("'")[1].split("|")))
    return []


def parse_annotation_entry(entry: str):
    return list(map(str.strip, entry.split("|")))


def filter_annotation_entries(entries: list, ann_filter_expression: str):
    for entry in entries:
        env = dict(zip(annotation_keys, parse_annotation_entry(entry)))
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

        if not filter_expression or eval(filter_expression, globals_whitelist, env):
            if ann_filter_expression:
                # filter annotation entries
                ann = record.info.get("ANN")
                if not ann:
                    filtered_ann = list(
                        filter_annotation_entries(ann, ann_filter_expression)
                    )
                    if not filtered_ann:
                        # skip this record if filter removed all annotations
                        continue
                    record.info["ANN"] = filtered_ann

            yield record


def check_filter_expression(expression):
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
