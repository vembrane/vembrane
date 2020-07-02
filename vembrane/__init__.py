__version__ = "0.1.0"

import argparse
import math
import re
from itertools import chain
from sys import stderr
from typing import Iterator, List

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


def eval_expression(
    expression: str, annotation: str, annotation_keys: List[str], env: dict
) -> bool:
    env["ANN"] = dict(zip(annotation_keys, parse_annotation_entry(annotation)))
    try:
        return eval(expression, globals_whitelist, env)
    except KeyError as ke:
        print(f"Unknown annotation {ke}, skipping", file=stderr)
        return False
    except NameError as ne:
        print(f"{ne}, skipping", file=stderr)
        return False


def filter_vcf(vcf: VariantFile, expression: str) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    for name in header.info:
        env[name] = None

    annotation_keys = get_annotation_keys(header)

    for record in vcf:
        # setup filter expression env
        env.clear()
        env["CHROM"] = record.chrom
        env["POS"] = record.pos
        env["REF"], env["ALT"] = chain(record.alleles)
        for key in record.info:
            if key != "ANN":
                env[key] = record.info[key]

        annotations = record.info.get("ANN", [])
        filtered_annotations = [
            annotation
            for annotation in annotations
            if eval_expression(expression, annotation, annotation_keys, env)
        ]

        if annotations and not filtered_annotations:
            # skip this record if filter removed all annotations
            continue
        elif len(annotations) != len(filtered_annotations):
            # update annotations if they have been actually filtered
            record.info["ANN"] = filtered_annotations
        yield record


def check_filter_expression(expression: str):
    if ".__" in expression:
        raise ValueError("basic sanity check failed")  # TODO: better error message


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="The file containing the variants.")
    parser.add_argument(
        "expression",
        help="Filter variants and annotations. If this removes all annotations, "
        "the variant is removed as well.",
    )
    args = parser.parse_args()

    with VariantFile(args.vcf) as vcf:
        with VariantFile("-", "w", header=vcf.header) as out:
            for record in filter_vcf(vcf, args.expression):
                out.write(record)
