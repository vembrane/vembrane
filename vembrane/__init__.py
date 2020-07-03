__version__ = "0.1.0"

import argparse
import ast
import math
import re
import yaml

from itertools import chain
from sys import stderr
from typing import Iterator, List

from pysam import VariantFile, VariantRecord, VariantHeader
from functools import lru_cache
from collections import defaultdict

from vembrane.errors import UnknownAnnotation, UnknownInfoField, InvalidExpression

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


def get_annotation_keys(header: VariantHeader,) -> List[str]:
    for rec in header.records:
        if rec.get("ID") == "ANN":
            return list(map(str.strip, rec.get("Description").split("'")[1].split("|")))
    return []


@lru_cache
def parse_annotation_entry(entry: str,) -> List[str]:
    return list(map(str.strip, entry.split("|")))


def eval_expression(
    expression: str, idx: int, annotation: str, annotation_keys: List[str], env: dict,
) -> bool:
    env["ANN"] = dict(zip(annotation_keys, parse_annotation_entry(annotation)))
    try:
        return eval(expression, globals_whitelist, env)
    except KeyError as ke:
        raise UnknownAnnotation(idx, ke)
    except NameError as ne:
        raise UnknownInfoField(idx, str(ne).split("'")[1])


def filter_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str = "ANN",
    keep_unmatched: bool = False,
) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    annotation_keys = get_annotation_keys(header)

    for idx, record in enumerate(vcf):
        # setup filter expression env
        env.clear()
        for name in header.info:
            env[name] = None

        env["CHROM"] = record.chrom
        env["POS"] = record.pos
        (env["REF"], env["ALT"]) = chain(record.alleles)
        for key in record.info:
            if key != ann_key:
                env[key] = record.info[key]

        annotations = dict(record.info).get(ann_key, [""])
        filtered_annotations = [
            annotation
            for annotation in annotations
            if eval_expression(expression, idx, annotation, annotation_keys, env,)
        ]

        if not filtered_annotations:
            # skip this record if filter removed all annotations
            continue
        elif not keep_unmatched and (len(annotations) != len(filtered_annotations)):
            # update annotations if they have been actually filtered
            record.info[ann_key] = filtered_annotations
        yield record


def statistics(
    records: Iterator[VariantRecord], vcf: VariantFile, filename: str
) -> Iterator[VariantRecord]:
    annotation_keys = get_annotation_keys(vcf.header)
    counter = defaultdict(lambda: defaultdict(lambda: 0))
    for record in records:
        for annotation in record.info["ANN"]:
            for key, value in zip(annotation_keys, parse_annotation_entry(annotation)):
                if value:
                    counter[key][value] += 1
        yield record

    # reduce dicts with many items, to just one counter
    for key, subdict in counter.items():
        if len(subdict) > 10:
            counter[key] = f"#{len(subdict)}"

    yaml.add_representer(defaultdict, yaml.representer.Representer.represent_dict)
    with open(filename, "w") as out:
        yaml.dump(dict(counter), out)


def check_filter_expression(expression: str,) -> str:
    if ".__" in expression:
        raise InvalidExpression(expression, "The expression must not contain '.__'")
    try:
        tree = ast.parse(expression, mode="eval")
        if isinstance(tree.body, (ast.BoolOp, ast.Compare)):
            return expression
        else:
            # TODO possibly check for ast.Call, func return type
            return expression
    except SyntaxError:
        raise InvalidExpression(
            expression, "The expression has to be syntactically correct."
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="The file containing the variants.")
    parser.add_argument(
        "expression",
        type=check_filter_expression,
        help="Filter variants and annotations. If this removes all annotations, "
        "the variant is removed as well.",
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
        "--statistics",
        "-s",
        metavar="FILE",
        default=None,
        help="Write statistics to this file.",
    )
    parser.add_argument(
        "--keep-unmatched",
        default=False,
        action="store_true",
        help="Keep all annotations of a variant if at least one of them passes "
        "the expression.",
    )
    args = parser.parse_args()

    with VariantFile(args.vcf) as vcf:
        fmt = ""
        if args.output_fmt == "bcf":
            fmt = "b"
        elif args.output_fmt == "uncompressed-bcf":
            fmt = "u"
        with VariantFile(args.output, "w" + fmt, header=vcf.header,) as out:
            records = filter_vcf(
                vcf,
                args.expression,
                keep_unmatched=args.keep_unmatched,
                ann_key=args.annotation_key,
            )
            if args.statistics is not None:
                records = statistics(records, vcf, args.statistics)

            try:
                for record in records:
                    out.write(record)
            except UnknownAnnotation as ua:
                print(ua, file=stderr)
                exit(1)
            except UnknownInfoField as ua:
                print(ua, file=stderr)
                exit(1)
