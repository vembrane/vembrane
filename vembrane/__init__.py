try:
    from importlib.metadata import version, PackageNotFoundError  # type: ignore
except ImportError:  # pragma: no cover
    from importlib_metadata import version, PackageNotFoundError  # type: ignore
try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"

import argparse
import ast
import math
import re
import sys
from collections import defaultdict
from functools import lru_cache
from itertools import chain
from sys import stderr
from typing import Any, Callable, Dict, Iterator, List, Tuple

import yaml
from pysam import VariantFile, VariantRecord, VariantHeader

from vembrane.ann_types import NA, type_info, ANN_TYPER
from vembrane.errors import (
    UnknownAnnotation,
    UnknownInfoField,
    InvalidExpression,
    UnknownFormatField,
    UnknownSample,
    VembraneError,
)

globals_whitelist = {
    **{
        "__builtins__": {},
        "__builtin__": {},
        "__file__": None,
        "__name__": None,
        "__doc__": None,
        "__package__": None,
    },
    **{
        mod.__name__: mod
        for mod in [any, all, min, max, re, list, dict, set, tuple, zip, map, sum]
    },
    **{name: mod for name, mod in vars(math).items() if not name.startswith("__")},
}


class Sample:
    def __init__(self, record_idx: int, sample: str, format_data: Dict[str, Any]):
        self._record_idx = record_idx
        self._sample = sample
        self._data = format_data

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError as ke:
            raise UnknownFormatField(self._record_idx, self._sample, ke)


class Format:
    def __init__(self, record_idx: int, sample_formats: Dict[Sample, Dict[str, Any]]):
        self._record_idx = record_idx
        self._sample_formats = sample_formats

    def __getitem__(self, item):
        try:
            return self._sample_formats[item]
        except KeyError as ke:
            raise UnknownSample(self._record_idx, ke)


class Info:
    def __init__(self, record_idx: int, info_dict: Dict[str, Dict[str, Any]]):
        self._record_idx = record_idx
        self._info_dict = info_dict
        self._typed_dict = {}

    def __getitem__(self, item):
        try:
            return self._typed_dict[item]
        except KeyError:
            try:
                untyped_value = self._info_dict[item]
            except KeyError as ke:
                raise UnknownInfoField(self._record_idx, ke)
            value = self._typed_dict[item] = type_info(untyped_value)
            return value


class Annotation:
    def __init__(self, record_idx: int, annotation_data: List[str], ann_conv):
        self._record_idx = record_idx
        self._data = {}
        self._annotation_data = annotation_data
        self._ann_conv = ann_conv

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            try:
                ann_idx, convert = self._ann_conv[item]
            except KeyError as ke:
                raise UnknownAnnotation(self._record_idx, ke)
            value = self._data[item] = convert(self._annotation_data[ann_idx])
            return value


def get_annotation_keys(header: VariantHeader, ann_key: str) -> List[str]:
    separator = "'"
    for rec in header.records:
        if rec.key == "VEP":
            separator = ":"
            continue
        if rec.get("ID") == ann_key:
            return list(
                map(str.strip, rec.get("Description").split(separator)[1].split("|"))
            )
    return []


@lru_cache(maxsize=32)
def parse_annotation_entry(entry: str,) -> List[str]:
    return list(map(str.strip, entry.split("|")))


class Expression:
    def __init__(self, expression: str, ann_key: str = "ANN"):
        # We use self._globals + self.func as a closure.
        # Do not reassign self._globals, but use .update() on it!
        self._globals = globals_whitelist.copy()
        self._func = eval(f"lambda: {expression}", self._globals, {})
        self._ann_key = ann_key
        self._has_ann = any(
            hasattr(node, "id") and isinstance(node, ast.Name) and node.id == ann_key
            for node in ast.walk(ast.parse(expression))
        )

    def annotation_key(self):
        return self._ann_key

    def filters_annotations(self):
        return self._has_ann

    def evaluate(
        self,
        idx: int,
        annotation: str,
        ann_conv: Dict[str, Tuple[int, Callable[[str], Any]]],
        env: dict,
    ) -> bool:
        if self._has_ann:
            env[self._ann_key] = Annotation(
                idx, parse_annotation_entry(annotation), ann_conv,
            )
        self._globals.update(env)
        return self._func()


class Environment:
    def __init__(self, expression: Expression, header: VariantHeader):
        self.expression = expression
        self._empty_env = {name: NA for name in header.info}
        self._env = {}
        annotation_keys = get_annotation_keys(header, expression.annotation_key())
        self.ann_conv = {
            entry.name: (ann_idx, entry.convert)
            for ann_idx, entry in enumerate(map(ANN_TYPER.get_entry, annotation_keys))
        }

    def update(self, idx, record):
        self.idx = idx
        env = self._env
        env.clear()
        env.update(self._empty_env)

        env["CHROM"] = record.chrom
        env["POS"] = record.pos
        env["ID"] = record.id
        (env["REF"], *env["ALT"]) = chain(record.alleles)
        env["QUAL"] = type_info(record.qual)
        env["FILTER"] = record.filter
        ann_key = self.expression.annotation_key()
        env["INFO"] = Info(
            idx, {key: record.info[key] for key in record.info if key != ann_key},
        )

        formats = {
            sample: Sample(
                idx, sample, {fmt: record.samples[sample][fmt] for fmt in record.format}
            )
            for sample in record.samples
        }

        env["FORMAT"] = Format(idx, formats)
        env["SAMPLES"] = list(record.samples)

    def evaluate(self, annotation) -> bool:
        return self.expression.evaluate(self.idx, annotation, self.ann_conv, self._env)


def filter_vcf(
    vcf: VariantFile, expression: Expression, keep_unmatched: bool = False,
) -> Iterator[VariantRecord]:

    ann_key = expression.annotation_key()
    env = Environment(expression, vcf.header)

    for idx, record in enumerate(vcf):
        env.update(idx, record)
        if expression.filters_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            annotations = record.info.get(ann_key, [""])
            #  … and only keep the annotations where the expression evaluates to true
            filtered_annotations = [
                annotation for annotation in annotations if env.evaluate(annotation)
            ]
            if not filtered_annotations:
                # skip this record if filter removed all annotations
                continue
            elif not keep_unmatched and (len(annotations) != len(filtered_annotations)):
                # update annotations if they have actually been filtered
                record.info[ann_key] = filtered_annotations
            yield record
        else:
            # otherwise, the annotations are irrelevant w.r.t. the expression,
            # so we can omit them
            if env.evaluate(""):
                yield record
            else:
                continue


def statistics(
    records: Iterator[VariantRecord], vcf: VariantFile, filename: str, ann_key: str
) -> Iterator[VariantRecord]:
    annotation_keys = get_annotation_keys(vcf.header, ann_key)
    counter = defaultdict(lambda: defaultdict(lambda: 0))
    for record in records:
        for annotation in record.info[ann_key]:
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
    expression = Expression(args.expression, ann_key=args.annotation_key)

    with VariantFile(args.vcf) as vcf:
        fmt = ""
        if args.output_fmt == "bcf":
            fmt = "b"
        elif args.output_fmt == "uncompressed-bcf":
            fmt = "u"

        header: VariantHeader = vcf.header
        header.add_meta("vembraneVersion", __version__)
        header.add_meta(
            "vembraneCmd",
            "vembrane "
            + " ".join(
                "'" + arg.replace("'", '"') + '"' if " " in arg else arg
                for arg in sys.argv[1:]
            ),
        )

        with VariantFile(args.output, "w" + fmt, header=header,) as out:
            records = filter_vcf(vcf, expression, keep_unmatched=args.keep_unmatched,)
            if args.statistics is not None:
                records = statistics(records, vcf, args.statistics, args.annotation_key)

            try:
                for record in records:
                    out.write(record)
            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
