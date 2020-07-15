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
from sys import stderr
from typing import Iterator, List, Dict, Any, Set, Tuple

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

    def update(self, items: Dict[str, Any]):
        self._data.update(items)

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError as ke:
            raise UnknownFormatField(self._record_idx, self._sample, ke)


class Format:
    def __init__(
        self,
        record_idx: int,
        sample_formats: Dict[Sample, Dict[str, Any]],
        format_keys: Set[str],
    ):
        self._record_idx = record_idx
        self._sample_formats = sample_formats
        self._samples = sample_formats.keys()
        self._format_keys = format_keys

    def update(self, idx: int, record: VariantRecord):
        self._record_idx = idx
        for sample_name in self._samples:
            self._sample_formats[sample_name].update(
                {
                    fmt: record.samples[sample_name].get(fmt, NA)
                    for fmt in self._format_keys
                }
            )

    def __getitem__(self, item):
        try:
            return self._sample_formats[item]
        except KeyError as ke:
            raise UnknownSample(self._record_idx, ke)


class Info:
    def __init__(self, record_idx: int, info_dict: Dict[str, Dict[str, Any]]):
        self._record_idx = record_idx
        self._info_dict = info_dict
        self._keys = info_dict.keys()

    def update(self, idx: int, record: VariantRecord):
        self._record_idx = idx
        self._info_dict.update(
            {key: type_info(record.info.get(key, NA)) for key in self._keys}
        )

    def __getitem__(self, item):
        try:
            return self._info_dict[item]
        except KeyError as ke:
            raise UnknownInfoField(self._record_idx, ke)


class Annotation:
    def __init__(
        self,
        record_idx: int,
        annotation_data: Dict[str, Any],
        available_keys: List[Tuple[int, str]],
    ):
        self._record_idx = record_idx
        self._data = annotation_data
        self._keys = available_keys

    def update(self, idx: int, annotation: str):
        self._record_idx = idx
        split = annotation.split("|")

        self._data.update(
            dict(
                map(
                    lambda v: ANN_TYPER.convert(v[0], v[1]),
                    ((k, split[i].strip()) for i, k in self._keys),
                )
            ),
        )

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError as ke:
            raise UnknownAnnotation(self._record_idx, ke)


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
def split_annotation_entry(entry: str,) -> List[str]:
    return list(map(str.strip, entry.split("|")))


class Expression:
    def __init__(self, expression: str, ann_key: str = "ANN"):
        self._ann_key = ann_key
        self._expression = expression
        self._has_ann = any(
            hasattr(node, "id") and isinstance(node, ast.Name) and node.id == ann_key
            for node in ast.walk(ast.parse(expression))
        )

    def annotation_key(self) -> str:
        return self._ann_key

    def filters_annotations(self) -> bool:
        return self._has_ann

    @property
    def expression(self) -> str:
        return self._expression


class Environment:
    def __init__(self, expression: Expression, vcf_header: VariantHeader):
        tree = ast.parse(expression.expression)

        # This set contains *all* nodes that refer to a Constant/Name/String/Function/…
        # (note that we have to check for Str explicitly for python 3.7 support).
        names = (
            set(node.id for node in ast.walk(tree) if hasattr(node, "id"))
            | set(
                node.value for node in ast.walk(tree) if isinstance(node, ast.Constant)
            )
            | set(node.s for node in ast.walk(tree) if isinstance(node, ast.Str))
        )

        # restrict names/strings/identifiers/functions/symbols to the ones actually
        # seen in the expression
        available_info_fields = set(vcf_header.info) & names
        available_sample_names = set(vcf_header.samples) & names
        available_format_keys = set(vcf_header.formats) & names
        avaible_vcf_fields = {
            "QUAL",
            "FILTER",
            "ID",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "SAMPLES",
        } & names
        available_symbols = globals_whitelist.keys() & names

        ann_field_name = expression.annotation_key()
        ann_keys = get_annotation_keys(vcf_header, ann_field_name)
        available_annotation_keys = [
            (i, k) for i, k in enumerate(ann_keys) if k in names
        ]

        self.idx = -1
        self.field_lookup = {
            "CHROM": lambda record: record.chrom,
            "POS": lambda record: record.pos,
            "ID": lambda record: record.id,
            "QUAL": lambda record: type_info(record.qual),
            "FILTER": lambda record: record.filter,
            "REF": lambda record: record.alleles[0],
            "ALT": lambda record: record.alleles[1:],
            "SAMPLES": lambda record: list(record.samples),
        }
        self.info = Info(
            record_idx=self.idx,
            info_dict={
                info_key: {}
                for info_key in available_info_fields
                if info_key != ann_field_name
            },
        )
        self.annotation = Annotation(
            record_idx=self.idx,
            annotation_data={field: None for (_, field) in available_annotation_keys},
            available_keys=available_annotation_keys,
        )
        self.sample_names = available_sample_names
        self.formats = Format(
            record_idx=self.idx,
            sample_formats={
                sample: Sample(record_idx=self.idx, sample=sample, format_data={})
                for sample in self.sample_names
            },
            format_keys=available_format_keys,
        )
        self.env = {
            "INFO": self.info,
            "FORMAT": self.formats,
            ann_field_name: self.annotation,
        }
        self.globals = {key: globals_whitelist[key] for key in available_symbols}

        # remove unneeded fields
        for field in self.field_lookup.keys() - set(avaible_vcf_fields):
            self.field_lookup.pop(field, None)

    def update(self, idx: int, record: VariantRecord):
        self.idx = idx
        for field, get_field_value in self.field_lookup.items():
            self.env[field] = get_field_value(record)
        self.info.update(idx, record)
        self.formats.update(idx, record)

    def update_annotation(self, annotation: str) -> bool:
        # update(…) has to be called once for the current record first
        # then update_annotation(…) for each annotation
        self.annotation.update(self.idx, annotation)
        return True

    def evaluate(self, expression: Expression) -> bool:
        return eval(expression.expression, self.globals, self.env)

    def __str__(self):
        return (
            f"\tINFO: {self.info._info_dict}\n\tANN: {self.annotation._data}\n"
            + "\n\t".join(
                f"{field}: {self.env[field]}" for field in self.field_lookup.keys()
            )
        )

    def __repr__(self):
        self.__str__()


def filter_vcf(
    vcf: VariantFile, expression: Expression, keep_unmatched: bool = False,
) -> Iterator[VariantRecord]:

    ann_key = expression.annotation_key()
    env = Environment(expression, vcf.header)

    for idx, record in enumerate(vcf):
        # update the environment with data from the record
        env.update(idx, record)

        # TODO move this block to `Environment.evaluate(…)`
        if expression.filters_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            annotations = dict(record.info).get(ann_key, [""])
            #  … and only keep the annotations where the expression evaluates to true
            filtered_annotations = [
                annotation
                for annotation in annotations
                if env.update_annotation(annotation) and env.evaluate(expression)
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
            if env.evaluate(expression):
                yield record
            else:
                continue


def statistics(
    records: Iterator[VariantRecord], vcf: VariantFile, filename: str, ann_key: str
) -> Iterator[VariantRecord]:
    annotation_keys = get_annotation_keys(vcf.header, ann_key)
    counter = defaultdict(lambda: defaultdict(lambda: 0))
    for record in records:
        for annotation in record.info["ANN"]:
            for key, value in zip(annotation_keys, split_annotation_entry(annotation)):
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
