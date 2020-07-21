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
from itertools import chain, islice
from sys import stderr
from typing import Iterator, List, Tuple

import yaml
from pysam import VariantFile, VariantRecord, VariantHeader
from pysam.libcbcf import (
    VariantRecordInfo,
    VariantRecordSample,
    VariantRecordFormat,
    VariantRecordSamples,
)
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


class Format:
    def __init__(
        self, record_idx: int, name: str, record_samples: VariantRecordSamples
    ):
        self._record_idx = record_idx
        self._name = name
        self._record_samples = record_samples
        self._sample_values = {}

    def __getitem__(self, sample):
        try:
            return self._sample_values[sample]
        except KeyError:
            try:
                record_sample = self._record_samples[sample]
            except KeyError:
                raise UnknownSample(self._record_idx, sample)
            value = record_sample[self._name]
            self._sample_values[sample] = value
            return value


class Formats:
    def __init__(
        self,
        record_idx: int,
        record_format: VariantRecordFormat,
        record_samples: VariantRecordSamples,
    ):
        self._record_idx = record_idx
        self._record_format = record_format
        self._record_samples = record_samples
        self._formats = {}

    def __getitem__(self, item):
        try:
            return self._formats[item]
        except KeyError:
            try:
                self._record_format[item]
            except KeyError:
                raise UnknownFormatField(self._record_idx, item)
            format_field = Format(self._record_idx, item, self._record_samples)
            self._formats[item] = format_field
            return format_field


class Info:
    def __init__(self, record_idx: int, record_info: VariantRecordInfo, ann_key: str):
        self._record_idx = record_idx
        self._record_info = record_info
        self._ann_key = ann_key
        self._info_dict = {}

    def __getitem__(self, item):
        try:
            return self._info_dict[item]
        except KeyError:
            try:
                if item == self._ann_key:
                    raise KeyError(item)
                untyped_value = self._record_info[item]
            except KeyError as ke:
                raise UnknownInfoField(self._record_idx, ke)
            value = self._info_dict[item] = type_info(untyped_value)
            return value


class Annotation:
    def __init__(self, ann_key: str, header: VariantHeader):
        self._record_idx = -1
        self._annotation_data = {}
        self._data = {}
        annotation_keys = get_annotation_keys(header, ann_key)
        self._ann_conv = {
            entry.name: (ann_idx, entry.convert)
            for ann_idx, entry in enumerate(map(ANN_TYPER.get_entry, annotation_keys))
        }

    def update(self, record_idx: int, annotation: str):
        self._record_idx = record_idx
        self._data.clear()
        self._annotation_data = split_annotation_entry(annotation)

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            try:
                ann_idx, convert = self._ann_conv[item]
            except KeyError as ke:
                raise UnknownAnnotation(self._record_idx, ke)
            raw_value = self._annotation_data[ann_idx].strip()
            value = self._data[item] = convert(raw_value)
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


def split_annotation_entry(entry: str,) -> List[str]:
    return entry.split("|")


UNSET = object()


class Environment(dict):
    def __init__(self, expression: str, ann_key: str, header: VariantHeader):
        self._ann_key: str = ann_key
        self._has_ann: bool = any(
            hasattr(node, "id") and isinstance(node, ast.Name) and node.id == ann_key
            for node in ast.walk(ast.parse(expression))
        )
        self._annotation: Annotation = Annotation(ann_key, header)
        self._globals = {}
        # We use self + self.func as a closure.
        self._globals = globals_whitelist.copy()
        self._func = eval(f"lambda: {expression}", self, {})

        self._getters = {
            "CHROM": self._get_chrom,
            "POS": self._get_pos,
            "ID": self._get_id,
            "ALT": self._get_alt,
            "REF": self._get_ref,
            "QUAL": self._get_qual,
            "FILTER": self._get_filter,
            "INFO": self._get_info,
            "FORMAT": self._get_format,
            "SAMPLES": self._get_samples,
        }
        self._empty_globals = {name: NA for name in header.info}
        self._empty_globals.update({name: UNSET for name in self._getters})
        self.record: VariantRecord = None
        self.idx: int = -1

    def filters_annotations(self):
        return self._has_ann

    def update_from_record(self, idx: int, record: VariantRecord):
        self.idx = idx
        self.record = record
        self._globals.update(self._empty_globals)

    def _get_chrom(self) -> str:
        value = self.record.chrom
        self._globals["CHROM"] = value
        return value

    def _get_pos(self) -> int:
        value = self.record.pos
        self._globals["POS"] = value
        return value

    def _get_id(self) -> str:
        value = self.record.id
        self._globals["ID"] = value
        return value

    def _get_ref_alt(self) -> Tuple[str, List[str]]:
        ref, *alt = chain(self.record.alleles)
        self._globals["REF"], self._globals["ALT"] = ref, alt
        return ref, alt

    def _get_ref(self) -> str:
        return self._get_ref_alt()[0]

    def _get_alt(self) -> List[str]:
        return self._get_ref_alt()[1]

    def _get_qual(self) -> float:
        value = type_info(self.record.qual)
        self._globals["QUAL"] = value
        return value

    def _get_filter(self) -> str:
        value = self.record.filter
        self._globals["FILTER"] = value
        return value

    def _get_info(self) -> Info:
        value = Info(self.idx, self.record.info, self._ann_key)
        self._globals["INFO"] = value
        return value

    def _get_format(self) -> Formats:
        value = Formats(self.idx, self.record.format, self.record.samples)
        self._globals["FORMAT"] = value
        return value

    def _get_samples(self) -> List[VariantRecordSample]:
        value = list(self.record.samples)
        self._globals["SAMPLES"] = value
        return value

    def __getitem__(self, item):
        if item == self._ann_key:
            return self._annotation
        value = self._globals[item]
        if value is UNSET:
            value = self._getters[item]()
        return value

    def evaluate(self, annotation: str = "") -> bool:
        if self._has_ann:
            self._annotation.update(self.idx, annotation)
        return self._func()


def filter_vcf(
    vcf: VariantFile, expression: str, ann_key: str, keep_unmatched: bool = False,
) -> Iterator[VariantRecord]:

    env = Environment(expression, ann_key, vcf.header)

    record: VariantRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.filters_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            try:
                annotations = record.info[ann_key]
            except KeyError:
                annotations = [""]
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
            if env.evaluate():
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
            for key, raw_value in zip(
                annotation_keys, split_annotation_entry(annotation)
            ):
                value = raw_value.strip()
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

        records = filter_vcf(
            vcf,
            args.expression,
            args.annotation_key,
            keep_unmatched=args.keep_unmatched,
        )

        try:
            first_record = list(islice(records, 1))
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)

        records = chain(first_record, records)

        with VariantFile(args.output, "w" + fmt, header=header,) as out:
            if args.statistics is not None:
                records = statistics(records, vcf, args.statistics, args.annotation_key)

            try:
                for record in records:
                    out.write(record)
            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
