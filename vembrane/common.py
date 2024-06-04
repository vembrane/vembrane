from __future__ import annotations

import argparse
import ast
import shlex
from collections import defaultdict
from typing import Iterable, Iterator

import pandas as pd
from intervaltree import Interval, IntervalTree

from .backend.backend_cyvcf2 import Cyvcf2Reader, Cyvcf2Writer
from .backend.backend_pysam import PysamReader, PysamWriter
from .backend.base import Backend, VCFHeader, VCFReader, VCFRecord
from .errors import InvalidExpressionError


def add_common_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--overwrite-number-info",
        nargs=1,
        action=AppendKeyValuePair,
        metavar="FIELD=NUMBER",
        default={},
        help="Overwrite the number specification for INFO fields "
        "given in the VCF header. "
        "Example: `--overwrite-number cosmic_CNT=.`",
    )
    parser.add_argument(
        "--overwrite-number-format",
        nargs=1,
        action=AppendKeyValuePair,
        metavar="FIELD=NUMBER",
        default={},
        help="Overwrite the number specification for FORMAT fields "
        "given in the VCF header. "
        "Example: `--overwrite-number-format DP=2`",
    )
    parser.add_argument(
        "--backend",
        "-b",
        default="cyvcf2",
        type=Backend.from_string,
        choices=[Backend.cyvcf2, Backend.pysam],
        help="Set the backend library.",
    )
    parser.add_argument(
        "--definitions",
        "-d",
        metavar="FILENAME",
        default=None,
        help="A file containing additional info and annotation definitions.",
    )


def check_expression(expression: str) -> str:
    if ".__" in expression:
        raise InvalidExpressionError(
            expression,
            "The expression must not contain '.__'",
        )
    try:
        tree = ast.parse(expression, mode="eval")
        if isinstance(tree.body, ast.BoolOp | ast.Compare):
            return expression
        else:
            # TODO possibly check for ast.Call, func return type
            return expression
    except SyntaxError as se:
        raise InvalidExpressionError(
            expression,
            "The expression has to be syntactically correct.",
        ) from se


def swap_quotes(s: str) -> str:
    return s.replace('"', '\\"').replace("'", '"').replace('\\"', "'")


def single_outer(s: str) -> bool:
    if '"' in s and "'" in s:
        return s.index('"') > s.index("'")
    elif '"' in s:
        return True
    elif "'" in s:
        return False
    return True


def normalize(s: str) -> str:
    return shlex.quote(swap_quotes(s) if not single_outer(s) else s)


def get_annotation_description_and_keys(header: VCFHeader, ann_key: str) -> list[str]:
    separator = "'"
    if header.contains_generic("VEP"):
        separator = ":"
    if h := header.infos.get(ann_key):
        split = h.get("Description").split(separator, 3)
        if len(split) == 3:
            prefix, key_string, suffix = h.get("Description").split(separator)
        else:
            prefix, key_string = h.get("Description").split(separator)
            suffix = ""
        description_string = separator.join([prefix, "{keys}", suffix])

        return description_string, list(map(str.strip, key_string.split("|")))
    return "", []


def get_annotation_keys(header: VCFHeader, ann_key: str) -> list[str]:
    return get_annotation_description_and_keys(header, ann_key)[1]


def split_annotation_entry(entry: str) -> list[str]:
    return entry.split("|")


class Auxiliary(dict):
    def __init__(self, path: str = None, df=None, environment=None):
        if path:
            self._df = pd.read_csv(path, sep="\t")
        elif df is not None:
            self._df = df
        else:
            raise AssertionError()

        for c in self._df.columns:
            self[c] = None
        self._tree: dict[IntervalTree] = None
        self._sets = dict()
        self.environment = environment

    def __getitem__(self, arg):
        # lazy set building
        if arg not in self._sets:
            self._sets[arg] = set(self._df[arg])
        return self._sets[arg]

    def overlap(
        self,
        chrom_column="chrom",
        start_column="start",
        end_column="end",
    ) -> Auxiliary:
        # lazy tree building
        # TODO: build lazy only chromosome-wise?
        if self._tree is None:
            self._tree = dict()
            available_chromsomes: set[str] = set(self._df[chrom_column])
            for chrom in available_chromsomes:
                d = self._df[self._df[chrom_column] == chrom]
                self._tree[chrom] = IntervalTree(
                    Interval(d[start_column], d[end_column], i) for i, d in d.iterrows()
                )

        t = self._tree[self.environment._get_chrom()]
        start, end = self.environment._get_pos(), self.environment._get_end() + 1
        indices = [i for _, _, i in t.overlap(start, end)]
        test = self._df.loc[indices]
        return Auxiliary(df=test, environment=self.environment)


class BreakendEvent:
    __slots__ = ["name", "keep", "records", "keep_records", "mate_pair"]

    def __init__(self, name: str, mate_pair: bool = False) -> None:
        self.name = name
        self.records: list[VCFRecord] = []
        self.keep_records: list[bool] = []
        self.keep = False
        self.mate_pair = mate_pair

    def add(self, record: VCFRecord, keep_record: bool):
        self.records.append(record)
        self.keep_records.append(keep_record)
        self.keep |= keep_record

    def emit(self) -> Iterator[VCFRecord]:
        assert self.keep
        yield from self.records
        self.records = []
        self.keep_records = []
        # do not reset self.keep!

    def is_mate_pair(self) -> bool:
        return self.mate_pair

    def __str__(self) -> str:
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name


def mate_key(mates: Iterable[str | None]) -> str:
    return "__MATES: " + ",".join(sorted(m for m in mates if m is not None))


class AppendTagExpression(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        assert len(values) == 1
        if not hasattr(namespace, self.dest) or getattr(namespace, self.dest) is None:
            setattr(namespace, self.dest, {})
        value = values[0].strip()
        key, value = value.strip().split("=", 1)
        expr = check_expression(value)
        getattr(namespace, self.dest)[key.strip()] = expr


class AppendKeyValuePair(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        assert len(values) == 1
        if not hasattr(namespace, self.dest) or getattr(namespace, self.dest) is None:
            setattr(namespace, self.dest, {})
        value = values[0].strip()
        key, value = value.strip().split("=", 1)
        getattr(namespace, self.dest)[key.strip()] = value


def read_auxiliary(aux: dict[str, str]) -> dict[Auxiliary]:
    # read auxiliary files
    return {name: Auxiliary(contents) for name, contents in aux.items()}


def create_reader(
    filename: str,
    backend: Backend = Backend.pysam,
    overwrite_number=None,
):
    # GT should always be "."
    if overwrite_number is None:
        overwrite_number = defaultdict(dict)
    if backend == Backend.pysam:
        return PysamReader(filename, overwrite_number)
    elif backend == Backend.cyvcf2:
        return Cyvcf2Reader(filename, overwrite_number)
    else:
        raise ValueError(f"{backend} is not a known backend.")


def create_writer(
    filename: str,
    fmt: str,
    template: VCFReader,
    backend: Backend = Backend.pysam,
):
    if backend == Backend.pysam:
        return PysamWriter(filename, fmt, template)
    elif backend == Backend.cyvcf2:
        return Cyvcf2Writer(filename, fmt, template)
    else:
        raise ValueError(f"{backend} is not a known backend.")
