import argparse
import ast
import contextlib
import shlex
import sys
from collections import defaultdict
from typing import Iterable, Iterator

from .backend.backend_cyvcf2 import Cyvcf2Reader, Cyvcf2Writer
from .backend.backend_pysam import PysamReader, PysamWriter
from .backend.base import Backend, VCFHeader, VCFReader, VCFRecord
from .errors import InvalidExpressionError
from .sequence_ontology import SequenceOntology


def add_common_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--aux",
        "-a",
        nargs=1,
        action=AppendKeyValuePair,
        metavar="NAME=PATH",
        default={},
        help="Path to an auxiliary file containing a set of symbols",
    )
    parser.add_argument(
        "--ontology",
        nargs=1,
        default=None,
        metavar="PATH",
        help="Path to an ontology in OBO format. "
        "May be compressed with gzip, bzip2 and xz. "
        "Defaults to built-in ontology (from sequenceontology.org).",
    )
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


def get_annotation_keys(header: VCFHeader, ann_key: str) -> list[str]:
    separator = "'"
    if header.contains_generic("VEP"):
        separator = ":"
    if h := header.infos.get(ann_key):
        return list(
            map(
                str.strip,
                h.get("Description").strip('"').split(separator)[1].split("|"),
            ),
        )
    return []


def split_annotation_entry(entry: str) -> list[str]:
    return entry.split("|")


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


def read_auxiliary(aux: dict[str, str]) -> dict[str, set[str]]:
    # read auxiliary files, split at any whitespace and store contents in a set
    def read_set(path: str) -> set[str]:
        with open(path) as f:
            return {line.rstrip() for line in f}

    return {name: read_set(contents) for name, contents in aux.items()}


def read_ontology(path: str | None) -> SequenceOntology | None:
    return SequenceOntology.from_obo(path) if path else None


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


@contextlib.contextmanager
def smart_open(filename=None, *args, **kwargs):
    fh = open(filename, *args, **kwargs) if filename and filename != "-" else sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()
