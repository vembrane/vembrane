import argparse
import ast
import shlex
from collections.abc import Iterable, Iterator

from pysam.libcbcf import VariantHeader, VariantRecord

from .errors import InvalidExpression


def check_expression(expression: str) -> str:
    if ".__" in expression:
        raise InvalidExpression(expression, "The expression must not contain '.__'")
    try:
        tree = ast.parse(expression, mode="eval")
        if isinstance(tree.body, ast.BoolOp | ast.Compare):
            return expression
        else:
            # TODO possibly check for ast.Call, func return type
            return expression
    except SyntaxError as se:
        raise InvalidExpression(
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


def get_annotation_keys(header: VariantHeader, ann_key: str) -> list[str]:
    separator = "'"
    for rec in header.records:
        if rec.key == "VEP":
            separator = ":"
            continue
        if rec.get("ID") == ann_key:
            return list(
                map(
                    str.strip,
                    rec.get("Description").strip('"').split(separator)[1].split("|"),
                ),
            )
    return []


def split_annotation_entry(entry: str) -> list[str]:
    return entry.split("|")


def is_bnd_record(record: VariantRecord) -> bool:
    return "SVTYPE" in record.info and record.info.get("SVTYPE", None) == "BND"


class BreakendEvent:
    __slots__ = ["name", "keep", "records", "keep_records", "mate_pair"]

    def __init__(self, name: str, mate_pair: bool = False) -> None:
        self.name = name
        self.records: list[VariantRecord] = []
        self.keep_records: list[bool] = []
        self.keep = False
        self.mate_pair = mate_pair

    def add(self, record: VariantRecord, keep_record: bool):
        self.records.append(record)
        self.keep_records.append(keep_record)
        self.keep |= keep_record

    def emit(self) -> Iterator[VariantRecord]:
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
