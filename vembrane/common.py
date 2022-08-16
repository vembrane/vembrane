import ast
from typing import Iterable, Iterator, List

from pysam.libcbcf import VariantHeader, VariantRecord

from .errors import InvalidExpression


def check_expression(expression: str) -> str:
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


def get_annotation_keys(header: VariantHeader, ann_key: str) -> List[str]:
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
                )
            )
    return []


def split_annotation_entry(entry: str) -> List[str]:
    return entry.split("|")


class BreakendEvent(object):
    __slots__ = ["name", "keep", "records", "keep_records", "mate_pair"]

    def __init__(self, name: str, mate_pair: bool = False):
        self.name = name
        self.records = []
        self.keep_records = []
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

    def __str__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name


def mate_key(mates: Iterable[str]) -> str:
    return "__MATES: " + ",".join(sorted(mates))
