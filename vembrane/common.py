import argparse
import ast
import shlex
from typing import Dict, Iterable, Iterator, List, Optional, Set
from sys import stderr

from cyvcf2 import VCF
from cyvcf2.cyvcf2 import Variant

from .errors import InvalidExpression


def header_keys(vcf, header_type):
    return set(x.info().get("ID") for x in vcf.header_iter() if x.type == header_type)


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


def get_annotation_keys(vcf: VCF, ann_key: str) -> List[str]:
    separator = "'"
    if vcf.contains("VEP"):
        separator = ":"
    if vcf.contains(ann_key):
        return list(
            map(
                str.strip,
                vcf.get_header_type(ann_key)["Description"]
                .strip('"')
                .split(separator)[1]
                .split("|"),
            )
        )
    return []


def split_annotation_entry(entry: str) -> List[str]:
    return entry.split("|")


def is_bnd_record(record: Variant) -> bool:
    return "SVTYPE" in record.INFO and record.INFO.get("SVTYPE", None) == "BND"


class BreakendEvent(object):
    __slots__ = ["name", "keep", "records", "keep_records", "mate_pair"]

    def __init__(self, name: str, mate_pair: bool = False):
        self.name = name
        self.records: List[Variant] = []
        self.keep_records: List[bool] = []
        self.keep = False
        self.mate_pair = mate_pair

    def add(self, record: Variant, keep_record: bool):
        self.records.append(record)
        self.keep_records.append(keep_record)
        self.keep |= keep_record

    def emit(self) -> Iterator[Variant]:
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


def mate_key(mates: Iterable[Optional[str]]) -> str:
    return "__MATES: " + ",".join(sorted(m for m in mates if m is not None))


class AppendTagExpression(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        assert len(values) == 1
        if not hasattr(namespace, self.dest) or getattr(namespace, self.dest) is None:
            setattr(namespace, self.dest, dict())
        value = values[0].strip()
        key, value = value.strip().split("=", 1)
        expr = check_expression(value)
        getattr(namespace, self.dest)[key.strip()] = expr


class AppendKeyValuePair(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        assert len(values) == 1
        if not hasattr(namespace, self.dest) or getattr(namespace, self.dest) is None:
            setattr(namespace, self.dest, dict())
        value = values[0].strip()
        key, value = value.strip().split("=", 1)
        getattr(namespace, self.dest)[key.strip()] = value


def read_auxiliary(aux: Dict[str, str]) -> Dict[str, Set[str]]:
    # read auxiliary files, split at any whitespace and store contents in a set
    def read_set(path: str) -> Set[str]:
        with open(path, "rt") as f:
            return set(line.rstrip() for line in f)

    return {name: read_set(contents) for name, contents in aux.items()}
