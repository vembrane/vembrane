import argparse
import ast
import contextlib
import shlex
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, TextIO, Type

from .backend.backend_cyvcf2 import Cyvcf2Reader, Cyvcf2Writer
from .backend.backend_pysam import PysamReader, PysamWriter
from .backend.base import Backend, VCFHeader, VCFReader, VCFRecord
from .errors import InvalidExpressionError, VembraneError
from .globals import _explicit_clear
from .sequence_ontology import SequenceOntology


def add_common_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field. "
        "This defaults to 'ANN', but tools might use other field names. "
        "For example, default VEP annotations can be parsed by setting 'CSQ' here.",
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
        "--context",
        help="Python statement defining a context for given Python expressions. "
        "Extends eventual definitions given via --context-file. "
        "Any global variables (or functions) become available in the Python "
        "expressions. Note that the code you pass here is not sandboxed and should "
        "be trusted. Carefully review any code you get from the internet or AI.",
    )
    parser.add_argument(
        "--context-file",
        help="Path to Python script defining a context for given Python expressions. "
        "Any global variables (or functions) become available in the Python "
        "expressions. Note that the code you pass here is not sandboxed and should "
        "be trusted. Carefully review any code you get from the internet or AI.",
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
        try:
            return list(
                map(
                    str.strip,
                    h.get("Description").strip('"').split(separator)[1].split("|"),
                ),
            )
        except Exception as e:
            raise VembraneError(
                "Could not parse annotation keys from header. Has your VCF been "
                "properly annotated (e.g. with VEP or SnpEff)?"
            ) from e
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


@dataclass
class Context:
    context_file: str | None = None
    context_file_path: Path | None = None
    context_statement: str | None = None
    internal_context: str | None = None

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "Context":
        context_file = None
        context_statement = None
        internal_context = None
        context_file_path = None
        if path := args.context_file:
            context_file_path = Path(path).absolute()
            parent = context_file_path.parent
            internal_context = f"import sys\nsys.path.insert(0, '{parent}')\n"
            with open(path) as f:
                context_file = f.read()
        if statement := args.context:
            context_statement = statement
        return cls(
            context_file=context_file,
            context_file_path=context_file_path,
            context_statement=context_statement,
            internal_context=internal_context,
        )

    def get_globals(self) -> dict[str, Any]:
        globals_dict: dict[str, Any] = {}
        try:
            if self.internal_context:
                exec(self.internal_context, globals_dict)
            if self.context_file:
                exec(self.context_file, globals_dict)
            if self.context_statement:
                exec(self.context_statement, globals_dict)
        except Exception as e:
            raise VembraneError("Error while executing context: {e}") from e
        for name in _explicit_clear:
            if name in globals_dict:
                del globals_dict[name]
        return globals_dict


def create_reader(
    filename: str | Path,
    backend: Backend = Backend.pysam,
    overwrite_number=None,
) -> VCFReader:
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
    filename: str | Path,
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
def smart_open(filename=None, *args, **kwargs) -> Iterator[TextIO]:
    fh = open(filename, *args, **kwargs) if filename and filename != "-" else sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


type Primitive = str | int | float | bool | None


# Inspired by https://stackoverflow.com/a/13793033
class Singleton(type):
    """A singleton metaclass: ensures that only one instance of a class exists."""

    _instances: dict[Type, Any] = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]


class HumanReadableDefaultsFormatter(argparse.HelpFormatter):
    """A custom argparse formatter that displays default values in a more
    human-readable way than the argparse's own defaults formatter which
    exposes Python types like dicts and None."""

    def _get_help_string(self, action):
        help_text = action.help
        if help_text is None:
            help_text = ""

        if action.default and action.default is not argparse.SUPPRESS:
            if isinstance(action.default, list):
                default_value = " ".join(map(str, action.default))
            elif isinstance(action.default, dict):
                default_value = " ".join(f"{k}={v}" for k, v in action.default.items())
            else:
                default_value = action.default
            help_text += f" (default: {default_value})"
        return help_text
