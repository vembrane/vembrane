import contextlib
import csv
import sys
from collections.abc import Iterator
from sys import stderr
from types import MappingProxyType
from typing import Any

import asttokens

from ..ann_types import NA
from ..backend.base import VCFReader, VCFRecord
from ..common import (
    AppendKeyValuePair,
    add_common_arguments,
    check_expression,
    create_reader,
    read_auxiliary,
)
from ..errors import HeaderWrongColumnNumberError, VembraneError
from ..globals import allowed_globals
from ..representations import Environment
from .filter import DeprecatedAction


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("table")
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "expression",
        type=check_expression,
        help="The expression for the output.",
    )
    parser.add_argument(
        "vcf",
        help="The file containing the variants.",
        nargs="?",
        default="-",
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--separator",
        "-s",
        default="\t",
        metavar="CHAR",
        help="Define the field separator (default: \\t).",
    )
    parser.add_argument(
        "--header",
        default="auto",
        metavar="TEXT",
        help='Override the automatically generated header. Provide "auto" (default) \
              to automatically generate the header from the expression. Provide a \
              comma separated string to manually set the header. Provide "none" to \
              disable any header output.',
    )
    parser.add_argument(
        "--long",
        help="Instead of using `for_each_sample` to generate multiple columns "
        "(wide format), "
        "use only one sample column but one row for each sample (long format).",
        action="store_true",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
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
    add_common_arguments(parser)


def tableize_vcf(
    vcf: VCFReader,
    expression: str,
    ann_key: str,
    overwrite_number: dict[str, dict[str, str]] = MappingProxyType({}),
    long: bool = False,
    auxiliary: dict[str, set[str]] = MappingProxyType({}),
) -> Iterator[tuple]:
    kwargs: dict[str, Any] = dict(auxiliary=auxiliary)
    if long:
        kwargs["evaluation_function_template"] = (
            "lambda: (({expression}) for SAMPLE in SAMPLES)"
        )
    else:
        expression = f"({expression})"
    env = Environment(expression, ann_key, vcf.header, **kwargs)

    record: VCFRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.expression_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            annotations = record.info[ann_key]
            if annotations is NA:
                num_ann_entries = len(env._annotation._ann_conv.keys())
                empty = "|" * num_ann_entries
                print(
                    f"No ANN field found in record {idx}, "
                    f"replacing with NAs (i.e. 'ANN={empty}')",
                    file=sys.stderr,
                )
                annotations = [empty]
            for annotation in annotations:
                if long:
                    yield from env.table(annotation)
                else:
                    yield env.table(annotation)
        else:
            if long:
                yield from env.table()
            else:
                yield env.table()


def generate_for_each_sample_expressions(s: str, vcf: VCFReader) -> list[str]:
    from asttokens.util import replace

    # parse the `for_each_sample(lambda var: inner) expression
    var, inner = _var_and_body(s)

    # here we parse the inner part again, so that position offsets start at 0
    tok = asttokens.ASTTokens(inner, parse=True)

    # *replace* `var` with each value from the list of samples
    samples = list(vcf.header.samples)
    expanded = []
    for sample in samples:
        replacements = []
        for t in tok.tokens:
            if t.type == 1 and t.string == var:
                pos = (t.startpos, t.endpos)
                replacements.append((*pos, f"'{sample}'"))
        expanded.append(replace(inner, replacements))
    return expanded


def generate_for_each_sample_column_names(s: str, vcf: VCFReader) -> list[str]:
    # parse the `for_each_sample(lambda var: inner) expression
    var, inner = _var_and_body(s)

    samples = list(vcf.header.samples)
    __globals = allowed_globals.copy()

    column_names = []
    for sample in samples:
        __globals[var] = sample
        column_name = eval(inner, __globals, {})
        if not isinstance(column_name, str | bytes):
            if hasattr(column_name, "__str__"):
                column_name = str(column_name)
            else:
                raise ValueError(
                    "The specified header expression does not evaluate to a string."
                    "Consider using `str(expression)` instead.",
                )
        column_names.append(column_name)
    return column_names


def _var_and_body(s):
    import ast

    tok = asttokens.ASTTokens(s, parse=True)
    tok.mark_tokens(tok.tree)
    var = None
    inner = ""

    # for_each_sample must be the top level function
    # since that is the only place it is allowed in order to split the header expression
    # such that there is one column per sample (instead of *one* list-valued column)

    # walk the resulting AST, find the "for_each_sample" ast.Call node,
    # and extract the variable name and lambda body code from that
    for node in asttokens.util.walk(tok.tree):
        if isinstance(node, ast.Call) and hasattr(node.func, "id"):
            if node.func.id == "for_each_sample":
                for arg in node.args:
                    if isinstance(arg, ast.Lambda):
                        var = tok.get_text(arg.args)
                        inner = tok.get_text(arg.body)
    return var, inner


def preprocess_expression(
    header: str,
    vcf: VCFReader,
    make_expression: bool = True,
) -> str:
    """
    Split the header expression at toplevel commas into parts.
    Then, if one of these parts starts with 'for_each_sample',
    that part is expanded for each sample in vcf.header.samples
    """
    parts: list[str] = get_toplevel(header)
    to_expand = list(
        filter(lambda x: x[1].startswith("for_each_sample"), enumerate(parts)),
    )
    if len(to_expand) > 0 and vcf is None:
        raise ValueError("If FORMAT is to be expanded, the VCF kwarg must not be none.")
    parts_expanded: list[list[str]] = [[p] for p in parts]
    func = (
        generate_for_each_sample_expressions
        if make_expression
        else generate_for_each_sample_column_names
    )
    for i, p in to_expand:
        expanded = func(p, vcf)
        parts_expanded[i] = expanded
    parts_flattened = [p for pp in parts_expanded for p in pp]
    return ", ".join(parts_flattened)


def get_header(args, vcf: VCFReader) -> list[str]:
    header = args.expression if args.header == "auto" else args.header
    if args.long:
        header = f"SAMPLE, {header}"
    return get_toplevel(preprocess_expression(header, vcf, args.header == "auto"))


def get_toplevel(header: str) -> list[str]:
    splitpos = [0]
    level = 0
    stack = []
    for i, c in enumerate(header):
        if c == "," and level == 0:
            splitpos.append(i + 1)
        elif c == "(":
            stack.append("(")
            level += 1
        elif c == ")":
            previous = stack[-1]
            if previous == "(":
                # found matching bracket pair
                stack.pop()
                level -= 1
            else:
                raise SyntaxError("No matching ( found.")
    if len(stack) != 0:
        raise SyntaxError("Imbalanced number of brackets.")
    splitpos.append(len(header) + 1)
    parts = []
    for start, end in zip(
        splitpos,
        splitpos[1:],
        strict=False,
    ):
        # remove leading + trailing whitespace
        parts.append(header[start : end - 1].strip())
    return parts


def get_row(row):
    if not isinstance(row, tuple):
        row = (row,)
    return row


@contextlib.contextmanager
def smart_open(filename=None, *args, **kwargs):
    fh = open(filename, *args, **kwargs) if filename and filename != "-" else sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def execute(args):
    aux = read_auxiliary(args.aux)
    overwrite_number = {
        "INFO": dict(args.overwrite_number_info),
        "FORMAT": dict(args.overwrite_number_format),
    }
    with create_reader(
        args.vcf,
        backend=args.backend,
        overwrite_number=overwrite_number,
    ) as vcf:
        expression = preprocess_expression(args.expression, vcf, True)
        if args.long:
            expression = f"SAMPLE, {expression}"
        rows = tableize_vcf(
            vcf,
            expression,
            args.annotation_key,
            overwrite_number=overwrite_number,
            long=args.long,
            auxiliary=aux,
        )

        try:
            with smart_open(args.output, "wt", newline="") as csvfile:
                writer = csv.writer(
                    csvfile,
                    delimiter=args.separator,
                    quoting=csv.QUOTE_MINIMAL,
                )
                if args.header != "none":
                    header = get_header(args, vcf)
                    n_header_cols = len(header)
                    expr_cols = get_toplevel(expression)
                    n_expr_cols = len(expr_cols)
                    if n_header_cols != n_expr_cols:
                        raise HeaderWrongColumnNumberError(
                            n_expr_cols,
                            expr_cols,
                            n_header_cols,
                            header,
                        )
                    writer.writerow(header)
                writer.writerows(get_row(row) for row in rows)
        except VembraneError as ve:
            print(ve, file=stderr)
            sys.exit(1)
