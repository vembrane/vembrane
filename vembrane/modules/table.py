import contextlib
import csv
import sys
from sys import stderr
from typing import Iterator, List, Optional

import asttokens
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import check_expression
from ..errors import VembraneError, HeaderWrongColumnNumber
from ..representations import Environment


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("table")
    parser.add_argument(
        "expression",
        type=check_expression,
        help="The expression for the output.",
    )
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
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
        "--all",
        "-a",
        default=False,
        action="store_true",
        help="Do not filter duplicate entries.",
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
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )


def tableize_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str,
) -> Iterator[tuple]:
    expression = f"({expression})"
    env = Environment(expression, ann_key, vcf.header)

    record: VariantRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.expression_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            try:
                annotations = record.info[ann_key]
            except KeyError:
                annotations = [""]
            for annotation in annotations:
                yield env.table(annotation)
        else:
            yield env.table()


def generate_for_each_sample_expressions(s: str, vcf: VariantFile) -> List[str]:
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


def generate_for_each_sample_column_names(s: str, vcf: VariantFile) -> List[str]:
    # parse the `for_each_sample(lambda var: inner) expression
    var, inner = _var_and_body(s)

    samples = list(vcf.header.samples)
    from vembrane.globals import allowed_globals

    __globals = allowed_globals.copy()

    column_names = []
    for sample in samples:
        __globals[var] = sample
        column_name = eval(inner, __globals, {})
        if not isinstance(column_name, (str, bytes)):
            if hasattr(column_name, "__str__"):
                column_name = str(column_name)
            else:
                raise ValueError(
                    "The specified header expression does not evaluate to a string."
                    "Consider using `str(expression)` instead."
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
        if isinstance(node, ast.Call):
            if hasattr(node.func, "id"):
                if node.func.id == "for_each_sample":
                    for arg in node.args:
                        if isinstance(arg, ast.Lambda):
                            var = tok.get_text(arg.args)
                            inner = tok.get_text(arg.body)
    return var, inner


def preprocess_header_expression(
    header: str, vcf: Optional[VariantFile] = None, make_expression: bool = True
) -> str:
    """
    Split the header expression at toplevel commas into parts.
    Then, if one of these parts starts with 'for_each_sample',
    that part is expanded for each sample in vcf.header.samples
    """
    parts = get_toplevel(header)
    to_expand = list(
        filter(lambda x: x[1].startswith("for_each_sample"), enumerate(parts))
    )
    if len(to_expand) > 0 and vcf is None:
        raise ValueError("If FORMAT is to be expanded, the VCF kwarg must not be none.")
    parts = [[p] for p in parts]
    func = (
        generate_for_each_sample_expressions
        if make_expression
        else generate_for_each_sample_column_names
    )
    for i, p in to_expand:
        expanded = func(p, vcf)
        parts[i] = expanded
    parts = [p for pp in parts for p in pp]
    return ", ".join(parts)


def get_header(args, vcf: Optional[VariantFile] = None) -> List[str]:
    if args.header == "auto":
        header = args.expression
    else:
        header = args.header
    return get_toplevel(
        preprocess_header_expression(header, vcf, args.header == "auto")
    )


def get_toplevel(header: str) -> List[str]:
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
    for start, end in zip(splitpos, splitpos[1:]):
        # remove leading + trailing whitespace
        parts.append(header[start : end - 1].strip())
    return parts


def get_row(row):
    if not isinstance(row, tuple):
        row = (row,)
    return row


@contextlib.contextmanager
def smart_open(filename=None, *args, **kwargs):
    if filename and filename != "-":
        fh = open(filename, *args, **kwargs)
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def execute(args):
    with VariantFile(args.vcf) as vcf:
        expression = preprocess_header_expression(args.expression, vcf, True)
        rows = tableize_vcf(
            vcf,
            expression,
            args.annotation_key,
        )

        try:
            with smart_open(args.output, "wt", newline="") as csvfile:
                writer = csv.writer(
                    csvfile, delimiter=args.separator, quoting=csv.QUOTE_MINIMAL
                )
                if args.header != "none":
                    header = get_header(args, vcf)
                    n_header_cols = len(header)
                    expr_cols = expression.split(", ")
                    n_expr_cols = len(expr_cols)
                    if n_header_cols != n_expr_cols:
                        raise HeaderWrongColumnNumber(
                            n_expr_cols, expr_cols, n_header_cols, header
                        )
                    writer.writerow(header)
                writer.writerows(get_row(row) for row in rows)
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)
