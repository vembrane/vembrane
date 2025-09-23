import csv
from collections.abc import Iterator
from enum import Enum
from typing import Any

import asttokens

from ..backend.base import VCFHeader, VCFReader, VCFRecord
from ..common import (
    Context,
    HumanReadableDefaultsFormatter,
    add_common_arguments,
    check_expression,
    create_reader,
    get_annotation_keys,
    read_auxiliary,
    read_ontology,
    smart_open,
)
from ..errors import HeaderWrongColumnNumberError, handle_vembrane_error
from ..globals import default_allowed_globals
from ..representations import FuncWrappedExpressionEnvironment
from ..sequence_ontology import SequenceOntology
from .filter import DeprecatedAction

ALL_EXPRESSION = "ALL"


def add_subcommmand(subparsers):
    parser = subparsers.add_parser(
        "table",
        help="Convert VCF/BCF records to tabular format.",
        description="Convert VCF/BCF records to tabular format.",
        formatter_class=HumanReadableDefaultsFormatter,
    )
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "expression",
        type=check_expression,
        help="A comma-separated tuple of expressions "
        "that define the table column contents. "
        f"Use {ALL_EXPRESSION} to output all fields.",
    )
    parser.add_argument(
        "vcf",
        help="Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin.",
        nargs="?",
        default="-",
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
        "--naming-convention",
        metavar="CONVENTION",
        default="dictionary",
        type=NamingConvention,
        choices=list(NamingConvention),
        help="The naming convention to use for column names "
        "when generating the header for the ALL expression.",
    )
    parser.add_argument(
        "--wide",
        help="Instead of using long format with a special SAMPLE column, "
        "generate multiple columns per sample "
        "with the `for_each_sample` utility function.",
        action="store_true",
    )
    parser.add_argument(
        "--long",
        action=DeprecatedAction,
        help="Long format is now the default. For wide format, use `--wide` instead.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )
    add_common_arguments(parser)


def tableize_vcf(
    vcf: VCFReader,
    expression: str,
    ann_key: str,
    overwrite_number: dict[str, dict[str, str]] | None = None,
    wide: bool = False,
    auxiliary: dict[str, set[str]] | None = None,
    auxiliary_globals: dict[str, Any] | None = None,
    ontology: SequenceOntology | None = None,
) -> Iterator[tuple]:
    if overwrite_number is None:
        overwrite_number = {}
    if auxiliary is None:
        auxiliary = {}

    kwargs: dict[str, Any] = dict(
        auxiliary_globals=auxiliary_globals, auxiliary=auxiliary, ontology=ontology
    )

    long_with_samples = not wide and list(vcf.header.samples)
    if long_with_samples:
        kwargs["evaluation_function_template"] = (
            "lambda: (({expression}) for SAMPLE in SAMPLES)"
        )
    else:
        expression = f"({expression})"
    env = FuncWrappedExpressionEnvironment(expression, ann_key, vcf.header, **kwargs)

    record: VCFRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.expression_annotations():
            annotations = env.get_record_annotations(idx, record)
            for annotation in annotations:
                if long_with_samples:
                    yield from env.table_row(annotation)
                else:
                    yield env.table_row(annotation)
        else:
            if long_with_samples:
                yield from env.table_row()
            else:
                yield env.table_row()


class NamingConvention(Enum):
    DICTIONARY = "dictionary"
    UNDERSCORE = "underscore"
    SLASH = "slash"


def construct_all_expression_and_header(
    header: VCFHeader, ann_key: str, naming_convention: NamingConvention
) -> tuple[str, list[str]]:
    all_header = ["SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "ID"]
    all_expression = ", ".join(all_header)

    def format_key(field: str, key: str) -> str:
        match naming_convention:
            case NamingConvention.DICTIONARY:
                return f"{field}['{key}']"
            case NamingConvention.UNDERSCORE:
                return f"{field}_{key}"
            case NamingConvention.SLASH:
                return f"{field}/{key}"

    info_keys = list(header.infos.keys())
    if info_keys:
        info_keys_without_ann = [key for key in info_keys if key != ann_key]
        info_expr = ", ".join(f'INFO["{key}"]' for key in info_keys_without_ann)
        all_header += [format_key("INFO", key) for key in info_keys_without_ann]
        all_expression += ", " + info_expr

    format_keys = list(header.formats.keys())
    if format_keys:
        all_header += [format_key("FORMAT", key) for key in format_keys]
        format_expr = ", ".join(
            f'(FORMAT["{key}"] or {{}}).get(SAMPLE, NA)' for key in format_keys
        )
        all_expression += ", " + format_expr

    annotation_keys = get_annotation_keys(header, ann_key)
    if annotation_keys:
        annotation_expr = ", ".join(f'ANN["{key}"]' for key in annotation_keys)
        all_header += [format_key(ann_key, key) for key in annotation_keys]
        all_expression += ", " + annotation_expr

    return all_expression, all_header


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
    __globals = default_allowed_globals.copy()

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
    expression: str,
    vcf: VCFReader,
    make_expression: bool = True,
) -> str:
    """
    Split the header expression at toplevel commas into parts.
    Then, if one of these parts starts with 'for_each_sample',
    that part is expanded for each sample in vcf.header.samples
    """

    if expression == ALL_EXPRESSION:
        return expression

    parts: list[str] = get_toplevel(expression)
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
    if args.expression == ALL_EXPRESSION:
        _, all_header = construct_all_expression_and_header(
            vcf.header, args.annotation_key, args.naming_convention
        )
        return all_header
    header = args.expression if args.header == "auto" else args.header
    if not args.wide:
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


@handle_vembrane_error
def execute(args):
    aux = read_auxiliary(args.aux)
    ontology = read_ontology(args.ontology)
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
        if expression == ALL_EXPRESSION:
            args.wide = False
            expression, _header = construct_all_expression_and_header(
                vcf.header, args.annotation_key, args.naming_convention
            )
        else:
            if not args.wide:
                if list(vcf.header.samples):
                    expression = f"SAMPLE, {expression}"
                else:
                    # if there are no samples, SAMPLE will be undefined;
                    # add an empty string instead
                    expression = f"'', {expression}"
        rows = tableize_vcf(
            vcf,
            expression,
            args.annotation_key,
            overwrite_number=overwrite_number,
            wide=args.wide,
            auxiliary=aux,
            auxiliary_globals=Context.from_args(args).get_globals(),
            ontology=ontology,
        )

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
