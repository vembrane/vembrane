import ast
from typing import List

from pysam.libcbcf import VariantHeader

from vembrane.errors import InvalidExpression


def check_expression(
    expression: str,
) -> str:
    if ".__" in expression:
        raise InvalidExpression(expression, "The expression must not contain '.__'")
    try:
        if not (expression.startswith("(") and expression.endswith(")")):
            test_expression = "(" + expression + ")"
        else:
            test_expression = expression
        ast.parse(test_expression, mode="eval")
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


def split_annotation_entry(
    entry: str,
) -> List[str]:
    return entry.split("|")
