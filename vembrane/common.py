import ast
from typing import List

from vembrane.errors import InvalidExpression

from cyvcf2.cyvcf2 import HREC


def check_expression(
    expression: str,
) -> str:
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


def get_annotation_keys(header: List[HREC], ann_key: str) -> List[str]:
    separator = "'"
    for rec in header:
        # if rec.key == "VEP":
        if True:
            separator = ":"
            continue
        if rec.info().get("ID") == ann_key:
            return list(
                map(
                    str.strip,
                    rec.info()
                    .get("Description")
                    .strip('"')
                    .split(separator)[1]
                    .split("|"),
                )
            )
    return []


def split_annotation_entry(
    entry: str,
) -> List[str]:
    return entry.split("|")
