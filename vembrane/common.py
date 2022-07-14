import ast
from typing import List

from vembrane.errors import InvalidExpression

from cyvcf2 import VCF


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


def get_annotation_keys(vcf: VCF, ann_key: str) -> List[str]:
    has_snpeff_annotations = "##SnpEffCmd=" in vcf.raw_header
    has_vep_annotations = "##VEP=" in vcf.raw_header

    # TODO use match statement when minimum required python version is 3.10
    if has_vep_annotations and not has_snpeff_annotations:
        separator = ":"
    elif has_snpeff_annotations and not has_vep_annotations:
        separator = "'"
    elif not has_vep_annotations and not has_vep_annotations:
        separator = ":"
    else:
        raise ValueError(
            "Ambiguous annotation, could be either of SnpEff or VEP annotations"
        )

    for rec in vcf.header_iter():
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
