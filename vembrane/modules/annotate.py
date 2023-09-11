from collections.abc import Iterable
from sys import stderr
from typing import Any, Callable, Dict, Iterator
import ast

import numpy as np
import yaml
from intervaltree import Interval, IntervalTree
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import check_expression
from ..representations import Environment

from functools import partial


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("annotate")
    parser.add_argument(
        "target",
        # type=check_expression,
        help="Assign the expression values to these targets.",
    )
    parser.add_argument(
        "expression",
        type=check_expression,
        help="Annotate variants and annotations.",
    )
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
    )
    parser.add_argument(
        "--definitions",
        "-d",
        metavar="FILENAME",
        default=None,
        help="A file containing additional info and annotation definitions.",
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )
    parser.add_argument(
        "--output-fmt",
        "-O",
        default="vcf",
        choices=["vcf", "bcf", "uncompressed-bcf"],
        help="Output format.",
    )


typeparser: Dict[str, Callable[[str], Any]] = {
    "Float": float,
    "Integer": int,
    "Character": lambda x: str(x)[0],
    "String": str,
    "Flag": bool,
}


def annotate_vcf(
    vcf: VariantFile,
    target: str,
    expression: str,
    ann_key: str,
) -> Iterator[VariantRecord]:
    env = Environment(expression, ann_key, vcf.header)

    def set_info(record, value, t):
        getattr(record, "info")[t] = value

    target_func = []
    # print(target, file=stderr)
    for e in ast.parse(target).body[0].value.elts:
        if e.value.id == "INFO":
            target_func.append(partial(set_info, t=e.slice.value))
        elif e.value.id == "ANN":
            pass
        else:
            pass  # throw error

    record: VariantRecord
    for idx, record in enumerate(vcf):
        ann_values = env.table()
        for t, v in zip(target_func, ann_values):
            t(record, v)

        # print(getattr(record, "info")["SNVHPOL"], file=stderr)
        # print(ast.parse(target), file=stderr)

        # print(dir(list(ast.iter_child_nodes(ast.parse(target, mode="eval")))[0]), file=stderr)
        # print("ast", dir(ast.parse(target, mode="eval")), file=stderr)
        # print("ann_values 2", target, ann_values, file=stderr)
        # for t, v in zip(ann_values):
        #     print(t, v)
        # for v, expression_value in zip(
        #     map(lambda x: x["info"], config["annotation"]["values"]),
        #     ann_values,
        # ):
        #     if not v["number"] == ".":
        #         number = int(v["number"])
        #     else:
        #         number = -1

        #     parse = typeparser[v["type"]]
        #     if number == -1:
        #         expression_value = list(map(parse, expression_value))
        #     elif number > 1:
        #         expression_value = list(map(parse, expression_value))
        #         assert len(expression_value) == number
        #     else:
        #         # number == 1
        #         assert isinstance(expression_value, str) or not isinstance(
        #             expression_value, Iterable
        #         )
        #         expression_value = parse(expression_value)
        #     record.info[v["vcf_name"]] = expression_value

        yield record


def execute(args):
    with VariantFile(args.vcf) as vcf:
        if args.definitions:
            # add or modify meta defitions
            with open(args.definitions, "r") as stream:
                try:
                    definitions = yaml.safe_load(stream)
                except yaml.YAMLError as e:
                    print(e, file=stderr)
                    exit(1)
            for name in (infos := definitions.get("info", {})):
                if name in vcf.header.info.keys():
                    continue
                values = infos[name]
                vcf.header.add_meta(
                    "INFO",
                    items=[
                        ("ID", name),
                        ("Number", values["number"]),
                        ("Type", values["type"]),
                        ("Description", values["description"]),
                    ],
                )

        expression = f"({args.expression})"
        target = args.target
        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
        with VariantFile(
            args.output,
            f"w{fmt}",
            header=vcf.header,
        ) as out:
            pass
            variants = annotate_vcf(
                vcf,
                target,
                expression,
                args.annotation_key,
            )
            for v in variants:
                out.write(v)

    # expression = []

    # with VariantFile(args.vcf) as vcf:
    #     # add new info
    #     for value in config["annotation"]["values"]:
    #         if (v := value.get("info", None)):
    #             vcf.header.add_meta(
    #                 "INFO",
    #                 items=[
    #                     ("ID", v["vcf_name"]),
    #                     ("Number", v["number"]),
    #                     ("Type", v["type"]),
    #                     ("Description", v["description"]),
    #                 ],
    #             )
    #             expression.append(v["expression"])

    #     expression = ",".join(expression)
    #     expression = f"({expression})"

    #     fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
    #     with VariantFile(
    #         args.output,
    #         f"w{fmt}",
    #         header=vcf.header,
    #     ) as out:
    #         variants = annotate_vcf(
    #             vcf,
    #             expression,
    #             args.annotation_key,
    #             ann_data=ann_data,
    #             config=config,
    #         )
    #         for v in variants:
    #             out.write(v)
