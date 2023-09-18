from collections.abc import Iterable
from sys import stderr
from typing import Any, Callable, Dict, Iterator
import ast

import numpy as np
import yaml
from intervaltree import Interval, IntervalTree
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import check_expression, add_annotations
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
    env = Environment(expression, ann_key, vcf.header, target=target)

    def set_info(record, value, t):
        getattr(record, "info")[t] = value

    def set_ann(record, value, t):
        # print(env._annotation, file=stderr)
        # getattr(record, "info")["ANN"] = str(env._annotation)
        pass

    target_func = []
    for e in ast.parse(target).body[0].value.elts:
        if e.value.id == "INFO":
            target_func.append(partial(set_info, t=e.slice.value))
        elif e.value.id == "ANN":
            # add new annotations definitions
            env._annotation.add(e.slice.value)
            target_func.append(partial(set_ann, t=e.slice.value))
        else:
            pass  # throw error

    record: VariantRecord
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        if env.expression_annotations():
            try:
                annotations = record.info[ann_key]
            except KeyError:
                num_ann_entries = len(env._annotation._ann_conv.keys())
                empty = "|" * num_ann_entries
                print(
                    f"No ANN field found in record {idx}, "
                    f"replacing with NAs (i.e. 'ANN={empty}')",
                    file=sys.stderr,
                )
                annotations = [empty]
            new_annotations = []
            for annotation in annotations:
                ann_values = env.table(annotation)
                for t, v in zip(target_func, ann_values):
                    t(record, v)
                    new_annotations.append(env._annotation)
            print(new_annotations, file=stderr)
            # record.info["ANN"] = new_annotations
        yield record


def execute(args):
    with VariantFile(args.vcf) as vcf:
        header = vcf.header
        if args.definitions:
            # add or modify meta defitions
            with open(args.definitions, "r") as stream:
                try:
                    definitions = yaml.safe_load(stream)
                except yaml.YAMLError as e:
                    print(e, file=stderr)
                    exit(1)
            for name in (infos := definitions.get("info", {})):
                if name in header.info.keys():
                    continue
                values = infos[name]
                header.add_meta(
                    "INFO",
                    items=[
                        ("ID", name),
                        ("Number", values["number"]),
                        ("Type", values["type"]),
                        ("Description", values["description"]),
                    ],
                )

            if (annotations := definitions.get("info", {})):
                header = add_annotations(vcf.header, args.annotation_key, annotations)
                # h = vcf.header.copy()
                # h.info.remove_header(args.annotation_key)
                # h = h.copy()
                # h.info.add(args.annotation_key,".","String",'whatever')
                # print("foooo", h, file=sys.stderr)
                # exit()
                # header=copy(vcf.header)
            # for name in (ann := definitions.get("ann", {})):
            #     add_annotation_key(vcf.header, args.annotation_key, name)

        expression = f"({args.expression})"
        target = args.target
        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]

        with VariantFile(
            args.output,
            f"w{fmt}",
            header=vcf.header,

        ) as out:
            print(dir(out.header), file=stderr)
            for v in vcf:
                v.info["ANN"] = "test"
                out.write(v)
                # print(v, file=stderr)
            # variants = annotate_vcf(
            #     vcf,
            #     target,
            #     expression,
            #     args.annotation_key,
            # )
            # for v in variants:
            #     out.write(v)
