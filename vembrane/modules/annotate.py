import re
import sys

# from intervaltree import Interval, IntervalTree
from collections.abc import Callable, Iterator
from sys import exit, stderr
from types import MappingProxyType
from typing import Any, List

# import sys
# import numpy as np
import yaml

from .. import __version__
from ..backend.base import VCFHeader, VCFReader, VCFRecord
from ..common import (
    AppendKeyValuePair,
    AppendTagExpression,
    Auxiliary,
    add_common_arguments,
    create_reader,
    create_writer,
    get_annotation_description_and_keys,
    read_auxiliary,
)
from ..errors import FilterTagNameInvalidError, VembraneError
from ..representations import Environment


def add_subcommand(subparsers):
    parser = subparsers.add_parser("annotate")
    parser.add_argument(
        "vcf",
        help="The file containing the variants.",
        nargs="?",
        default="-",
    )
    parser.add_argument(
        "--tag",
        "-t",
        help="Tag records using the FILTER field.\n"
        "Note: tag names cannot contain `;` or whitespace and must not be '0'.\n"
        'Example: `--tag q_above_30="not (QUAL<=30)"`',
        nargs=1,
        metavar="TAG=EXPRESSION",
        action=AppendTagExpression,
    )
    parser.add_argument(
        "--info",
        "-i",
        help="Add info to each record.\n"
        "Note: info names cannot contain `;` or whitespace and must not be '0'.\n"
        'Example: `--info log_qual="log(QUAL)"`',
        nargs=1,
        metavar="KEY=EXPRESSION",
        action=AppendTagExpression,
    )
    parser.add_argument(
        "--ann",
        "-a",
        help="Add an annotation to each record.\n"
        "Note: annotation names cannot contain `;` or whitespace and must not be '0'.\n"
        'Example: `--ann log_qual="log(QUAL)"`',
        nargs=1,
        metavar="KEY=EXPRESSION",
        action=AppendTagExpression,
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--aux",
        nargs=1,
        action=AppendKeyValuePair,
        metavar="NAME=PATH",
        default={},
        help="Path to an auxiliary file containing a set of symbols",
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
    parser.add_argument(
        "--tag-mode",
        "-m",
        default="pass",
        choices=["pass", "fail"],
        help="Set, whether to tag records that pass the tag expression(s), "
        "or records that fail them."
        "By default, vembrane tags records for which the tag expression(s) pass. "
        "This allows for descriptive tag names such as `q_at_least_30`, "
        "which would correspond to an expression `QUAL >= 30`. "
        "However, the VCF specification (`v4.4`) defines tags to be set when a "
        "filter expression is failed, so vembrane also offers the `fail` mode.",
    )
    add_common_arguments(parser)


typeparser: dict[str, Callable[[str], Any]] = {
    "Float": float,
    "Integer": int,
    "Character": lambda x: str(x)[0],
    "String": str,
    "Flag": bool,
}


def annotate_vcf(
    reader: VCFReader,
    tag_targets: List[str],
    info_targets: List[str],
    ann_targets: List[str],
    tag_expressions: List[str],
    info_expressions: List[str],
    ann_expressions: List[str],
    ann_key: str,
    ann_keys: List[str],
    invert_tag_expression: bool = False,
    auxiliary: dict[str, Auxiliary] = MappingProxyType({}),
) -> Iterator[VCFRecord]:
    tag_expressions = "".join(f"({e})," for e in tag_expressions)
    info_expressions = "".join(f"({e})," for e in info_expressions)
    ann_expressions = "".join(f"({e})," for e in ann_expressions)
    if not ann_expressions:
        ann_expressions = "()"

    tag_info_expressions = (
        f"(({tag_expressions}), ({info_expressions}))"  # ({ann_expressions})
    )
    ann_expressions = f"({ann_expressions})"

    env = Environment(
        [tag_info_expressions, ann_expressions],
        ann_key,
        reader.header,
        auxiliary=auxiliary,
    )

    if ann_keys:
        ann_keys_dict = {key: i for i, key in enumerate(ann_keys)}

    for idx, record in enumerate(reader):
        env.update_from_record(idx, record)
        tag_values, info_values = env.table()

        for target, value in zip(info_targets, info_values, strict=True):
            record.info[target] = value

        for target, value in zip(tag_targets, tag_values, strict=True):
            if invert_tag_expression:
                value = not value
            if value:
                record.filter.add(target)

        if ann_targets:
            if ann_key not in record.info:
                print(
                    f"No ANN field found in record {idx}",
                    file=sys.stderr,
                )
            else:
                new_annotations = []
                for a in record.info[ann_key]:
                    ann_values = env.table(a, n=1)
                    # print(ann_values, file=sys.stderr)
                    a_values = a.split("|")
                    a_values.extend((len(ann_keys) - len(a_values)) * [""])
                    for target, value in zip(ann_targets, ann_values, strict=True):
                        a_values[ann_keys_dict[target]] = str(value)
                    new_annotations.append("|".join(a_values))
                record.info[ann_key] = new_annotations
        yield record

    # config_chrom_column: str = config["annotation"]["columns"]["chrom"]
    # config_start_column: str = config["annotation"]["columns"]["start"]
    # config_stop_column: str = config["annotation"]["columns"]["stop"]

    # available_chromsomes: set[str] = set(np.unique(ann_data[config_chrom_column]))

    # tree: dict[str, IntervalTree] = {}
    # chrom_ann_data: dict[str, Any] = {}
    # for chrom in available_chromsomes:
    #     d = ann_data[ann_data[config_chrom_column] == chrom]
    #     chrom_ann_data[chrom] = d
    #     tree[chrom] = IntervalTree(
    #         Interval(d[config_start_column], d[config_stop_column], i)
    #         for i, d in enumerate(d)
    #     )

    # current_chrom = None
    # current_ann_data = None

    # record: VCFRecord
    # for idx, record in enumerate(vcf):
    #     chrom = None
    #     if current_chrom != record.chrom:
    #         current_chrom = record.chrom

    #         # find the correct chrom name
    #         tmp = current_chrom
    #         if tmp.lower().startswith("chr"):
    #             tmp = tmp[3:]
    #         for prefix in ["", "chr", "Chr", "CHR"]:
    #             if prefix + tmp in available_chromsomes:
    #                 chrom = prefix + tmp
    #                 current_ann_data = chrom_ann_data[chrom]
    #                 t = tree[chrom]

    #     if chrom:
    #         indices = np.fromiter((i for _, _, i in t[record.start]), dtype=int)
    #         if len(indices):
    #             env.update_data(current_ann_data[(np.array(indices))])
    #             env.update_from_record(idx, record)
    #             ann_values = env.table()

    #             for v, expression_value in zip(
    #                 (x["value"] for x in config["annotation"]["values"]),
    #                 ann_values,
    #                 strict=True,
    #             ):
    #                 number = int(v["number"]) if v["number"] != "." else -1

    #                 parse = typeparser[v["type"]]
    #                 if number == -1:
    #                     ev = list(map(parse, expression_value))
    #                 elif number > 1:
    #                     ev = list(map(parse, expression_value))
    #                     assert len(ev) == number
    #                 else:
    #                     # number == 1
    #                     assert isinstance(expression_value, str) or not isinstance(
    #                         expression_value,
    #                         Iterable,
    #                     )
    #                     ev = parse(expression_value)
    #                 record.info[v["vcf_name"]] = ev

    #     yield record


def check_tag(tag: str):
    if re.search(r"^0$|[\s;]", tag):
        raise FilterTagNameInvalidError(tag)


def execute(args):
    aux = read_auxiliary(args.aux)
    # build expression
    # expression = ",".join(
    #     f'{value["expression"]}'
    #     for value in (x["value"] for x in config["annotation"]["values"])
    # )
    # expression = f"({expression})"
    tag_targets, tag_expressions = (
        zip(*args.tag.items(), strict=True) if args.tag else ([], [])
    )
    info_targets, info_expressions = (
        zip(*args.info.items(), strict=True) if args.info else ([], [])
    )
    ann_targets, ann_expressions = (
        zip(*args.ann.items(), strict=True) if args.ann else ([], [])
    )

    with create_reader(
        args.vcf,
        backend=args.backend,
        # overwrite_number=overwrite_number,
    ) as reader:
        header: VCFHeader = reader.header
        reader.header.add_generic("vembraneVersion", __version__)

        definitions = None
        if args.definitions:
            # add or modify meta defitions
            with open(args.definitions, "r") as stream:
                try:
                    definitions = yaml.safe_load(stream)
                except yaml.YAMLError as e:
                    print(e, file=stderr)
                    exit(1)

        # add tag definitions
        for target in tag_targets:
            # we just go on if the definition already exists
            # maybe we should throw an error?
            if target in header.filters:
                continue

            if not definitions or not (
                definition := definitions.get("tag", {}).get(target, None)
            ):
                header.add_filter(target, "")  # ADD definition without description
                # TODO: show a warning
                continue

            # # definition does not exist in header ...
            # if not definitions:
            #     # ... and no definitions are provided
            #     raise VembraneError(
            #         f"""Error: please provide definition file by
            #         -d including a definition for tag {target}""",
            #     )

            # if not (definition := definitions.get("tag", {}).get(target, None)):
            #     # ... and no definitions for tag in file
            #     raise VembraneError(
            #         f"Error: No definition found in definition file for tag {target}",
            #     )

            # all correct, add filter
            header.add_filter(target, definition["description"])

        # add info definitions
        for target in info_targets:
            # we just go on if the definition already exists
            # maybe we should throw an error?
            if target in header.infos:
                continue

            # definition does not exist in header ...
            if not definitions:
                # ... and no definitions are provided
                raise VembraneError(
                    f"""Error: please provide definition file by
                    -d including a definition for info {target}""",
                )

            if not (definition := definitions.get("info", {}).get(target, None)):
                # ... and no definitions for tag in file
                raise VembraneError(
                    f"Error: No definition found in definition file for info {target}",
                )

            # all correct, add info
            header.add_info(
                target,
                definition["number"],
                definition["type"],
                definition["description"],
            )

        # add annotation definitions
        new_ann_keys = None
        if ann_targets:
            ann_key = args.annotation_key
            ann_description, existing_ann_keys = get_annotation_description_and_keys(
                header,
                ann_key,
            )
            if ann_key not in header.infos:
                raise VembraneError(
                    f"""Error: Currently Vembrane needs an existing
                    annotation field with key {target} to add an annotation""",
                )
            new_ann_keys = existing_ann_keys
            for target in ann_targets:
                # we just go on if the definition already exists
                # maybe we should throw an error?
                if target in header.infos:
                    continue

                # definition does not exist in header ...
                if not definitions:
                    # ... and no definitions are provided
                    raise VembraneError(
                        f"""Error: please provide definition file by
                        -d including a definition for ann {target}""",
                    )

                if not (definition := definitions.get("ann", {}).get(target, None)):
                    # ... and no definitions for tag in file
                    raise VembraneError(
                        f"""
                        Error: No definition found in
                        definition file for ann {target}
                        """,
                    )

                new_ann_keys.append(target)
                ann_meta = header.infos[ann_key]
                header.update_info(
                    ann_key,
                    ann_meta["Number"],
                    ann_meta["Type"],
                    ann_description.format(keys=" | ".join(new_ann_keys)),
                )

        # #TODO?
        # FilterAlreadyDefinedError(tag_key)
        # check_expression(expression)
        # expression = swap_quotes(expr) if single_outer(expr) else expr

        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
        with create_writer(args.output, fmt, reader, backend=args.backend) as writer:
            records = annotate_vcf(
                reader,
                tag_targets,
                info_targets,
                ann_targets,
                tag_expressions,
                info_expressions,
                ann_expressions,
                args.annotation_key,
                ann_keys=new_ann_keys,
                invert_tag_expression=(args.tag_mode == "fail"),
                auxiliary=aux,
            )
            for record in records:
                writer.write(record)
