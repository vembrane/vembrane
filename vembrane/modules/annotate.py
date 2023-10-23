import re

# from intervaltree import Interval, IntervalTree
from collections.abc import Callable, Iterator
from sys import exit, stderr
from typing import Any

# import sys
# import numpy as np
import yaml

from ..backend.base import VCFHeader, VCFReader, VCFRecord
from ..common import (
    AppendTagExpression,
    add_common_arguments,
    create_reader,
    create_writer,
    get_annotation_description_and_keys,
)
from ..errors import FilterTagNameInvalidError, VembraneError
from ..representations import Environment


def add_subcommmand(subparsers):
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
    add_common_arguments(parser)


typeparser: dict[str, Callable[[str], Any]] = {
    "Float": float,
    "Integer": int,
    "Character": lambda x: str(x)[0],
    "String": str,
    "Flag": bool,
}


def annotate_vcf(
    vcf: VCFReader,
    targets: str,
    expression: str,
    ann_key: str,
) -> Iterator[VCFRecord]:
    env = Environment(expression, ann_key, vcf.header)
    tag_targets, info_targets, ann_targets = targets

    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        tag_values, info_values, ann_values = env.table()
        for target, value in zip(info_targets, info_values, strict=True):
            record.info[target] = value
        # if env.expression_annotations():
        #     try:
        #         annotations = record.info[ann_key]
        #     except KeyError:
        #         num_ann_entries = len(env._annotation._ann_conv.keys())
        #         empty = "|" * num_ann_entries
        #         print(
        #             f"No ANN field found in record {idx}, "
        #             f"replacing with NAs (i.e. 'ANN={empty}')",
        #             file=sys.stderr,
        #         )
        #         annotations = [empty]
        #     new_annotations = []
        #     for annotation in annotations:
        #         ann_values = env.table(annotation)
        #         for t, v in zip(target_func, ann_values):
        #             t(record, v)
        #             new_annotations.append(env._annotation)
        #     print(new_annotations, file=stderr)
        #     # record.info["ANN"] = new_annotations
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
    # build expression
    # expression = ",".join(
    #     f'{value["expression"]}'
    #     for value in (x["value"] for x in config["annotation"]["values"])
    # )
    # expression = f"({expression})"
    tag_keys, tag_expressions = (
        zip(*args.tag.items(), strict=True) if args.tag else ([], [])
    )
    info_keys, info_expressions = (
        zip(*args.info.items(), strict=True) if args.info else ([], [])
    )
    ann_keys, ann_expressions = (
        zip(*args.ann.items(), strict=True) if args.ann else ([], [])
    )

    tag_expressions = "".join(f"({e})," for e in tag_expressions)
    info_expressions = "".join(f"({e})," for e in info_expressions)
    ann_expressions = "".join(f"({e})," for e in ann_expressions)

    expression = f"(({tag_expressions}), ({info_expressions}), ({ann_expressions}))"
    targets = [tag_keys, info_keys, ann_keys]

    with create_reader(
        args.vcf,
        backend=args.backend,
        # overwrite_number=overwrite_number,
    ) as reader:
        header: VCFHeader = reader.header

        definitions = None
        if args.definitions:
            # add or modify meta defitions
            with open(args.definitions, "r") as stream:
                try:
                    definitions = yaml.safe_load(stream)
                except yaml.YAMLError as e:
                    print(e, file=stderr)
                    exit(1)

        for key in tag_keys:
            # we just go on if the definition already exists
            # maybe we should throw an error?
            if key in header.filters:
                continue

            # definition does not exist in header ...
            if not definitions:
                # ... and no definitions are provided
                raise VembraneError(
                    f"""Error: please provide definition file by
                    -d including a definition for tag {key}""",
                )

            if not (definition := definitions.get("tag", {}).get(key, None)):
                # ... and no definitions for tag in file
                raise VembraneError(
                    f"Error: No definition found in definition file for tag {key}",
                )

            # all correct, add filter
            header.add_filter(key, definition["description"])

        for key in info_keys:
            # we just go on if the definition already exists
            # maybe we should throw an error?
            if key in header.infos:
                continue

            # definition does not exist in header ...
            if not definitions:
                # ... and no definitions are provided
                raise VembraneError(
                    f"""Error: please provide definition file by
                    -d including a definition for info {key}""",
                )

            if not (definition := definitions.get("info", {}).get(key, None)):
                # ... and no definitions for tag in file
                raise VembraneError(
                    f"Error: No definition found in definition file for tag {key}",
                )

            # all correct, add info
            header.add_info(
                key,
                definition["number"],
                definition["type"],
                definition["description"],
            )

        for key in ann_keys:
            ann_key = args.annotation_key
            if ann_key not in header.infos:
                raise VembraneError(
                    f"""Error: Currently Vembrane needs an existing
                    annotation field with key {key} to add an annotation""",
                )

            ann_description, ann_keys = get_annotation_description_and_keys(
                header,
                ann_key,
            )

            new_ann_keys = []
            # add annotation definitions
            # here things become complicated
            for key, _ in definitions.get("ann", {}).items():
                new_ann_keys.append(key)

            new_ann_keys = ann_keys + new_ann_keys
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
                targets,
                expression,
                args.annotation_key,
            )
            for record in records:
                writer.write(record)
