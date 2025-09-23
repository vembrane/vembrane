from collections.abc import Callable, Iterable, Iterator
from sys import exit, stderr
from typing import Any

import numpy as np
import yaml
from intervaltree import Interval, IntervalTree  # type: ignore

from vembrane.errors import handle_vembrane_error  # type: ignore

from ..backend.base import VCFReader, VCFRecord
from ..common import (
    Context,
    HumanReadableDefaultsFormatter,
    add_common_arguments,
    check_expression,
    create_reader,
    create_writer,
)
from ..representations import FuncWrappedExpressionEnvironment


def add_subcommmand(subparsers):
    parser = subparsers.add_parser(
        "annotate",
        help="Add new INFO field annotations to a VCF/BCF from other data sources, "
        "using a configuration file.",
        description="Add new INFO field annotations to a VCF/BCF from other data "
        "sources, using a configuration file.",
        formatter_class=HumanReadableDefaultsFormatter,
    )
    parser.add_argument(
        "config",
        type=check_expression,
        help="The configuration file.",
    )
    parser.add_argument(
        "vcf",
        help="Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin.",
        nargs="?",
        default="-",
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
    expression: str,
    ann_key: str,
    ann_data,
    config: dict,
    auxiliary_globals: dict[str, Any] | None = None,
) -> Iterator[VCFRecord]:
    env = FuncWrappedExpressionEnvironment(
        expression, ann_key, vcf.header, auxiliary_globals=auxiliary_globals
    )

    config_chrom_column: str = config["annotation"]["columns"]["chrom"]
    config_start_column: str = config["annotation"]["columns"]["start"]
    config_stop_column: str = config["annotation"]["columns"]["stop"]

    available_chromsomes: set[str] = set(np.unique(ann_data[config_chrom_column]))

    tree: dict[str, IntervalTree] = {}
    chrom_ann_data: dict[str, Any] = {}
    for available_chrom in available_chromsomes:
        d = ann_data[ann_data[config_chrom_column] == available_chrom]
        chrom_ann_data[available_chrom] = d
        tree[available_chrom] = IntervalTree(
            Interval(d[config_start_column], d[config_stop_column], i)
            for i, d in enumerate(d)
        )

    current_chrom = None
    current_ann_data = None

    record: VCFRecord
    for idx, record in enumerate(vcf):
        chrom: str | None = None
        if current_chrom != record.contig:
            current_chrom = record.contig

            # find the correct chrom name
            tmp = current_chrom
            if tmp.lower().startswith("chr"):
                tmp = tmp[3:]
            for prefix in ["", "chr", "Chr", "CHR"]:
                if prefix + tmp in available_chromsomes:
                    chrom = prefix + tmp
                    current_ann_data = chrom_ann_data[chrom]
                    t = tree[chrom]

        if chrom:
            indices = np.fromiter((i for _, _, i in t[record.position]), dtype=int)
            if len(indices):
                env.update_data(current_ann_data[(np.array(indices))])  # type: ignore
                env.update_from_record(idx, record)
                ann_values = env.table_row()

                for v, expression_value in zip(
                    (x["value"] for x in config["annotation"]["values"]),
                    ann_values,
                    strict=True,
                ):
                    number = int(v["number"]) if v["number"] != "." else -1

                    parse = typeparser[v["type"]]
                    if number == -1:
                        ev = list(map(parse, expression_value))
                    elif number > 1:
                        ev = list(map(parse, expression_value))
                        assert len(ev) == number
                    else:
                        # number == 1
                        assert isinstance(expression_value, str) or not isinstance(
                            expression_value,
                            Iterable,
                        )
                        ev = parse(expression_value)
                    record.info[v["vcf_name"]] = ev

        yield record


@handle_vembrane_error
def execute(args):
    with open(args.config) as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as e:
            print(e, file=stderr)
            exit(1)

    # load annotation data
    ann_data: np.ndarray = np.genfromtxt(
        config["annotation"]["file"],
        delimiter=config["annotation"].get("delimiter", "\t"),
        names=True,
        dtype=None,
        encoding=None,
    )

    # build expression
    expression = ",".join(
        f"{value['expression']}"
        for value in (x["value"] for x in config["annotation"]["values"])
    )
    expression = f"({expression})"

    with create_reader(args.vcf) as vcf:
        # add new info
        for value in config["annotation"]["values"]:
            v = value["value"]
            vcf.header.add_meta(
                "INFO",
                items=[
                    ("ID", v["vcf_name"]),
                    ("Number", v["number"]),
                    ("Type", v["type"]),
                    ("Description", v["description"]),
                ],
            )

        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
        with create_writer(
            args.output,
            f"w{fmt}",
            header=vcf.header,
        ) as out:
            variants = annotate_vcf(
                vcf,
                expression,
                args.annotation_key,
                ann_data=ann_data,
                config=config,
                auxiliary_globals=Context.from_args(args).get_globals(),
            )
            for v in variants:
                out.write(v)
