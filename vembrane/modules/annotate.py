from typing import Iterator

import numpy as np
import yaml
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import check_expression
from ..representations import Environment

from intervaltree import Interval, IntervalTree
from collections.abc import Iterable

from sys import stderr


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("annotate")
    parser.add_argument(
        "config",
        type=check_expression,
        help="The configuration file.",
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


typeparser = {
    "Float": float,
    "Integer": int,
    "Character": lambda x: str(x)[0],
    "String": str,
    "Flag": bool,
}


def annotate_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str,
    ann_data,
    config: dict,
) -> Iterator[tuple]:
    env = Environment(expression, ann_key, vcf.header)

    config_chrom_column = config["annotation"]["columns"]["chrom"]
    config_start_column = config["annotation"]["columns"]["start"]
    config_stop_column = config["annotation"]["columns"]["stop"]

    available_chromsomes = set(np.unique(ann_data[config_chrom_column]))

    tree = dict()
    chrom_ann_data = dict()
    for chrom in available_chromsomes:
        d = ann_data[ann_data[config_chrom_column] == chrom]
        chrom_ann_data[chrom] = d
        tree[chrom] = IntervalTree(
            Interval(d[config_start_column], d[config_stop_column], i)
            for i, d in enumerate(d)
        )

        # tree = IntervalTree.from_tuples(interval_tuples))

    current_chrom = None
    current_ann_data = None

    record: VariantRecord
    for idx, record in enumerate(vcf):
        chrom = None
        if current_chrom != record.chrom:
            current_chrom = record.chrom

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
            indices = np.fromiter((i for _, _, i in t[record.start]), dtype=np.int)
            if len(indices):
                env.update_data(current_ann_data[(np.array(indices))])
                env.update_from_record(idx, record)
                ann_values = env.table()

                for v, expression_value in zip(
                    map(lambda x: x["value"], config["annotation"]["values"]),
                    ann_values,
                ):
                    if not v["number"] == ".":
                        number = int(v["number"])
                    else:
                        number = -1

                    parse = typeparser[v["type"]]
                    if number == -1:
                        expression_value = list(map(parse, expression_value))
                    elif number > 1:
                        expression_value = list(map(parse, expression_value))
                        assert len(expression_value) == number
                    else:
                        # number == 1
                        assert isinstance(expression_value, str) or not isinstance(
                            expression_value, Iterable
                        )
                        expression_value = parse(expression_value)
                    record.info[v["vcf_name"]] = expression_value

        yield record


def execute(args):
    with open(args.config, "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as e:
            print(e, file=stderr)
            exit(1)

    # load annotation data
    ann_data = np.genfromtxt(
        config["annotation"]["file"],
        delimiter=config["annotation"].get("delimiter", "\t"),
        names=True,
        dtype=None,
        encoding=None,
    )
    # ann_data = pd.read_csv(config["annotation"]["file"], sep="\t", header=0)
    # ann_data = dict(tuple(ann_data.groupby("chrom")))

    # build expression
    expression = ",".join(
        f'{value["expression"]}'
        for value in map(lambda x: x["value"], config["annotation"]["values"])
    )
    expression = f"({expression})"

    with VariantFile(args.vcf) as vcf:
        # add new info
        for value in config["annotation"]["values"]:
            value = value["value"]
            vcf.header.add_meta(
                "INFO",
                items=[
                    ("ID", value["vcf_name"]),
                    ("Number", value["number"]),
                    ("Type", value["type"]),
                    ("Description", value["description"]),
                ],
            )

        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
        with VariantFile(
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
            )
            for v in variants:
                out.write(v)
