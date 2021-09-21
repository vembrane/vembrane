import contextlib
import csv
import sys
import yaml
import pandas as pd

from sys import stderr
from typing import Iterator, List, Optional

import asttokens
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import check_expression
from ..errors import VembraneError, HeaderWrongColumnNumber
from ..representations import Environment

from itertools import chain
import numpy as np


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


def annotate_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str,
    ann_data: dict,
    config: dict,
) -> Iterator[tuple]:

    env = Environment(expression, ann_key, vcf.header)

    current_chrom = None
    current_index = None
    current_data = None
    indices = []

    record: VariantRecord
    for idx, record in enumerate(vcf):
        if not current_chrom == record.chrom:
            current_chrom = record.chrom
            chrom = None

            # find the correct chrom name
            tmp = current_chrom
            if tmp.lower().startswith("chr"):
                tmp = tmp[3:]
            for prefix in ["", "chr", "Chr", "CHR"]:
                if prefix + tmp in ann_data:
                    chrom = prefix + tmp
                    current_data = ann_data[chrom].to_records()
                    current_index = 0
                    indices = []
        
        if chrom:
            # append possible intervals
            while current_index < len(current_data) and (
                current_data[current_index]["chromStart"] < record.start
            ):
                indices.append(current_index)
                current_index += 1

            # copy only overlapping intervals
            valid_indices = []
            for index in indices:
                if current_data[index]["chromEnd"] > record.start:
                    valid_indices.append(index)

            indices = valid_indices
            if len(indices):
                env.update_data(current_data[(np.array(indices))])
                env.update_from_record(idx, record)
                ann_values = env.table()

                for name, value in zip(
                    map(lambda x: x["value"]["vcf_name"], config["annotation"]["values"]),
                    ann_values,
                ):
                    record.info[name] = float(value)

        yield record


def execute(args):
    with open(args.config, "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # load annotation data
    ann_data = pd.read_csv(config["annotation"]["file"], sep="\t", header=0)
    ann_data = dict(tuple(ann_data.groupby("chrom")))

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

        with VariantFile("out.vcf", "w", header=vcf.header) as o:
            variants = annotate_vcf(
                vcf,
                expression,
                "ann",
                ann_data=ann_data,
                config=config,
            )
            for v in variants:
                o.write(v)
