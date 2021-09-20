import contextlib
import csv
import sys
import yaml
import pandas

from sys import stderr
from typing import Iterator, List, Optional

import asttokens
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import check_expression
from ..errors import VembraneError, HeaderWrongColumnNumber
from ..representations import Environment


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
) -> Iterator[tuple]:

    env = Environment(expression, ann_key, vcf.header)

    current_chrom = None

    record: VariantRecord
    for idx, record in enumerate(vcf):
        if not current_chrom == record.chrom:
            current_chrom = record.chrom
            data_iter = ann_data["chr" + record.chrom]
            env.update_data(data_iter)

        env.update_from_record(idx, record)
    #     if env.expression_annotations():
    #         # if the expression contains a reference to the ANN field
    #         # get all annotations from the record.info field
    #         # (or supply an empty ANN value if the record has no ANN field)
    #         try:
    #             annotations = record.info[ann_key]
    #         except KeyError:
    #             annotations = [""]
    #         for annotation in annotations:
    #             yield env.table(annotation)
    #     else:
    #         yield env.table()


def execute(args):
    with open(args.config, "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # load annotation data    
    ann_data = pandas.read_csv(config["annotation"]["file"], sep="\t", header=0)
    ann_data = dict(tuple(ann_data.groupby('chrom')))

    # build expression
    expression = "\n".join(f'INFO["{value["vcf_name"]}"]={value["expression"]}' for value in map(lambda x: x["value"], config["annotation"]["values"]))
    print(expression)

    with VariantFile(args.vcf) as vcf:
        # add new info
        for value in config["annotation"]["values"]:
            value = value["value"]
            vcf.header.add_meta('INFO', items=[('ID',value["vcf_name"]), ('Number',value["number"]), ('Type','Float'), ('Description', value["description"])])

        with VariantFile("out.vcf", "w", header=vcf.header) as o:
            variants = annotate_vcf(
                vcf,
                expression,
                "ann",
                ann_data=ann_data,
            ) 