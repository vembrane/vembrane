import sys
import yaml

from sys import stderr
from collections import defaultdict
from typing import Iterator
from itertools import islice, chain
from pysam.libcbcf import VariantFile, VariantHeader, VariantRecord

from ..common import (
    check_expression,
    get_annotation_keys,
    split_annotation_entry,
)
from ..errors import VembraneError
from ..representations import Environment
from .. import __version__


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("filter")
    parser.add_argument(
        "expression",
        type=check_expression,
        help="Filter variants and annotations. If this removes all annotations, "
        "the variant is removed as well.",
    )
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
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
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--statistics",
        "-s",
        metavar="FILE",
        default=None,
        help="Write statistics to this file.",
    )
    parser.add_argument(
        "--keep-unmatched",
        default=False,
        action="store_true",
        help="Keep all annotations of a variant if at least one of them passes "
        "the expression.",
    )
    parser.add_argument(
        "--events",
        default=False,
        action="store_true",
        help="The variant file contains breakend events. This requires two passes, so streaming is not possible.",
    )


def filter_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str,
    keep_unmatched: bool = False,
    bnd_only: bool = False,
    events: set = None,
) -> Iterator[VariantRecord]:

    env = Environment(expression, ann_key, vcf.header)

    record: VariantRecord
    for idx, record in enumerate(vcf):
        svtype = record.info.get("SVTYPE", None)
        event = record.info.get("EVENT", None)

        force_record_output = False

        if svtype == "BND":
            if events is not None:
                if event in events:
                    force_record_output = True
                else:
                    continue
        elif bnd_only:
            continue

        env.update_from_record(idx, record)
        if env.expression_annotations():
            # if the expression contains a reference to the ANN field
            # get all annotations from the record.info field
            # (or supply an empty ANN value if the record has no ANN field)
            try:
                annotations = record.info[ann_key]
            except KeyError:
                annotations = [""]
            #  … and only keep the annotations where the expression evaluates to true
            filtered_annotations = [
                annotation for annotation in annotations if env.evaluate(annotation)
            ]

            if not filtered_annotations:
                # skip this record if filter removed all annotations
                continue
            elif (
                not keep_unmatched and (len(annotations) != len(filtered_annotations))
            ) or force_record_output:
                # update annotations if they have actually been filtered
                record.info[ann_key] = filtered_annotations
            yield record
        else:
            # otherwise, the annotations are irrelevant w.r.t. the expression,
            # so we can omit them
            if force_record_output or env.evaluate():
                yield record
            else:
                continue


def statistics(
    records: Iterator[VariantRecord], vcf: VariantFile, filename: str, ann_key: str
) -> Iterator[VariantRecord]:
    annotation_keys = get_annotation_keys(vcf.header, ann_key)
    counter = defaultdict(lambda: defaultdict(lambda: 0))
    for record in records:
        for annotation in record.info[ann_key]:
            for key, raw_value in zip(
                annotation_keys, split_annotation_entry(annotation)
            ):
                value = raw_value.strip()
                if value:
                    counter[key][value] += 1
        yield record

    # reduce dicts with many items, to just one counter
    for key, subdict in counter.items():
        if len(subdict) > 10:
            counter[key] = f"#{len(subdict)}"

    yaml.add_representer(defaultdict, yaml.representer.Representer.represent_dict)
    with open(filename, "w") as out:
        yaml.dump(dict(counter), out)


def execute(args):
    events = None
    if args.events:
        # first pass
        with VariantFile(args.vcf) as vcf:
            bnds = filter_vcf(
                vcf,
                args.expression,
                args.annotation_key,
                keep_unmatched=args.keep_unmatched,
            )
            events = set(record.info["EVENT"] for record in bnds)

    with VariantFile(args.vcf) as vcf:
        fmt = ""
        if args.output_fmt == "bcf":
            fmt = "b"
        elif args.output_fmt == "uncompressed-bcf":
            fmt = "u"

        header: VariantHeader = vcf.header
        header.add_meta("vembraneVersion", __version__)
        header.add_meta(
            "vembraneCmd",
            "vembrane "
            + " ".join(
                "'" + arg.replace("'", '"') + '"' if " " in arg else arg
                for arg in sys.argv[1:]
            ),
        )

        records = filter_vcf(
            vcf,
            args.expression,
            args.annotation_key,
            keep_unmatched=args.keep_unmatched,
            events=events,
        )

        try:
            first_record = list(islice(records, 1))
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)

        records = chain(first_record, records)

        with VariantFile(
            args.output,
            "w" + fmt,
            header=header,
        ) as out:
            if args.statistics is not None:
                records = statistics(records, vcf, args.statistics, args.annotation_key)

            try:
                for record in records:
                    out.write(record)
            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
