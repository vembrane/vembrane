import argparse
import sys
from collections import defaultdict
from itertools import chain, islice
from sys import stderr
from typing import Dict, Iterator, Set, Tuple

import yaml
from pysam.libcbcf import VariantFile, VariantHeader, VariantRecord

from .. import __version__
from ..common import (
    BreakendEvent,
    check_expression,
    get_annotation_keys,
    mate_key,
    split_annotation_entry,
)
from ..errors import VembraneError
from ..representations import Environment


class DeprecatedAction(argparse.Action):
    def __call__(self, *args, **kwargs):
        print(
            f"{'/'.join(self.option_strings)} is deprecated.\n{self.help}",
            file=sys.stderr,
        )
        exit(1)


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("filter")
    parser.register("action", "deprecated", DeprecatedAction)
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
        "--aux",
        "-a",
        nargs=2,
        action="append",
        metavar=("NAME", "PATH"),
        default=[],
        help="Path to an auxiliary file containing a set of symbols",
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
        "--preserve-order",
        default=False,
        action="store_true",
        help="Ensures that the order of the output matches that of the \
              input. This is only useful if the input contains breakends (BNDs) \
              since the order of all other variants is preserved anyway.",
    )
    parser.add_argument(
        "--overwrite-number",
        help="Deprecated. "
        "Use --overwrite-number-info or --overwrite-number-format instead.",
        action="deprecated",
    )
    parser.add_argument(
        "--overwrite-number-info",
        nargs=2,
        action="append",
        metavar=("FIELD", "NUMBER"),
        default=[],
        help="Overwrite the number specification for INFO fields "
        "given in the VCF header. "
        "Example: `--overwrite-number cosmic_CNT .`",
    )
    parser.add_argument(
        "--overwrite-number-format",
        nargs=2,
        action="append",
        metavar=("FIELD", "NUMBER"),
        default=[],
        help="Overwrite the number specification for FORMAT fields "
        "given in the VCF header. "
        "Example: `--overwrite-number-format DP 2`",
    )


def test_and_update_record(
    env: Environment,
    idx: int,
    record: VariantRecord,
    ann_key: str,
    keep_unmatched: bool,
) -> (VariantRecord, bool):
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
        if keep_unmatched:
            filtered = any(map(env.evaluate, annotations))
            return record, filtered
        else:
            filtered_annotations = [
                annotation for annotation in annotations if env.evaluate(annotation)
            ]

            if len(annotations) != len(filtered_annotations):
                # update annotations if they have actually been filtered
                record.info[ann_key] = filtered_annotations

            return record, len(filtered_annotations) > 0
    else:
        # otherwise, the annotations are irrelevant w.r.t. the expression,
        # so we can omit them
        return record, env.evaluate()


def filter_vcf(
    vcf: VariantFile,
    expression: str,
    ann_key: str,
    keep_unmatched: bool = False,
    preserve_order: bool = False,
    auxiliary: Dict[str, Set[str]] = {},
    overwrite_number: Dict[str, Dict[str, str]] = {},
) -> Iterator[VariantRecord]:
    env = Environment(expression, ann_key, vcf.header, auxiliary, overwrite_number)

    info_keys = set(vcf.header.info.keys())
    has_svtype = "SVTYPE" in info_keys

    def is_bnd_record(record: VariantRecord) -> bool:
        return has_svtype and record.info.get("SVTYPE", None) == "BND"

    def get_event_name(record: VariantRecord) -> Tuple[str, str]:
        mate_ids = record.info.get("MATEID", [])
        event_name = record.info.get("EVENT", None)

        if len(mate_ids) > 1 and not event_name:
            raise ValueError(
                f"Filtering of BND records with multiple mates is unsupported "
                f"(see VCF 4.3, section 5.4.3 'Multiple mates'):\n{str(record)}"
            )

        mate_id = mate_ids[0] if len(mate_ids) == 1 else None
        mate_pair = mate_key([record.id, mate_id]) if mate_ids else None

        return event_name, mate_pair

    record: VariantRecord
    if not preserve_order:
        # If the order is not important, emit records that pass the filter expression
        # as we encounter them
        # However, breakends have to be considered jointly, so keep track of the
        # respective events.
        events: Dict[str, BreakendEvent] = dict()
        for idx, record in enumerate(vcf):
            record, keep = test_and_update_record(
                env, idx, record, ann_key, keep_unmatched
            )

            # Breakend records *may* have the "EVENT" specified, but don't have to.
            # In that case only the MATEID *may* bee available
            # (which may contain more than one ID)
            if is_bnd_record(record):
                event_name, mate_pair_name = get_event_name(record)

                # if EVENT is set, it has priority over MATEID.
                event_name = event_name or mate_pair_name

                # if both EVENT and MATEID are not set
                # we can only treat this BND record as a regular one,
                # (or somehow try to guess its correct mate)
                if not event_name:
                    print(
                        f"Warning: Encountered breakend record at index {idx} "
                        f"without either of MATEID or EVENT specified, "
                        f"treating it as a regular record:\n{str(record)}",
                        file=stderr,
                    )
                    if keep:
                        yield record
                    continue

                event = events.get(event_name, None)

                # if there's already an associated event
                if event:
                    # add this record to the event
                    event.add(record, keep)

                    # if we already know that the event is a "PASS"…
                    if event.keep:
                        # … emit all records associated with it
                        yield from event.emit()

                        # in the case of a simple mate pair, we can delete the event
                        # at this point, because no more records will be added to it
                        if event.is_mate_pair():
                            del events[mate_pair_name]
                else:
                    # if there's no entry for the event or mate pair yet, create one
                    is_mate_pair = mate_pair_name and mate_pair_name == event_name
                    event = BreakendEvent(event_name, is_mate_pair)
                    event.add(record, keep)
                    events[event_name] = event
            elif keep:
                yield record

        if len(events) > 0:
            # output BNDs if any are left unprocessed
            for event_name, event in events.items():
                if event_name and event.keep:
                    yield from event.emit()
    elif preserve_order:
        # If order *is* important, the first pass cannot emit any records but only
        # keep track of breakend events. The records will only be emitted during the
        # second pass.
        events: Set[str] = set()

        def fallback_name(record: VariantRecord) -> str:
            event_name, mate_pair_name = get_event_name(record)

            # There may be BND record which have neither of
            # INFO['EVENT'], INFO['MATEID'] or even ID specified.
            # In that case, we still need a pseudo-event name,
            # which we construct from the record's index in the file.
            return event_name or mate_pair_name or record.id or f"DUMMY: {idx}"

        for idx, record in enumerate(vcf):
            if is_bnd_record(record):
                record, keep = test_and_update_record(
                    env, idx, record, ann_key, keep_unmatched
                )
                if keep:
                    event_name = fallback_name(record)
                    events.add(event_name)

        # The second pass can now yield records in the correct order
        vcf.reset()
        for idx, record in enumerate(vcf):
            if is_bnd_record(record):
                event_name = fallback_name(record)
                if event_name in events:
                    yield record
            else:
                record, keep = test_and_update_record(
                    env, idx, record, ann_key, keep_unmatched
                )
                if keep:
                    yield record


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


def read_auxiliary(aux) -> Dict[str, Set[str]]:
    # read auxiliary files, split at any whitespace and store contents in a set
    def read_set(path: str) -> Set[str]:
        with open(path, "rt") as f:
            return set(line.rstrip() for line in f)

    return {name: read_set(contents) for name, contents in aux}


def execute(args):
    aux = read_auxiliary(args.aux)
    with VariantFile(args.vcf) as vcf:
        header: VariantHeader = vcf.header
        header.add_meta("vembraneVersion", __version__)
        # NOTE: If .modules.filter.execute might be used as a library function
        #       in the future, we should not record sys.argv directly below.
        header.add_meta(
            "vembraneCmd",
            "vembrane "
            + " ".join(
                "'" + arg.replace("'", '"') + '"' if " " in arg else arg
                for arg in sys.argv[1:]
            ),
        )

        overwrite_number = {
            "INFO": dict(args.overwrite_number_info),
            "FORMAT": dict(args.overwrite_number_format),
        }

        records = filter_vcf(
            vcf,
            args.expression,
            args.annotation_key,
            keep_unmatched=args.keep_unmatched,
            preserve_order=args.preserve_order,
            auxiliary=aux,
            overwrite_number=overwrite_number,
        )

        try:
            first_record = list(islice(records, 1))
        except VembraneError as ve:
            print(ve, file=stderr)
            exit(1)

        records = chain(first_record, records)

        if args.statistics is not None:
            records = statistics(records, vcf, args.statistics, args.annotation_key)

        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]
        with VariantFile(
            args.output,
            f"w{fmt}",
            header=header,
        ) as out:
            try:
                for record in records:
                    out.write(record)

            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
