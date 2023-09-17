import argparse
import sys
from collections import defaultdict
from itertools import chain, islice
from sys import stderr
from typing import Dict, Iterator, List, Optional, Set, Tuple

import yaml

from .. import __version__
from ..backend.base import Backend, VCFReader, VCFRecord
from ..common import (
    AppendKeyValuePair,
    BreakendEvent,
    check_expression,
    create_reader,
    create_writer,
    get_annotation_keys,
    is_bnd_record,
    mate_key,
    normalize,
    read_auxiliary,
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
        nargs=1,
        action=AppendKeyValuePair,
        default={},
        metavar="NAME=PATH",
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
        "--overwrite-number-info",
        nargs=1,
        action=AppendKeyValuePair,
        metavar="FIELD=NUMBER",
        default={},
        help="Overwrite the number specification for INFO fields "
        "given in the VCF header. "
        "Example: `--overwrite-number cosmic_CNT=.`",
    )
    parser.add_argument(
        "--overwrite-number-format",
        nargs=1,
        action=AppendKeyValuePair,
        metavar="FIELD=NUMBER",
        default={},
        help="Overwrite the number specification for FORMAT fields "
        "given in the VCF header. "
        "Example: `--overwrite-number-format DP=2`",
    )
    parser.add_argument(
        "--backend",
        "-b",
        default="pysam",
        type=Backend.from_string,
        choices=[Backend.pysam, Backend.cyvcf2],
        help="Set the backend library.",
    )


def test_and_update_record(
    env: Environment,
    idx: int,
    record: VCFRecord,
    ann_key: str,
    keep_unmatched: bool,
) -> Tuple[VCFRecord, bool]:
    try:
        return _test_and_update_record(env, idx, record, ann_key, keep_unmatched)
    except VembraneError as e:
        raise e
    except Exception as e:
        print(f"Encountered an error while processing record {idx}", file=stderr)
        print(str(record), file=stderr)
        raise e


def _test_and_update_record(
    env: Environment,
    idx: int,
    record: VCFRecord,
    ann_key: str,
    keep_unmatched: bool,
) -> Tuple[VCFRecord, bool]:
    env.update_from_record(idx, record)
    if env.expression_annotations():
        # if the expression contains a reference to the ANN field
        # get all annotations from the record.info field
        # (or supply an empty ANN value if the record has no ANN field)
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
    reader: VCFReader,
    expression: str,
    ann_key: str,
    keep_unmatched: bool = False,
    preserve_order: bool = False,
    auxiliary: Dict[str, Set[str]] = {},
    overwrite_number: Dict[str, Dict[str, str]] = {},
) -> Iterator[VCFRecord]:
    env = Environment(expression, ann_key, reader.header, auxiliary, overwrite_number)
    has_mateid_key = reader.header.infos.get("MATEID", None) is not None
    has_event_key = reader.header.infos.get("EVENT", None) is not None

    def get_event_name(
        record: VCFRecord,
        has_mateid_key=has_mateid_key,
        has_event_key=has_event_key,
    ) -> Tuple[Optional[str], Optional[str]]:
        print("fooo1", file=sys.stderr)
        mate_ids: List[str] = record.info.get("MATEID", []) if has_mateid_key else []
        event_name: Optional[str] = (
            record.info.get("EVENT", None) if has_event_key else None
        )

        if len(mate_ids) > 1 and not event_name:
            raise ValueError(
                f"Filtering of BND records with multiple mates is unsupported "
                f"(see VCF 4.3, section 5.4.3 'Multiple mates'):\n{str(record)}"
            )

        mate_id = mate_ids[0] if len(mate_ids) == 1 else None
        mate_pair = mate_key([record.id, mate_id]) if mate_ids else None

        return event_name, mate_pair

    record: VCFRecord
    if not preserve_order:
        # If the order is not important, emit records that pass the filter expression
        # as we encounter them
        # However, breakends have to be considered jointly, so keep track of the
        # respective events.
        event_dict: Dict[str, BreakendEvent] = dict()
        for idx, record in enumerate(reader):
            record, keep = test_and_update_record(
                env, idx, record, ann_key, keep_unmatched
            )

            # Breakend records *may* have the "EVENT" specified, but don't have to.
            # In that case only the MATEID *may* be available
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

                event = event_dict.get(event_name, None)

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
                            del event_dict[mate_pair_name]
                else:
                    # if there's no entry for the event or mate pair yet, create one
                    is_mate_pair = mate_pair_name and mate_pair_name == event_name
                    event = BreakendEvent(event_name, is_mate_pair)
                    event.add(record, keep)
                    event_dict[event_name] = event
            elif keep:
                yield record

        if len(event_dict) > 0:
            # output BNDs if any are left unprocessed
            for event_name, event in event_dict.items():
                if event_name and event.keep:
                    yield from event.emit()
    elif preserve_order:
        # If order *is* important, the first pass cannot emit any records but only
        # keep track of breakend events. The records will only be emitted during the
        # second pass.
        event_set: Set[str] = set()

        def fallback_name(record: VCFRecord) -> str:
            event_name, mate_pair_name = get_event_name(record)

            # There may be BND record which have neither of
            # INFO['EVENT'], INFO['MATEID'] or even ID specified.
            # In that case, we still need a pseudo-event name,
            # which we construct from the record's index in the file.
            return event_name or mate_pair_name or record.id or f"DUMMY: {idx}"

        for idx, record in enumerate(reader):
            if is_bnd_record(record):
                record, keep = test_and_update_record(
                    env, idx, record, ann_key, keep_unmatched
                )
                if keep:
                    event_name = fallback_name(record)
                    event_set.add(event_name)

        # The second pass can now yield records in the correct order
        reader.reset()
        for idx, record in enumerate(reader):
            if is_bnd_record(record):
                event_name = fallback_name(record)
                if event_name in event_set:
                    yield record
            else:
                record, keep = test_and_update_record(
                    env, idx, record, ann_key, keep_unmatched
                )
                if keep:
                    yield record


def statistics(
    records: Iterator[VCFRecord], vcf: VCFReader, filename: str, ann_key: str
) -> Iterator[VCFRecord]:
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


def execute(args) -> None:
    aux = read_auxiliary(args.aux)
    with create_reader(args.vcf, backend=args.backend) as reader:
        # header: dict = vcf.header
        reader.header.add_generic("vembraneVersion", __version__)
        # NOTE: If .modules.filter.execute might be used as a library function
        #       in the future, we should not record sys.argv directly below.
        cmd_parts = [normalize(arg) if " " in arg else arg for arg in sys.argv[3:]]
        expr = " ".join(a.strip() for a in args.expression.split("\n"))
        expr = normalize(expr)

        reader.header.add_generic(
            "vembraneCmd",
            "vembrane " + expr + " " + " ".join(cmd_parts),
        )

        overwrite_number = {
            "INFO": dict(args.overwrite_number_info),
            "FORMAT": dict(args.overwrite_number_format),
        }

        records = filter_vcf(
            reader,
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
            records = statistics(records, reader, args.statistics, args.annotation_key)

        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]

        with create_writer(args.output, fmt, reader, backend=args.backend) as writer:
            try:
                for record in records:
                    writer.write(record)

            except VembraneError as ve:
                print(ve, file=stderr)
                exit(1)
