import argparse
import sys
from collections import defaultdict
from collections.abc import Iterator
from itertools import chain, islice
from sys import stderr
from typing import Any

import yaml

from .. import __version__
from ..backend.base import VCFReader, VCFRecord
from ..common import (
    BreakendEvent,
    Context,
    HumanReadableDefaultsFormatter,
    add_common_arguments,
    check_expression,
    create_reader,
    create_writer,
    get_annotation_keys,
    mate_key,
    normalize,
    read_auxiliary,
    read_ontology,
    split_annotation_entry,
)
from ..errors import VembraneError, handle_vembrane_error
from ..representations import FuncWrappedExpressionEnvironment
from ..sequence_ontology import SequenceOntology


class DeprecatedAction(argparse.Action):
    def __call__(self, *args, **kwargs):
        print(
            f"{'/'.join(self.option_strings)} is deprecated.\n{self.help}",
            file=sys.stderr,
        )
        sys.exit(1)


def add_subcommmand(subparsers):
    parser = subparsers.add_parser(
        "filter",
        help="Filter VCF/BCF records and annotations using a python expression.",
        description="Filter VCF/BCF records and annotations "
        "based on a user-defined Python expression."
        "Only records for which the expression evaluates to True are kept.",
        formatter_class=HumanReadableDefaultsFormatter,
    )
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "expression",
        type=check_expression,
        help="A python expression to filter variants. "
        "The expression must evaluate to bool. "
        "All VCF/BCF fields can be accessed. "
        "Additionally, annotation fields can be accessed, see `--annotation-key`. "
        "If all annotations of a record are filtered out, "
        "the entire record is removed.",
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
    add_common_arguments(parser)


def test_and_update_record(
    env: FuncWrappedExpressionEnvironment,
    idx: int,
    record: VCFRecord,
    ann_key: str,
    keep_unmatched: bool,
) -> tuple[VCFRecord, bool]:
    try:
        return _test_and_update_record(env, idx, record, ann_key, keep_unmatched)
    except VembraneError as e:
        raise e
    except Exception as e:
        raise VembraneError.from_record_and_exception(idx, record, e) from e


def _test_and_update_record(
    env: FuncWrappedExpressionEnvironment,
    idx: int,
    record: VCFRecord,
    ann_key: str,
    keep_unmatched: bool,
) -> tuple[VCFRecord, bool]:
    env.update_from_record(idx, record)
    if env.expression_annotations():
        annotations = env.get_record_annotations(idx, record)

        #  … and only keep the annotations where the expression evaluates to true
        if keep_unmatched:
            filtered = any(map(env.is_true, annotations))
            return record, filtered
        else:
            filtered_annotations = [
                annotation for annotation in annotations if env.is_true(annotation)
            ]

            if len(annotations) != len(filtered_annotations):
                # update annotations if they have actually been filtered
                record.info[ann_key] = filtered_annotations

            return record, len(filtered_annotations) > 0
    else:
        # otherwise, the annotations are irrelevant w.r.t. the expression,
        # so we can omit them
        return record, env.is_true()


def filter_vcf(
    reader: VCFReader,
    expression: str,
    ann_key: str,
    keep_unmatched: bool = False,
    preserve_order: bool = False,
    auxiliary: dict[str, set[str]] | None = None,
    auxiliary_globals: dict[str, Any] | None = None,
    ontology: SequenceOntology | None = None,
) -> Iterator[VCFRecord]:
    if auxiliary is None:
        auxiliary = {}
    env = FuncWrappedExpressionEnvironment(
        expression,
        ann_key,
        reader.header,
        auxiliary,
        ontology,
        auxiliary_globals=auxiliary_globals,
    )
    has_mateid_key = reader.header.infos.get("MATEID", None) is not None
    has_event_key = reader.header.infos.get("EVENT", None) is not None

    def get_event_name(
        record: VCFRecord,
        has_mateid_key: bool = has_mateid_key,
        has_event_key: bool = has_event_key,
    ) -> tuple[str | None, str | None]:
        mate_ids: list[str] | str = (
            record.info.get("MATEID", []) if has_mateid_key else []
        )
        if isinstance(mate_ids, str):
            # some callers annotate MATEID as single string
            mate_ids = [mate_ids]
        event_name: str | None = (
            record.info.get("EVENT", None) if has_event_key else None
        )

        if len(mate_ids) > 1 and not event_name:
            raise ValueError(
                f"Filtering of BND records with multiple mates is unsupported "
                f"(see VCF 4.3, section 5.4.3 'Multiple mates'):\n{str(record)}",
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
        event_dict: dict[str, BreakendEvent] = {}
        for idx, record in enumerate(reader):
            record, keep = test_and_update_record(
                env,
                idx,
                record,
                ann_key,
                keep_unmatched,
            )

            # Breakend records *may* have the "EVENT" specified, but don't have to.
            # In that case only the MATEID *may* be available
            # (which may contain more than one ID)
            if record.is_bnd_record:
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
                        if event.is_mate_pair() and mate_pair_name:
                            del event_dict[mate_pair_name]
                else:
                    # if there's no entry for the event or mate pair yet, create one
                    is_mate_pair = bool(mate_pair_name) and mate_pair_name == event_name
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
        event_set: set[str] = set()

        def fallback_name(record: VCFRecord) -> str:
            event_name, mate_pair_name = get_event_name(record)

            # There may be BND record which have neither of
            # INFO['EVENT'], INFO['MATEID'] or even ID specified.
            # In that case, we still need a pseudo-event name,
            # which we construct from the record's index in the file.
            return event_name or mate_pair_name or record.id or f"DUMMY: {idx}"

        for idx, record in enumerate(reader):
            if record.is_bnd_record:
                record, keep = test_and_update_record(
                    env,
                    idx,
                    record,
                    ann_key,
                    keep_unmatched,
                )
                if keep:
                    event_name = fallback_name(record)
                    event_set.add(event_name)

        # The second pass can now yield records in the correct order
        reader.reset()
        for idx, record in enumerate(reader):
            if record.is_bnd_record:
                event_name = fallback_name(record)
                if event_name in event_set:
                    yield record
            else:
                record, keep = test_and_update_record(
                    env,
                    idx,
                    record,
                    ann_key,
                    keep_unmatched,
                )
                if keep:
                    yield record


def statistics(
    records: Iterator[VCFRecord],
    vcf: VCFReader,
    filename: str,
    ann_key: str,
) -> Iterator[VCFRecord]:
    annotation_keys = get_annotation_keys(vcf.header, ann_key)
    counter: defaultdict[str, defaultdict[str, Any]] = defaultdict(
        lambda: defaultdict(lambda: 0)
    )
    for record in records:
        for annotation in record.info[ann_key]:
            for key, raw_value in zip(
                annotation_keys,
                split_annotation_entry(annotation),
                strict=True,
            ):
                value = raw_value.strip()
                if value:
                    counter[key][value] += 1
        yield record

    # reduce dicts with many items, to just one counter
    for key, subdict in counter.items():
        if len(subdict) > 10:
            counter[key] = f"#{len(subdict)}"  # type: ignore

    yaml.add_representer(defaultdict, yaml.representer.Representer.represent_dict)
    with open(filename, "w") as out:
        yaml.dump(dict(counter), out)


@handle_vembrane_error
def execute(args) -> None:
    aux = read_auxiliary(args.aux)
    ontology = read_ontology(args.ontology)

    overwrite_number = {
        "INFO": dict(args.overwrite_number_info),
        "FORMAT": dict(args.overwrite_number_format),
    }
    with create_reader(
        args.vcf,
        backend=args.backend,
        overwrite_number=overwrite_number,
    ) as reader:
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

        records = filter_vcf(
            reader,
            args.expression,
            args.annotation_key,
            keep_unmatched=args.keep_unmatched,
            preserve_order=args.preserve_order,
            auxiliary=aux,
            auxiliary_globals=Context.from_args(args).get_globals(),
            ontology=ontology,
        )

        first_record = list(islice(records, 1))

        records = chain(first_record, records)

        if args.statistics is not None:
            records = statistics(records, reader, args.statistics, args.annotation_key)

        fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]

        with create_writer(args.output, fmt, reader, backend=args.backend) as writer:
            for record in records:
                writer.write(record)
