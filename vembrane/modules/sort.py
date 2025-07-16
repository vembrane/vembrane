import heapq
import math
import sys
import tempfile
from functools import total_ordering
from itertools import batched
from pathlib import Path
from typing import Any, Callable, Iterable, List, Self, Sequence, Tuple, TypeVar

from vembrane import __version__
from vembrane.ann_types import NA
from vembrane.backend.base import Backend, VCFReader, VCFRecord, VCFWriter
from vembrane.common import (
    add_common_arguments,
    check_expression,
    create_reader,
    create_writer,
    normalize,
)
from vembrane.errors import VembraneError
from vembrane.representations import Annotation, SourceEnvironment


def add_subcommand(subparsers):
    parser = subparsers.add_parser(
        "sort",
        description="Sort VCF records by one or multiple Python expressions that "
        "encode keys for the desired order. This feature loads the entire VCF file "
        "into memory in order to maximize performance. It is thus meant to sort small, "
        "already filtered VCF files, e.g. for prioritizing records for the human eye. "
        "For large VCF files, the only relevant sorting is usually by position, "
        "which is better done with e.g. bcftools (and usually the sorting that variant "
        "callers output).",
    )
    parser.add_argument(
        "vcf",
        help="The VCF/BCF file containing the variants. If not specified, "
        "reads from STDIN.",
        nargs="?",
        default="-",
    )
    parser.add_argument(
        "expression",
        type=check_expression,
        help="Python expression (or tuple of expressions) returning orderable values "
        "to sort the VCF records by (ascending, smallest values coming first). "
        "If multiple expressions are provided as a tuple, they are prioritized from "
        "left to right with lowest priority on the right. "
        "NA/NaN values are sorted to the end.",
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
        "--chunk-size",
        type=int,
        default=1000,
        help="Number of VCF records to sort in memory. If the VCF file exceeds this "
        "number of records, external sorting is used.",
    )
    add_common_arguments(parser)


def is_na(value: Any) -> bool:
    return (
        value is NA or value is None or (isinstance(value, float) and math.isnan(value))
    )


T = TypeVar("T")

type NestedSequence[T] = Sequence[T | NestedSequence[T]]


# Key wrapper classes for sorting.
# sorted() uses <, max and min use < and >, so we need to implement both.
# For maximizing speed, we avoid delegation, saving an additional method call.
class KeyBase:
    @classmethod
    def wrap(cls, key: Any) -> Self | NestedSequence[Self]:
        if isinstance(key, Iterable) and not isinstance(key, str):
            # apply recursively to tuples, lists, etc.
            return tuple(cls.wrap(k) for k in key)
        else:
            return cls(key)

    def __init__(self, key: Any):
        self.key = key

    # rationale for the type ignore below: the code ensures that other is of same type
    # hence we avoid additional checks here for performance reasons.
    def __eq__(self, other: Self) -> bool:  # type: ignore
        if is_na(self.key) and is_na(other.key):
            return True
        return self.key == other.key


@total_ordering
class KeyAscending(KeyBase):
    def __lt__(self, other: Self) -> bool:
        if is_na(self.key):
            return False
        if is_na(other.key):
            return True
        return self.key < other.key

    def __gt__(self, other: Self) -> bool:
        if is_na(self.key):
            return True
        if is_na(other.key):
            return False
        return self.key > other.key


@total_ordering
class KeyDescending(KeyBase):
    def __lt__(self, other: Self) -> bool:
        if is_na(self.key):
            return True
        if is_na(other.key):
            return False
        return self.key > other.key

    def __gt__(self, other: Self) -> bool:
        if is_na(self.key):
            return False
        if is_na(other.key):
            return True
        return self.key < other.key


type SortKey = (
    KeyDescending | KeyAscending | NestedSequence[KeyDescending | KeyAscending]
)


def apply_default_key(
    key: Any,
) -> SortKey:
    if isinstance(key, (KeyDescending, KeyAscending)):
        return key
    else:
        return KeyAscending.wrap(key)


auxiliary_globals: dict[str, Any] = {
    "asc": KeyAscending.wrap,
    "desc": KeyDescending.wrap,
}


def execute(args) -> None:
    overwrite_number = {
        "INFO": dict(args.overwrite_number_info),
        "FORMAT": dict(args.overwrite_number_format),
    }
    fmt = {"vcf": "", "bcf": "b", "uncompressed-bcf": "u"}[args.output_fmt]

    try:
        with create_reader(
            args.vcf,
            backend=args.backend,
            overwrite_number=overwrite_number,
        ) as reader:
            # TODO: unify the reader/writer setup across subcommands
            reader.header.add_generic("vembraneVersion", __version__)
            reader.header.add_generic(
                "vembraneCmd",
                "vembrane "
                + " ".join(
                    normalize(arg) if " " in arg else arg for arg in sys.argv[1:]
                ),
            )

            annotation = Annotation(args.annotation_key, reader.header)

            sort_keys = SourceEnvironment(
                args.expression,
                args.annotation_key,
                reader.header,
                auxiliary_globals=auxiliary_globals,
            )

            def get_sort_key(item: Tuple[int, VCFRecord]) -> SortKey:
                idx, record = item
                sort_keys.update_from_record(idx, record)

                if sort_keys.expression_annotations():
                    # if the expression refers to annotations, we obtain the
                    # minimal keys across all annotations
                    annotations = annotation.get_record_annotations(idx, record)

                    def eval_with_ann(annotation):
                        sort_keys.update_annotation(annotation)
                        return apply_default_key(
                            eval(
                                sort_keys.compiled,
                                sort_keys,
                            )
                        )

                    return min(map(eval_with_ann, annotations))
                else:
                    return apply_default_key(eval(sort_keys.compiled, sort_keys))

            with create_writer(
                args.output, fmt, reader, backend=args.backend
            ) as writer:
                external_sort(
                    reader,
                    writer,
                    key=get_sort_key,
                    chunk_size=args.chunk_size,
                    backend=args.backend,
                    overwrite_number=overwrite_number,
                )
    except VembraneError as ve:
        print(ve, file=sys.stderr)
        sys.exit(1)


def external_sort(
    reader: VCFReader,
    writer: VCFWriter,
    key: Callable[[Tuple[int, VCFRecord]], SortKey],
    chunk_size: int,
    backend: Backend,
    overwrite_number,
):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        chunk_files: List[Path] = []

        last_chunk = None

        def store_chunk(chunk):
            path = temp_path / f"chunk_{len(chunk_files)}.bcf"
            with create_writer(path, "", reader, backend=backend) as writer:
                for _, record in chunk:
                    writer.write(record)
            chunk_files.append(path)

        for chunk in batched(enumerate(reader), n=chunk_size):
            if last_chunk is not None:
                store_chunk(last_chunk)
            chunk = sorted(chunk, key=key)

            last_chunk = chunk
        if last_chunk:
            if not chunk_files:
                # only one chunk, nothing stored yet
                # directly write the sorted chunk
                for _, record in last_chunk:
                    writer.write(record)
                return
            else:
                # store last chunk
                store_chunk(last_chunk)
        else:
            # empty input
            return

        # now merge the chunks
        readers = {}
        for i, chunk_file in enumerate(chunk_files):
            reader = create_reader(
                chunk_file, backend=backend, overwrite_number=overwrite_number
            )
            readers[i] = enumerate(reader), reader

        def next_record_to_heap(reader_idx: int) -> None:
            reader_iter, reader = readers[reader_idx]
            try:
                idx, record = next(reader_iter)
            except StopIteration:
                reader.close()
                del readers[reader_idx]
                return
            heapq.heappush(heap, (key((idx, record)), record, reader_idx))

        heap: List[Tuple[SortKey, VCFRecord, int]] = []
        for reader_idx in range(len(readers)):
            next_record_to_heap(reader_idx)

        while True:
            # pop the smallest item from the heap
            try:
                _, record, reader_idx = heapq.heappop(heap)
            except IndexError:
                # heap empty, everything has been processed
                assert not readers
                return

            writer.write(record)
            next_record_to_heap(reader_idx)
