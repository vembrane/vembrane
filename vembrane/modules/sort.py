import math
import sys
from typing import Any

from vembrane import __version__
from vembrane.ann_types import NA
from vembrane.common import (
    add_common_arguments,
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
        "sort_keys",
        nargs="+",
        help="Python expressions returning orderable values to sort the VCF records "
        "by (ascending, smallest values coming first). If multiple expressions are "
        "provided, they are prioritized from left to "
        "right with lowest priority on the right. NA/NaN values are sorted to the end.",
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
            reader.header.add_generic("vembraneVersion", __version__)
            reader.header.add_generic(
                "vembraneCmd",
                "vembrane "
                + " ".join(
                    normalize(arg) if " " in arg else arg for arg in sys.argv[1:]
                ),
            )

            annotation = Annotation(args.annotation_key, reader.header)

            # load records into memory
            # since sorting by something else than the position is only reasonable
            # for already small, filtered VCF files, this is not an issue in practice.
            records = list(enumerate(reader))

            sort_keys = [
                SourceEnvironment(expr, args.annotation_key, reader.header)
                for expr in args.sort_keys
            ]

        def get_sort_key(item):
            idx, record = item
            for env in sort_keys:
                env.update_from_record(idx, record)

            def eval_env(env: SourceEnvironment) -> Any:
                def handle_na(value) -> Any:
                    if (
                        value is NA
                        or value is None
                        or (isinstance(value, float) and math.isnan(value))
                    ):
                        return NAToEnd()
                    return value

                if env.expression_annotations():
                    annotations = annotation.get_record_annotations(idx, record)

                    def eval_with_ann(annotation):
                        env.update_annotation(annotation)
                        return handle_na(eval(env.compiled, env))

                    return min(map(eval_with_ann, annotations))
                else:
                    return handle_na(eval(env.compiled, env))

            return tuple(eval_env(env) for env in sort_keys)

        records.sort(key=get_sort_key)

        with create_writer(args.output, fmt, reader, backend=args.backend) as writer:
            for _, record in records:
                writer.write(record)
    except VembraneError as ve:
        print(ve, file=sys.stderr)
        sys.exit(1)


class NAToEnd:
    def __lt__(self, other: Any) -> bool:
        return False
