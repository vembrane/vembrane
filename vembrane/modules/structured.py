from io import TextIOWrapper
import json
import textwrap
from typing import Any, Iterator
import yaml

from vembrane.backend.base import VCFReader
from vembrane.common import add_common_arguments, create_reader, smart_open
from vembrane.errors import VembraneError

from vembrane.representations import ModifiableEnvironment
import yte


def add_subcommand(subparsers):
    parser = subparsers.add_parser("structured")
    parser.add_argument(
        "template",
        help="File containing a YTE template with the desired structure per record and "
        "expressions that retrieve data from the VCF record.",
    )
    parser.add_argument(
        "vcf",
        help="The file containing the variants.",
        nargs="?",
        default="-",
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--output-fmt",
        choices=["json", "jsonl", "yaml"],
        help="Output format. If not specified, can be automatically determined from "
        "the --output file extension.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )
    add_common_arguments(parser)


ConvertedRecords = Iterator[list[Any] | dict[Any] | Any]


def process_vcf(
    vcf: VCFReader,
    template: str,
    ann_key: str,

) -> ConvertedRecords:
    env = ModifiableEnvironment(ann_key, vcf.header)
    for idx, record in enumerate(vcf):
        env.update_from_record(idx, record)
        # TODO: analogous to EvalEnvironment, we should implement a way to
        # work around the float32 issue, see the AST parsing in
        # EvalEnvironment.__init__
        # This needs however an extension of YTE internals.
        converted = yte.process_yaml(template, variables=env)
        yield converted


def write_records_jsonl(output_file: TextIOWrapper, records: ConvertedRecords):
    for record in records:
        print(json.dumps(record), file=output_file)


def write_records_json(output_file: TextIOWrapper, records: ConvertedRecords):
    print("[", file=output_file)
    for idx, record in enumerate(records):
        if idx > 0:
            print(",", file=output_file)
        print(json.dumps(record), file=output_file)
    print("]", file=output_file)


def write_records_yaml(output_file: TextIOWrapper, records: ConvertedRecords):
    for record in records:
        head, tail = yaml.dump(record).split("\n", 1)
        head = f"- {head}"
        tail = textwrap.indent(tail, "  ")
        
        print(head, file=output_file)
        print(tail, file=output_file)


def execute(args):
    overwrite_number = {
        "INFO": dict(args.overwrite_number_info),
        "FORMAT": dict(args.overwrite_number_format),
    }

    with open(args.template, "r") as template:
        template = template.read()

    if args.output_fmt is None:
        if args.output.endswith(".json"):
            args.output_fmt = "json"
        elif args.output.endswith(".jsonl"):
            args.output_fmt = "jsonl"
        elif args.output.endswith(".yaml") or args.output.endswith(".yml"):
            args.output_fmt = "yaml"
        elif args.output == "-":
            raise VembraneError(
                "--output-fmt must be specified when writing to STDOUT."
            )
        else:
            raise VembraneError(
                f"Unsupported file format: {args.output}, only .json, .jsonl and .yaml are supported."
            )
        
    if args.output_fmt == "jsonl":
        write_records = write_records_jsonl
    elif args.output_fmt == "json":
        write_records = write_records_json
    elif args.output_fmt == "yaml":
        write_records = write_records_yaml

    with create_reader(
        args.vcf,
        backend=args.backend,
        overwrite_number=overwrite_number,
    ) as reader, smart_open(args.output, mode="w") as writer:
        write_records(writer, process_vcf(reader, template, args.annotation_key))