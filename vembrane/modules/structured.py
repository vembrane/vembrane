import json
import textwrap
from io import TextIOWrapper
from typing import Any, Dict, Iterator

import yaml
import yte

from vembrane.ann_types import NA
from vembrane.backend.base import VCFHeader, VCFReader
from vembrane.common import add_common_arguments, create_reader, smart_open
from vembrane.errors import VembraneError
from vembrane.representations import (
    Annotation,
    SourceEnvironment,
)


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


class CodeHandler(yte.CodeHandler):
    def __init__(self, ann_key: str, header: VCFHeader) -> None:
        self.ann_key = ann_key
        self.header = header
        self._envs: Dict[str, SourceEnvironment] = {}
        self._record = None
        self._record_idx = None
        self._annotation = None

    def update_from_record(self, idx: int, record: VCFReader) -> None:
        self._record = record
        self._record_idx = idx
        for env in self._envs.values():
            env.update_from_record(idx, record)

    def update_from_annotation(self, annotation: str | None) -> None:
        self._annotation = annotation
        if annotation is not None:
            for env in self._envs.values():
                if env.expression_annotations():
                    env.update_annotation(annotation)

    def _env(self, source: str) -> SourceEnvironment:
        env = self._envs.get(source)
        if env is None:
            env = SourceEnvironment(source, self.ann_key, self.header)
            env.update_from_record(self._record_idx, self._record)
            if self._annotation is not None:
                env.update_annotation(self._annotation)
            self._envs[source] = env
        return env

    def eval(self, expr: str, variables: Dict[str, Any]) -> Any:
        env = self._env(expr)
        return eval(env.compiled, variables, env)

    def exec(self, source: str, variables: Dict[str, Any]) -> Any:
        env = self._env(source)
        return exec(env.compiled, variables, env)


class ValueHandler(yte.ValueHandler):
    def postprocess_atomic_value(self, value: Any) -> Any:
        processed = super().postprocess_atomic_value(value)
        if processed is NA:
            return None
        return processed


def process_vcf(
    vcf: VCFReader,
    template: str,
    ann_key: str,
    variables: Dict[str, Any] | None = None,
) -> ConvertedRecords:
    code_handler = CodeHandler(ann_key, vcf.header)
    value_handler = ValueHandler()
    if variables is None:
        variables = {}

    annotation = Annotation(ann_key, vcf.header)

    for idx, record in enumerate(vcf):
        code_handler.update_from_record(idx, record)
        annotations = annotation.get_record_annotations(idx, record)

        for ann in annotations:
            code_handler.update_from_annotation(ann)
            converted = yte.process_yaml(
                template,
                code_handler=code_handler,
                value_handler=value_handler,
                variables=variables,
            )
            converted["component"] = [
                entry for entry in converted["component"] if is_valid(entry)
            ]
            yield converted


def is_valid(record):
    if type(record) == dict:
        return all([is_valid(value) for value in record.values()])
    elif type(record) == list:
        return all([is_valid(entry) for entry in record])
    else:
        return True if record else False


def write_records_jsonl(output_file: TextIOWrapper, records: ConvertedRecords) -> None:
    for record in records:
        print(json.dumps(record), file=output_file)


def write_records_json(output_file: TextIOWrapper, records: ConvertedRecords) -> None:
    print("[", file=output_file)
    for idx, record in enumerate(records):
        if idx > 0:
            print(",", file=output_file)
        print(json.dumps(record), file=output_file)
    print("]", file=output_file)


def write_records_yaml(output_file: TextIOWrapper, records: ConvertedRecords) -> None:
    for record in records:
        head, tail = yaml.dump(record).split("\n", 1)
        head = f"- {head}"
        tail = textwrap.indent(tail, "  ")

        print(head, file=output_file)
        print(tail, file=output_file)


def process(
    template: str,
    output_fmt: str | None,
    output: str | None,
    vcf: str,
    annotation_key: str,
    overwrite_number_info: dict[str, str],
    overwrite_number_format,
    backend,
    variables: Dict[str, Any] | None = None,
) -> None:
    overwrite_number = {
        "INFO": dict(overwrite_number_info),
        "FORMAT": dict(overwrite_number_format),
    }

    if output_fmt is None:
        if output.endswith(".json"):
            output_fmt = "json"
        elif output.endswith(".jsonl"):
            output_fmt = "jsonl"
        elif output.endswith(".yaml") or output.endswith(".yml"):
            output_fmt = "yaml"
        elif output == "-":
            raise VembraneError(
                "--output-fmt must be specified when writing to STDOUT."
            )
        else:
            raise VembraneError(
                f"Unsupported file format: {output}, only .json, .jsonl and "
                ".yaml are supported."
            )

    if output_fmt == "jsonl":
        write_records = write_records_jsonl
    elif output_fmt == "json":
        write_records = write_records_json
    elif output_fmt == "yaml":
        write_records = write_records_yaml

    with (
        create_reader(
            vcf,
            backend=backend,
            overwrite_number=overwrite_number,
        ) as reader,
        smart_open(output, mode="w") as writer,
    ):
        write_records(
            writer, process_vcf(reader, template, annotation_key, variables=variables)
        )


def execute(args):
    with open(args.template, "r") as template:
        template = template.read()

    process(
        template=template,
        output_fmt=args.output_fmt,
        output=args.output,
        vcf=args.vcf,
        annotation_key=args.annotation_key,
        overwrite_number_info=args.overwrite_number_info,
        overwrite_number_format=args.overwrite_number_format,
        backend=args.backend,
    )
