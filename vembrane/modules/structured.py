import json
import sys
import textwrap
from typing import Any, Callable, Dict, Iterator, TextIO

import yaml
import yte  # type: ignore

from vembrane.ann_types import NA
from vembrane.backend.base import VCFHeader, VCFReader, VCFRecord
from vembrane.common import (
    Context,
    HumanReadableDefaultsFormatter,
    Primitive,
    add_common_arguments,
    create_reader,
    smart_open,
)
from vembrane.errors import VembraneError, handle_vembrane_error
from vembrane.representations import (
    Annotation,
    SourceEnvironment,
)


def add_subcommand(subparsers):
    parser = subparsers.add_parser(
        "structured",
        help="Create structured output (e.g., JSON/YAML) "
        "from a VCF/BCF file using a YTE template.",
        description="Create structured output from a VCF/BCF and a YTE template.",
        formatter_class=HumanReadableDefaultsFormatter,
    )
    parser.add_argument(
        "template",
        help="File containing a YTE template with the desired structure per record "
        "and expressions that retrieve data from the VCF/BCF record.",
    )
    parser.add_argument(
        "vcf",
        help="Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin.",
        nargs="?",
        default="-",
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


ConvertedRecords = Iterator[list[Any] | dict[Any, Any] | Any]


class CodeHandler(yte.CodeHandler):
    def __init__(
        self,
        ann_key: str,
        header: VCFHeader,
        allowed_globals: dict[str, Any] | None,
        auxiliary_globals: dict[str, Any] | None,
    ) -> None:
        self.ann_key = ann_key
        self.header = header
        self._envs: Dict[str, SourceEnvironment] = {}
        self._record: VCFRecord | None = None
        self._record_idx: int | None = None
        self._annotation: str | None = None
        self._auxiliary_globals = auxiliary_globals
        self._allowed_globals = allowed_globals

    def update_from_record(self, idx: int, record: VCFRecord) -> None:
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
            env = SourceEnvironment(
                source,
                self.ann_key,
                self.header,
                allowed_globals=self._allowed_globals,
                auxiliary_globals=self._auxiliary_globals,
            )
            env.update_from_record(self._record_idx, self._record)  # type: ignore
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

    def has_annotation(self) -> bool:
        return any(env.expression_annotations() for env in self._envs.values())


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
    allowed_globals: dict[str, Any] | None = None,
    auxiliary_globals: dict[str, Any] | None = None,
    postprocess: Callable[[Any], Any] | None = None,
) -> ConvertedRecords:
    code_handler = CodeHandler(ann_key, vcf.header, allowed_globals, auxiliary_globals)
    value_handler = ValueHandler()

    annotation = Annotation(ann_key, vcf.header)

    for idx, record in enumerate(vcf):
        try:
            # TODO: For now structural variants are not supported
            # but should be added later.
            if record.is_sv_record:
                print(
                    f"Warning: Record is a structural variant which are currently "
                    f"not supported and will be skipped.\n"
                    f"Record: {record}"
                )
                continue
            code_handler.update_from_record(idx, record)

            annotations: list[str] | None

            if annotation:
                annotations = annotation.get_record_annotations(idx, record)
            else:
                annotations = None

            # type ignore below because mypy is too stupid!
            for ann_entry in annotations or [None]:  # type: ignore
                if ann_entry is not None:
                    code_handler.update_from_annotation(ann_entry)
                converted = yte.process_yaml(
                    template,
                    code_handler=code_handler,
                    value_handler=value_handler,
                )
                if postprocess is not None:
                    converted = postprocess(converted)
                yield converted
                if not code_handler.has_annotation():
                    # If annotation is not used in template, only process first entry.
                    # There can be no subsequent entries leading to different results.
                    break

        except Exception as e:
            raise VembraneError.from_record_and_exception(idx, record, e) from e


def write_records_jsonl(output_file: TextIO, records: ConvertedRecords) -> None:
    for record in records:
        print(json.dumps(record), file=output_file)


def write_records_json(output_file: TextIO, records: ConvertedRecords) -> None:
    print("[", file=output_file)
    for idx, record in enumerate(records):
        if idx > 0:
            print(",", file=output_file)
        print(json.dumps(record), file=output_file)
    print("]", file=output_file)


def write_records_yaml(output_file: TextIO, records: ConvertedRecords) -> None:
    for record in records:
        head, tail = yaml.dump(record).split("\n", 1)
        head = f"- {head}"
        tail = textwrap.indent(tail, "  ")

        print(head, file=output_file)
        print(tail, file=output_file)


@handle_vembrane_error
def process(
    template: str,
    output_fmt: str | None,
    output: str,
    vcf: str,
    annotation_key: str,
    overwrite_number_info: dict[str, str],
    overwrite_number_format,
    backend,
    allowed_globals: Dict[str, Any] | None = None,
    auxiliary_globals: Dict[str, Any] | None = None,
    postprocess: (
        Callable[[Primitive | dict | list], Primitive | dict | list] | None
    ) = None,
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
    else:
        raise VembraneError(
            f"Unsupported output format: {output_fmt}, only json, jsonl and yaml are "
            "supported."
        )

    with (
        create_reader(
            vcf,
            backend=backend,
            overwrite_number=overwrite_number,
        ) as reader,
        smart_open(output, mode="w") as writer,
    ):
        write_records(
            writer,
            process_vcf(
                reader,
                template,
                annotation_key,
                allowed_globals,
                auxiliary_globals,
                postprocess=postprocess,
            ),
        )


def execute(args):
    try:
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
            auxiliary_globals=Context.from_args(args).get_globals(),
        )
    except VembraneError as ve:
        print(ve, file=sys.stderr)
        sys.exit(1)
