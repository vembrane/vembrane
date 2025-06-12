import csv
import importlib.resources
from collections import defaultdict
from typing import Any

# intervaltree is untyped, so we use type: ignore to suppress type checking errors
from intervaltree import IntervalTree  # type: ignore

from vembrane.common import Primitive, add_common_arguments
from vembrane.modules import structured

PROFILE_DIR = (
    importlib.resources.files("vembrane.modules") / "assets" / "fhir" / "profiles"
)


def add_subcommand(subparsers):
    parser = subparsers.add_parser("fhir")
    parser.add_argument(
        "vcf",
        help="The file containing the variants.",
        nargs="?",
        default="-",
    )
    parser.add_argument(
        "sample",
        help="The sample to use for generating FHIR output.",
    )
    parser.add_argument(
        "assembly",
        help="The reference assembly used for read mapping.",
        choices=["GRCh38", "GRCh37"],
    )
    parser.add_argument(
        "--url",
        "-u",
        help="Generic url used as identifier by FHIR e.g. http://<institute>/<department>/VCF",
        default=None,
    )
    parser.add_argument(
        "--status",
        "-s",
        help="Status of findings. E.g. final, preliminary, ...",
        default=None,
    )
    parser.add_argument(
        "--coordinates",
        "-c",
        help="Coordinate system of base positions.",
        default=None,
        choices=["0", "1"],
    )
    parser.add_argument(
        "--profile",
        required=True,
        choices=[path.with_suffix("").name for path in PROFILE_DIR.iterdir()],
        help="The FHIR profile to use for generating the output.",
    )
    parser.add_argument(
        "--sample-allelic-frequency",
        default=None,
        help="Python expression calculating the the samples allelic frequency"
        "as percentage. E.g. \"FORMAT['AF'][sample][0] * 100\" ",
    )
    parser.add_argument(
        "--sample-allelic-read-depth",
        default="FORMAT['AD'][sample][1]",
        help="Python expression accessing the the samples allelic read depth."
        "Default is: \"FORMAT['AD'][sample][1]\"",
    )
    parser.add_argument(
        "--confidence-status",
        default=None,
        help="Python expression for calculating the variants confidence status being "
        "High, Intermediate or Low. E.g. \"'High' if QUAL >= 20 else ('Intermediate' "
        "if QUAL >= 10 else 'Low')\"",
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


def load_tsv(path, skip_comments=True):
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        if skip_comments:
            return [line for line in reader if not line[0].startswith("#")]
        else:
            return list(reader)


class Cytobands:
    def __init__(self):
        self._data = defaultdict(IntervalTree)
        records = load_tsv(
            importlib.resources.files("vembrane.modules")
            / "assets"
            / "fhir"
            / "cytobands.txt"
        )
        for chrom, start, end, band, _ in records:
            self._data[chrom].addi(int(start), int(end), band)

    def get(self, chrom: str, pos: int) -> str | None:
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        intervals = self._data.get(chrom)
        if intervals:
            hits = intervals[pos]
            if hits:
                return next(iter(hits)).data
        return None


class Assemblies:
    def __init__(self):
        self._data = {
            assembly: loinc
            for assembly, loinc in load_tsv(
                importlib.resources.files("vembrane.modules")
                / "assets"
                / "fhir"
                / "assemblies.txt",
                skip_comments=False,
            )
        }

    def get(self, assembly: str) -> str | None:
        return self._data.get(assembly)


class Chromosomes:
    def __init__(self):
        self._data = {
            chrom: (chrom_display, loinc)
            for chrom, chrom_display, loinc in load_tsv(
                importlib.resources.files("vembrane.modules")
                / "assets"
                / "fhir"
                / "chromosomes.txt",
                skip_comments=False,
            )
        }

    def get(self, chrom: str) -> tuple[str, str] | None:
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        return self._data.get(chrom)


def execute(args):
    with open(PROFILE_DIR / f"{args.profile}.yaml", mode="r") as infile:
        template = infile.read()
    structured.process(
        template=template,
        output_fmt=args.output_fmt,
        output=args.output,
        vcf=args.vcf,
        annotation_key=args.annotation_key,
        overwrite_number_info=args.overwrite_number_info,
        overwrite_number_format=args.overwrite_number_format,
        backend=args.backend,
        variables={
            "sample": args.sample,
            "url": args.url,
            "status": args.status,
            "assembly": args.assembly,
            "coordinates": args.coordinates,
            "sample_allelic_frequency": args.sample_allelic_frequency,
            "sample_allelic_read_depth": args.sample_allelic_read_depth,
            "confidence_status": args.confidence_status,
            "cytobands": Cytobands(),
            "assemblies": Assemblies(),
            "chromosomes": Chromosomes(),
        },
        postprocess=postprocess_fhir_record,
    )


def postprocess_fhir_record(record: Primitive | dict | list) -> Primitive | dict | list:
    assert isinstance(record, dict), "bug: FHIR record is expected to be a dictionary"

    def is_valid(entry: Any) -> bool:
        if isinstance(entry, dict):
            return all(is_valid(value) for value in entry.values())
        elif isinstance(entry, list):
            return all(is_valid(item) for item in entry)
        else:
            return entry is not None

    record["component"] = list(filter(is_valid, record["component"]))
    return record
