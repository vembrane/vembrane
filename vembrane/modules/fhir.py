import csv
import importlib
from collections import defaultdict

from intervaltree import IntervalTree

from vembrane.common import add_common_arguments
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
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
    )
    parser.add_argument(
        "--profile",
        required=True,
        choices=[path.with_suffix("").name for path in PROFILE_DIR.iterdir()],
        help="The FHIR profile to use for generating the output.",
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


def load_tsv(path, skip_comment=True):
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        if skip_comment:
            next(reader)
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
        for chrom, start, end, band in records:
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
                skip_comment=False,
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
                skip_comment=False,
            )
        }

    def get(self, chrom: str) -> tuple[str, str] | None:
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        return self._data.get(chrom)


class CodingChangeType:
    def __init__(self):
        self._data = {
            type_short: (type_display, loinc)
            for type_short, type_display, loinc in load_tsv(
                importlib.resources.files("vembrane.modules")
                / "assets"
                / "fhir"
                / "coding_change_type.txt",
                skip_comment=False,
            )
        }

    def get(self, type_short: str) -> tuple[str, str] | None:
        return self._data.get(type_short)


def execute(args):
    with open(PROFILE_DIR / f"{args.profile}.yaml", "r") as infile:
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
            "cytobands": Cytobands(),
            "assemblies": Assemblies(),
            "chromosomes": Chromosomes(),
            "coding_change_type": CodingChangeType(),
        },
    )
