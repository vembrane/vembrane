import csv
import importlib.resources
import sys
from collections import defaultdict
from typing import Any

# intervaltree is untyped, so we use type: ignore to suppress type checking errors
from intervaltree import IntervalTree  # type: ignore

from vembrane.common import (
    Context,
    HumanReadableDefaultsFormatter,
    Primitive,
    Singleton,
    add_common_arguments,
)
from vembrane.errors import VembraneError
from vembrane.globals import default_allowed_globals
from vembrane.modules import structured

PROFILE_DIR = (
    importlib.resources.files("vembrane.modules") / "assets" / "fhir" / "profiles"
)


def add_subcommand(subparsers):
    parser = subparsers.add_parser(
        "fhir",
        help="Generate FHIR records from VCF/BCF files.",
        description="Generate FHIR records from VCF/BCF files.",
        formatter_class=HumanReadableDefaultsFormatter,
    )
    parser.add_argument(
        "vcf",
        help="Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin.",
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
        choices=Assemblies().names(),
    )
    parser.add_argument(
        "--url",
        "-u",
        help="Generic url used as identifier by FHIR e.g. "
        "http://<institute>/<department>/VCF",
        default=None,
    )
    parser.add_argument(
        "--status",
        "-s",
        help="Status of findings. E.g. final, preliminary, ...",
        default=None,
    )
    parser.add_argument(
        "--profile",
        required=True,
        choices=[path.with_suffix("").name for path in PROFILE_DIR.iterdir()],
        help="The FHIR profile to use for generating the output, see "
        "https://github.com/vembrane/vembrane/tree/main/vembrane/modules/assets/"
        "fhir/profiles for available profiles and the degree of support.",
    )
    parser.add_argument(
        "--id-source",
        help="URL to the source of IDs found in the ID column of the VCF/BCF file. "
        "IDs are only used if this is given.",
    )
    parser.add_argument(
        "--genomic-source-class",
        help="The genomic source class of the given variants as defined by "
        "LOINC: https://loinc.org/48002-0. Either provide the name as a string or "
        "a Python expression that evaluates to the name, e.g., for Varlociraptor "
        '\'"Somatic" if INFO["PROB_SOMATIC"] > 0.95 else ...\'.',
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
        "--detection-limit",
        type=int,
        help="Detection limit / sensitivity of the analysis in percent (e.g. 95).",
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


class Cytobands(metaclass=Singleton):
    def __init__(self):
        self._data = defaultdict(IntervalTree)
        records = load_tsv(
            importlib.resources.files("vembrane.modules")
            / "assets"
            / "fhir"
            / "cytobands.tsv"
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


class Assemblies(metaclass=Singleton):
    def __init__(self):
        self._data = {
            assembly: loinc
            for assembly, loinc in load_tsv(
                importlib.resources.files("vembrane.modules")
                / "assets"
                / "fhir"
                / "assemblies.tsv",
                skip_comments=False,
            )
        }

    def get(self, assembly: str) -> str | None:
        return self._data.get(assembly)

    def names(self) -> list[str]:
        """
        Get a list of all assembly names.
        """
        return list(self._data.keys())


class Chromosomes(metaclass=Singleton):
    def __init__(self):
        self._data = {
            chrom: (chrom_display, loinc)
            for chrom, chrom_display, loinc in load_tsv(
                importlib.resources.files("vembrane.modules")
                / "assets"
                / "fhir"
                / "chromosomes.tsv",
                skip_comments=False,
            )
        }

    def get(self, chrom: str) -> tuple[str, str] | None:
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        return self._data.get(chrom)


class GenomicSourceClasses(metaclass=Singleton):
    def __init__(self):
        self._data = {
            name: code
            for name, code in load_tsv(
                importlib.resources.files("vembrane.modules")
                / "assets"
                / "fhir"
                / "genomic_source_classes.tsv"
            )
        }

    def code(self, name: str | None) -> str | None:
        """
        Get the LOINC code for the given genomic source class name.
        """
        if name is None:
            return None
        return self._data[name]

    def is_name(self, expr: str | None) -> bool:
        """
        Check if the given expression is a valid genomic source class name.
        """
        if expr is None:
            return False
        return expr in self._data

    def preprocess_expression(self, expr: str | None) -> str | None:
        """
        Preprocess the expression to ensure it is a valid genomic source class name.
        If the expression is a string, it will be wrapped in quotes.
        """
        if expr is None:
            return None
        if self.is_name(expr):
            return f'"{expr}"'
        return expr


def execute(args):
    try:
        if args.detection_limit is not None:
            if not (0 < args.detection_limit <= 100):
                raise ValueError("Detection limit must be between 1 and 100.")

        genomic_source_classes = GenomicSourceClasses()

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
            allowed_globals=default_allowed_globals
            | {
                "sample": args.sample,
                "url": args.url,
                "status": args.status,
                "assembly": args.assembly,
                "sample_allelic_frequency": args.sample_allelic_frequency,
                "sample_allelic_read_depth": args.sample_allelic_read_depth,
                "confidence_status": args.confidence_status,
                "cytobands": Cytobands(),
                "assemblies": Assemblies(),
                "chromosomes": Chromosomes(),
                "genomic_source_classes": genomic_source_classes,
                "genomic_source_class": genomic_source_classes.preprocess_expression(
                    args.genomic_source_class
                ),
                "id_source": args.id_source,
                "detection_limit": args.detection_limit,
            },
            auxiliary_globals=Context.from_args(args).get_globals(),
            postprocess=postprocess_fhir_record,
        )
    except VembraneError as ve:
        print(ve, file=sys.stderr)
        sys.exit(1)


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
