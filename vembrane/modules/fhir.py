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


class Cytobands:
    def __init__(self):
        self._data = defaultdict(IntervalTree)

        with open(
            importlib.resources.files("vembrane.modules")
            / "assets"
            / "fhir"
            / "cytobands.txt",
            "r",
        ) as asset:
            reader = csv.reader(asset, delimiter="\t")
            next(reader)  # skip leading comment
            for record in reader:
                chrom, start, end, band = record[:4]
                # cytoband records are 0-based, end-exclusive
                self._data[chrom].addi(int(start), int(end), band)

    def get(self, chrom: str, pos: int) -> str | None:
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        try:
            return self._data[chrom][pos].pop().data
        except KeyError:
            return None


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
        variables={"sample": args.sample, "cytobands": Cytobands()},
    )
