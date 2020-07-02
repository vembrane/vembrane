__version__ = "0.1.0"

from typing import Iterator

from pysam import VariantFile, VariantRecord
from sys import argv

# import stuff we want to be available in eval by default:
import re, argparse


def filter_vcf(vcf: VariantFile, expression: str) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    for name in header.info:
        env[name] = None

    annotation_keys = []
    for rec in header.records:
        if rec.get("ID") == "ANN":
            annotation_keys = list(
                map(str.strip, rec.get("Description").split("'")[1].split("|"))
            )
            break

    for record in vcf:
        for key in record.info:
            env[key] = record.info[key]
        ann = env.get("ANN", [])
        env["ANNO"] = dict(
            zip(
                annotation_keys, zip(*(list(map(str.strip, a.split("|"))) for a in ann))
            )
        )
        if eval(expression, dict(), env):
            yield record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="The file containing the variants.")
    parser.add_argument("expression", help="An expression to filter the variants.")
    args = parser.parse_args()

    with VariantFile(args.vcf) as vcf:
        with VariantFile("-", "w", header=vcf.header) as out:
            for record in filter_vcf(vcf, args.expression):
                out.write(record)


if __name__ == "__main__":
    main()
