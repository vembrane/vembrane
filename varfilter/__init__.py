__version__ = "0.1.0"

from typing import Iterator

from pysam import VariantFile, VariantRecord
from sys import argv

# import stuff we want to be available in eval by default:
import re


def filter_vcf(vcf: VariantFile, expression: str) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    for name in header.info:
        env[name] = None

    for record in vcf:
        for key in record.info:
            env[key] = record.info[key]
        if eval(expression, locals=env):
            yield record


def main():
    expression = argv[2]
    with VariantFile(argv[1]) as vcf:
        with VariantFile('-', 'w', header=vcf.header) as out:
            for record in filter_vcf(vcf, expression):
                out.write(record)
