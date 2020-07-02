__version__ = "0.1.0"

from typing import Iterator

from pysam import VariantFile, VariantRecord
from sys import argv

# import stuff we want to be available in eval by default:
import re
import math
from math import log, log2, log10, log1p

globals_whitelist = {
    **{
        "__builtins__": None,
        "__builtin__": None,
        "__file__": None,
        "__name__": None,
        "__doc__": None,
        "__package__": None,
    },
    **{
        mod.__name__: mod
        for mod in [any, all, min, max, re, log, log2, log10, log1p, list, dict, zip]
    },
}


def filter_vcf(vcf: VariantFile, expression: str) -> Iterator[VariantRecord]:
    header = vcf.header

    env = dict()

    for name in header.info:
        env[name] = None

    annotation_keys = []
    for rec in header.records:
        if rec.get("ID") == "ANN":
            annotation_keys = list(map(str.strip, rec.get("Description").split("'")[1].split("|")))
            break

    for record in vcf:
        for key in record.info:
            env[key] = record.info[key]
        ann = env.get("ANN", [])
        env["ANNO"] = dict(zip(annotation_keys, zip(*(list(map(str.strip, a.split('|'))) for a in ann))))
        if eval(expression, locals=env):
            yield record


def main():
    expression = argv[2]
    with VariantFile(argv[1]) as vcf:
        with VariantFile("-", "w", header=vcf.header) as out:
            for record in filter_vcf(vcf, expression):
                out.write(record)
