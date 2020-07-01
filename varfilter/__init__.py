__version__ = "0.1.0"

from typing import Iterator

from pysam import VariantFile, VariantRecord
from sys import argv

env = {
    "__builtins__": None,
    "__file__": None,
    "__name__": None,
    "globals": None,
    "locals": None
}


def filter_vcf(vcf: VariantFile, expression: str) -> Iterator[VariantRecord]:
    header = vcf.header
    for name in header.info:
        vars()[name] = None

    for record in vcf:
        for key in record.info:
            vars()[key] = record.info[key]
        # TODO properly restrict env and locals
        available_vars = locals()
        if eval(expression, env, available_vars):
            yield record


if __name__ == "__main__":
    expression = argv[2]
    with VariantFile(argv[1]) as vcf:
        with VariantFile('-', 'w', header=vcf.header) as out:
            for record in filter_vcf(vcf, expression):
                out.write(record)
