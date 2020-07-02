from pathlib import Path
import os

from pysam import VariantFile
import pytest

from varfilter import __version__, filter_vcf

CASES = Path(__file__).parent.joinpath("testcases")


def test_version():
    assert __version__ == "0.1.0"


@pytest.mark.parametrize(
    "testcase", [d for d in os.listdir(CASES) if not d.startswith(".")]
)
def test_filter(testcase):
    path = CASES.joinpath(testcase)

    expression = open(path.joinpath("expression.txt")).read().strip()
    vcf = VariantFile(path.joinpath("test.vcf"))

    expected = list(VariantFile(path.joinpath("expected.vcf")))
    result = list(filter_vcf(vcf, expression))

    assert result == expected
