from pathlib import Path
import os

from pysam import VariantFile
import pytest
import yaml
from vembrane import errors
from vembrane import __version__, filter_vcf

CASES = Path(__file__).parent.joinpath("testcases")


def test_version():
    assert __version__ == "0.1.0"


@pytest.mark.parametrize(
    "testcase", [d for d in os.listdir(CASES) if not d.startswith(".")]
)
def test_filter(testcase):
    path = CASES.joinpath(testcase)

    with open(path.joinpath("config.yaml")) as config_fp:
        config = yaml.load(config_fp, Loader=yaml.FullLoader)

    vcf = VariantFile(path.joinpath("test.vcf"))
    if "raises" in config:
        exception = getattr(errors, config["raises"])

        from vembrane import check_filter_expression

        with pytest.raises(exception):
            # FIXME we have to explicitly check the filter expression here
            # until we change from calling filter_vcf
            # to actually invoking vembrane.main
            check_filter_expression(config.get("filter_expression"))
            list(
                filter_vcf(
                    vcf,
                    config.get("filter_expression"),
                    config.get("ann_key", "ANN"),
                    config.get("keep_unmatched", False),
                )
            )
    else:
        expected = list(VariantFile(path.joinpath("expected.vcf")))
        result = list(
            filter_vcf(
                vcf,
                config.get("filter_expression"),
                config.get("ann_key", "ANN"),
                config.get("keep_unmatched", False),
            )
        )
        assert result == expected
