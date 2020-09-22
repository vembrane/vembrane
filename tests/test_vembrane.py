import os
from pathlib import Path

import pytest
import yaml
from pysam import VariantFile

from vembrane import errors, __version__
from vembrane.common import check_expression
from vembrane.modules.filter import filter_vcf
from vembrane.modules.table import tableize_vcf

CASES = Path(__file__).parent.joinpath("testcases")


def test_version():
    assert __version__ != "unknown"


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

        with pytest.raises(exception):
            # FIXME we have to explicitly check the filter expression here
            # until we change from calling filter_vcf
            # to actually invoking vembrane.main
            check_expression(config.get("expression"))
            if config["function"] == "filter":
                list(
                    filter_vcf(
                        vcf,
                        config.get("expression"),
                        config.get("ann_key", "ANN"),
                        config.get("keep_unmatched", False),
                    )
                )
            elif config["function"] == "table":
                list(
                    tableize_vcf(
                        vcf,
                        config.get("expression"),
                        config.get("ann_key", "ANN"),
                    )
                )
            else:
                assert False
    else:
        if config["function"] == "filter":
            expected = list(VariantFile(path.joinpath("expected.vcf")))
            result = list(
                filter_vcf(
                    vcf,
                    config.get("expression"),
                    config.get("ann_key", "ANN"),
                    config.get("keep_unmatched", False),
                )
            )
            assert result == expected
        elif config["function"] == "table":
            separator = config.get("separator", "\t")
            expected = list(
                map(
                    lambda x: x.strip("\n"),
                    open(path.joinpath("expected.tsv"), "r").readlines(),
                )
            )
            result = list(
                map(
                    lambda x: separator.join(map(str, x)),
                    tableize_vcf(
                        vcf,
                        config.get("expression"),
                        config.get("ann_key", "ANN"),
                    ),
                )
            )
            assert result == expected[1:]  # avoid header check by now
        else:
            assert config["function"] in {"filter", "table"}, "Unknown subcommand"
