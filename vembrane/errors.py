from pysam import VariantRecord


class VembraneError(Exception):
    """Basic exception for errors raised by vembrane"""

    def __str__(self) -> str:
        return self.args[0]


class UnknownAnnotation(VembraneError):
    """Unknown annotation entry"""

    def __init__(self, record_idx: int, record: VariantRecord, key: str) -> None:
        super().__init__(
            f"No ANN entry for '{key}' in record {record_idx}:\n{str(record)}"
        )
        self.record_idx = record_idx
        self.record = record
        self.key = key


class UnknownSample(VembraneError, KeyError):
    """Unknown Sample"""

    def __init__(self, record_idx: int, record: VariantRecord, sample: str) -> None:
        super().__init__(
            f"No sample with name '{sample}' in record {record_idx}:\n{str(record)}"
        )
        self.record_idx = record_idx
        self.record = record
        self.field = sample


class UnknownFormatField(VembraneError, KeyError):
    """Unknown FORMAT key"""

    def __init__(self, record_idx: int, record: VariantRecord, field: str) -> None:
        super().__init__(
            f"No FORMAT field '{field}' in record {record_idx}:\n{str(record)}"
        )
        self.record_idx = record_idx
        self.record = record_idx
        self.field = field


class UnknownInfoField(VembraneError):
    """Unknown INFO key"""

    def __init__(self, record_idx: int, record: VariantRecord, field: str) -> None:
        super().__init__(
            f"No INFO field '{field}' in record {record_idx}:\n{str(record)}"
        )
        self.record = record_idx
        self.field = field


class InvalidExpression(VembraneError):
    """Filter expression is invalid"""

    def __init__(self, expression: str, reason: str) -> None:
        super().__init__(
            f"The provided expression '{expression}' is invalid. Reason: {reason}"
        )
        self.expression = expression


class MoreThanOneAltAllele(VembraneError):
    """vembrane only supports one ALT allele per record"""

    def __init__(self) -> None:
        msg = (
            "vembrane only supports records with one alternative allele.\n"
            "Please split multi-allelic records first, for example with "
            "`bcftools norm -m-any […]` or "
            "`gatk LeftAlignAndTrimVariants […] --split-multi-allelics` or "
            "`vcfmulti2oneallele […]`"
        )
        super().__init__(msg)


class NotExactlyOneValue(VembraneError):
    """There may only be one value in VCF fields with `Number=1`"""

    def __init__(self, field, nvalues, record_idx) -> None:
        msg = (
            f"record {record_idx} has {nvalues} values in {field}, "
            f"but defined `Number=1`.\n"
            "vembrane enforces fields with `Number=1` to actually only have "
            "one number, "
            "(see VCF specification v4.3 section 1.4.2).\n"
            "To override this behaviour, "
            f"use `--overwrite-number-[info|format] {field} NUM`, "
            "where `NUM` is one of: \n"
            "- A: the field has one value per alternate allele\n"
            "- R: the field has one value for each possible allele, "
            "including the reference.\n"
            "- .: the number of possible values varies, is unknown or unbounded\n"
            "- G (FORMAT only): the field has one value for each possible genotype\n"
            "- {number}: the field has exactly {NUMBER} values\n"
        )
        super().__init__(msg)


class HeaderWrongColumnNumber(VembraneError):
    """
    table --header expression generates different number of columns than main expression
    """

    def __init__(self, n_expr_cols, expr_cols, n_header_cols, header_cols) -> None:
        msg = (
            "The provided --header expression generates a different\n"
            "number of columns than the main expression:\n"
            f"{n_expr_cols} main expression columns: {expr_cols}\n"
            f"{n_header_cols} --header expression columns: {header_cols}\n"
        )
        super().__init__(msg)


class FilterAlreadyDefined(VembraneError):
    def __init__(self, tag: str) -> None:
        super().__init__(
            f"Filter {tag} already defined in the header. Choose a different name."
        )


class FilterTagNameInvalid(VembraneError):
    def __init__(self, tag: str) -> None:
        super().__init__(
            f"Filter '{tag}' contains invalid characters (whitespace or semicolon) "
            f"or is '0'."
        )
