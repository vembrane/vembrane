import sys
from typing import TYPE_CHECKING, Any, List, Self

if TYPE_CHECKING:
    from vembrane.backend.base import VCFRecord


class VembraneError(Exception):
    """Basic exception for errors raised by vembrane"""

    @classmethod
    def from_record_and_exception(
        cls, idx: int, record: "VCFRecord", e: Exception
    ) -> Self:
        """
        Create a VembraneError from a VCFRecord and an Exception.
        This is useful to preserve the record context in the error message.
        """
        return cls(
            f"Error processing record {idx}: {e}\nRecord: {record}",
        )

    def __str__(self) -> str:
        return self.args[0]


class UnknownAnnotationError(VembraneError):
    """Unknown annotation entry"""

    def __init__(self, record, key: str) -> None:
        super().__init__(
            f"No ANN entry for '{key}' in record {record.record_idx}:\n{str(record)}",
        )
        self.record_idx = record.record_idx
        self.record = record
        self.key = key


class MalformedAnnotationError(VembraneError):
    """Malformed annotation entry"""

    def __init__(self, record, key: str, ann_idx: int) -> None:
        super().__init__(
            f"The annotation index {ann_idx} ('{key}') "
            f"in record {record.record_idx} is out of bounds. "
            f"The ANN field might be malformed for this record:\n{str(record)}",
        )
        self.record_idx = record.record_idx
        self.record = record
        self.key = key
        self.ann_idx = ann_idx


class UnknownSampleError(VembraneError, KeyError):
    """Unknown Sample"""

    def __init__(self, record, sample: str) -> None:
        super().__init__(
            f"No sample with name '{sample}' in record \
                {record.record_idx}:\n{str(record)}",
        )
        self.record_idx = record.record_idx
        self.record = record
        self.field = sample


class UnknownFormatFieldError(VembraneError, KeyError):
    """Unknown FORMAT key"""

    def __init__(self, record, field: str) -> None:
        super().__init__(
            f"No FORMAT field '{field}' in record \
                {record.record_idx}:\n{str(record)}",
        )
        self.record_idx = record.record_idx
        self.record = record
        self.field = field


class UnknownInfoFieldError(VembraneError):
    """Unknown INFO key"""

    def __init__(self, record, field: str) -> None:
        super().__init__(
            f"No INFO field '{field}' in record {record.record_idx}:\n{str(record)}",
        )
        self.record_idx = record.record_idx
        self.record = record
        self.field = field


class InvalidExpressionError(VembraneError):
    """Filter expression is invalid"""

    def __init__(self, expression: str, reason: str) -> None:
        super().__init__(
            f"The provided expression '{expression}' is invalid. Reason: {reason}",
        )
        self.expression = expression


class MoreThanOneAltAlleleError(VembraneError):
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


class NotExactlyOneValueError(VembraneError):
    """There may only be one value in VCF fields with `Number=1`"""

    def __init__(self, field: str, nvalues: int, record_idx: int) -> None:
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


class HeaderWrongColumnNumberError(VembraneError):
    """
    table --header expression generates different number of columns than main expression
    """

    def __init__(
        self,
        n_expr_cols: int,
        expr_cols: List[str],
        n_header_cols: int,
        header_cols: List[str],
    ) -> None:
        msg = (
            "The provided --header expression generates a different\n"
            "number of columns than the main expression:\n"
            f"{n_expr_cols} main expression columns: {expr_cols}\n"
            f"{n_header_cols} --header expression columns: {header_cols}\n"
        )
        super().__init__(msg)


class FilterAlreadyDefinedError(VembraneError):
    def __init__(self, tag: str) -> None:
        super().__init__(
            f"Filter {tag} already defined in the header. Choose a different name.",
        )


class FilterTagNameInvalidError(VembraneError):
    def __init__(self, tag: str) -> None:
        super().__init__(
            f"Filter '{tag}' contains invalid characters (whitespace or semicolon) "
            f"or is '0'.",
        )


class NonBoolTypeError(VembraneError):
    def __init__(self, value: Any):
        super().__init__(
            "The expression does not evaluate to bool, "
            f"but to {type(value)} ({value}).\n"
            "If you wish to use truthy values, "
            "explicitly wrap the expression in `bool(…)`, "
            "or aggregate multiple values via `any(…)` or `all(…)`.",
        )


class UnsupportedChromName(VembraneError):
    def __init__(self, chrom: str):
        super().__init__(f"Unsupported chromosome name: {chrom}.")


def handle_vembrane_error(func):
    """
    Decorator to handle VembraneError exceptions and print a user-friendly message.
    """

    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except VembraneError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    return wrapper
