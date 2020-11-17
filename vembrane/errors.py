class VembraneError(Exception):
    """Basic exception for errors raised by vembrane"""


class UnknownAnnotation(VembraneError):
    """Unknown annotation entry"""

    def __init__(self, record, key, msg=None):
        if msg is None:
            msg = f"No ANN entry for '{key}' in record {record}"
        super(UnknownAnnotation, self).__init__(msg)
        self.record = record
        self.key = key


class UnknownSample(VembraneError, KeyError):
    """Unknown Sample"""

    def __init__(self, record, sample, msg=None):
        if msg is None:
            msg = f"No sample with name '{sample}' in record {record}"
        super(UnknownSample, self).__init__(msg)
        self.record = record
        self.field = sample


class UnknownFormatField(VembraneError, KeyError):
    """Unknown FORMAT key"""

    def __init__(self, record, field, msg=None):
        if msg is None:
            msg = f"No FORMAT field '{field}' in record {record}"
        super(UnknownFormatField, self).__init__(msg)
        self.record = record
        self.field = field


class UnknownInfoField(VembraneError):
    """Unknown INFO key"""

    def __init__(self, record, field, msg=None):
        if msg is None:
            msg = f"No INFO field '{field}' in record {record}"
        super(UnknownInfoField, self).__init__(msg)
        self.record = record
        self.field = field


class InvalidExpression(VembraneError):
    """Filter expression is invalid"""

    def __init__(self, expression, reason, msg=None):
        if msg is None:
            msg = f"The provided expression '{expression}' is invalid. Reason: {reason}"
        super(InvalidExpression, self).__init__(msg)
        self.expression = expression


class MoreThanOneAltAllele(VembraneError):
    """vembrane only supports one ALT allele per record"""

    def __init__(self):
        msg = (
            "vembrane only supports records with one alternative allele.\n"
            "Please split multi-allelic records first, for example with "
            "`bcftools norm -m-any […]` or "
            "`gatk LeftAlignAndTrimVariants […] --split-multi-allelics` or "
            "`vcfmulti2oneallele […]`"
        )
        super(MoreThanOneAltAllele, self).__init__(msg)
