from __future__ import annotations

import re
from collections import defaultdict
from ctypes import c_float
from sys import stderr
from typing import Any, Callable, Iterable, Union

import numpy as np

from .errors import MoreThanOneAltAlleleError, NotExactlyOneValueError
from .sequence_ontology import Consequences, Term


def float32(val: str) -> float:
    return c_float(float(val)).value


# If NoValue inherits from str, re.search("something", NoValue()) does not error
# but just comes up empty-handed, which is convenient behaviour.
# This way, we do not have to special case / monkey patch / wrap the regex module.
class NoValue(str):
    warnings: set[str] = set()

    def __lt__(self, other) -> bool:
        return False

    def __gt__(self, other) -> bool:
        return False

    def __le__(self, other) -> bool:
        return False

    def __ge__(self, other) -> bool:
        return False

    def __eq__(self, other) -> bool:
        return False

    def __ne__(self, other) -> bool:
        return True

    def __bool__(self) -> bool:
        return False

    def __str__(self) -> str:
        # nonexistent fields will result in an empty string
        # should be configurable in upcoming versions
        return ""

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"

    def __hash__(self) -> int:
        return super().__hash__()

    def __getattr__(self, item):
        if item not in self.warnings:
            self.warnings.add(item)
            print(
                f"Warning: Trying to access non-existent attribute '{item}' of NoValue."
                " Returning NA instead.\n"
                "This warning will only be printed once per attribute.\n"
                "It either indicates a typo in the attribute name or "
                "the access of an attribute of a field with a custom type "
                "which is empty/has no value.",
                file=stderr,
            )
        return self


NA = NoValue()


class InfoTuple:
    """A container that lazily evaluates None to NA in case of access."""

    def __init__(self, values) -> None:
        self.values = values

    def __getitem__(self, spec):
        values = self.values.__getitem__(spec)
        if isinstance(values, tuple) and any(v is None for v in values):
            return tuple((NA if v is None else v) for v in values)
        if values is None:
            return NA
        return values

    def __str__(self) -> str:
        return self.values.__str__()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(values={self.values!r})"

    def __len__(self) -> int:
        return len(self.values)

    def __eq__(self, other):
        value_type = type(self.values)
        if isinstance(other, InfoTuple) and isinstance(other.values, value_type):
            return self.values == other.values
        elif isinstance(other, value_type):
            return self.values == other
        else:
            raise TypeError(
                f"Incomparable: "
                f"{type(self)} (value type: {type(self.values)}), "
                f"{type(other)}"
                + (
                    f" (value type: {type(other.values)})"
                    if isinstance(other, InfoTuple)
                    else ""
                ),
            )

    def __hash__(self):
        return hash(self.values)


IntFloatStr = Union[int, float, str]
NvIntFloatStr = Union[IntFloatStr, NoValue]
NvInt = Union[int, NoValue]
NvFloat = Union[float, NoValue]


def type_info(
    value,
    number=".",
    field=None,
    record_idx=None,
) -> NvIntFloatStr | InfoTuple:
    if value is None:
        return NA
    if value is NA:
        if number == "0":
            return False
        return NA
    if number == "0":
        return False if value is NA else value
    if isinstance(value, np.ndarray):
        value = tuple(value.tolist())
    elif isinstance(value, list):
        value = tuple(value)
    if number == "A":
        if isinstance(value, tuple):
            if len(value) > 1:
                raise MoreThanOneAltAlleleError()
            return value[0] if value[0] is not None else NA
        return value
    if number == "R":
        if len(value) > 2:
            raise MoreThanOneAltAlleleError()
        if len(value) == 1:
            return InfoTuple(tuple(value) + (None,))
        else:
            return InfoTuple(tuple(value))
    if number == "1":
        if isinstance(value, tuple):
            if len(value) != 1:
                raise NotExactlyOneValueError(field, len(value), record_idx)
            return value[0] if value[0] is not None else NA
        return value if value is not None else NA
    if isinstance(value, tuple):
        return InfoTuple(value)
    else:
        return InfoTuple((value,))


class PosRange:
    def __init__(self, start: NvInt, end: NvInt, raw: str) -> None:
        self.start = start
        self.end = end
        self._raw = raw

    raw = property(lambda self: self._raw, None, None)

    @classmethod
    def from_snpeff_str(cls, value: str) -> PosRange:
        pos, length = (int(v.strip()) for v in value.split("/"))
        return cls(pos, pos + length, value)

    @classmethod
    def from_vep_str(cls, value: str) -> PosRange:
        if "-" in value:
            s, e = map(str.strip, value.split("-", 1))
            start: NvInt = NA if s == "?" else int(s)
            end: NvInt = NA if e == "?" else int(e)
            return cls(start, end, value)
        elif len(value) > 0:
            start = int(value.strip())
            return cls(start, start + 1, value)
        else:
            return cls(NA, NA, value)

    @property
    def length(self):
        if self.end is not NA and self.start is not NA:
            return self.end - self.start
        else:
            return NA

    def __str__(self) -> str:
        return f"(start: {self.start}, end: {self.end}, length: {self.length})"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(start={self.start!r}, end={self.end!r})"


def na_func() -> NoValue:
    return NA


AnnotationType = Union[IntFloatStr, Iterable[IntFloatStr], NoValue]


class AnnotationEntry:
    def __init__(
        self,
        name: str,
        typefunc: Callable[[str], Any] = str,
        nafunc: Callable[[], Any] = na_func,
        description: str | None = None,
    ) -> None:
        self._name = name
        self._typefunc = typefunc
        self._nafunc = nafunc
        self._description = description

    @property
    def name(self):
        return self._name

    def convert(self, value: str) -> Any:
        if value:
            return self._typefunc(value)
        else:
            return self._nafunc()

    def description(self):
        return self._description

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"name={self._name!r}, "
            f"typefunc={self._typefunc!r}, "
            f"nafunc={self._nafunc!r}, "
            f"description={self._description!r}"
            f")"
        )


class AnnotationListEntry(AnnotationEntry):
    def __init__(
        self,
        name: str,
        sep: str,
        typefunc: Callable[[str], Any] | None = None,
        inner_typefunc: Callable[[str], Any] = lambda x: x,
        **kwargs,
    ) -> None:
        typefunc = (
            typefunc
            if typefunc
            else lambda v: [inner_typefunc(x.strip()) for x in v.split(sep)]
        )
        super().__init__(name, typefunc, nafunc=lambda: [], **kwargs)


class AnnotationSetEntry(AnnotationEntry):
    def __init__(
        self,
        name: str,
        sep: str,
        typefunc: Callable[[str], Any] | None = None,
        inner_typefunc: Callable[[str], Any] = lambda x: x,
        nafunc: Callable[[], Any] = lambda: set(),
        **kwargs,
    ) -> None:
        typefunc = (
            typefunc
            if typefunc
            else lambda v: {inner_typefunc(x.strip()) for x in v.split(sep)}
        )
        super().__init__(name, typefunc, nafunc=nafunc, **kwargs)


class ConsequenceEntry(AnnotationSetEntry):
    def __init__(self, name: str, **kwargs) -> None:
        super().__init__(
            name,
            sep="&",
            typefunc=lambda v: Consequences(Term(t.strip()) for t in v.split("&")),
            nafunc=lambda: Consequences(),
            **kwargs,
        )


class AnnotationListDictEntry(AnnotationEntry):
    def __init__(self, name: str, nafunc=lambda: None, **kwargs) -> None:
        def typefunc(x):
            key_value_tuples = (v.split(":", maxsplit=1) for v in x.split("&"))
            mapping = defaultdict(list)
            for key, value in key_value_tuples:
                mapping[key].append(value)
            # make sure non-existent keys do not return list() on __getitem__ calls
            # in the filter expression (without actually building a new dict from
            # the defaultdict)
            mapping.default_factory = nafunc
            return mapping

        super().__init__(name=name, typefunc=typefunc, **kwargs)


class RangeTotal:
    def __init__(self, r: range, total: int, raw: str) -> None:
        self.range = r
        self.total = total
        self._raw = raw

    raw = property(lambda self: self._raw, None, None)

    @classmethod
    def from_vep_string(cls, value: str) -> RangeTotal:
        v = value.split("/")
        r = [int(s) for s in v[0].split("-")]
        if len(r) == 1:
            return cls(range(r[0], r[0] + 1), int(v[1]), value)
        elif len(r) == 2:
            return cls(range(r[0], r[1] + 1), int(v[1]), value)
        else:
            raise ValueError(
                "Found more than two values separated by '-', "
                "expected only a single int, or two ints separated by '-'.",
            )

    def __str__(self) -> str:
        if len(self.range) == 1:
            return f"number / total: {self.range.start} / {self.total}"
        else:
            return (
                f"range / total: "
                f"{self.range.start} - {self.range.stop - 1} / {self.total}"
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}"
            f"(range=range({self.range.start!r}, {self.range.stop!r}), "
            f"total={self.total!r})"
        )


class AnnotationRangeTotalEntry(AnnotationEntry):
    def __init__(self, name: str, **kwargs) -> None:
        super().__init__(name, typefunc=RangeTotal.from_vep_string, **kwargs)


class NumberTotal:
    def __init__(self, number: int, total: int, raw: str) -> None:
        self.number = number
        self.total = total
        self._raw = raw

    raw = property(lambda self: self._raw, None, None)

    @classmethod
    def from_vep_string(cls, value: str) -> NumberTotal:
        v = value.split("/")
        return cls(int(v[0]), int(v[1]), value)

    def __str__(self) -> str:
        return f"number / total: {self.number} / {self.total}"

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(number={self.number!r}, total={self.total!r})"
        )


class AnnotationNumberTotalEntry(AnnotationEntry):
    def __init__(self, name: str, **kwargs) -> None:
        super().__init__(name, typefunc=NumberTotal.from_vep_string, **kwargs)


PREDICTION_SCORE_REGEXP: re.Pattern = re.compile(r"(.*)\((.*)\)")


class AnnotationPredictionScoreEntry(AnnotationEntry):
    def __init__(self, name: str, **kwargs) -> None:
        def typefunc(x: str) -> dict[str, float]:
            match = PREDICTION_SCORE_REGEXP.findall(x)[0]
            return {match[0]: float32(match[1])}

        super().__init__(name, typefunc=typefunc, nafunc=lambda: {}, **kwargs)


class DefaultAnnotationEntry(AnnotationEntry):
    def __init__(self, name: str) -> None:
        super().__init__(name)


class AnnotationTyper:
    def __init__(self, mapping: dict[str, AnnotationEntry]) -> None:
        self._mapping = mapping

    def get_entry(self, key: str) -> AnnotationEntry:
        entry = self._mapping.get(key)
        if not entry:
            print(
                f"No type information available for '{key}', defaulting to `str`. "
                f"If you would like to have a custom type for this, "
                f"please consider filing an issue at "
                f"https://github.com/vembrane/vembrane/issues",
                file=stderr,
            )
            self._mapping[key] = DefaultAnnotationEntry(key)
            entry = self._mapping[key]
        return entry

    def convert(self, key: str, value: str) -> tuple[str, AnnotationType]:
        entry = self.get_entry(key)
        return entry.name, entry.convert(value)


KNOWN_ANN_TYPE_MAP_SNPEFF = {
    "Allele": AnnotationEntry("Allele"),
    "Annotation": ConsequenceEntry("Annotation"),
    "Annotation_Impact": AnnotationEntry("Annotation_Impact"),
    "Gene_Name": AnnotationEntry("Gene_Name"),
    "Gene_ID": AnnotationEntry("Gene_ID"),
    "Feature_Type": AnnotationEntry("Feature_Type"),
    "Feature_ID": AnnotationEntry("Feature_ID"),
    "Transcript_BioType": AnnotationEntry("Transcript_BioType"),
    "Rank": AnnotationEntry("Rank"),
    "HGVS.c": AnnotationEntry("HGVS.c"),
    "HGVS.p": AnnotationEntry("HGVS.p"),
    "cDNA.pos / cDNA.length": AnnotationEntry(
        "cDNA.pos / cDNA.length",
        PosRange.from_snpeff_str,
    ),
    "CDS.pos / CDS.length": AnnotationEntry(
        "CDS.pos / CDS.length",
        PosRange.from_snpeff_str,
    ),
    "AA.pos / AA.length": AnnotationEntry(
        "AA.pos / AA.length",
        PosRange.from_snpeff_str,
    ),
    "Distance": AnnotationEntry("Distance", str),
    "ERRORS / WARNINGS / INFO": AnnotationListEntry(
        "ERRORS / WARNINGS / INFO",
        "/",
    ),
}

# fields, types and description taken from:
# https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#other_fields
KNOWN_ANN_TYPE_MAP_VEP = {
    "Location": AnnotationEntry(
        "Location",
        description="In standard coordinate format (chr:start or chr:start-end)",
    ),
    "Allele": AnnotationEntry(
        "Allele",
        description="The variant allele used to calculate the consequence",
    ),
    "Gene": AnnotationEntry("Gene", description="Ensembl stable ID of affected gene"),
    "Feature": AnnotationEntry("Feature", description="Ensembl stable ID of feature"),
    "Feature_type": AnnotationEntry(
        "Feature_type",
        description="Type of feature. "
        "Currently one of Transcript, RegulatoryFeature, MotifFeature.",
    ),
    "Consequence": ConsequenceEntry(
        "Consequence",
        description="Consequence type of this variant",
    ),
    "cDNA_position": AnnotationEntry("cDNA_position", PosRange.from_vep_str),
    "CDS_position": AnnotationEntry("CDS_position", PosRange.from_vep_str),
    "Protein_position": AnnotationEntry("Protein_position", PosRange.from_vep_str),
    "HGSVc": AnnotationEntry("HGSVc"),
    "HGSVp": AnnotationEntry("HGSVp"),
    "REF_ALLELE": AnnotationEntry("REF_ALLELE", description="The reference allele"),
    "IMPACT": AnnotationEntry(
        "IMPACT",
        description="The impact modifier for the consequence type",
    ),
    "SYMBOL": AnnotationEntry("SYMBOL", description="The gene symbol"),
    "VARIANT_CLASS": AnnotationEntry(
        "VARIANT_CLASS",
        description="Sequence Ontology variant class",
    ),
    "SYMBOL_SOURCE": AnnotationEntry(
        "SYMBOL_SOURCE",
        description="The source of the gene symbol",
    ),
    "STRAND": AnnotationEntry(
        "STRAND",
        int,
        description="The DNA strand (1 or -1) on which the transcript/feature lies",
    ),
    "ENSP": AnnotationEntry(
        "ENSP",
        description="The Ensembl protein identifier of the affected transcript",
    ),
    "FLAGS": AnnotationListEntry(
        "FLAGS",
        sep="&",
        description="Transcript quality flags: "
        "cds_start_NF: CDS 5' incomplete, cds_end_NF: CDS 3' incomplete",
    ),
    "SWISSPROT": AnnotationEntry(
        "SWISSPROT",
        description="Best match UniProtKB/Swiss-Prot accession of protein product",
    ),
    "TREMBL": AnnotationEntry(
        "TREMBL",
        description="Best match UniProtKB/TrEMBL accession of protein product",
    ),
    "UNIPARC": AnnotationEntry(
        "UNIPARC",
        description="Best match UniParc accession of protein product",
    ),
    "HGVSc": AnnotationEntry("HGVSc", description="The HGVS coding sequence name"),
    "HGVSp": AnnotationEntry("HGVSp", description="The HGVS protein sequence name"),
    "HGVSg": AnnotationEntry("HGVSg", description="The HGVS genomic sequence name"),
    "HGVS_OFFSET": AnnotationEntry(
        "HGVS_OFFSET",
        int,
        description="Indicates by how many bases the HGVS notations for this variant "
        "have been shifted",
    ),
    # "NEAREST": AnnotationEntry("NEAREST",
    #               description="Identifier(s) of nearest transcription start site"),
    "SIFT": AnnotationPredictionScoreEntry(
        "SIFT",
        description="The SIFT prediction and/or score,"
        " with both given as prediction(score)",
    ),
    "PolyPhen": AnnotationPredictionScoreEntry(
        "PolyPhen",
        description="The PolyPhen prediction and/or score",
    ),
    "MOTIF_NAME": AnnotationEntry(
        "MOTIF_NAME",
        description="The source and identifier of a transcription factor binding "
        "profile aligned at this position",
    ),
    "MOTIF_POS": AnnotationEntry(
        "MOTIF_POS",
        int,
        description="The relative position of the variation in the aligned TFBP",
    ),
    "HIGH_INF_POS": AnnotationEntry(
        "HIGH_INF_POS",
        typefunc="Y".__eq__,
        description="A flag indicating if the variant falls in a high information "
        "position of a transcription factor binding profile (TFBP)",
    ),
    "MOTIF_SCORE_CHANGE": AnnotationEntry(
        "MOTIF_SCORE_CHANGE",
        float32,
        description="The difference in motif score of the reference and variant "
        "sequences for the TFBP",
    ),
    "CELL_TYPE": AnnotationListEntry(
        "CELL_TYPE",
        ",",
        description="List of cell types and classifications for regulatory feature",
    ),
    "CANONICAL": AnnotationEntry(
        "CANONICAL",
        typefunc="YES".__eq__,
        description="A flag indicating if the transcript is denoted as the canonical "
        "transcript for this gene",
    ),
    "CCDS": AnnotationEntry(
        "CCDS",
        description="The CCDS identifer for this transcript, where applicable",
    ),
    "INTRON": AnnotationRangeTotalEntry(
        "INTRON", description="The intron number (out of total number)"
    ),
    "EXON": AnnotationRangeTotalEntry(
        "EXON",
        description="The exon index range (out of total number of exons)",
    ),
    "DOMAINS": AnnotationListDictEntry(
        "DOMAINS",
        description="The source and identifer of any overlapping protein domains",
    ),
    "DISTANCE": AnnotationEntry(
        "DISTANCE",
        int,
        description="Shortest distance from variant to transcript",
    ),
    "IND": AnnotationEntry("IND", description="Individual name"),
    # "ZYG": AnnotationEntry(
    #     "ZYG", description="Zygosity of individual genotype at this locus"
    # ),
    # "SV": AnnotationEntry("SV", description="IDs of overlapping structural variants"),
    # "FREQS": AnnotationEntry(
    #     "FREQS", description="Frequencies of overlapping variants used in filtering"
    # ),
    "AF": AnnotationEntry(
        "AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes",
    ),
    "AFR_AF": AnnotationEntry(
        "AFR_AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes combined "
        "African population",
    ),
    "AMR_AF": AnnotationEntry(
        "AMR_AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes combined "
        "American population",
    ),
    "ASN_AF": AnnotationEntry(
        "ASN_AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes combined "
        "Asian population",
    ),
    "EUR_AF": AnnotationEntry(
        "EUR_AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes combined "
        "European population",
    ),
    "EAS_AF": AnnotationEntry(
        "EAS_AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes combined "
        "East Asian population",
    ),
    "SAS_AF": AnnotationEntry(
        "SAS_AF",
        float32,
        description="Frequency of existing variant in 1000 Genomes combined "
        "South Asian population",
    ),
    "AA_AF": AnnotationEntry(
        "AA_AF",
        float32,
        description="Frequency of existing variant in NHLBI-ESP "
        "African American population",
    ),
    "EA_AF": AnnotationEntry(
        "EA_AF",
        float32,
        description="Frequency of existing variant in NHLBI-ESP "
        "European American population",
    ),
    "gnomAD_AF": AnnotationEntry(
        "gnomAD_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "combined population",
    ),
    "gnomAD_AFR_AF": AnnotationEntry(
        "gnomAD_AFR_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "African/American population",
    ),
    "gnomAD_AMR_AF": AnnotationEntry(
        "gnomAD_AMR_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "American population",
    ),
    "gnomAD_ASJ_AF": AnnotationEntry(
        "gnomAD_ASJ_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "Ashkenazi Jewish population",
    ),
    "gnomAD_EAS_AF": AnnotationEntry(
        "gnomAD_EAS_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "East Asian population",
    ),
    "gnomAD_FIN_AF": AnnotationEntry(
        "gnomAD_FIN_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes Finnish population",
    ),
    "gnomAD_NFE_AF": AnnotationEntry(
        "gnomAD_NFE_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "Non-Finnish European population",
    ),
    "gnomAD_OTH_AF": AnnotationEntry(
        "gnomAD_OTH_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "combined other combined populations",
    ),
    "gnomAD_SAS_AF": AnnotationEntry(
        "gnomAD_SAS_AF",
        float32,
        description="Frequency of existing variant in gnomAD exomes "
        "South Asian population",
    ),
    "MAX_AF": AnnotationEntry(
        "MAX_AF",
        float32,
        description="Maximum observed allele frequency in 1000 Genomes, ESP and gnomAD",
    ),
    "MAX_AF_POPS": AnnotationListEntry(
        "MAX_AF_POPS",
        sep="&",
        description="Populations in which maximum allele frequency was observed",
    ),
    "CLIN_SIG": AnnotationListEntry(
        "CLIN_SIG",
        "&",
        description="ClinVar clinical significance of the dbSNP variant",
    ),
    "BIOTYPE": AnnotationEntry(
        "BIOTYPE",
        description="Biotype of transcript or regulatory feature",
    ),
    "APPRIS": AnnotationEntry(
        "APPRIS",
        description="Annotates alternatively spliced transcripts as primary or "
        "alternate based on a range of computational methods. "
        "NB: not available for GRCh37",
    ),
    "TSL": AnnotationEntry(
        "TSL",
        description="Transcript support level. NB: not available for GRCh37",
    ),
    "PUBMED": AnnotationListEntry(
        "PUBMED",
        description="Pubmed ID(s) of publications that cite existing variant",
        sep="&",
    ),
    "SOMATIC": AnnotationListEntry(
        "SOMATIC",
        description="Somatic status of existing variant(s); "
        "multiple values correspond to multiple values "
        "in the Existing_variation field",
        sep="&",
    ),
    "PHENO": AnnotationListEntry(
        "PHENO",
        sep="&",
        description="Indicates if existing variant is associated with a phenotype, "
        "disease or trait; "
        "multiple values correspond to multiple values "
        "in the Existing_variation field",
    ),
    # TODO list type with unconfirmed sep
    "GENE_PHENO": AnnotationListEntry(
        "GENE_PHENO",
        sep="&",
        description="Indicates if overlapped gene is associated with a phenotype, "
        "disease or trait",
    ),
    "ALLELE_NUM": AnnotationEntry(
        "ALLELE_NUM",
        int,
        description="Allele number from input; "
        "0 is reference, 1 is first alternate etc",
    ),
    # "MINIMISED": AnnotationEntry(
    #     "MINIMISED",
    #     description="Alleles in this variant have been converted to minimal "
    #                 "representation before consequence calculation",
    # ),
    # TODO parse flag:
    #  "indicates if this block of consequence data was picked by
    #  --flag_pick or --flag_pick_allele"
    # "PICK": AnnotationEntry(
    #     "PICK",
    #     description="Indicates if this block of consequence data was picked by "
    #                 "--flag_pick or --flag_pick_allele",
    # ),
    # "BAM_EDIT": AnnotationEntry(
    #     "BAM_EDIT", description="Indicates success or failure of edit using BAM file"
    # ),
    "GIVEN_REF": AnnotationEntry(
        "GIVEN_REF",
        description="Reference allele from input",
    ),
    "USED_REF": AnnotationEntry(
        "USED_REF",
        description="Reference allele as used to get consequences",
    ),
    # TODO parse flag:
    #  see https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#other_fields
    # "REFSEQ_MATCH": AnnotationEntry(
    #     "REFSEQ_MATCH",
    #     description="The RefSeq transcript match status; contains a number of flags "
    #                 "indicating whether this RefSeq transcript matches "
    #                 "the underlying reference sequence and/or an Ensembl transcript",
    # ),
    "OverlapBP": AnnotationEntry(
        "OverlapBP",
        int,
        description="Number of base pairs overlapping with the corresponding "
        "structural variation feature",
    ),
    "OverlapPC": AnnotationEntry(
        "OverlapPC",
        float32,
        description="Percentage of corresponding structural variation feature "
        "overlapped by the given input",
    ),
    # "CHECK_REF": AnnotationEntry(
    #     "CHECK_REF",
    #     description="Reports variants where the input reference does not match "
    #                 "the expected reference",
    # ),
    "AMBIGUITY": AnnotationEntry(
        "AMBIGUITY",
        description="IUPAC allele ambiguity code",
    ),
    "Amino_acids": AnnotationListEntry(
        "Amino_acids",
        description="Reference and variant amino acids",
        sep="/",
    ),
    "Codons": AnnotationListEntry(
        "Codons",
        description="Reference and variant codon sequence",
        sep="/",
    ),
    # TODO HGNC_ID description
    "HGNC_ID": AnnotationEntry("HGNC_ID"),
    "MANE": AnnotationEntry(
        "MANE",
        description="Matched Annotation from NCBI and EMBL-EBI (MANE).",
    ),
    "MANE_SELECT": AnnotationEntry(
        "MANE_SELECT",
        description="Matched Annotation from NCBI and EMBL-EBI (MANE) "
        "canonical transcript, indicated by the respective RefSeq NM ID."
        "For more info, see: "
        "https://www.ensembl.org/info/genome/genebuild/mane.html",
    ),
    "MANE_PLUS_CLINICAL": AnnotationEntry(
        "MANE_PLUS_CLINICAL",
        description="Matched Annotation from NCBI and EMBL-EBI (MANE) "
        "transcripts beyond MANE_SELECT, that are important clinically. "
        "Indicated by their RefSeq NM IDs."
        "For more info, see: "
        "https://www.ensembl.org/info/genome/genebuild/mane.html",
    ),
    "GO": AnnotationEntry("GO", description="Gene ontology (GO) terms."),
    "miRNA": AnnotationEntry(
        "miRNA",
        description="Determines where in the secondary structure of a miRNA "
        "a variant falls",
    ),
    "Existing_variation": AnnotationListEntry(
        "Existing_variation",
        sep="&",
        description="Identifier(s) of co-located known variants",
    ),
    "LoFtool": AnnotationEntry(
        "LoFtool",
        typefunc=float32,
        description="Provides a rank of genic intolerance and "
        "consequent susceptibility to disease based on the ratio of Loss-of-function "
        "(LoF) to synonymous mutations.",
    ),
    "REVEL": AnnotationEntry(
        "REVEL",
        typefunc=float32,
        description="Estimate of the pathogenicity of missense variants. "
        "Please cite the REVEL publication alongside the VEP if you use this "
        "resource: https://www.ncbi.nlm.nih.gov/pubmed/27666373",
    ),
    "ExACpLI": AnnotationEntry(
        "ExACpLI",
        typefunc=float32,
        description="Probability of a gene being loss-of-function intolerant (pLI). "
        "Please cite the respective ExAC publication alongside the VEP if you use this "
        "resource: https://www.ncbi.nlm.nih.gov/pubmed/27535533",
    ),
    "am_class": AnnotationEntry(
        "am_class",
        description="AlphaMissense classification of variant "
        "into one of three discrete categories: "
        "'Likely pathogenic', 'Likely benign', or 'ambiguous'"
        "Please cite the AlphaMissense publication alongside Ensembl VEP "
        "if you use this resource: https://doi.org/10.1126/science.adg7492",
    ),
    "am_pathogenicity": AnnotationEntry(
        "am_pathogenicity",
        typefunc=float32,
        description="AlphaMissense score of a variant being "
        "likely pathogenic, likely benign, or ambiguous."
        "Please cite the AlphaMissense publication alongside Ensembl VEP "
        "if you use this resource: https://doi.org/10.1126/science.adg7492",
    ),
    "SpliceAI_pred_DS_AG": AnnotationEntry(
        "SpliceAI_pred_DS_AG",
        typefunc=float32,
        description="Probability of variant being splice-altering acceptor gain. "
        "Please cite the respective ExAC publication alongside the VEP if you use this "
        "resource: https://www.ncbi.nlm.nih.gov/pubmed/30661751",
    ),
    "SpliceAI_pred_DS_AL": AnnotationEntry(
        "SpliceAI_pred_DS_AL",
        typefunc=float32,
        description="Probability of variant being splice-altering acceptor loss. "
        "Please cite the respective ExAC publication alongside the VEP if you use this "
        "resource: https://www.ncbi.nlm.nih.gov/pubmed/30661751",
    ),
    "SpliceAI_pred_DS_DG": AnnotationEntry(
        "SpliceAI_pred_DS_DG",
        typefunc=float32,
        description="Probability of variant being splice-altering donor gain. "
        "Please cite the respective ExAC publication alongside the VEP if you use this "
        "resource: https://www.ncbi.nlm.nih.gov/pubmed/30661751",
    ),
    "SpliceAI_pred_DS_DL": AnnotationEntry(
        "SpliceAI_pred_DS_DL",
        typefunc=float32,
        description="Probability of variant being splice-altering donor loss. "
        "Please cite the respective ExAC publication alongside the VEP if you use this "
        "resource: https://www.ncbi.nlm.nih.gov/pubmed/30661751",
    ),
}

# VEP now uses gnomADe_* and gnomADg_* instead of gnomAD_*.
# Therefore, we need to add these new entries to the KNOWN_ANN_TYPE_MAP_VEP
NEW_GNOMAD_ENTRIES = {}
for key, typer in KNOWN_ANN_TYPE_MAP_VEP.items():
    if key.startswith("gnomAD_"):
        sub_key = key.split("_", 1)[1]
        for t in ("e", "g"):
            name = f"gnomAD{t}_{sub_key}"
            description = typer.description()
            if t == "g":
                description = description.replace("gnomAD exomes", "gnomAD genomes")
            entry = AnnotationEntry(
                name=name,
                typefunc=typer._typefunc,
                description=description,
            )
            NEW_GNOMAD_ENTRIES[name] = entry
KNOWN_ANN_TYPE_MAP_VEP.update(NEW_GNOMAD_ENTRIES)
KNOWN_ANN_TYPE_MAP = {**KNOWN_ANN_TYPE_MAP_SNPEFF, **KNOWN_ANN_TYPE_MAP_VEP}

ANN_TYPER = AnnotationTyper(KNOWN_ANN_TYPE_MAP)
