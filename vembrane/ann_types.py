from typing import Union, Iterable, Tuple, Dict, Callable, Any


class NoValue:
    def __lt__(self, other):
        return False

    def __gt__(self, other):
        return False

    def __le__(self, other):
        return False

    def __ge__(self, other):
        return False

    def __eq__(self, other):
        return False


NA = NoValue()


class InfoTuple:
    """A container that lazily evaluates None to NA in case of access."""

    def __init__(self, values):
        self.values = values

    def __getitem__(self, spec):
        values = self.values.__getitem__(spec)
        if isinstance(values, tuple) and any(v is None for v in values):
            return tuple((NA if v is None else v) for v in values)
        if values is None:
            return NA
        return values


IntFloatStr = Union[int, float, str]


def type_info(value):
    if value is None:
        return NA
    if isinstance(value, tuple):
        return InfoTuple(value)
    return value


class PosLength:
    def __init__(self, pos: int, length: int):
        self.pos = pos
        self.length = length

    @classmethod
    def from_snpeff_str(cls, value: str) -> "PosLength":
        pos, length = [int(v.strip()) for v in value.split("/")]
        return cls(pos, length)

    @classmethod
    def from_vep_str(cls, value: str) -> "PosLength":
        start, end = (
            #  the "-" is optional, so vep either has either start and end position
            [int(v.strip()) for v in value.split("-")]
            if "-" in value
            #  or start position only
            else [int(value.strip())] * 2
        )
        return cls(start, end - start)

    def __str__(self):
        return f"(pos: {self.pos}, length: {self.length})"

    def __repr__(self):
        return self.__str__()


def na_func() -> NoValue:
    return NA


AnnotationType = Union[IntFloatStr, Iterable["IntFloatStr"], NoValue]


class AnnotationEntry:
    def __init__(
        self,
        name: str,
        typefunc: Callable[[str], Any] = str,
        nafunc: Callable[[], Any] = na_func,
    ):
        self._name = name
        self._typefunc = typefunc
        self._nafunc = nafunc

    def convert(self, value: str) -> Tuple[str, Any]:
        if value:
            return self._name, self._typefunc(value)
        else:
            return self._name, self._nafunc()


class DefaultAnnotationEntry(AnnotationEntry):
    def __init__(self, name: str):
        # TODO issue INFO: No type conversion information about {name},
        #  please open an issue here: https://github.com/vembrane/vembrane/issues
        super().__init__(name)


class AnnotationTyper:
    def __init__(self, mapping: Dict[str, AnnotationEntry]):
        self._mapping = mapping

    def convert(self, key: str, value: str) -> Tuple[str, AnnotationType]:
        return self._mapping.get(key, DefaultAnnotationEntry(key)).convert(value)


KNOWN_ANN_TYPE_MAP_SNPEFF = {
    "Allele": AnnotationEntry("Allele"),
    "Annotation": AnnotationEntry("Annotation"),
    "Annotation_Impact": AnnotationEntry("Annotation_Impact"),
    "Gene_Name": AnnotationEntry("Gene_Name"),
    "Gene_ID": AnnotationEntry("Gene_ID"),
    "Feature_Type": AnnotationEntry("Feature_Type"),
    "Feature_ID": AnnotationEntry("Feature_ID"),
    "Transcript_BioType": AnnotationEntry("Transcript_BioType"),
    "Rank": AnnotationEntry("Rank"),
    "HGVS.c": AnnotationEntry("HGVS.c"),
    "HGVS.p": AnnotationEntry("HGVS.p"),
    "cDNA.pos / cDNA.length": AnnotationEntry("cDNA", PosLength.from_snpeff_str),
    "CDS.pos / CDS.length": AnnotationEntry("CDS", PosLength.from_snpeff_str),
    "AA.pos / AA.length": AnnotationEntry("AA", PosLength.from_snpeff_str),
    "Distance": AnnotationEntry("Distance", str),
    "ERRORS / WARNINGS / INFO": AnnotationEntry(
        "ERRORS / WARNINGS / INFO",
        lambda x: [v.strip() for v in x.split("/")],
        lambda: [],
    ),
}

KNOWN_ANN_TYPE_MAP_VEP = {
    "Allele": AnnotationEntry("Allele"),
    "Consequence": AnnotationEntry("Consequence"),
    "IMPACT": AnnotationEntry("IMPACT"),
    "SYMBOL": AnnotationEntry("SYMBOL"),
    "Gene": AnnotationEntry("Gene"),
    "Feature_type": AnnotationEntry("Feature_type"),
    "Feature": AnnotationEntry("Feature"),
    "EXON": AnnotationEntry("EXON"),
    "INTRON": AnnotationEntry("INTRON"),
    "HGSVc": AnnotationEntry("HGSVc"),
    "HGSVp": AnnotationEntry("HGSVp"),
    "cDNA_position": AnnotationEntry("cDNA", PosLength.from_vep_str),
    "CDS_position": AnnotationEntry("CDS", PosLength.from_vep_str),
    "Protein_position": AnnotationEntry("Protein", PosLength.from_vep_str),
    "CLIN_SIG": AnnotationEntry(
        "CLIN_SIG", lambda x: [v.strip() for v in x.split("&")], lambda: []
    ),
}

KNOWN_ANN_TYPE_MAP = {**KNOWN_ANN_TYPE_MAP_SNPEFF, **KNOWN_ANN_TYPE_MAP_VEP}

ANN_TYPER = AnnotationTyper(KNOWN_ANN_TYPE_MAP)
