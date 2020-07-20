from sys import stderr
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


class PosRange:
    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end

    @classmethod
    def from_snpeff_str(cls, value: str) -> "PosRange":
        pos, length = [int(v.strip()) for v in value.split("/")]
        return cls(pos, pos + length)

    @classmethod
    def from_vep_str(cls, value: str) -> "PosRange":
        start, end = (
            #  the "-" is optional, so vep either has either start and end position
            [NA if v == "?" else int(v) for v in map(str.strip, value.split("-"))]
            if "-" in value
            #  or start position only
            else [int(value.strip())] * 2
        )
        return cls(start, end)

    @property
    def length(self):
        if self.end != NA and self.start != NA:
            return self.end - self.start
        else:
            return NA

    def __str__(self):
        return f"(start: {self.start}, end: {self.end}, length: {self.length})"

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
        description: str = None,
    ):
        self._name = name
        self._typefunc = typefunc
        self._nafunc = nafunc
        self._description = description

    @property
    def name(self):
        return self._name

    def convert(self, value: str) -> Tuple[str, Any]:
        if value:
            return self._typefunc(value)
        else:
            return self._nafunc()

    def description(self):
        return self._description


class AnnotationListEntry(AnnotationEntry):
    def __init__(
        self, name: str, sep: str, typefunc: Callable[[str], Any] = None, **kwargs
    ):
        typefunc = typefunc if typefunc else lambda v: [x.strip() for x in v.split(sep)]
        super().__init__(name, typefunc, nafunc=lambda: [], **kwargs)


class DefaultAnnotationEntry(AnnotationEntry):
    def __init__(self, name: str):
        super().__init__(name)


class AnnotationTyper:
    def __init__(self, mapping: Dict[str, AnnotationEntry]):
        self._mapping = mapping

    def get_entry(self, key: str) -> Tuple[str, AnnotationType]:
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

    def convert(self, key: str, value: str) -> Tuple[str, AnnotationType]:
        entry = self.get_entry(key)
        return entry.name, entry.convert(value)


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
    "cDNA.pos / cDNA.length": AnnotationEntry("cDNA", PosRange.from_snpeff_str),
    "CDS.pos / CDS.length": AnnotationEntry("CDS", PosRange.from_snpeff_str),
    "AA.pos / AA.length": AnnotationEntry("AA", PosRange.from_snpeff_str),
    "Distance": AnnotationEntry("Distance", str),
    "ERRORS / WARNINGS / INFO": AnnotationListEntry("ERRORS / WARNINGS / INFO", "/",),
}

# fields, types and description taken from:
# https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#other_fields
KNOWN_ANN_TYPE_MAP_VEP = {
    "Location": AnnotationEntry(
        "Location",
        description="In standard coordinate format (chr:start or chr:start-end)",
    ),
    "Allele": AnnotationEntry(
        "Allele", description="The variant allele used to calculate the consequence"
    ),
    "Gene": AnnotationEntry("Gene", description="Ensembl stable ID of affected gene"),
    "Feature": AnnotationEntry("Feature", description="Ensembl stable ID of feature"),
    "Feature_type": AnnotationEntry(
        "Feature_type",
        description="Type of feature. "
        "Currently one of Transcript, RegulatoryFeature, MotifFeature.",
    ),
    "Consequence": AnnotationEntry(
        "Consequence", description="Consequence type of this variant"
    ),
    "cDNA_position": AnnotationEntry("cDNA", PosRange.from_vep_str),
    "CDS_position": AnnotationEntry("CDS", PosRange.from_vep_str),
    "Protein_position": AnnotationEntry("Protein", PosRange.from_vep_str),
    "HGSVc": AnnotationEntry("HGSVc"),
    "HGSVp": AnnotationEntry("HGSVp"),
    "REF_ALLELE": AnnotationEntry("REF_ALLELE", description="The reference allele"),
    "IMPACT": AnnotationEntry(
        "IMPACT", description="The impact modifier for the consequence type"
    ),
    "SYMBOL": AnnotationEntry("SYMBOL", description="The gene symbol"),
    "VARIANT_CLASS": AnnotationEntry(
        "VARIANT_CLASS", description="Sequence Ontology variant class"
    ),
    "SYMBOL_SOURCE": AnnotationEntry(
        "SYMBOL_SOURCE", description="The source of the gene symbol"
    ),
    "STRAND": AnnotationEntry(
        "STRAND",
        int,
        description="The DNA strand (1 or -1) on which the transcript/feature lies",
    ),
    "ENSP": AnnotationEntry(
        "ENSP", description="The Ensembl protein identifier of the affected transcript"
    ),
    # TODO parse flags: "cds_start_NF: CDS 5' incomplete, cds_end_NF: CDS 3' incomplete"
    # "FLAGS": AnnotationEntry(
    #     "FLAGS",
    #     description="Transcript quality flags:
    #     cds_start_NF: CDS 5' incomplete, cds_end_NF: CDS 3' incomplete",
    # ),
    "SWISSPROT": AnnotationEntry(
        "SWISSPROT",
        description="Best match UniProtKB/Swiss-Prot accession of protein product",
    ),
    "TREMBL": AnnotationEntry(
        "TREMBL", description="Best match UniProtKB/TrEMBL accession of protein product"
    ),
    "UNIPARC": AnnotationEntry(
        "UNIPARC", description="Best match UniParc accession of protein product"
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
    # TODO custom type:
    #  "the SIFT prediction and/or score, with both given as prediction(score)"
    # "SIFT": AnnotationEntry(
    #     "SIFT",
    #     description="The SIFT prediction and/or score,"
    #     " with both given as prediction(score)",
    # ),
    # TODO custom type:
    #  "the PolyPhen prediction and/or score"
    # "PolyPhen": AnnotationEntry(
    #     "PolyPhen", description="The PolyPhen prediction and/or score"
    # ),
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
    # TODO parse flag:
    #  "a flag indicating if the variant falls in a high information position of a
    #  transcription factor binding profile (TFBP)"
    # "HIGH_INF_POS": AnnotationEntry(
    #     "HIGH_INF_POS",
    #     description="A flag indicating if the variant falls in a high information "
    #                 "position of a transcription factor binding profile (TFBP)",
    # ),
    "MOTIF_SCORE_CHANGE": AnnotationEntry(
        "MOTIF_SCORE_CHANGE",
        float,
        description="The difference in motif score of the reference and variant "
        "sequences for the TFBP",
    ),
    "CELL_TYPE": AnnotationListEntry(
        "CELL_TYPE",
        ",",
        description="List of cell types and classifications for regulatory feature",
    ),
    # TODO parse flag:
    #  "a flag indicating if the transcript is denoted as the canonical transcript
    #  for this gene"
    # "CANONICAL": AnnotationEntry(
    #     "CANONICAL",
    #     description="A flag indicating if the transcript is denoted as the canonical "
    #                 "transcript for this gene",
    # ),
    "CCDS": AnnotationEntry(
        "CCDS", description="The CCDS identifer for this transcript, where applicable"
    ),
    # TODO custom type: "the intron number (out of total number)"
    # "INTRON": AnnotationEntry(
    #     "INTRON", description="The intron number (out of total number)"
    # ),
    # TODO custom type: "the exon number (out of total number)"
    # "EXON": AnnotationEntry(
    #     "EXON", description="The exon number (out of total number)"
    # ),
    # "DOMAINS": AnnotationEntry(
    #     "DOMAINS",
    #     description="The source and identifer of any overlapping protein domains",
    # ),
    "DISTANCE": AnnotationEntry(
        "DISTANCE", int, description="Shortest distance from variant to transcript"
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
        "AF", float, description="Frequency of existing variant in 1000 Genomes"
    ),
    "AFR_AF": AnnotationEntry(
        "AFR_AF",
        float,
        description="Frequency of existing variant in 1000 Genomes combined "
        "African population",
    ),
    "AMR_AF": AnnotationEntry(
        "AMR_AF",
        float,
        description="Frequency of existing variant in 1000 Genomes combined "
        "American population",
    ),
    "ASN_AF": AnnotationEntry(
        "ASN_AF",
        float,
        description="Frequency of existing variant in 1000 Genomes combined "
        "Asian population",
    ),
    "EUR_AF": AnnotationEntry(
        "EUR_AF",
        float,
        description="Frequency of existing variant in 1000 Genomes combined "
        "European population",
    ),
    "EAS_AF": AnnotationEntry(
        "EAS_AF",
        float,
        description="Frequency of existing variant in 1000 Genomes combined "
        "East Asian population",
    ),
    "SAS_AF": AnnotationEntry(
        "SAS_AF",
        float,
        description="Frequency of existing variant in 1000 Genomes combined "
        "South Asian population",
    ),
    "AA_AF": AnnotationEntry(
        "AA_AF",
        float,
        description="Frequency of existing variant in NHLBI-ESP "
        "African American population",
    ),
    "EA_AF": AnnotationEntry(
        "EA_AF",
        float,
        description="Frequency of existing variant in NHLBI-ESP "
        "European American population",
    ),
    "gnomAD_AF": AnnotationEntry(
        "gnomAD_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "combined population",
    ),
    "gnomAD_AFR_AF": AnnotationEntry(
        "gnomAD_AFR_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "African/American population",
    ),
    "gnomAD_AMR_AF": AnnotationEntry(
        "gnomAD_AMR_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "American population",
    ),
    "gnomAD_ASJ_AF": AnnotationEntry(
        "gnomAD_ASJ_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "Ashkenazi Jewish population",
    ),
    "gnomAD_EAS_AF": AnnotationEntry(
        "gnomAD_EAS_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "East Asian population",
    ),
    "gnomAD_FIN_AF": AnnotationEntry(
        "gnomAD_FIN_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "Finnish population",
    ),
    "gnomAD_NFE_AF": AnnotationEntry(
        "gnomAD_NFE_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "Non-Finnish European population",
    ),
    "gnomAD_OTH_AF": AnnotationEntry(
        "gnomAD_OTH_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "combined other combined populations",
    ),
    "gnomAD_SAS_AF": AnnotationEntry(
        "gnomAD_SAS_AF",
        float,
        description="Frequency of existing variant in gnomAD exomes "
        "South Asian population",
    ),
    "MAX_AF": AnnotationEntry(
        "MAX_AF",
        float,
        description="Maximum observed allele frequency in "
        "1000 Genomes, ESP and gnomAD",
    ),
    # "MAX_AF_POPS": AnnotationEntry(
    #     "MAX_AF_POPS",
    #     description="Populations in which maximum allele frequency was observed",
    # ),
    "CLIN_SIG": AnnotationListEntry(
        "CLIN_SIG",
        "&",
        description="ClinVar clinical significance of the dbSNP variant",
    ),
    "BIOTYPE": AnnotationEntry(
        "BIOTYPE", description="Biotype of transcript or regulatory feature"
    ),
    "APPRIS": AnnotationEntry(
        "APPRIS",
        description="Annotates alternatively spliced transcripts as primary or "
        "alternate based on a range of computational methods. "
        "NB: not available for GRCh37",
    ),
    "TSL": AnnotationEntry(
        "TSL", description="Transcript support level. " "NB: not available for GRCh37"
    ),
    "PUBMED": AnnotationEntry(
        "PUBMED", description="Pubmed ID(s) of publications that cite existing variant"
    ),
    # TODO list type with unknown sep
    # "SOMATIC": AnnotationEntry(
    #     "SOMATIC",
    #     description="Somatic status of existing variant(s); "
    #                 "multiple values correspond to multiple values "
    #                 "in the Existing_variation field",
    # ),
    # TODO list type with unknown sep
    # "PHENO": AnnotationEntry(
    #     "PHENO",
    #     description="Indicates if existing variant is associated with a phenotype, "
    #                 "disease or trait; "
    #                 "multiple values correspond to multiple values "
    #                 "in the Existing_variation field",
    # ),
    # TODO list type with unknown sep
    # "GENE_PHENO": AnnotationEntry(
    #     "GENE_PHENO",
    #     description="Indicates if overlapped gene is associated with a phenotype, "
    #                 "disease or trait",
    # ),
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
        "GIVEN_REF", description="Reference allele from input"
    ),
    "USED_REF": AnnotationEntry(
        "USED_REF", description="Reference allele as used to get consequences"
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
        float,
        description="Percentage of corresponding structural variation feature "
        "overlapped by the given input",
    ),
    # "CHECK_REF": AnnotationEntry(
    #     "CHECK_REF",
    #     description="Reports variants where the input reference does not match "
    #                 "the expected reference",
    # ),
    "AMBIGUITY": AnnotationEntry(
        "AMBIGUITY", description="IUPAC allele ambiguity code"
    ),
}

KNOWN_ANN_TYPE_MAP = {**KNOWN_ANN_TYPE_MAP_SNPEFF, **KNOWN_ANN_TYPE_MAP_VEP}

ANN_TYPER = AnnotationTyper(KNOWN_ANN_TYPE_MAP)
