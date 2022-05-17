# Custom ANN types

Since arbitrary information may be encoded in the annotation field, `vembrane` has custom parsers for most of them.
The tables below give an overview of annotations that are recognized by vembrane; any unrecognized annotation is left unchanged (i.e. is of type `str`).
The descriptions are mostly taken from the `snpEff` and `vep` website/documentation, if available.

## vep

Annotations with custom types:

|Name in vep `ANN`|Name in `vembrane`|Type|Description|Example expression|
|---|---|---|---|---|
|`cDNA_position`|`cDNA`|`PosRange` with properties `start`, `end` and `length`| |`'ANN["cDNA"].start < 42'`|
|`CDS_position`|`CDS`|`PosRange` with properties `start`, `end` and `length`| |`'ANN["CDS"].end > 42'`|
|`Protein_position`|`AA`|`PosRange` with properties `start`, `end` and `length`| |`'ANN["AA"].length == 42'`|
|`STRAND`|`STRAND`|`int`|The DNA strand (1 or -1) on which the transcript/feature lies| |
|`FLAGS`|`FLAGS`|`List[str]`|Transcript quality flags: `cds_start_NF`: CDS 5' incomplete, `cds_end_NF`: CDS 3' incomplete| |
|`HGVS_OFFSET`|`HGVS_OFFSET`|`int`|Indicates by how many bases the HGVS notations for this variant have been shifted| |
|`SIFT`|`SIFT`|`Dict[str, float]`|The SIFT prediction and/or score, with both given as prediction(score)|`'ANN["SIFT"]["tolerated"] > 0.05'`|
|`PolyPhen`|`PolyPhen`|`Dict[str, float]`|The PolyPhen prediction and/or score|`'ANN["PolyPhen"]["probably_damaging"] > 0.9'`|
|`MOTIF_POS`|`MOTIF_POS`|`int`|The relative position of the variation in the aligned TFBP| |
|`HIGH_INF_POS`|`HIGH_INF_POS`|`bool`|A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)| |
|`MOTIF_SCORE_CHANGE`|`MOTIF_SCORE_CHANGE`|`float`|The difference in motif score of the reference and variant sequences for the TFBP| |
|`CELL_TYPE`|`CELL_TYPE`|`List[str]`|List of cell types and classifications for regulatory feature| |
|`CANONICAL`|`CANONICAL`|`bool`|A flag indicating if the transcript is denoted as the canonical transcript for this gene| |
|`INTRON`|`INTRON`|`NumberTotal` with properties `number` and `total`|The intron number (out of total number)|`'ANN["INTRON"].number == 2'`|
|`EXON`|`EXON`|`RangeTotal` with properties `range` and `total`|The known range of the exon indices (out of the total number of exons)|`'ANN["EXON"].total >= 2'`|
|`DOMAINS`|`DOMAINS`|`List[Dict[str, Any]]`|The source and identifer of any overlapping protein domains|`'ANN["DOMAINS"] is not NA and ("PTHR16515" in ANN["DOMAINS"]["PANTHER"])'`|
|`DISTANCE`|`DISTANCE`|`int`|Shortest distance from variant to transcript| |
|`AF`|`AF`|`float`|Frequency of existing variant in 1000 Genomes| |
|`AFR_AF`|`AFR_AF`|`float`|Frequency of existing variant in 1000 Genomes combined African population| |
|`AMR_AF`|`AMR_AF`|`float`|Frequency of existing variant in 1000 Genomes combined American population| |
|`ASN_AF`|`ASN_AF`|`float`|Frequency of existing variant in 1000 Genomes combined Asian population| |
|`EUR_AF`|`EUR_AF`|`float`|Frequency of existing variant in 1000 Genomes combined European population| |
|`EAS_AF`|`EAS_AF`|`float`|Frequency of existing variant in 1000 Genomes combined East Asian population| |
|`SAS_AF`|`SAS_AF`|`float`|Frequency of existing variant in 1000 Genomes combined South Asian population| |
|`AA_AF`|`AA_AF`|`float`|Frequency of existing variant in NHLBI-ESP African American population| |
|`EA_AF`|`EA_AF`|`float`|Frequency of existing variant in NHLBI-ESP European American population| |
|`gnomAD_AF`|`gnomAD_AF`|`float`|Frequency of existing variant in gnomAD exomes combined population| |
|`gnomAD_AFR_AF`|`gnomAD_AFR_AF`|`float`|Frequency of existing variant in gnomAD exomes African/American population| |
|`gnomAD_AMR_AF`|`gnomAD_AMR_AF`|`float`|Frequency of existing variant in gnomAD exomes American population| |
|`gnomAD_ASJ_AF`|`gnomAD_ASJ_AF`|`float`|Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population| |
|`gnomAD_EAS_AF`|`gnomAD_EAS_AF`|`float`|Frequency of existing variant in gnomAD exomes East Asian population| |
|`gnomAD_FIN_AF`|`gnomAD_FIN_AF`|`float`|Frequency of existing variant in gnomAD exomes Finnish population| |
|`gnomAD_NFE_AF`|`gnomAD_NFE_AF`|`float`|Frequency of existing variant in gnomAD exomes Non-Finnish European population| |
|`gnomAD_OTH_AF`|`gnomAD_OTH_AF`|`float`|Frequency of existing variant in gnomAD exomes combined other combined populations| |
|`gnomAD_SAS_AF`|`gnomAD_SAS_AF`|`float`|Frequency of existing variant in gnomAD exomes South Asian population| |
|`MAX_AF`|`MAX_AF`|`float`|Maximum observed allele frequency in 1000 Genomes, ESP and gnomAD| |
|`MAX_AF_POPS`|`MAX_AF_POPS`|`List[str]`|Populations in which maximum allele frequency was observed| |
|`CLIN_SIG`|`CLIN_SIG`|`List[str]`|ClinVar clinical significance of the dbSNP variant|`'"uncertain_significance" in ANN["CLIN_SIG"]'`|
|`PUBMED`|`PUBMED`|`List[str]`|Pubmed ID(s) of publications that cite existing variant| |
|`SOMATIC`|`SOMATIC`|`List[str]`|Somatic status of existing variant(s); multiple values correspond to multiple values in the Existing_variation field| |
|`PHENO`|`PHENO`|`List[str]`|Indicates if existing variant is associated with a phenotype, disease or trait; multiple values correspond to multiple values in the Existing_variation field| |
|`GENE_PHENO`|`GENE_PHENO`|`List[str]`|Indicates if overlapped gene is associated with a phenotype, disease or trait| |
|`ALLELE_NUM`|`ALLELE_NUM`|`int`|Allele number from input; 0 is reference, 1 is first alternate etc| |
|`OverlapBP`|`OverlapBP`|`int`|Number of base pairs overlapping with the corresponding structural variation feature| |
|`OverlapPC`|`OverlapPC`|`float`|Percentage of corresponding structural variation feature overlapped by the given input| |
|`Amino_acids`|`Amino_acids`|`List[str]`|Reference and variant amino acids| |
|`Codons`|`Codons`|`List[str]`|Reference and variant codon sequence| |
|`Existing_variation`|`Existing_variation`|`List[str]`|Identifier(s) of co-located known variants| |
|`LoFtool`|`LoFtool`|`float`|Provides a rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations."| |
|`REVEL`|`REVEL`|`float`|Estimate of the pathogenicity of missense variants.| |
|`ExACpLI`|`ExACpLI`|`float`|Probabililty of a gene being loss-of-function intolerant (pLI).| |


Annotations with type `str`:

|Name in vep `ANN`|Name in `vembrane`|Type|Description|Example expression|
|---|---|---|---|---|
|`Location`|`Location`|`str`|In standard coordinate format (chr:start or chr:start-end)| |
|`Allele`|`Allele`|`str`|The variant allele used to calculate the consequence| |
|`Gene`|`Gene`|`str`|Ensembl stable ID of affected gene| |
|`Feature`|`Feature`|`str`|Ensembl stable ID of feature| |
|`Feature_type`|`Feature_type`|`str`|Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature.| |
|`Consequence`|`Consequence`|`str`|Consequence type of this variant| |
|`HGSVc`|`HGSVc`|`str`| | |
|`HGSVp`|`HGSVp`|`str`| | |
|`HGVSc`|`HGVSc`|`str`|The HGVS coding sequence name| |
|`HGVSp`|`HGVSp`|`str`|The HGVS protein sequence name| |
|`HGVSg`|`HGVSg`|`str`|The HGVS genomic sequence name| |
|`REF_ALLELE`|`REF_ALLELE`|`str`|The reference allele| |
|`IMPACT`|`IMPACT`|`str`|The impact modifier for the consequence type|`'"HIGH" in ANN["IMPACT"]'`|
|`SYMBOL`|`SYMBOL`|`str`|The gene symbol| |
|`VARIANT_CLASS`|`VARIANT_CLASS`|`str`|Sequence Ontology variant class| |
|`SYMBOL_SOURCE`|`SYMBOL_SOURCE`|`str`|The source of the gene symbol| |
|`ENSP`|`ENSP`|`str`|The Ensembl protein identifier of the affected transcript| |
|`SWISSPROT`|`SWISSPROT`|`str`|Best match UniProtKB/Swiss-Prot accession of protein product| |
|`TREMBL`|`TREMBL`|`str`|Best match UniProtKB/TrEMBL accession of protein product| |
|`UNIPARC`|`UNIPARC`|`str`|Best match UniParc accession of protein product| |
|`MOTIF_NAME`|`MOTIF_NAME`|`str`|The source and identifier of a transcription factor binding profile aligned at this position| |
|`CCDS`|`CCDS`|`str`|The CCDS identifer for this transcript, where applicable| |
|`IND`|`IND`|`str`|Individual name| |
|`BIOTYPE`|`BIOTYPE`|`str`|Biotype of transcript or regulatory feature| |
|`APPRIS`|`APPRIS`|`str`|Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods. NB: not available for GRCh37| |
|`TSL`|`TSL`|`str`|Transcript support level. NB: not available for GRCh37| |
|`GIVEN_REF`|`GIVEN_REF`|`str`|Reference allele from input| |
|`USED_REF`|`USED_REF`|`str`|Reference allele as used to get consequences| |
|`AMBIGUITY`|`AMBIGUITY`|`str`|IUPAC allele ambiguity code| |
|`HGNC_ID`|`HGNC_ID`|`str`| | |
|`MANE`|`MANE`|`str`|Matched Annotation from NCBI and EMBL-EBI (MANE).| |
|`MANE_SELECT`|`MANE_SELECT`|`str`|Matched Annotation from NCBI and EMBL-EBI (MANE) canonical transcript.| |
|`MANE_PLUS_CLINICAL`|`MANE_PLUS_CLINICAL`|`str`|MANE transcripts beyond MANE_SELECT that are clinically relevant.| |
|`GO`|`GO`|`str`|Gene ontology (GO) terms.| |
|`miRNA`|`miRNA`|`str`|Determines where in the secondary structure of a miRNA a variant falls| |


## snpEff
Annotations with custom types:

|Name in snpEff `ANN`|Name in `vembrane`|Type|Description|Example expression|
|---|---|---|---|---|
|`cDNA.pos / cDNA.length`|`cDNA`|`PosRange` with properties `start`, `end` and `length`| |`vembrane filter 'ANN["cDNA"].start < 42'`|
|`CDS.pos / CDS.length`|`CDS`|`PosRange` with properties `start`, `end` and `length`| |`vembrane filter 'ANN["CDS"].end > 42'`|
|`AA.pos / AA.length`|`AA`|`PosRange` with properties `start`, `end` and `length`| |`vembrane filter 'ANN["AA"].length == 42'`|
|`ERRORS / WARNINGS / INFO`|`ERRORS / WARNINGS / INFO`|`List[str]`| | |


Annotations with type `str`:

|Name in snpEff `ANN`|Name in `vembrane`|Type|Description|Example expression|
|---|---|---|---|---|
|`Allele`|`Allele`|`str`| | |
|`Annotation`|`Annotation`|`str`| | |
|`Annotation_Impact`|`Annotation_Impact`|`str`| | |
|`Gene_Name`|`Gene_Name`|`str`| | |
|`Gene_ID`|`Gene_ID`|`str`| | |
|`Feature_Type`|`Feature_Type`|`str`| | |
|`Feature_ID`|`Feature_ID`|`str`| | |
|`Transcript_BioType`|`Transcript_BioType`|`str`| | |
|`Rank`|`Rank`|`str`| | |
|`HGVS.c`|`HGVS.c`|`str`| | |
|`HGVS.p`|`HGVS.p`|`str`| | |
|`Distance`|`Distance`|`str`| | |
