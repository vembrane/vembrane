import os
import contextlib
import csv
import sys
from sys import stderr
from typing import Any, Dict, Iterator, List, Set

import asttokens
from pysam.libcbcf import VariantFile, VariantRecord

from ..common import AppendKeyValuePair, check_expression, read_auxiliary
from ..errors import HeaderWrongColumnNumber, VembraneError
from ..globals import allowed_globals
from ..representations import Environment
from .filter import DeprecatedAction

from sqlalchemy import create_engine, ForeignKey
from sqlalchemy import Column
from sqlalchemy import Table
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy import Float

from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship


def add_subcommmand(subparsers):
    parser = subparsers.add_parser("database")
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="Output file, if not specified, output is written to STDOUT.",
    )


consequences = {
    key.strip(): value
    for value, key in enumerate(open("consequences.txt", "r").readlines())
}


def consequence_bitvector(consequence):
    return sum(2 ** consequences[c] for c in consequence.split("&"))


impacts = {
    "MODIFIER": 0,
    "LOW": 1,
    "MODERATE": 2,
    "HIGH": 3,
}


def execute(args):
    if os.path.exists("test.sqlite"):
        os.remove("test.sqlite")
    engine = create_engine("sqlite:///test.sqlite")

    db_session = scoped_session(sessionmaker(bind=engine, autoflush=True))

    Base = declarative_base()
    with VariantFile(args.vcf) as vcf:

        class Field(Base):
            __tablename__ = "fields"
            name = Column(String, primary_key=True)
            table = Column(String, primary_key=True)
            type = Column(String)
            description = Column(String)

        class Sample_Has_Variant(Base):
            __tablename__ = "sample_has_variant"
            variant_id = Column(Integer, ForeignKey("variants.id"), primary_key=True)
            sample_id = Column(Integer, ForeignKey("samples.id"), primary_key=True)
            genotype = Column(Integer, index=True)
            dp = Column(Integer, index=True)
            sample = relationship("Sample", back_populates="variants")
            variant = relationship("Variant", back_populates="samples")

        for format in vcf.header.formats.values():
            if format.name in ["DP", "Genotype"]:
                continue
            setattr(
                Sample_Has_Variant,
                format.name,
                Column(getattr(sys.modules[__name__], format.type)),
            )

        class Variant(Base):
            __tablename__ = "variants"
            id = Column(Integer, primary_key=True)
            chr = Column(String, index=True)
            pos = Column(Integer, index=True)
            ref = Column(String, index=True)
            alt = Column(String, index=True)
            qual = Column(Float, index=True)
            samples = relationship("Sample_Has_Variant", back_populates="variant")
            ann = relationship("Annotation", backref="variant")
            # for info in vcf.header.info.values():
            # mq = Column(Integer, index=True)
            # ac = Column(Integer, index=True)
            # af = Column(Float, index=True)
            # an = Column(Integer, index=True)
            # nhomalt = Column(Integer, index=True)
            # consequences = Column(Integer)
            # impacts = Column(Integer)

        # add variables (info fields) dynamically to the variant class
        for info in vcf.header.info.values():
            setattr(
                Variant, info.name, Column(getattr(sys.modules[__name__], info.type))
            )

        class Annotation(Base):
            __tablename__ = "annotations"
            id = Column(Integer, primary_key=True)
            variant_id = Column(Integer, ForeignKey("variants.id"))
            transcript = Column(Integer, index=True)
            gene = Column(Integer, index=True)
            consequence = Column(String, index=True)

        # this is the vep annotation TODO: support SNPEff annotation
        vep = (
            ("impact", "String"),
            ("symbol", "String"),
            ("feature_type", "String"),
            ("biotype", "String"),
            ("exon", "Integer"),
            ("intron", "Integer"),
            ("hgvsc", "String"),
            ("hgvsp", "String"),
            ("cdna_position", "Integer"),
            ("cds_position", "Integer"),
            ("protein_position", "Integer"),
            ("amino_acids", "String"),
            ("codons", "String"),
            ("existing_variation", "String"),
            ("distance", "Integer"),
            ("strand", "Integer"),
            ("flags", "String"),
            ("symbol_source", "String"),
            ("hgnc_id", "Integer"),
        )

        for annotation_name, annotation_type in vep:
            setattr(
                Annotation,
                annotation_name,
                Column(getattr(sys.modules[__name__], annotation_type)),
            )

        class Sample(Base):
            __tablename__ = "samples"
            id = Column(Integer, primary_key=True)
            name = Column(String, index=True)
            variants = relationship("Sample_Has_Variant", back_populates="sample")

        Base.metadata.create_all(bind=engine)

        ann_format = list(
            map(
                str.lower,
                vcf.header.info["CSQ"].description.rsplit(" ", 1)[-1].split("|"),
            )
        )

        print("insert info fields", file=sys.stderr)
        objects = [
            Field(
                name=info.name,
                table="variants",
                type=info.type,
                description=info.description,
            )
            for info in vcf.header.info.values()
        ]
        db_session.bulk_save_objects(objects)

        print("insert annotation fields", file=sys.stderr)
        objects = [
            Field(
                name=annotation_name,
                table="annotations",
                type=annotation_type,
                description="",
            )
            for annotation_name, annotation_type in vep
        ]
        db_session.bulk_save_objects(objects)

        print("insert format fields", file=sys.stderr)
        objects = [
            Field(
                name=format.name,
                table="sample_has_variant",
                type=format.type,
                description=format.description,
            )
            for format in vcf.header.formats.values()
            if format.name not in ["DP", "Genotype"]
        ]
        db_session.bulk_save_objects(objects)

        print("insert samples", file=sys.stderr)
        objects = [
            Sample(id=sample_id, name=str(s))
            for sample_id, s in enumerate(vcf.header.samples)
        ]
        db_session.bulk_save_objects(objects)

        print("insert variants", file=sys.stderr)
        samples = [str(s) for s in vcf.header.samples]
        objects = []
        annotation_id = 0
        for variant_id, variant in enumerate(vcf):
            # if variant_id == 25000:
            #     break

            # create variant
            info = variant.info
            if variant_id % 10000 == 0:
                print(variant_id)
                db_session.bulk_save_objects(objects)
                objects = []

            # create relation of sample - variant
            for sample_id, s in enumerate(samples):
                format = variant.samples[s]
                gt = format["GT"]
                if gt == (None,):  # dont write sample, if it doesn't own the variant
                    continue
                genotype = sum(2**i * x for i, x in enumerate(gt))
                if genotype == 0:  # dont write sample, if it doesn't own the variant
                    continue
                objects.append(
                    Sample_Has_Variant(
                        sample_id=sample_id,
                        variant_id=variant_id,
                        genotype=genotype,
                        dp=format.get("DP", None) or format.get("DPI", None),
                    ),
                )

            all_consequences = 0
            all_impacts = 0

            # add annotations
            for a in info["CSQ"]:
                ann = {key: value for key, value in zip(ann_format, a.split("|"))}
                consequence_bits = consequence_bitvector(ann["consequence"])
                all_consequences |= consequence_bits
                impact_bit = 2 ** impacts[ann["impact"]]
                all_impacts |= impact_bit
                hgnc_id = (
                    int(ann["hgnc_id"].removeprefix("HGNC:"))
                    if ann["hgnc_id"]
                    else None
                )

                if ann["feature"]:
                    transcript = int(ann["feature"].removeprefix("ENST"))
                else:
                    transcript = None

                if ann["gene"]:
                    gene = int(ann["gene"].removeprefix("ENSG"))
                else:
                    gene = None
                # assert ann["feature_type"]=="Transcript"
                # print(ann["biotype"])
                objects.append(
                    Annotation(
                        id=annotation_id,
                        variant_id=variant_id,
                        transcript=transcript,
                        consequence=consequence_bits,
                        impact=ann["impact"],
                        symbol=ann["symbol"],
                        gene=gene,
                        # feature_type=ann["feature_type"],
                        biotype=ann["biotype"],
                        exon=ann["exon"],
                        intron=ann["intron"],
                        hgvsc=ann["hgvsc"],
                        hgvsp=ann["hgvsp"],
                        cdna_position=ann["cdna_position"],
                        cds_position=ann["cds_position"],
                        protein_position=ann["protein_position"],
                        amino_acids=ann["amino_acids"],
                        codons=ann["codons"],
                        existing_variation=ann["existing_variation"],
                        distance=ann["distance"],
                        strand=ann["strand"],
                        flags=ann["flags"],
                        symbol_source=ann["symbol_source"],
                        hgnc_id=hgnc_id,
                    ),
                )
                annotation_id += 1

            # add variant
            objects.append(
                Variant(
                    id=variant_id,
                    chr=variant.chrom,
                    pos=variant.pos,
                    ref=variant.alts[0],
                    alt=variant.ref,
                    qual=variant.qual,
                    # consequences=all_consequences,
                    # impacts=all_impacts,
                    MQ=info.get("MQ", None),
                    AC=info.get("AC", (None,))[0],
                    AF=info.get("AF", (None,))[0],
                    AN=info.get("AN", None),
                    nhomalt=info.get("nhomalt", (None,))[0],
                ),
            )
        db_session.bulk_save_objects(objects)

    db_session.commit()
    db_session.flush()
    db_session.close()
