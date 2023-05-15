import os
import contextlib
import csv
import sys
from sys import stderr
from typing import Any, Dict, Iterator, List, Set

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

from vembrane.ann_types import KNOWN_ANN_TYPE_MAP_VEP
from vembrane.common import get_annotation_keys

def add_subcommmand(subparsers):
    parser = subparsers.add_parser("database")
    parser.register("action", "deprecated", DeprecatedAction)
    parser.add_argument(
        "vcf", help="The file containing the variants.", nargs="?", default="-"
    )
    parser.add_argument(
        "output",
        help="Output database file.",
    )
    parser.add_argument(
        "--annotation-key",
        "-k",
        metavar="FIELDNAME",
        default="ANN",
        help="The INFO key for the annotation field.",
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
    if os.path.exists(args.output):
        os.remove(args.output)
    engine = create_engine(f"sqlite:///{args.output}")

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

        for key in get_annotation_keys(vcf.header, args.annotation_key):
            print(key, KNOWN_ANN_TYPE_MAP_VEP[key].database_type, KNOWN_ANN_TYPE_MAP_VEP[key].database_name)

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
                vcf.header.info[args.annotation_key]
                .description.rsplit(" ", 1)[-1]
                .split("|"),
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

        def detuple(x, k=None):
            if isinstance(x, tuple):
                if len(x) == 1:
                    return x[0]
                if len(x) == 2:
                    return x[1]
                else:
                    # print("unexpected tuple size:", x, k)
                    return str(x)
            return x

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

                formats = {
                    k: detuple(v, k)
                    for k, v in format.items()
                    if k not in ["GT", "DP", "DPI"]
                }
                objects.append(
                    Sample_Has_Variant(
                        sample_id=sample_id,
                        variant_id=variant_id,
                        genotype=genotype,
                        dp=format.get("DP", None) or format.get("DPI", None),
                        **formats,
                    ),
                )

            all_consequences = 0
            all_impacts = 0

            # add annotations
            for a in info[args.annotation_key]:
                ann = {key: value for key, value in zip(ann_format, a.split("|"))}
                consequence_bits = consequence_bitvector(ann["consequence"])
                all_consequences |= consequence_bits
                impact_bit = 2 ** impacts[ann["impact"]]
                all_impacts |= impact_bit
                # hgnc_id = (
                #     int(ann["hgnc_id"].removeprefix("HGNC:"))
                #     if ann["hgnc_id"]
                #     else None
                # )

                if ann["feature"]:
                    transcript = int(ann["feature"].removeprefix("ENST"))
                else:
                    transcript = None

                del ann["allele"]
                del ann["feature"]
                objects.append(
                    Annotation(
                        id=annotation_id,
                        variant_id=variant_id,
                        transcript=transcript,
                        **ann,
                    ),
                )
                annotation_id += 1

            # get infos but exclude the annotations
            infos = {
                k: detuple(v) for k, v in info.items() if k not in [args.annotation_key]
            }

            # add variant
            objects.append(
                Variant(
                    id=variant_id,
                    chr=variant.chrom,
                    pos=variant.pos,
                    ref=variant.alts[0],
                    alt=variant.ref,
                    qual=variant.qual,
                    **infos,
                ),
            )
        db_session.bulk_save_objects(objects)

    db_session.commit()
    db_session.flush()
    db_session.close()
