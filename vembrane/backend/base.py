from abc import abstractmethod, abstractproperty
from enum import Enum
from typing import List, Optional, Tuple

from ..ann_types import NA
from ..errors import UnknownSampleError


class Backend(Enum):
    pysam = 0
    cyvcf2 = 1

    def __str__(self):
        return self.name.lower()

    def __repr__(self):
        return str(self)

    @staticmethod
    def from_string(s):
        try:
            return Backend[s]
        except KeyError:
            return s


class VCFRecordInfo:
    __slots__ = ()

    @abstractmethod
    def __contains__(self, item):
        raise NotImplementedError

    @abstractmethod
    def __getitem__(self, item):
        raise NotImplementedError

    @abstractmethod
    def __setitem__(self, key, value):
        raise NotImplementedError

    @abstractmethod
    def keys(self):
        raise NotImplementedError

    def get(self, key, default=None):
        if key in self:
            return self[key]
        return default


class VCFRecordSamples:
    pass


class NoValueDict:
    def __contains__(self, item):
        try:
            value = self[item]
        except KeyError:
            return False
        return value is not NA


class DefaultGet:
    def get(self, item, default=NA):
        v = self[item]
        if v is not NA:
            return v
        else:
            return default


class VCFRecordFormat(NoValueDict):
    @abstractmethod
    def __setitem__(self, key, value):
        raise NotImplementedError

    def get(self, sample, default=None):
        try:
            return self[sample]
        except UnknownSampleError:
            return default

    @abstractmethod
    def keys(self):
        raise NotImplementedError

    def __repr__(self):
        return str(
            {
                sample: self.__getitem__(sample)
                for i, sample in enumerate(self._header.samples)
            },
        )


class VCFRecordFilter:
    __slots__ = ()

    @abstractmethod
    def add(self, tag: str):
        raise NotImplementedError

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError


class VCFRecord:
    __slots__ = ("_raw_record", "record_idx", "_header")

    @abstractmethod
    def __init__(self, record, record_idx: int, header):
        self._raw_record = record
        self.record_idx = record_idx
        self._header = header

    @abstractproperty
    def contig(self) -> str:
        raise NotImplementedError

    def chrom(self) -> str:
        return self.contig

    @abstractproperty
    def position(self) -> int:
        raise NotImplementedError

    def start(self) -> int:
        return self.position

    @abstractproperty
    def stop(self) -> int:
        raise NotImplementedError

    @abstractproperty
    def id(self) -> str:
        raise NotImplementedError

    @property
    def alleles(self) -> Tuple[str, ...]:
        return self.reference_allele, *self.alt_alleles

    @abstractproperty
    def reference_allele(self) -> str:
        raise NotImplementedError

    @abstractproperty
    def alt_alleles(self) -> Tuple[str]:
        raise NotImplementedError

    @abstractproperty
    def quality(self) -> float:
        raise NotImplementedError

    @abstractproperty
    def filter(self) -> VCFRecordFilter:
        raise NotImplementedError

    @abstractproperty
    def info(self) -> VCFRecordInfo:
        raise NotImplementedError

    @abstractproperty
    def formats(self) -> "VCFRecordFormats":
        raise NotImplementedError

    @abstractproperty
    def header(self) -> "VCFHeader":
        raise NotImplementedError

    @property
    def is_bnd_record(self) -> bool:
        return "SVTYPE" in self.info and self.info.get("SVTYPE", None) == "BND"

    @property
    def is_sv_record(self) -> bool:
        is_sv = self.info.get("SVTYPE", "") in [
            "DEL",
            "INS",
            "DUP",
            "CNV",
            "INV",
            "INDEL",
            "BND",
            "TRA",
        ] or any(
            alt
            in [
                "<DEL>",
                "<DEL:ME>",
                "<INS>",
                "<INS:ME>",
                "<DUP>",
                "<DUP:TANDEM>",
                "<CNV>",
                "<CNV:TR>",
                "<INV>",
            ]
            for alt in self.alt_alleles
        )
        return is_sv

    @property
    def end(self):
        if self.is_bnd_record:
            return NA
        else:
            # record.stop is pysams unified approach to get the end position of
            # a variant. It either considers the alt allele or the END field,
            # depending on the record. Stop is 0-based, but for consistency
            # with POS we convert into 1-based. Since for 1-based coordinates,
            # the expectation is that END is inclusive instead of exclusive
            # (as it is with 0-based), we do not need to add 1.
            return self.stop

    @abstractmethod
    def __repr__(self):
        return self._raw_record.__str__()

    @abstractmethod
    def __str__(self):
        return self._raw_record.__str__()

    def __eq__(self, other: object):
        if not isinstance(other, VCFRecord):
            return NotImplemented
        return all(
            (
                self.contig == other.contig,
                self.id == other.id,
                self.alleles == other.alleles,
                self.position == other.position,
                self.quality == other.quality,
                set(self.filter) == set(other.filter),
                # because we might end up comparing NA == NA which is False by design,
                # we replace those with None here for the equality check,
                # i.e. None == None is True
                all(
                    (self.info.get(key) or None) == (other.info.get(key) or None)
                    for key in self.header.infos  # type: ignore
                ),
                all(
                    (self.formats.get(key) or None) == (other.formats.get(key) or None)
                    for key in self.header.formats  # type: ignore
                ),
            ),
        )


class VCFRecordFormats(NoValueDict):
    @abstractmethod
    def __init__(
        self,
        record: VCFRecord,
    ):
        raise NotImplementedError

    @abstractmethod
    def keys(self):
        raise NotImplementedError

    def get(self, key: str, default=None):
        if key not in self:
            return default


class VCFReader:
    __slots__ = ("_file",)

    @abstractmethod
    def __init__(self, filename: str):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()

    def __iter__(self):
        return self

    @property
    def file(self):
        return self._file

    @abstractmethod
    def __next__(self) -> VCFRecord:
        raise NotImplementedError

    @property
    def header(self):
        return self._header

    @abstractmethod
    def reset(self):
        raise NotImplementedError


class VCFHeader:
    @abstractmethod
    def __init__(self, reader: VCFReader):
        raise NotImplementedError

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError

    @abstractmethod
    def __next__(self):
        raise NotImplementedError

    @abstractproperty
    def samples(self) -> List[str]:
        raise NotImplementedError

    @abstractproperty
    def filters(self) -> List[str]:
        raise NotImplementedError

    @abstractproperty
    def records(self):
        raise NotImplementedError

    @abstractproperty
    def infos(self) -> VCFRecordInfo:
        raise NotImplementedError

    @abstractproperty
    def formats(self) -> VCFRecordFormats:
        raise NotImplementedError

    @abstractmethod
    def add_meta(
        self,
        key: str,
        value: Optional[str] = None,
        items: Optional[List[Tuple[str, str]]] = None,
    ):
        raise NotImplementedError

    @abstractmethod
    def add_generic(self, key: str, value: str):
        raise NotImplementedError

    @abstractmethod
    def add_filter(self, id: str, description: str):
        raise NotImplementedError

    @abstractmethod
    def contains_generic(self, key: str):
        raise NotImplementedError


class VCFWriter:
    __slots__ = ("_file",)

    @abstractmethod
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        raise NotImplementedError

    @abstractmethod
    def write(self, record: VCFRecord):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()
