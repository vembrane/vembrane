from abc import abstractmethod, abstractproperty
from enum import Enum
from typing import List, Optional, Tuple

from ..ann_types import NA


class Backend(Enum):
    pysam = (0,)
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
    def get(self, item, default=None):
        raise NotImplementedError


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


class VCFRecordFormat(NoValueDict, DefaultGet):
    @abstractmethod
    def __setitem__(self, key, value):
        raise NotImplementedError

    def __repr__(self):
        return str(
            {
                sample: self.__getitem__(sample)
                for i, sample in enumerate(self._header.samples)
            }
        )


class VCFRecordFilter:
    @abstractmethod
    def add(self, tag: str):
        raise NotImplementedError

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError


class VCFRecord:
    @abstractmethod
    def __init__(self, filename: str):
        raise NotImplementedError

    @abstractproperty
    def contig(self) -> str:
        raise NotImplementedError

    @abstractproperty
    def position(self) -> int:
        raise NotImplementedError

    @abstractproperty
    def stop(self) -> int:
        raise NotImplementedError

    @abstractproperty
    def id(self) -> str:
        raise NotImplementedError

    @property
    def alleles(self) -> Tuple[str]:
        return (self.reference_allele, *self.alt_alleles)

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
    def format(self) -> VCFRecordFormat:
        raise NotImplementedError

    @abstractproperty
    def formats(self) -> VCFRecordFormat:
        raise NotImplementedError

    @abstractmethod
    def __repr__(self):
        raise NotImplementedError

    @abstractmethod
    def __str__(self):
        raise NotImplementedError

    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError


class VCFRecordFormats(NoValueDict):
    @abstractmethod
    def __init__(
        self,
        record: VCFRecord,
    ):
        raise NotImplementedError


class VCFReader:
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

    @abstractmethod
    def __next__(self) -> VCFRecord:
        raise NotImplementedError

    @property
    def header(self):
        return self._header

    @abstractproperty
    def records(self):
        raise NotImplementedError

    @abstractmethod
    def reset(self):
        pass


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


class VCFWriter:
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


# class VCFWriter:
#     @abstractmethod
#     def __init__(self, filename: str):
#         raise NotImplementedError

#     def __exit__(self, exc_type, exc_val, exc_tb):
#         print("exit", file=stderr)
#         self._file.close()
