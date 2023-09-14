from abc import abstractmethod, abstractproperty
from typing import Any, Dict, List, Tuple


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
    def filter(self) -> List[str]:
        raise NotImplementedError

    @abstractproperty
    def info(self) -> Dict[str, Any]:
        raise NotImplementedError

    @abstractproperty
    def format(self) -> Dict[str, Dict[str, Any]]:
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

    @abstractmethod
    def add_meta(self, name: str, value: str):
        raise NotImplementedError

    def __iter__(self):
        return self

    @abstractmethod
    def __next__(self) -> VCFRecord:
        raise NotImplementedError

    @property
    def header(self):
        return self._header


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
