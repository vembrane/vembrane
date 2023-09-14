from abc import abstractmethod, abstractproperty
from sys import stderr


class VCFRecord:
    @abstractmethod
    def __init__(self, filename: str):
        raise NotImplementedError

    @abstractproperty
    def info(self):
        raise NotImplementedError


class VCFHeader:
    @abstractmethod
    def __init__(self, filename: str):
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
