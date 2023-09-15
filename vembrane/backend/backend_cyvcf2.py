import cyvcf2
from cyvcf2.cyvcf2 import Variant
from typing import Tuple
from sys import stderr

from collections import OrderedDict, defaultdict

from vembrane.backend.base import (
    VCFHeader,
    VCFReader,
    VCFRecord,
    VCFRecordFilter,
    VCFRecordFormat,
    VCFRecordInfo,
    VCFWriter,
)

class Cyvcf2VCFRecord(VCFRecord):
    def __init__(self, record: Variant):
        self._record = record

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()

    def __eq__(self, other):
        return self._record.__repr__() == other._record.__repr__() # TODO: implement real check for values

class Cyvcf2VCFReader(VCFReader):
    def __init__(self, filename: str):
        self.filename = filename
        self._file = cyvcf2.VCF(self.filename)
        self._header = Cyvcf2VCFHeader(self)

    # @abstractmethod
    # def add_meta(self, key: str, value: str, items: Tuple[str, str]):
    #     raise NotImplementedError

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        return Cyvcf2VCFRecord(self._iter_file.__next__())

    # @property
    # def header(self):
    #     return self._header

    # @abstractmethod
    # def reset(self):
        # pass

class Cyvcf2VCFHeader(VCFHeader):
    def __init__(self, reader: Cyvcf2VCFReader):
        self._reader = reader
        self._data = []
        self._data_category = defaultdict(OrderedDict)
        
        # print(reader._file.get_header_type("vembraneVersion"))
        # print(dir(reader._file), file=stderr)
        # print(reader._file.hdr, file=stderr)

        for r in reader._file.header_iter():
            print(reader._file.get_header_type("fileformat", order=""), r.type, r.__dir__(), file=stderr)
            d = r.info()
            d["key"] = r.key
            d["value"] = r.value
            d["type"] = r.type
            self._data.append(d)
            if r.type not in self._data_category:
                self._data_category[r.type] = dict()
            if "ID" in d:
                self._data_category[r.type][d["ID"]] = d

        @abstractproperty
        def records(self):
            raise NotImplementedError

    

class Cyvcf2VCFWriter(VCFWriter):
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        pass
        # self.filename = filename
        # self._file = pysam.VariantFile(
        #     self.filename, f"w{fmt}", header=template.header._header
        # )
        # self._header = Cyvcf2VCFHeader(self)

    def write(self, record: Cyvcf2VCFRecord):
        self._file.write(record._record)


# class PysamVCFRecord(VCFRecord):
#     def __init__(self, record: VariantRecord):
#         self._record = record