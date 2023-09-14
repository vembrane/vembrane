import pysam
from vembrane.backend.base import VCFReader, VCFWriter, VCFHeader, VCFRecord
from collections import defaultdict, OrderedDict

from sys import stderr


class PysamVCFRecord(VCFRecord):
    def __init__(self, record):
        self._record = record

    @property
    def info(self):
        return self._record.info


class PysamVCFReader(VCFReader):
    def __init__(self, filename: str):
        self.filename = filename
        self._file = pysam.VariantFile(self.filename)
        self._header = PysamVCFHeader(self)

    def add_meta(self, name, value):
        self._file.header.add_meta(name, value)

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        return PysamVCFRecord(self._iter_file.__next__())


class PysamVCFHeader(VCFHeader):
    def __init__(self, reader: VCFReader):
        self._reader = reader
        self._header = reader._file.header
        self._data = []
        self._data_category = defaultdict(OrderedDict)
        for r in self._header.records:
            d = dict(r)
            d["key"] = r.key
            d["value"] = r.value
            d["type"] = r.type
            self._data.append(d)
            if not r.type in self._data_category:
                self._data_category[r.type] = dict()
            if "ID" in d:
                self._data_category[r.type][d["ID"]] = d

    @property
    def info(self):
        return self._data_category["INFO"]

    def __iter__(self):
        self._iter_index = 0
        return self

    def __next__(self):
        if self._iter_index >= len(self._data):
            raise StopIteration
        ret = self._data[self._iter_index]
        self._iter_index += 1
        return ret

    @property
    def samples(self):
        return self._header.samples


class PysamVCFWriter(VCFWriter):
    def __init__(self, filename: str, fmt: str, reader: VCFReader):
        self.filename = filename
        self._file = pysam.VariantFile(
            self.filename, f"w{fmt}", header=reader.header._header
        )
        self._header = PysamVCFHeader(self)

    def write(self, record: PysamVCFRecord):
        self._file.write(record._record)
