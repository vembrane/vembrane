from collections import OrderedDict, defaultdict
from typing import Tuple

import pysam
from pysam import VariantRecord

from vembrane.backend.base import (
    VCFHeader,
    VCFReader,
    VCFRecord,
    VCFRecordFilter,
    VCFRecordFormat,
    VCFRecordFormats,
    VCFRecordInfo,
    VCFWriter,
)

from ..errors import UnknownSample


class PysamRecord(VCFRecord):
    def __init__(self, record: VariantRecord):
        self._record = record

    @property
    def contig(self) -> str:
        return self._record.chrom

    @property
    def position(self) -> int:
        return self._record.pos

    @property
    def stop(self) -> int:
        return self._record.stop

    @property
    def id(self) -> str:
        return self._record.id

    @property
    def reference_allele(self) -> str:
        return self._record.ref

    @property
    def alt_alleles(self) -> Tuple[str]:
        return tuple(self._record.alts)

    @property
    def quality(self) -> float:
        return self._record.qual

    @property
    def filter(self) -> VCFRecordFilter:
        return PysamRecordFilter(self._record)

    @property
    def info(self) -> VCFRecordInfo:
        return PysamRecordInfo(self._record)

    @property
    def format(self) -> VCFRecordFormat:
        return PysamRecordFormat(self._record)

    @property
    def formats(self) -> VCFRecordFormats:
        return PysamRecordFormats(self._record)

    @property
    def samples(self):
        return self._record.samples

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()

    def __eq__(self, other):
        return self._record == other._record


class PysamRecordFormats(VCFRecordFormats):
    def __init__(self, record: VariantRecord):
        self._record = record

    def __getitem__(self, key):
        return PysamRecordFormat(key, self._record)


class PysamRecordFormat(VCFRecordFormat):
    def __init__(
        self,
        format_key: str,
        record: VariantRecord,
    ):
        self._format_key = format_key
        self._record = record

    def __getitem__(self, sample):
        try:
            self._record.samples.__contains__(sample)
        except KeyError:
            raise UnknownSample(self._record_idx, self._record, sample)
        return self._record.samples[sample][self._format_key]

    def __setitem__(self, key, value):
        self._record.format[key] = value

    def __contains__(self, item):
        return item in self._record.format

    def get(self, item, default=None):
        return self._record.format.get(item, default)


class PysamRecordInfo(VCFRecordInfo):
    def __init__(
        self,
        record: VariantRecord,
    ):
        self._record = record

    def __getitem__(self, item):
        return self._record.info[item]

    def __setitem__(self, key, value):
        self._record.info[key] = value

    def __contains__(self, item):
        return item in self._record.info

    def get(self, item, default=None):
        return self._record.info.get(item, default)


class PysamRecordFilter(VCFRecordFilter):
    def __init__(self, record: VariantRecord):
        self._record = record

    def __iter__(self):
        yield from self._record.filter

    def add(self, tag: str):
        self._record.filter.add(tag)


class PysamReader(VCFReader):
    def __init__(self, filename: str):
        self.filename = filename
        self._file = pysam.VariantFile(self.filename)
        self._header = PysamHeader(self)

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        return PysamRecord(self._iter_file.__next__())

    def reset(self):
        self._file.reset()


class PysamHeader(VCFHeader):
    def __init__(self, reader: PysamReader):
        self._reader = reader
        self._header = reader._file.header
        self._data = []
        self._data_category = defaultdict(OrderedDict)
        self._data_generic = dict()
        for r in self._header.records:
            if r.type == "GENERIC":
                self._data_generic[r.key] = r.value
                continue
            d = dict(r)
            self._data.append(d)
            if "ID" in d:
                self._data_category[r.type][d["ID"]] = d

    def contains_generic(self, key: str):
        return key in self._data_generic

    @property
    def infos(self):
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

    @property
    def filters(self):
        return self._header.filters

    @property
    def records(self):
        return self._data_category

    # def add_meta(
    #     self,
    #     key: str,
    #     value: Optional[str] = None,
    #     items: Optional[List[Tuple[str, str]]] = None,
    # ):
    #     self._file.header.add_meta(key, value, items)

    def add_generic(self, key: str, value: str):
        self._header.add_meta(key, value)

    def add_filter(self, id: str, description: str):
        self._header.add_meta(
            key="FILTER", items=[("ID", id), ("Description", description)]
        )


class PysamWriter(VCFWriter):
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        self.filename = filename
        self._file = pysam.VariantFile(
            self.filename, f"w{fmt}", header=template.header._header
        )
        self._header = PysamHeader(self)

    def write(self, record: PysamRecord):
        self._file.write(record._record)
