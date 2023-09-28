from collections import OrderedDict, defaultdict
from typing import Dict, List, Optional, Tuple

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

from ..ann_types import NA, type_info
from ..errors import UnknownInfoField, UnknownSample


class PysamRecord(VCFRecord):
    __slots__ = ("_record", "_header")

    def __init__(self, record: VariantRecord, header: VCFHeader):
        self._record = record
        self._header = header

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
        return PysamRecordInfo(self._record, self._header)

    @property
    def format(self) -> VCFRecordFormat:
        return PysamRecordFormat(self._record)

    @property
    def formats(self) -> VCFRecordFormats:
        return PysamRecordFormats(self._record, self._header)

    @property
    def samples(self):
        return self._record.samples

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()


class PysamRecordFormats(VCFRecordFormats):
    __slots__ = ("_record", "_header")

    def __init__(self, record: VariantRecord, header: VCFHeader):
        self._header = header
        self._record = record

    def __getitem__(self, key):
        return PysamRecordFormat(key, self._record, self._header)


class PysamRecordFormat(VCFRecordFormat):
    __slots__ = ("_record", "_header")

    def __init__(
        self,
        format_key: str,
        record: VariantRecord,
        header: VCFHeader,
    ):
        self._format_key = format_key
        self._record = record
        self._header = header

    def __getitem__(self, sample):
        if not self.__contains__(sample):
            raise UnknownSample(self._record_idx, self._record, sample)
        meta = self._header.formats[self._format_key]
        return type_info(self._record.samples[sample][self._format_key], meta["Number"])

    def __setitem__(self, key, value):
        self._record.format[key] = value

    def __contains__(self, sample):
        return sample in self._record.samples.keys()


class PysamRecordInfo(VCFRecordInfo):
    __slots__ = ("_record", "_header")

    def __init__(
        self,
        record: VariantRecord,
        header: VCFHeader,
    ):
        self._record = record
        self._header = header

    def __getitem__(self, key):
        if key == "END":
            return get_end(self._record)
        if key not in self._header.infos.keys():
            raise UnknownInfoField(0, self._record, key)  # TODO: self._record_idx
        meta = self._header.infos[key]
        if not self.__contains__(key):
            return type_info(NA, meta["Number"])
        return type_info(self._record.info[key], meta["Number"])

    def __setitem__(self, key, value):
        self._record.info[key] = value

    def __contains__(self, key):
        return key in self._record.info.keys()


class PysamRecordFilter(VCFRecordFilter):
    __slots__ = "_record"

    def __init__(self, record: VariantRecord):
        self._record = record

    def __iter__(self):
        yield from self._record.filter

    def add(self, tag: str):
        self._record.filter.add(tag)


class PysamReader(VCFReader):
    __slots__ = (
        "filename",
        "_iter_file",
        "_header",
    )

    def __init__(
        self,
        filename: str,
        overwrite_number: Dict[str, Dict[str, str]] = {},
    ):
        self.filename = filename
        self._file = pysam.VariantFile(self.filename)
        self._header = PysamHeader(self, overwrite_number)

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        return PysamRecord(self._iter_file.__next__(), self._header)

    def reset(self):
        self._file.reset()


class PysamHeader(VCFHeader):
    __slots__ = ("_record", "_header", "_data", "_data_category", "_data_generic")

    def __init__(self, reader: PysamReader, overwrite_number={}):
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

        # override numbers
        for category, items in overwrite_number.items():
            for key, value in items.items():
                self._data_category[category][key]["Number"] = value

    def contains_generic(self, key: str):
        return key in self._data_generic

    @property
    def infos(self):
        return self._data_category["INFO"]

    @property
    def formats(self):
        return self._data_category["FORMAT"]

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

    def add_generic(self, key: str, value: str):
        self._header.add_meta(key, value)

    def add_filter(self, id: str, description: str):
        self._header.add_meta(
            key="FILTER", items=[("ID", id), ("Description", description)]
        )

    def add_meta(
        self,
        key: str,
        value: Optional[str] = None,
        items: Optional[List[Tuple[str, str]]] = None,
    ):
        raise NotImplementedError


class PysamWriter(VCFWriter):
    __slots__ = ("filename", "_header", "_file")

    def __init__(self, filename: str, fmt: str, template: VCFReader):
        self.filename = filename
        self._file = pysam.VariantFile(
            self.filename, f"w{fmt}", header=template.header._header
        )
        self._header = PysamHeader(self)

    def write(self, record: PysamRecord):
        self._file.write(record._record)


def get_end(record: VariantRecord):
    if is_bnd_record(record):
        return NA
    else:
        # record.stop is pysams unified approach to get the end position of a variant.
        # It either considers the alt allele or the END field, depending on the record.
        # Stop is 0-based, but for consistency with POS we convert into 1-based.
        # Since for 1-based coordinates, the expectation is that END is inclusive
        # instead of exclusive (as it is with 0-based), we do not need to add 1.
        return record.stop


def is_bnd_record(record: VariantRecord) -> bool:
    return "SVTYPE" in record.info and record.info.get("SVTYPE", None) == "BND"
