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
    __slots__ = ("_record", "_header", "record_idx")

    def __init__(self, record: VariantRecord, header: VCFHeader, record_idx: int):
        self._record = record
        self._header = header
        self.record_idx = record_idx

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
        return PysamRecordInfo(self)

    @property
    def formats(self) -> VCFRecordFormats:
        return PysamRecordFormats(self)

    @property
    def samples(self):
        return self._header.samples

    @property
    def header(self) -> VCFHeader:
        return self._header

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()


class PysamRecordFormats(VCFRecordFormats):
    __slots__ = "_record"

    def __init__(self, record: PysamRecord):
        self._record = record

    def __getitem__(self, key):
        return PysamRecordFormat(key, self._record)


class PysamRecordFormat(VCFRecordFormat):
    __slots__ = ("_record", "_header", "_raw_record")

    def __init__(
        self,
        format_key: str,
        record: PysamRecord,
    ):
        self._format_key = format_key
        self._record = record
        self._header = record._header
        self._raw_record: VariantRecord = record._record

    def __getitem__(self, sample):
        if not self.__contains__(sample):
            raise UnknownSample(self.record, sample)
        meta = self._header.formats[self._format_key]
        return type_info(
            self._raw_record.samples[sample][self._format_key], meta["Number"]
        )

    def __setitem__(self, key, value):
        self._raw_record.format[key] = value

    def __contains__(self, sample):
        return sample in self._header.samples


class PysamRecordInfo(VCFRecordInfo):
    __slots__ = ("_record", "_raw_record")

    def __init__(
        self,
        record: PysamRecord,
    ):
        self._record = record
        self._raw_record: VariantRecord = record._record

    def __getitem__(self, key):
        if key == "END":
            return self._record.end
        if key not in (infos := self._record._header.infos).keys():
            raise UnknownInfoField(self._record, key)
        meta = infos[key]
        if not self.__contains__(key):
            return type_info(NA, meta["Number"])
        return type_info(self._raw_record.info[key], meta["Number"])

    def __setitem__(self, key, value):
        self._raw_record.info[key] = value

    def __contains__(self, key):
        return key in self._raw_record.info.keys()


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
        "_current_record_idx",
    )

    def __init__(
        self,
        filename: str,
        overwrite_number: Dict[str, Dict[str, str]] = {},
    ):
        self.filename = filename
        self._file = pysam.VariantFile(self.filename)
        self._header = PysamHeader(self, overwrite_number)
        self._current_record_idx = 0

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        self._current_record_idx += 1
        return PysamRecord(
            self._iter_file.__next__(), self._header, self._current_record_idx
        )

    def reset(self):
        self._current_record_idx = 0
        self._file.reset()


class PysamHeader(VCFHeader):
    __slots__ = (
        "_record",
        "_header",
        "_metadata",
        "_metadata_category",
        "_metadata_generic",
        "_samples",
    )

    def __init__(self, reader: PysamReader, overwrite_number={}):
        self._reader = reader
        self._header = reader._file.header
        self._metadata = []
        self._metadata_category = defaultdict(OrderedDict)
        self._metadata_generic = dict()
        self._samples = OrderedDict((s, None) for s in self._header.samples)

        for r in self._header.records:
            if r.type == "GENERIC":
                self._metadata_generic[r.key] = r.value
                continue
            d = dict(r)
            self._metadata.append(d)
            if "ID" in d:
                self._metadata_category[r.type][d["ID"]] = d

        # override numbers
        for category, items in overwrite_number.items():
            for key, value in items.items():
                if key in self._metadata_category[category]:
                    self._metadata_category[category][key]["Number"] = value

    def contains_generic(self, key: str):
        return key in self._metadata_generic

    @property
    def infos(self):
        return self._metadata_category["INFO"]

    @property
    def formats(self):
        return self._metadata_category["FORMAT"]

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
        return self._samples

    @property
    def filters(self):
        return self._header.filters

    @property
    def records(self):
        return self._metadata_category

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
