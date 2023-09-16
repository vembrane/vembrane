from collections import OrderedDict, defaultdict
from typing import Tuple

from cyvcf2.cyvcf2 import VCF, Variant, Writer

from vembrane.backend.base import (
    VCFHeader,
    VCFReader,
    VCFRecord,
    VCFRecordFilter,
    VCFRecordFormat,
    VCFRecordInfo,
    VCFWriter,
)


class Cyvcf2VCFReader(VCFReader):
    def __init__(self, filename: str):
        self.filename = filename
        self._file = VCF(self.filename)
        self._header = Cyvcf2VCFHeader(self)

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        return Cyvcf2VCFRecord(self._iter_file.__next__(), self._header)

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

        for r in reader._file.header_iter():
            if r.type == "GENERIC":
                continue
            d = r.info()
            self._data.append(d)
            if "ID" in d:
                self._data_category[r.type][d["ID"]] = d

    @property
    def records(self):
        return self._data_category

    def contains_generic(self, key: str):
        return self._reader._file.contains(key)

    @property
    def infos(self):
        return self._data_category["INFO"]

    @property
    def samples(self):
        return self._reader._file.samples

    def add_generic(
        self,
        key: str,
        value: str,
    ):
        self._reader._file.add_to_header(f"##{key}={value}")

    def add_filter(self, id: str, description: str):
        self._reader._file.add_filter_to_header({"ID": id, "Description": description})

    @property
    def filters(self):
        return self._data_category["FILTER"]


class Cyvcf2VCFRecord(VCFRecord):
    def __init__(self, record: Variant, header: Cyvcf2VCFHeader):
        self._record = record
        self._header = header

    @property
    def contig(self) -> str:
        return self._record.CHROM

    @property
    def position(self) -> int:
        return self._record.POS

    @property
    def id(self) -> str:
        return self._record.ID

    @property
    def reference_allele(self) -> str:
        return self._record.REF

    @property
    def alt_alleles(self) -> Tuple[str]:
        return tuple(self._record.ALT)

    @property
    def quality(self) -> float:
        return self._record.QUAL

    @property
    def filter(self) -> VCFRecordFilter:
        return Cyvcf2RecordFilter(self._record)

    @property
    def info(self) -> VCFRecordInfo:
        return Cyvcf2RecordInfo(self._record, self._header)

    @property
    def format(self) -> VCFRecordFormat:
        return Cyvcf2RecordInfo(self._record)

    # @property
    # def samples(self):
    #     return self._record.samples

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()

    def __eq__(self, other):
        return (
            self._record.__repr__() == other._record.__repr__()
        )  # TODO: implement real check for values


class Cyvcf2RecordFilter(VCFRecordFilter):
    def __init__(self, record: Variant):
        self._record = record

    def __iter__(self):
        yield from self._record.FILTERS

    def add(self, tag: str):
        self._record.FILTERS.append(tag)


class Cyvcf2RecordInfo(VCFRecordInfo):
    def __init__(
        self,
        record: Variant,
        header: Cyvcf2VCFHeader,
    ):
        self._record = record
        self._header = header

    def __getitem__(self, key):
        value = self._record.INFO[key]
        # for some reasons cyvcf2 doesn't split String lists, a known circumstance
        meta = self._header.infos[key]
        number, typ = meta["Number"], meta["Type"]
        if typ == "String" and number == ".":
            value = value.split(",")
        return value

    def __setitem__(self, key, value):
        # for some reasons cyvcf2 doesn't split String lists, a known circumstance
        meta = self._header.infos[key]
        number, typ = meta["Number"], meta["Type"]
        if typ == "String" and number == ".":
            value = ",".join(value)

        self._record.INFO[key] = value

    def __contains__(self, item):
        return item in self._record.INFO

    def get(self, item, default=None):
        return self._record.INFO.get(item, default)


class Cyvcf2VCFWriter(VCFWriter):
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        self._file = Writer(filename, template._file, mode=f"w{fmt}")

    def write(self, record: Cyvcf2VCFRecord):
        self._file.write_record(record._record)
