from collections import OrderedDict, defaultdict
from typing import Tuple

from cyvcf2.cyvcf2 import VCF, Variant, Writer

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

from ..ann_types import NA
from ..errors import UnknownSample


class Cyvcf2Reader(VCFReader):
    def __init__(self, filename: str):
        self.filename = filename
        self._file = VCF(self.filename)
        self._header = Cyvcf2Header(self)

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        return Cyvcf2Record(self._iter_file.__next__(), self._header, self._file)

    # @property
    # def header(self):
    #     return self._header

    # @abstractmethod
    # def reset(self):
    # pass


class Cyvcf2Header(VCFHeader):
    def __init__(self, reader: Cyvcf2Reader):
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
    def formats(self):
        return self._data_category["FORMAT"]

    @property
    def filters(self):
        return self._data_category["FILTER"]

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


class Cyvcf2Record(VCFRecord):
    def __init__(self, record: Variant, header: Cyvcf2Header, file: VCF):
        self._record = record
        self._header = header
        self._file = file

    @property
    def contig(self) -> str:
        return self._record.CHROM

    @property
    def position(self) -> int:
        return self._record.POS

    @property
    def stop(self) -> int:
        return self._record.POS + len(self._record.REF)

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
        return Cyvcf2RecordFormat(self._record, self._file, self._header)

    @property
    def formats(self) -> VCFRecordFormats:
        return Cyvcf2RecordFormats(self._record, self._header)

    @property
    def samples(self):
        return self._file.samples

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()

    def __eq__(self, other):
        return (
            self._record.__repr__() == other._record.__repr__()
        )  # TODO: implement real check for values


class Cyvcf2RecordFormats(VCFRecordFormats):
    def __init__(self, record: Variant, header: Cyvcf2Header):
        self._record = record
        self._header = header

    def __getitem__(self, item):
        return Cyvcf2RecordFormat(item, self._record, self._header)


class Cyvcf2RecordFormat(VCFRecordFormat):
    def __init__(self, format_key: str, record, header: Cyvcf2Header):
        self._header = header
        self._record = record
        self._format_key = format_key

    def __getitem__(self, sample):
        value = self._record.format(self._format_key)
        i = self._header.samples.index(sample)
        if i == -1:
            raise UnknownSample(self._record_idx, self._record, sample)
        if self._format_key == "GT":  # genotype
            return tuple(
                NA if gt == -1 else gt for gt in self._record.genotypes[i][:-1]
            )
        if self._header.formats[self._format_key]["Number"] == "1":
            return value[i].tolist()[0]
        return tuple(value[i].tolist())


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
        header: Cyvcf2Header,
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
            if len(value) == 1:
                value = value[0].split("/")
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


class Cyvcf2Writer(VCFWriter):
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        self._file = Writer(filename, template._file, mode=f"w{fmt}")

    def write(self, record: Cyvcf2Record):
        self._file.write_record(record._record)
