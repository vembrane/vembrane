from collections import OrderedDict, defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
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

from ..ann_types import NA, type_info
from ..errors import UnknownInfoField, UnknownSample


class Cyvcf2Reader(VCFReader):
    __slots__ = (
        "filename",
        "_iter_file",
        "_header",
        "_overwrite_number",
        "_current_record_idx",
    )

    def __init__(
        self,
        filename: str,
        overwrite_number: Dict[str, Dict[str, str]] = {},
    ):
        self.filename = filename
        self._file = VCF(self.filename)
        self._header = Cyvcf2Header(self, overwrite_number)
        self._overwrite_number = overwrite_number
        self._current_record_idx = 0

    def __iter__(self):
        self._iter_file = self._file.__iter__()
        return self

    def __next__(self):
        self._current_record_idx += 1
        return Cyvcf2Record(
            self._iter_file.__next__(),
            self._header,
            self._file,
            self._current_record_idx,
        )

    def reset(self):
        # cyvcv2 doesnt have a reset function
        self._file.close()
        self._file = VCF(self.filename)
        self._header = Cyvcf2Header(self, self._overwrite_number)
        self._current_record_idx = 0
        # TODO: may this workaround lead to problems?


class Cyvcf2Header(VCFHeader):
    __slots__ = ("_reader", "_data", "_data_category")

    def __init__(
        self, reader: Cyvcf2Reader, overwrite_number: Dict[str, Dict[str, str]]
    ):
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

        # override numbers
        for category, items in overwrite_number.items():
            for key, value in items.items():
                self._data_category[category][key]["Number"] = value

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

    def __iter__(self):
        raise NotImplementedError

    def __next__(self):
        raise NotImplementedError

    def add_meta(
        self,
        key: str,
        value: Optional[str] = None,
        items: Optional[List[Tuple[str, str]]] = None,
    ):
        raise NotImplementedError


class Cyvcf2Record(VCFRecord):
    __slots__ = ("_record", "_header", "_file", "_record_idx")

    def __init__(
        self, record: Variant, header: Cyvcf2Header, file: VCF, record_idx: int
    ):
        self._record = record
        self._header = header
        self._file = file
        self._record_idx = record_idx

    @property
    def contig(self) -> str:
        return self._record.CHROM

    @property
    def position(self) -> int:
        return self._record.POS

    @property
    def stop(self) -> int:
        if "END" in self.info:
            return self.info["END"]
        return self._record.POS + len(self._record.REF) - 1

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
    def formats(self) -> VCFRecordFormats:
        return Cyvcf2RecordFormats(self._record, self._header)

    @property
    def samples(self):
        return self._file.samples

    @property
    def header(self) -> VCFHeader:
        return self._header

    def __repr__(self):
        return self._record.__repr__()

    def __str__(self):
        return self._record.__str__()


class Cyvcf2RecordFormats(VCFRecordFormats):
    __slots__ = ("_record", "_header")

    def __init__(self, record: Variant, header: Cyvcf2Header):
        self._record = record
        self._header = header

    def __getitem__(self, key: str):
        return Cyvcf2RecordFormat(key, self._record, self._header)

    def __contains__(self, key):
        return key in self._record.FORMAT


class Cyvcf2RecordFormat(VCFRecordFormat):
    __slots__ = ("_header", "_record", "_format_key")

    def __init__(self, format_key: str, record, header: Cyvcf2Header):
        self._header = header
        self._record = record
        self._format_key = format_key

    def __getitem__(self, sample):
        i = self._header.samples.index(sample)
        if i == -1:
            raise UnknownSample(self._record_idx, self._record, sample)
        if self._format_key == "GT":  # genotype
            value = tuple(
                None if gt == -1 else gt for gt in self._record.genotypes[i][:-1]
            )
            return type_info(value)

        value = self._record.format(self._format_key)[i]
        meta = self._header.formats[self._format_key]
        number = meta["Number"]
        if meta["Type"] == "String":
            if not number == "1":
                value = value.split(",")
        if number == "1":
            if isinstance(value, np.ndarray):
                value = value[0].item()
            # cyvcf2 gives min int for unknown integer values
            if meta["Type"] == "Integer" and value == np.iinfo(np.int32).min:
                return NA
        return type_info(value, number)

    def __setitem__(self, key, value):
        raise NotImplementedError

    def __contains__(self, sample):
        return sample in self._header.samples

    def __eq__(self, other):
        return all(self[sample] == other[sample] for sample in self._header.samples)


class Cyvcf2RecordFilter(VCFRecordFilter):
    __slots__ = "_record"

    def __init__(self, record: Variant):
        self._record = record

    def __iter__(self):
        yield from self._record.FILTERS

    def add(self, tag: str):
        # cyvcf needs a semicolon separated string
        filter = self._record.FILTER
        if filter:
            filter += f";{tag}"
        else:
            filter = tag
        self._record.FILTER = filter

    def __contains__(self, key):
        return key in self._record.FILTERS


class Cyvcf2RecordInfo(VCFRecordInfo):
    __slots__ = ("_record", "_header", "_record_idx")

    def __init__(
        self,
        record: Cyvcf2Record,
        header: Cyvcf2Header,
    ):
        self._record = record
        self._header = header

    def __getitem__(self, key):
        # if key == "END":
        #     return get_end(self._record)
        if key not in self._header.infos.keys():
            raise UnknownInfoField(self.record.record_idx, self._record, key)

        meta = self._header.infos[key]
        if not self.__contains__(key):
            return type_info(NA, meta["Number"])
        if meta["Type"] == "Flag":
            return key in self
        value = self._record.INFO[key]
        number, typ = meta["Number"], meta["Type"]

        # for some reasons cyvcf2 doesn't split String lists, a known circumstance
        if typ == "String" and not number == "1":
            value = value.split(",")

        return type_info(value, number)

    def __setitem__(self, key, value):
        # for some reasons cyvcf2 doesn't split String lists, a known circumstance
        meta = self._header.infos[key]
        number, typ = meta["Number"], meta["Type"]
        if typ == "String" and number == ".":
            value = ",".join(value)
        self._record.INFO[key] = value

    def __contains__(self, key):
        return self._record.INFO.get(key, None) is not None


class Cyvcf2Writer(VCFWriter):
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        self._file = Writer(filename, template._file, mode=f"w{fmt}")

    def write(self, record: Cyvcf2Record):
        self._file.write_record(record._record)
