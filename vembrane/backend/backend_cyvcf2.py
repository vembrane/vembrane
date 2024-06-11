from collections import OrderedDict, defaultdict
from typing import Dict, Tuple

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
from ..errors import UnknownInfoFieldError, UnknownSampleError


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
        overwrite_number: Dict[str, Dict[str, str]] = None,
    ):
        if overwrite_number is None:
            overwrite_number = {}
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
            self._current_record_idx,
            self._header,
            self._file,
        )

    def reset(self):
        # cyvcv2 doesnt have a reset function
        metadata_generic = self._header._metadata_generic.copy()
        self._file.close()
        self._file = VCF(self.filename)
        self._header = Cyvcf2Header(self, self._overwrite_number)
        for k, v in metadata_generic.items():
            self._header.add_generic(k, v)
        self._current_record_idx = 0
        # TODO: may this workaround lead to problems?


class Cyvcf2Header(VCFHeader):
    __slots__ = (
        "_reader",
        "_data",
        "_meta_category",
        "_metadata_generic",
        "_raw_header",
    )

    def __init__(
        self,
        reader: Cyvcf2Reader,
        overwrite_number: Dict[str, Dict[str, str]],
    ):
        self._reader = reader
        self._data = []
        self._meta_category = defaultdict(OrderedDict)
        self._metadata_generic = dict()
        self._raw_header = reader._file.raw_header.split("\n")

        for r in reader._file.header_iter():
            if r.type == "GENERIC":
                continue
            d = r.info()
            self._data.append(d)
            if "ID" in d:
                self._meta_category[r.type][d["ID"]] = d

        specific_keys = {r.type for r in reader._file.header_iter()} | {"contig"}
        generic_entries = [
            r.lstrip("#").split("=", 1)
            for r in reader._file.raw_header.split("\n")
            if r.startswith("##")
        ]
        generic_entries = {k: v for (k, v) in generic_entries if k not in specific_keys}
        for k, v in generic_entries.items():
            self._metadata_generic[k] = v

        # override numbers
        for category, items in overwrite_number.items():
            for key, value in items.items():
                if key in self._meta_category[category]:
                    self._meta_category[category][key]["Number"] = value

    @property
    def records(self):
        return self._meta_category

    def contains_generic(self, key: str):
        return key in self._metadata_generic

    def get_generic(self, key: str):
        return self._metadata_generic[key]

    @property
    def infos(self):
        return self._meta_category["INFO"]

    @property
    def formats(self):
        return self._meta_category["FORMAT"]

    @property
    def filters(self):
        return self._meta_category["FILTER"]

    @property
    def samples(self):
        return self._reader._file.samples

    def add_generic(
        self,
        key: str,
        value: str,
    ):
        self._metadata_generic[key] = value
        self._reader._file.add_to_header(f"##{key}={value}")
        self._raw_header.insert(-2, f"##{key}={value}")

    def add_filter(
        self,
        id: str,
        description: str,
    ):
        self._reader._file.add_filter_to_header(
            {
                "ID": id,
                "Description": description,
            },
        )
        self._meta_category["FILTER"][id] = {
            "ID": id,
            "Description": description,
            "HeaderType": "FILTER",
        }
        self._raw_header.insert(-2, f'##FILTER=<ID={id},Description="{description}">')

    def __iter__(self):
        raise NotImplementedError

    def __next__(self):
        raise NotImplementedError

    def add_info(
        self,
        id: str,
        number: str,
        type: str,
        description: str,
    ):
        self._reader._file.add_info_to_header(
            {
                "ID": id,
                "Type": type,
                "Description": description,
                "Number": number,
            },
        )
        self._meta_category["INFO"][id] = self._reader._file.get_header_type(id)
        self._raw_header.insert(
            -2,
            f'##INFO=<ID={id},Number={number},Type={type},Description="{description}">',
        )

    @property
    def raw(self):
        # print(self._raw_header)
        # exit()
        return "\n".join(self._raw_header)

    def update_info(
        self,
        id: str,
        number: str,
        type: str,
        description: str,
    ):
        # cyvcf2 "update" doesn't work correctly
        self._raw_header = [
            (
                f"##INFO=<ID={id},Number={number},Type={type},Description={description}>"
                if line.startswith(f"##INFO=<ID={id},")
                else line
            )
            for line in self._raw_header
        ]


class Cyvcf2Record(VCFRecord):
    __slots__ = ("_file",)

    def __init__(
        self,
        record: Variant,
        record_idx: int,
        header: Cyvcf2Header,
        file: VCF,
    ):
        super().__init__(record, record_idx, header)
        self._file = file

    @property
    def contig(self) -> str:
        return self._raw_record.CHROM

    @property
    def position(self) -> int:
        return self._raw_record.POS

    @property
    def stop(self) -> int:
        if "END" in self.info:
            return self.info["END"]
        return self._raw_record.POS + len(self._raw_record.REF) - 1

    @property
    def id(self) -> str:
        return self._raw_record.ID

    @property
    def reference_allele(self) -> str:
        return self._raw_record.REF

    @property
    def alt_alleles(self) -> Tuple[str]:
        return tuple(self._raw_record.ALT)

    @property
    def quality(self) -> float:
        return self._raw_record.QUAL

    @property
    def filter(self) -> VCFRecordFilter:
        return Cyvcf2RecordFilter(self)

    @property
    def info(self) -> VCFRecordInfo:
        return Cyvcf2RecordInfo(self)

    @property
    def formats(self) -> VCFRecordFormats:
        return Cyvcf2RecordFormats(self)

    @property
    def samples(self):
        return self._file.samples

    @property
    def header(self) -> VCFHeader:
        return self._header


class Cyvcf2RecordFormats(VCFRecordFormats):
    __slots__ = ("_record", "_raw_format", "_format")

    def __init__(self, record: Cyvcf2Record):
        self._record = record
        self._raw_format = record._raw_record.FORMAT
        self._format = {}

    def __getitem__(self, key: str):
        if key not in self._format:
            self._format[key] = Cyvcf2RecordFormat(key, self._record)
        return self._format[key]

    def __contains__(self, key):
        return key in self._raw_format


class Cyvcf2RecordFormat(VCFRecordFormat):
    __slots__ = ("_header", "_record", "_raw_record", "_format_key")

    def __init__(self, format_key: str, record: Cyvcf2Record):
        self._record = record
        self._raw_record = record._raw_record
        self._header = record._header
        self._format_key = format_key

    def __getitem__(self, sample):
        i = self._header.samples.index(sample)
        if i == -1:
            raise UnknownSampleError(self._record, sample)
        if self._format_key == "GT":  # genotype
            value = tuple(
                None if gt == -1 else gt for gt in self._raw_record.genotypes[i][:-1]
            )
            return type_info(value, ".")

        value = self._raw_record.format(self._format_key)[i]
        meta = self._header.formats[self._format_key]
        number = meta["Number"]
        if meta["Type"] == "String" and not number == "1":
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

    def __init__(self, record: Cyvcf2Record):
        self._record = record._raw_record

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
    __slots__ = ("_record", "_raw_record", "_header")

    def __init__(
        self,
        record: Cyvcf2Record,
    ):
        self._record = record
        self._raw_record = record._raw_record
        self._header = record._header

    def __getitem__(self, key):
        try:
            meta = self._header.infos[key]
            number = meta["Number"]
        except KeyError as ke:
            raise UnknownInfoFieldError(self._record, key) from ke

        try:
            value = self._raw_record.INFO[key]
            typ = meta["Type"]
        except KeyError:
            print(
                f"Warning: "
                f"record {self._record.record_idx} is missing a value for key {key}, "
                f"returning NA instead."
                f"\n{self._record}\n",
            )
            return type_info(NA, number)

        # for some reason cyvcf2 doesn't split String lists, a known circumstance
        if typ == "String" and not number == "1":
            value = value.split(",")

        return type_info(value, number)

    def __setitem__(self, key, value):
        # for some reason cyvcf2 doesn't split String lists, a known circumstance
        meta = self._header.infos[key]
        number, typ = meta["Number"], meta["Type"]
        if typ == "String" and number != "1":
            value = ",".join(value)
        self._raw_record.INFO[key] = value

    def __contains__(self, key):
        return self._raw_record.INFO.get(key, None) is not None


class Cyvcf2Writer(VCFWriter):
    def __init__(self, filename: str, fmt: str, template: VCFReader):
        self._file = Writer.from_string(filename, template.header.raw, mode=f"w{fmt}")
        # self._file = Writer(filename, template._file, mode=f"w{fmt}")

    def write(self, record: Cyvcf2Record):
        self._file.write_record(record._raw_record)
