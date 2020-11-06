import ast
from itertools import chain
from typing import Dict, List, Tuple

from pysam.libcbcf import (
    VariantRecordSamples,
    VariantRecordFormat,
    VariantRecordInfo,
    VariantHeader,
    VariantRecord,
    VariantRecordSample,
)

from .ann_types import (
    NA,
    type_info,
    ANN_TYPER,
    MoreThanOneAltAllele,
)
from .common import get_annotation_keys, split_annotation_entry
from .errors import (
    UnknownSample,
    UnknownFormatField,
    UnknownInfoField,
    UnknownAnnotation,
)
from .globals import allowed_globals, custom_functions


class NoValueDict:
    def __contains__(self, item):
        try:
            value = self[item]
        except KeyError:
            return False
        return value is not NA


class Format(NoValueDict):
    def __init__(
        self,
        record_idx: int,
        name: str,
        number: str,
        record_samples: VariantRecordSamples,
    ):
        self._record_idx = record_idx
        self._name = name
        self._number = number
        self._record_samples = record_samples
        self._sample_values = {}

    def __getitem__(self, sample):
        try:
            return self._sample_values[sample]
        except KeyError:
            try:
                record_sample = self._record_samples[sample]
            except KeyError:
                raise UnknownSample(self._record_idx, sample)
            value = type_info(record_sample[self._name], self._number)
            self._sample_values[sample] = value
            return value


class Formats(NoValueDict):
    def __init__(
        self,
        record_idx: int,
        header_format_fields: Dict[str, str],
        record_format: VariantRecordFormat,
        record_samples: VariantRecordSamples,
    ):
        self._record_idx = record_idx
        self._header_format_fields = header_format_fields
        self._record_format = record_format
        self._record_samples = record_samples
        self._formats = {}

    def __getitem__(self, item):
        try:
            return self._formats[item]
        except KeyError:
            try:
                self._record_format[item]
            except KeyError:
                raise UnknownFormatField(self._record_idx, item)
            number = self._header_format_fields[item]
            format_field = Format(self._record_idx, item, number, self._record_samples)
            self._formats[item] = format_field
            return format_field


class Info(NoValueDict):
    def __init__(
        self,
        record_idx: int,
        record_info: VariantRecordInfo,
        header_info_fields: Dict[str, str],
        ann_key: str,
    ):
        self._record_idx = record_idx
        self._record_info = record_info
        self._header_info_fields = header_info_fields
        self._ann_key = ann_key
        self._info_dict = {}

    def __getitem__(self, item):
        try:
            return self._info_dict[item]
        except KeyError:
            try:
                if item == self._ann_key:
                    raise KeyError(item)
                untyped_value = self._record_info[item]
            except KeyError as ke:
                if item in self._header_info_fields:
                    value = NA
                else:
                    raise UnknownInfoField(self._record_idx, ke)
            else:
                value = self._info_dict[item] = type_info(
                    untyped_value, self._header_info_fields[item]
                )
            return value


class Annotation(NoValueDict):
    def __init__(self, ann_key: str, header: VariantHeader):
        self._record_idx = -1
        self._annotation_data = {}
        self._data = {}
        annotation_keys = get_annotation_keys(header, ann_key)
        self._ann_conv = {
            entry.name: (ann_idx, entry.convert)
            for ann_idx, entry in enumerate(map(ANN_TYPER.get_entry, annotation_keys))
        }

    def update(self, record_idx: int, annotation: str):
        self._record_idx = record_idx
        self._data.clear()
        self._annotation_data = split_annotation_entry(annotation)

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            try:
                ann_idx, convert = self._ann_conv[item]
            except KeyError as ke:
                raise UnknownAnnotation(self._record_idx, ke)
            raw_value = self._annotation_data[ann_idx].strip()
            value = self._data[item] = convert(raw_value)
            return value


UNSET = object()


class Environment(dict):
    def __init__(self, expression: str, ann_key: str, header: VariantHeader):
        self._ann_key: str = ann_key
        self._has_ann: bool = any(
            hasattr(node, "id") and isinstance(node, ast.Name) and node.id == ann_key
            for node in ast.walk(ast.parse(expression))
        )
        self._annotation: Annotation = Annotation(ann_key, header)
        self._globals = {}
        # We use self + self.func as a closure.
        self._globals = allowed_globals.copy()
        self._globals.update(custom_functions(self))
        self._func = eval(f"lambda: {expression}", self, {})

        self._getters = {
            "CHROM": self._get_chrom,
            "POS": self._get_pos,
            "ID": self._get_id,
            "ALT": self._get_alt,
            "REF": self._get_ref,
            "QUAL": self._get_qual,
            "FILTER": self._get_filter,
            "INFO": self._get_info,
            "FORMAT": self._get_format,
            "SAMPLES": self._get_samples,
            "INDEX": self._get_index,
        }

        # vembrane only supports bi-allelic records (i.e. one REF, one ALT allele).
        # Hence, for all fields in the header whith `number in {"A", "R"}`
        # we check if there is indeed only 1 value ("A") or 2 values ("R")
        # and abort otherwise.
        # Then, in the case of `number == "A"`, the value tuples only have one entry,
        # so that the value can be accessed directly and need not be accessed via
        # an index operation.
        self._numbers = {
            kind: {
                record.get("ID"): record.get("Number")
                for record in header.records
                if record.type == kind
            }
            for kind in set(r.type for r in header.records)
        }
        # At the moment, only INFO and FORMAT records are checked
        self._header_info_fields = self._numbers.get("INFO", dict())
        self._header_format_fields = self._numbers.get("FORMAT", dict())
        self._empty_globals = {name: UNSET for name in self._getters}
        self.record: VariantRecord = None
        self.idx: int = -1

    def expression_annotations(self):
        return self._has_ann

    def update_from_record(self, idx: int, record: VariantRecord):
        self.idx = idx
        self.record = record
        self._globals.update(self._empty_globals)

    def _get_chrom(self) -> str:
        value = self.record.chrom
        self._globals["CHROM"] = value
        return value

    def _get_pos(self) -> int:
        value = self.record.pos
        self._globals["POS"] = value
        return value

    def _get_id(self) -> str:
        value = self.record.id
        self._globals["ID"] = value
        return value

    def _get_ref_alt(self) -> Tuple[str, List[str]]:
        ref, *alt = chain(self.record.alleles)
        self._globals["REF"], self._globals["ALT"] = ref, alt
        return ref, alt

    def _get_ref(self) -> str:
        return self._get_ref_alt()[0]

    def _get_alt(self) -> str:
        alts = self._get_ref_alt()[1]
        if len(alts) > 1:
            raise MoreThanOneAltAllele()
        return alts[0]

    def _get_qual(self) -> float:
        value = type_info(self.record.qual)
        self._globals["QUAL"] = value
        return value

    def _get_filter(self) -> List[str]:
        value = list(self.record.filter)
        self._globals["FILTER"] = value
        return value

    def _get_info(self) -> Info:
        value = Info(
            self.idx,
            self.record.info,
            self._header_info_fields,
            self._ann_key,
        )
        self._globals["INFO"] = value
        return value

    def _get_format(self) -> Formats:
        value = Formats(
            self.idx,
            self._header_format_fields,
            self.record.format,
            self.record.samples,
        )
        self._globals["FORMAT"] = value
        return value

    def _get_samples(self) -> List[VariantRecordSample]:
        value = list(self.record.samples)
        self._globals["SAMPLES"] = value
        return value

    def __getitem__(self, item):
        if item == self._ann_key:
            return self._annotation
        value = self._globals[item]
        if value is UNSET:
            value = self._getters[item]()
        return value

    def _get_index(self) -> int:
        self._globals["INDEX"] = self.idx
        return self.idx

    def evaluate(self, annotation: str = "") -> bool:
        if self._has_ann:
            self._annotation.update(self.idx, annotation)
        return self._func()

    def table(self, annotation: str = "") -> tuple:
        if self._has_ann:
            self._annotation.update(self.idx, annotation)
        return self._func()
