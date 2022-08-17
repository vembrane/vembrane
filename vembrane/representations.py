import ast
from itertools import chain
from typing import Dict, List, Set, Tuple

from pysam.libcbcf import VariantHeader, VariantRecord, VariantRecordSamples

from .ann_types import ANN_TYPER, NA, MoreThanOneAltAllele, type_info
from .common import get_annotation_keys, split_annotation_entry
from .errors import (
    UnknownAnnotation,
    UnknownFormatField,
    UnknownInfoField,
    UnknownSample,
)
from .globals import _explicit_clear, allowed_globals, custom_functions


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
        record: VariantRecord,
        name: str,
        number: str,
        record_samples: VariantRecordSamples,
    ):
        self._record_idx = record_idx
        self._record = record
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
                raise UnknownSample(self._record_idx, self._record, sample)
            value = type_info(
                record_sample[self._name], self._number, self._name, self._record_idx
            )
            self._sample_values[sample] = value
            return value


class Formats(NoValueDict):
    def __init__(
        self,
        record_idx: int,
        record: VariantRecord,
        header_format_fields: Dict[str, str],
    ):
        self._record = record
        self._record_idx = record_idx
        self._header_format_fields = header_format_fields
        self._record_format = record.format
        self._record_samples = record.samples
        self._formats = {}

    def __getitem__(self, item):
        try:
            return self._formats[item]
        except KeyError:
            try:
                self._record_format[item]
            except KeyError:
                raise UnknownFormatField(self._record_idx, self._record, item)
            number = self._header_format_fields[item]
            format_field = Format(
                self._record_idx, self._record, item, number, self._record_samples
            )
            self._formats[item] = format_field
            return format_field


class Info(NoValueDict):
    def __init__(
        self,
        record_idx: int,
        record: VariantRecord,
        header_info_fields: Dict[str, str],
        ann_key: str,
    ):
        self._record_idx = record_idx
        self._record = record
        self._record_info = record.info
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
            except KeyError:
                if item in self._header_info_fields:
                    value = NA
                else:
                    raise UnknownInfoField(self._record_idx, self._record, item)
            else:
                value = self._info_dict[item] = type_info(
                    untyped_value,
                    self._header_info_fields[item],
                    item,
                    self._record_idx,
                )
            return value


class Annotation(NoValueDict):
    def __init__(self, ann_key: str, header: VariantHeader):
        self._record_idx = -1
        self._record = None
        self._annotation_data = {}
        self._data = {}
        annotation_keys = get_annotation_keys(header, ann_key)
        self._ann_conv = {
            entry.name: (ann_idx, entry.convert)
            for ann_idx, entry in enumerate(map(ANN_TYPER.get_entry, annotation_keys))
        }

    def update(self, record_idx: int, record: VariantRecord, annotation: str):
        self._record_idx = record_idx
        self._record = record
        self._data.clear()
        self._annotation_data = split_annotation_entry(annotation)

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            try:
                ann_idx, convert = self._ann_conv[item]
            except KeyError:
                raise UnknownAnnotation(self._record_idx, self._record, item)
            raw_value = self._annotation_data[ann_idx].strip()
            value = self._data[item] = convert(raw_value)
            return value


UNSET = object()


class WrapFloat32Visitor(ast.NodeTransformer):
    def visit_Constant(self, node):
        from ctypes import c_float

        if not isinstance(node.value, float):
            return node

        return ast.Constant(c_float(node.value).value)


class Environment(dict):
    def __init__(
        self,
        expression: str,
        ann_key: str,
        header: VariantHeader,
        auxiliary: Dict[str, Set[str]] = {},
        overwrite_number: Dict[str, Dict[str, str]] = {},
        evaluation_function_template: str = "lambda: {expression}",
    ):
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
        self._globals["SAMPLES"] = list(header.samples)
        # REF/ALT alleles are cached separately to raise "MoreThanOneAltAllele"
        # only if ALT (but not REF) is accessed (and ALT has multiple entries).
        self._alleles = None

        func = evaluation_function_template.format(expression=expression)

        # The VCF specification only allows 32bit floats.
        # Comparisons such as `INFO["some_float"] > CONST` may yield unexpected results,
        # because `some_float` is a 32bit float while `CONST` is a python 64bit float.
        # Example: `c_float(0.6) = 0.6000000238418579 > 0.6`.
        # We can work around this by wrapping user provided floats in `ctypes.c_float`,
        # which should follow the same IEEE 754 specification as the VCF spec, as most C
        # compilers should follow this standard (https://stackoverflow.com/a/17904539):

        # parse the expression, obtaining an AST
        expression_ast = ast.parse(func, mode="eval")

        # wrap each float constant in numpy.float32
        expression_ast = WrapFloat32Visitor().visit(expression_ast)

        # housekeeping
        expression_ast = ast.fix_missing_locations(expression_ast)

        # compile the now-fixed code-tree
        func = compile(expression_ast, filename="<string>", mode="eval")

        self.update(**_explicit_clear)
        self._func = eval(func, self, {})

        self._getters = {
            "AUX": self._get_aux,
            "CHROM": self._get_chrom,
            "POS": self._get_pos,
            "ID": self._get_id,
            "ALT": self._get_alt,
            "REF": self._get_ref,
            "QUAL": self._get_qual,
            "FILTER": self._get_filter,
            "INFO": self._get_info,
            "FORMAT": self._get_format,
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
                record.get("ID"): overwrite_number.get(kind, {}).get(record.get("ID"))
                or record.get("Number")
                for record in header.records
                if record.type == kind
            }
            for kind in set(r.type for r in header.records)
        }

        # always explicitly set "Number" for certain fields
        # which get special pysam treatment:
        # - `FORMAT["GT"]` is always parsed as a list of integer values
        self._numbers.get("FORMAT", {})["GT"] = "."

        # At the moment, only INFO and FORMAT records are checked
        self._header_info_fields = self._numbers.get("INFO", dict())
        self._header_format_fields = self._numbers.get("FORMAT", dict())
        self._empty_globals = {name: UNSET for name in self._getters}
        self.record: VariantRecord = None
        self.idx: int = -1
        self.aux = auxiliary

    def expression_annotations(self):
        return self._has_ann

    def update_from_record(self, idx: int, record: VariantRecord):
        self.idx = idx
        self.record = record
        self._globals.update(self._empty_globals)
        self._alleles = None

    def update_data(self, data):
        self._globals["DATA"] = data

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

    def _get_alleles(self) -> Tuple[str, List[str]]:
        alleles = self._alleles
        if not alleles:
            alleles = self._alleles = tuple(chain(self.record.alleles))
        return alleles

    def _get_ref(self) -> str:
        alleles = self._get_alleles()
        value = alleles[0]
        self._globals["REF"] = value
        return value

    def _get_alt(self) -> str:
        alleles = self._get_alleles()
        if len(alleles) > 2:
            raise MoreThanOneAltAllele()
        value = alleles[1] if len(alleles) == 2 else NA
        self._globals["ALT"] = value
        return value

    def _get_qual(self) -> float:
        value = NA if self.record.qual is None else self.record.qual
        self._globals["QUAL"] = value
        return value

    def _get_filter(self) -> List[str]:
        value = list(self.record.filter)
        self._globals["FILTER"] = value
        return value

    def _get_info(self) -> Info:
        value = Info(
            self.idx,
            self.record,
            self._header_info_fields,
            self._ann_key,
        )
        self._globals["INFO"] = value
        return value

    def _get_format(self) -> Formats:
        value = Formats(
            self.idx,
            self.record,
            self._header_format_fields,
        )
        self._globals["FORMAT"] = value
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

    def _get_aux(self) -> Dict[str, Set[str]]:
        self._globals["AUX"] = self.aux
        return self.aux

    def evaluate(self, annotation: str = "") -> bool:
        if self._has_ann:
            self._annotation.update(self.idx, self.record, annotation)
        return self._func()

    def table(self, annotation: str = "") -> tuple:
        if self._has_ann:
            self._annotation.update(self.idx, self.record, annotation)
        return self._func()
