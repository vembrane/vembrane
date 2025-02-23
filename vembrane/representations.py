import ast
from types import CodeType, MappingProxyType
from typing import Any

from .ann_types import ANN_TYPER, NA, MoreThanOneAltAlleleError, NvFloat, NvInt
from .backend.base import VCFHeader, VCFRecord, VCFRecordFormats, VCFRecordInfo
from .common import get_annotation_keys, split_annotation_entry
from .errors import MalformedAnnotationError, NonBoolTypeError, UnknownAnnotationError
from .globals import _explicit_clear, allowed_globals, custom_functions


class NoValueDict:
    def __contains__(self, item) -> bool:
        try:
            value = self[item]
        except KeyError:
            return False
        return value is not NA


class DefaultGet:
    def get(self, item, default=NA):
        v = self[item]
        if v is not NA:
            return v
        else:
            return default


class Annotation(NoValueDict, DefaultGet):
    def __init__(self, ann_key: str, header: VCFHeader) -> None:
        self._record_idx = -1
        self._record: VCFRecord | None = None
        self._annotation_data: list[str] = []
        self._data: dict[str, Any] = {}
        annotation_keys = get_annotation_keys(header, ann_key)
        self._ann_conv = {
            entry.name: (ann_idx, entry.convert)
            for ann_idx, entry in enumerate(map(ANN_TYPER.get_entry, annotation_keys))
        }

    def update(self, record_idx: int, record: VCFRecord, annotation: str):
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
            except KeyError as ke2:
                raise UnknownAnnotationError(
                    self._record,
                    item,
                ) from ke2
            if ann_idx >= len(self._annotation_data):
                raise MalformedAnnotationError(
                    self._record_idx,
                    self._record,
                    item,
                    ann_idx,
                ) from None
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
        ann_key: str,
        header: VCFHeader,
        auxiliary: dict[str, set[str]] = MappingProxyType({}),
    ) -> None:
        self._ann_key: str = ann_key
        self._annotation: Annotation = Annotation(ann_key, header)
        self._globals: dict[str, Any] = {}
        # We use self + self.func as a closure.
        self._globals = dict(allowed_globals)
        self._globals.update(custom_functions(self))
        self._globals["SAMPLES"] = list(header.samples)
        # REF/ALT alleles are cached separately to raise "MoreThanOneAltAllele"
        # only if ALT (but not REF) is accessed (and ALT has multiple entries).
        self._alleles = None

        self.update(**_explicit_clear)

        self._getters = {
            "AUX": self._get_aux,
            "CHROM": self._get_chrom,
            "POS": self._get_pos,
            "END": self._get_end,
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

        # always explicitly set "Number" for certain fields
        # which get special pysam treatment:
        # - `FORMAT["GT"]` is always parsed as a list of integer values
        # self._numbers.get("FORMAT", {})["GT"] = "."

        # At the moment, only INFO and FORMAT records are checked
        self._empty_globals = {name: UNSET for name in self._getters}
        self.record: VCFRecord = None
        self.idx: int = -1
        self.aux = auxiliary

    def update_from_record(self, idx: int, record: VCFRecord):
        self.idx = idx
        self.record = record
        self._globals.update(self._empty_globals)
        self._alleles = None

    def update_data(self, data):
        self._globals["DATA"] = data

    def update_annotation(self, annotation):
        self._annotation.update(self.idx, self.record, annotation)

    def _get_chrom(self) -> str:
        value = self.record.contig
        self._globals["CHROM"] = value
        return value

    def _get_pos(self) -> int:
        value = self.record.position
        self._globals["POS"] = value
        return value

    def _get_end(self) -> NvInt:
        value = self.record.end
        self._globals["END"] = value
        return value

    def _get_id(self) -> str | None:
        value = self.record.id
        self._globals["ID"] = value
        return value

    def _get_alleles(self) -> tuple[str, ...]:
        alleles = self._alleles
        if not alleles:
            alleles = self._alleles = self.record.alleles
        return alleles

    def _get_ref(self) -> str:
        alleles = self._get_alleles()
        value = alleles[0]
        self._globals["REF"] = value
        return value

    def _get_alt(self) -> str:
        alleles = self._get_alleles()
        if len(alleles) > 2:
            raise MoreThanOneAltAlleleError
        value = alleles[1] if len(alleles) == 2 else NA
        self._globals["ALT"] = value
        return value

    def _get_qual(self) -> NvFloat:
        value: NvFloat = NA if self.record.quality is None else self.record.quality
        self._globals["QUAL"] = value
        return value

    def _get_filter(self) -> list[str]:
        value = list(self.record.filter)
        self._globals["FILTER"] = value
        return value

    def _get_info(self) -> VCFRecordInfo:
        value = self.record.info
        self._globals["INFO"] = value
        return value

    def _get_format(self) -> VCFRecordFormats:
        value = self.record.formats
        self._globals["FORMAT"] = value
        return value

    def __getitem__(self, item):
        print(item)
        if item == self._ann_key:
            return self._annotation
        value = self._globals[item]
        if value is UNSET:
            value = self._getters[item]()
        return value

    def _get_index(self) -> int:
        self._globals["INDEX"] = self.idx
        return self.idx

    def _get_aux(self) -> dict[str, set[str]]:
        self._globals["AUX"] = self.aux
        return self.aux


class EvalEnvironment(Environment):
    def __init__(
        self,
        expression: str,
        ann_key: str,
        header: VCFHeader,
        auxiliary: dict[str, set[str]] = MappingProxyType({}),
        evaluation_function_template: str = "lambda: {expression}",
    ) -> None:
        super().__init__(ann_key, header, auxiliary)

        self._has_ann: bool = any(
            hasattr(node, "id") and isinstance(node, ast.Name) and node.id == ann_key
            for node in ast.walk(ast.parse(expression))
        )

        func_str: str = evaluation_function_template.format(expression=expression)

        # The VCF specification only allows 32bit floats.
        # Comparisons such as `INFO["some_float"] > CONST` may yield unexpected results,
        # because `some_float` is a 32bit float while `CONST` is a python 64bit float.
        # Example: `c_float(0.6) = 0.6000000238418579 > 0.6`.
        # We can work around this by wrapping user provided floats in `ctypes.c_float`,
        # which should follow the same IEEE 754 specification as the VCF spec, as most C
        # compilers should follow this standard (https://stackoverflow.com/a/17904539):

        # parse the expression, obtaining an AST
        expression_ast = ast.parse(func_str, mode="eval")

        # wrap each float constant in numpy.float32
        expression_ast = WrapFloat32Visitor().visit(expression_ast)

        # housekeeping
        expression_ast = ast.fix_missing_locations(expression_ast)

        # compile the now-fixed code-tree
        func: CodeType = compile(expression_ast, filename="<string>", mode="eval")

        self._func = eval(func, self, {})

    def expression_annotations(self):
        return self._has_ann

    def evaluate(self, annotation: str = "") -> bool:
        if self._has_ann:
            self.update_annotation(annotation)
        keep = self._func()
        if not isinstance(keep, bool):
            raise NonBoolTypeError(keep)
        return keep

    def table(self, annotation: str = "") -> tuple:
        if self._has_ann:
            self.update_annotation(annotation)
        return self._func()


class ModifiableEnvironment(Environment):
    def __setitem__(self, key, value):
        self._globals[key] = value
        return value
