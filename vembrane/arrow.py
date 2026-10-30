import ast
import sys
from collections.abc import Mapping
from typing import Any

import pyarrow as pa

from vembrane.ann_types import ANN_TYPER, NA, InfoTuple, is_na
from vembrane.backend.base import VCFHeader
from vembrane.errors import VembraneError


class ArrowTypes:
    warnings: set[str] = set()

    mapping: dict[str, pa.DataType] = {
        "int": pa.int64(),
        "float": pa.float64(),
        "str": pa.string(),
        "bool": pa.bool_(),
        "list[int]": pa.list_(pa.int64()),
        "list[float]": pa.list_(pa.float64()),
        "list[str]": pa.list_(pa.string()),
        "list[bool]": pa.list_(pa.bool_()),
        "map[str,float]": pa.map_(pa.string(), pa.float64()),
        "map[str,list[str]]": pa.map_(pa.string(), pa.list_(pa.string())),
        "RangeTotal": pa.struct(
            [
                pa.field("start", pa.int64(), nullable=False),
                pa.field("stop", pa.int64(), nullable=False),
                pa.field("total", pa.int64(), nullable=False),
            ]
        ),
        "NumberTotal": pa.struct(
            [
                pa.field("number", pa.int64(), nullable=False),
                pa.field("total", pa.int64(), nullable=False),
            ]
        ),
        "PosRange": pa.struct(
            [
                pa.field("start", pa.int64(), nullable=True),
                pa.field("end", pa.int64(), nullable=True),
            ]
        ),
    }

    # types of the record fields that an expression can name directly
    field_types: dict[str, str] = {
        "CHROM": "str",
        "POS": "int",
        "END": "int",
        "ID": "str",
        "REF": "str",
        "ALT": "str",
        "QUAL": "float",
        "FILTER": "list[str]",
        "INDEX": "int",
        "SAMPLE": "str",
    }

    header_types: dict[str, str] = {
        "Integer": "int",
        "Float": "float",
        "String": "str",
        "Character": "str",
        "Flag": "bool",
    }

    def __init__(self):
        self.python_types = {}
        self.arrow_types = {}

    def declare(
        self, colname: str, expression: str, header: VCFHeader, ann_key: str
    ) -> bool:
        """
        Set the type of a column that just reads one field, as declared by the
        VCF header or the annotation types.
        Return False if the column expression does something else.
        """
        py_type = self.declared_type(expression, header, ann_key)
        if py_type is None:
            return False
        self.arrow_types[colname] = self.python_type_to_arrow_type(py_type)
        self.python_types[colname] = py_type
        return True

    @classmethod
    def declared_type(
        cls, expression: str, header: VCFHeader, ann_key: str
    ) -> str | None:
        try:
            node = ast.parse(expression, mode="eval").body
        except SyntaxError:
            return None

        match node:
            case ast.Name(id=name) if name in cls.field_types:
                return cls.field_types[name]
            case ast.Subscript(
                value=ast.Name(id="INFO"), slice=ast.Constant(value=str() as key)
            ) if key in header.infos:
                return cls.header_type(header.infos[key])
            case (
                # FORMAT["DP"][SAMPLE], or (FORMAT["DP"] or {}).get(SAMPLE, NA) for ALL
                ast.Subscript(
                    value=ast.Subscript(
                        value=ast.Name(id="FORMAT"),
                        slice=ast.Constant(value=str() as key),
                    )
                )
                | ast.Call(
                    func=ast.Attribute(
                        value=ast.BoolOp(
                            op=ast.Or(),
                            values=[
                                ast.Subscript(
                                    value=ast.Name(id="FORMAT"),
                                    slice=ast.Constant(value=str() as key),
                                ),
                                ast.Dict(keys=[]),
                            ],
                        ),
                        attr="get",
                    )
                )
            ) if key in header.formats:
                # both backends return genotypes as tuples of allele indices
                if key == "GT":
                    return "list[int]"
                return cls.header_type(header.formats[key])
            case ast.Subscript(
                value=ast.Name(id=name), slice=ast.Constant(value=str() as key)
            ) if name in ("ANN", ann_key):
                return ANN_TYPER.get_entry(key).value_type
        return None

    @classmethod
    def header_type(cls, meta: dict[str, str]) -> str | None:
        number = meta["Number"]
        if number == "0":
            return "bool"
        py_type = cls.header_types.get(meta["Type"])
        if py_type is None or number in ("1", "A"):
            return py_type
        return f"list[{py_type}]"

    def infer(self, colname: str, values: list[Any]) -> None:
        if all(value is NA for value in values):
            if colname not in self.warnings:
                self.warnings.add(colname)
                print(
                    f"Warning: Column '{colname}' is NA-only "
                    f"in the first {len(values)} rows or no rows present. "
                    "Assuming type 'string'."
                    "This warning will only be printed once per column.",
                    file=sys.stderr,
                )
            self.arrow_types[colname] = pa.string()
            self.python_types[colname] = "str"
            return

        py_type = self.python_type(values, colname, is_column_values=True)
        self.arrow_types[colname] = self.python_type_to_arrow_type(py_type)
        self.python_types[colname] = py_type

    @classmethod
    def python_type(
        cls, value: Any, colname: str, is_column_values: bool = False
    ) -> str:
        if isinstance(value, str):
            return "str"

        if isinstance(value, Mapping):
            # try infering key and value types from the first non-NA items,
            # defaulting to string
            k_type = "str"
            v_type = "str"
            for k, v in value.items():
                if not is_na(k) and not is_na(v):
                    k_type = cls.python_type(k, colname)
                    v_type = cls.python_type(v, colname)
                    break
            return f"map[{k_type},{v_type}]"

        if isinstance(value, (list, tuple, set, InfoTuple)):
            types = {
                cls.python_type(item, colname) for item in value if not is_na(item)
            }
            if len(types) > 1:
                raise VembraneError(
                    f"Column {colname} has inconsistent types: {','.join(types)}"
                )
            elif not types:
                # empty input
                py_type = "str"
            else:
                py_type = types.pop()

            if is_column_values:
                return py_type
            else:
                return f"list[{py_type}]"

        return type(value).__name__

    @classmethod
    def python_type_to_arrow_type(cls, py_type: str) -> pa.DataType:
        try:
            return cls.mapping[py_type]
        except KeyError as e:
            raise VembraneError(
                "Expression yields unsupported data type for parquet format: "
                f"{py_type}. Supported types are {', '.join(cls.mapping)}"
            ) from e

    def handle_values(self, colname: str, values: list[Any]) -> Any:
        try:
            return self._to_array(colname, values)
        except (pa.ArrowTypeError, pa.ArrowInvalid) as e:
            raise VembraneError(
                f"Column '{colname}' has a value that does not fit its parquet type "
                f"'{self.python_types[colname]}': {e}\n"
                "For columns that do more than read a single field, the type is "
                "inferred from the first row group, so a larger "
                "--parquet-row-group-size may help."
            ) from e

    def _to_array(self, colname: str, values: list[Any]) -> Any:
        py_type = self.python_types[colname]
        arrow_type = self.arrow_types[colname]
        match py_type:
            case "RangeTotal":
                return pa.array(
                    [
                        None
                        if is_na(value)
                        else (value.range.start, value.range.stop, value.total)
                        for value in values
                    ],
                    type=arrow_type,
                )
            case "NumberTotal":
                return pa.array(
                    [
                        None if is_na(value) else (value.number, value.total)
                        for value in values
                    ],
                    type=arrow_type,
                )
            case "PosRange":
                return pa.array(
                    [
                        None
                        if is_na(value)
                        else (value.start or None, value.end or None)
                        for value in values
                    ],
                    type=arrow_type,
                )

        return pa.array(
            [None if is_na(value) else value for value in values], type=arrow_type
        )

    @property
    def schema(self) -> pa.Schema:
        return pa.schema(
            [
                pa.field(colname, coltype, nullable=True)
                for colname, coltype in self.arrow_types.items()
            ]
        )
