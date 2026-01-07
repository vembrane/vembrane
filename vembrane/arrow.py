from typing import Any

import pyarrow as pa

from vembrane.ann_types import NA, is_na
from vembrane.errors import VembraneError


class ArrowTypes:
    mapping: dict[str, pa.DataType] = {
        "int": pa.int64(),
        "float": pa.float64(),
        "str": pa.string(),
        "bool": pa.bool_(),
        "list[int]": pa.list_(pa.int64()),
        "list[float]": pa.list_(pa.float64()),
        "list[str]": pa.list_(pa.string()),
        "list[bool]": pa.list_(pa.bool_()),
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
        "Term": pa.string(),
        "list[Term]": pa.list_(pa.string()),
    }

    def __init__(self):
        self.python_types = {}
        self.arrow_types = {}

    def infer(self, colname: str, values: list[Any]) -> None:
        if all(value is NA for value in values):
            raise VembraneError(
                f"Column {colname} is NA only in the first {len(values)} rows. "
                "This is unsupported."
            )

        py_type = self.python_type(values, colname, is_column_values=True)
        self.arrow_types[colname] = self.python_type_to_arrow_type(py_type)
        self.python_types[colname] = py_type

    @classmethod
    def python_type(
        cls, value: Any, colname: str, is_column_values: bool = False
    ) -> str:
        if isinstance(value, (list, tuple)):
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

        if isinstance(value, str):
            return "str"

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
