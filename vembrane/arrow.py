from typing import Any

import pyarrow as pa

from vembrane.ann_types import NA, NoValue
from vembrane.errors import VembraneError


class ArrowTypes:
    mapping = {
        "int": pa.int64(),
        "float": pa.float64(),
        "str": pa.string(),
        "bool": pa.bool_(),
        "RangeTotal": pa.struct(
            {"start": pa.int64(), "stop": pa.int64(), "total": pa.int64()}
        ),
        "NumberTotal": pa.struct({"number": pa.int64(), "total": pa.int64()}),
        "PosRange": pa.struct({"start": pa.int64(), "end": pa.int64()}),
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
        types = {type(value).__name__ for value in values if value is not NA}
        if len(types) > 1:
            raise VembraneError(
                f"Column {colname} has inconsistent types: {','.join(types)}"
            )
        elif not types:
            # empty input
            py_type = "str"
        else:
            py_type = types.pop()
        self.arrow_types[colname] = self.python_type_to_arrow_type(py_type)
        self.python_types[colname] = py_type

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
        def handle_na(value, *attrs):
            if value is None or isinstance(value, NoValue):
                return None
            else:
                for attr in attrs:
                    value = getattr(value, attr)
                return value

        py_type = self.python_types[colname]
        arrow_type = self.arrow_types[colname]
        match py_type:
            case "RangeTotal":
                return pa.array(
                    [
                        [handle_na(value, "range", "start") for value in values],
                        [handle_na(value, "range", "stop") for value in values],
                        [handle_na(value, "total") for value in values],
                    ],
                    type=arrow_type,
                )
            case "NumberTotal":
                return pa.array(
                    [
                        [handle_na(value, "number") for value in values],
                        [handle_na(value, "total") for value in values],
                    ],
                    type=arrow_type,
                )
            case "PosRange":
                return pa.array(
                    [
                        [handle_na(value, "start") for value in values],
                        [handle_na(value, "end") for value in values],
                    ],
                    type=arrow_type,
                )

        return pa.array(map(handle_na, values), type=arrow_type)

    @property
    def schema(self) -> pa.Schema:
        return pa.schema(
            [
                pa.field(colname, coltype, nullable=True)
                for colname, coltype in self.arrow_types.items()
            ]
        )
