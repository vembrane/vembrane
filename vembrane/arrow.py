from typing import Any

import pyarrow as pa

from vembrane.ann_types import RangeTotal
from vembrane.errors import VembraneError


class ArrowTypes:
    mapping = {
        "int": pa.int64(),
        "float": pa.float64(),
        "str": pa.string(),
        "bool": pa.bool_(),
        "RangeTotal": pa.string(),
    }

    @classmethod
    def python_type_to_arrow_type(cls, value: Any) -> pa.DataType:
        py_type = type(value).__name__
        try:
            return cls.mapping[py_type]
        except KeyError as e:
            raise VembraneError(
                "Expression yields unsupported data type for parquet format: "
                f"{py_type}. Supported types are {', '.join(cls.mapping)}"
            ) from e

    @classmethod
    def handle_value(cls, value: Any) -> Any:
        if isinstance(value, RangeTotal):
            return str(value)
        return value
