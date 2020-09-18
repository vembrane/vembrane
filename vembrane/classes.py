from typing import Iterator, List, Dict, Any

from vembrane.errors import (
    UnknownAnnotation,
    UnknownInfoField,
    InvalidExpression,
    UnknownFormatField,
    UnknownSample,
    VembraneError,
)


class Sample:
    def __init__(self, record_idx: int, sample: str, format_data: Dict[str, Any]):
        self._record_idx = record_idx
        self._sample = sample
        self._data = format_data

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError as ke:
            raise UnknownFormatField(self._record_idx, self._sample, ke)


class Format:
    def __init__(self, record_idx: int, sample_formats: Dict[Sample, Dict[str, Any]]):
        self._record_idx = record_idx
        self._sample_formats = sample_formats

    def __getitem__(self, item):
        try:
            return self._sample_formats[item]
        except KeyError as ke:
            raise UnknownSample(self._record_idx, ke)


class Info:
    def __init__(self, record_idx: int, info_dict: Dict[str, Dict[str, Any]]):
        self._record_idx = record_idx
        self._info_dict = info_dict

    def __getitem__(self, item):
        try:
            return self._info_dict[item]
        except KeyError as ke:
            raise UnknownInfoField(self._record_idx, ke)


class Annotation:
    def __init__(self, record_idx: int, annotation_data: Dict[str, Any]):
        self._record_idx = record_idx
        self._data = annotation_data

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError as ke:
            raise UnknownAnnotation(self._record_idx, ke)
