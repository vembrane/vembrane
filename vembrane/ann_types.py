from typing import Union, Iterable, List


class NoValue:
    def __lt__(self, other):
        return False

    def __gt__(self, other):
        return False

    def __le__(self, other):
        return False

    def __ge__(self, other):
        return False

    def __eq__(self, other):
        return False


NA = NoValue()


class InfoTuple:
    """A container that lazily evaluates None to NA in case of access."""

    def __init__(self, values):
        self.values = values

    def __getitem__(self, spec):
        values = self.values.__getitem__(spec)
        if isinstance(values, tuple) and any(v is None for v in values):
            return tuple((NA if v is None else v) for v in values)
        if values is None:
            return NA
        return values


IntFloatStr = Union[int, float, str]
AutoType = Union[IntFloatStr, "AutoType", Iterable["AutoType"]]


def try_auto_type(value: str) -> AutoType:
    if not value:
        return None
    for t in (int, float):  # try int and float first
        try:
            return t(value)
        except ValueError:
            continue
    for sep in ("/", "&"):  # if that didn't work, try splitting at "/" or "&"
        if sep in value:  # check if the value contains a separator
            # and re-try auto typing on every element of the value split by sep
            values = list(map(try_auto_type, map(str.strip, value.split(sep))))
            if values and all(values):  # if none of these valeus is None
                return values  # return the list of autotyped values
            else:
                return str(value)  # otherwise return the original string
        else:  # try the next separator
            continue
    return str(value)  # if all else fails, return the original string


def type_ann(
    key: str, value: str
) -> Union[IntFloatStr, Iterable["IntFloatStr"], NoValue]:
    if value:
        return KNOWN_ANN_TYPE_MAP.get(key, try_auto_type)(value)
    else:
        return NA


def type_int_pair(value: str) -> Union[List[int], NoValue]:
    if value:
        return [int(v.strip()) for v in value.split("/")]
    else:
        return NA


def type_info(value):
    if value is None:
        return NA
    if isinstance(value, tuple):
        return InfoTuple(value)
    return value


KNOWN_ANN_TYPE_MAP = {
    "Allele": str,
    "Annotation": str,
    "Annotation_Impact": str,
    "Gene_Name": str,
    "Gene_ID": str,
    "Feature_Type": str,
    "Feature_ID": str,
    "Transcript_BioType": str,
    "Rank": str,
    "HGVS.c": str,
    "HGVS.p": str,
    "cDNA.pos / cDNA.length": type_int_pair,
    "CDS.pos / CDS.length": type_int_pair,
    "AA.pos / AA.length": type_int_pair,
    "Distance": try_auto_type,
    "ERRORS / WARNINGS / INFO": lambda x: [v.strip() for v in x.split("/")],
    "CLIN_SIG": lambda x: [v.strip() for v in x.split("&")],
}
