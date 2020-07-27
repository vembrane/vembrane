import math
import re

# The builtins list below was generated with:
#    python -c 'print(
#        *sorted(o for o in dir(__builtins__) if o.islower() and not o.startswith("_")),
#        sep="\n",
#    )'
from .ann_types import NA

_builtins = {
    obj.__name__: obj
    for obj in (
        abs,
        all,
        any,
        # # ascii,
        # bin,
        bool,
        # # breakpoint,
        # # bytearray,
        # # bytes,
        # # callable,
        chr,
        # # classmethod,
        # # compile,
        # # complex,
        # # copyright,
        # # credits,
        # # delattr,
        dict,
        # # dir,
        # divmod,
        enumerate,
        # # eval,
        # # exec,
        # # exit,
        filter,
        float,
        # format,
        # frozenset,
        # # getattr,
        # # globals,
        # # hasattr,
        # # hash,
        # # help,
        # hex,
        # # id,
        # # input,
        int,
        # isinstance,
        # # issubclass,
        iter,
        len,
        # # license,
        list,
        # # locals,
        map,
        max,
        # # memoryview,
        min,
        next,
        # # object,
        # oct,
        # # open,
        ord,
        # pow,
        # # print,
        # # property,
        # # quit,
        range,
        # # repr,
        reversed,
        round,
        set,
        # # setattr,
        # # slice,
        sorted,
        # # staticmethod,
        str,
        sum,
        # # super,
        tuple,
        # # type,
        # # vars,
        zip,
    )
}

_modules = {mod.__name__: mod for mod in (re,)}

_math_exports = {
    name: mod for name, mod in vars(math).items() if not name.startswith("__")
}

_explicit_clear = {
    "__builtins__": {},
    "__builtin__": {},
    "__file__": None,
    "__name__": None,
    "__doc__": None,
    "__package__": None,
}

allowed_globals = {
    **_builtins,
    **_modules,
    **_math_exports,
    **_explicit_clear,
    "NA": NA,
}
