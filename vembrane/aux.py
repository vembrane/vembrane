import ast
from collections import defaultdict
from typing import Iterable, Any, Set


def dump(node, annotate_fields=True, include_attributes=False, indent="  "):
    """
    Return a formatted dump of the tree in *node*.  This is mainly useful for
    debugging purposes.  The returned string will show the names and the values
    for fields.  This makes the code impossible to evaluate, so if evaluation is
    wanted *annotate_fields* must be set to False.  Attributes such as line
    numbers and column offsets are not dumped by default.  If this is wanted,
    *include_attributes* can be set to True.
    """

    def _format(node, level=0):
        if isinstance(node, ast.AST):
            fields = [(a, _format(b, level)) for a, b in ast.iter_fields(node)]
            if include_attributes and node._attributes:
                fields.extend(
                    [(a, _format(getattr(node, a), level)) for a in node._attributes]
                )
            return "".join(
                [
                    node.__class__.__name__,
                    "(",
                    ", ".join(
                        ("%s=%s" % field for field in fields)
                        if annotate_fields
                        else (b for a, b in fields)
                    ),
                    ")",
                ]
            )
        elif isinstance(node, list):
            lines = ["["]
            lines.extend(
                (indent * (level + 2) + _format(x, level + 2) + "," for x in node)
            )
            if len(lines) > 1:
                lines.append(indent * (level + 1) + "]")
            else:
                lines[-1] += "]"
            return "\n".join(lines)
        return repr(node)

    if not isinstance(node, ast.AST):
        raise TypeError("expected AST, got %r" % node.__class__.__name__)
    return _format(node)


class NoopSet(set):
    def add(self, item):
        pass

    def clear(self, *args, **kwargs):
        pass

    def pop(self, *args, **kwargs):
        pass

    def remove(self, *args, **kwargs):
        pass

    def __contains__(self, *args, **kwargs):
        return True

    def copy(self) -> "NoopSet":
        return self

    def intersection(self, *s: Iterable[Any]) -> Set[Any]:
        return set(s)

    def __rand__(self, other) -> Set[Any]:
        return other


ALL = NoopSet()


class FieldLister(ast.NodeVisitor):
    def __init__(self):
        self.field_accesses = defaultdict(set)
        self.names = set()

    def visit(self, root: ast.AST):
        for node in ast.walk(root):
            if hasattr(node, "id"):
                self.names.add(node.id)
            if isinstance(node, ast.Constant):
                self.names.add(node.value)
            if isinstance(node, ast.Str):
                self.names.add(node.s)
        super().visit(root)

    def visit_Subscript(self, node: ast.Subscript) -> Any:
        # expressions of the form "INFO['DP']" or "ANN['Annotation_Impact']".
        # `node.slice` can be either Index, Slice or ExtSlice.
        # Anything that isn't Index makes this expression "complex"
        if isinstance(node.value, (ast.Name, ast.Attribute)) and isinstance(
            node.slice, ast.Index
        ):
            field = node.value.id
            key = node.slice.value
            if isinstance(key, (ast.Str, ast.Constant)):
                # python < 3.8
                if isinstance(key, ast.Str):
                    self.field_accesses[field].add(key.s)
                # python >= 3.8
                elif isinstance(key, ast.Constant):
                    self.field_accesses[field].add(key.value)
            else:
                # if the Index value is not a Str/Constant,
                # all keys of the field have to be available
                self.field_accesses[field] = ALL
        self.generic_visit(node)

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return self.__str__()
