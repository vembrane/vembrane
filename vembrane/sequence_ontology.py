from os import PathLike
from typing import NewType, Self, TextIO

import networkx
import obonet

DEFAULT_SEQUENCE_ONTOLOGY_FILE = "sequence_ontology.2024-06-06.obo.xz"

Id = NewType("Id", str)

_missing = object()


class CachedProperty(object):
    """
    A decorator for caching the result of a property.
    Replaces itself with the cached value after the first call.

    Inspired by https://github.com/pydanny/cached-property / django's cached_property
    """

    def __init__(self, func):
        self.__name__ = func.__name__
        self.func = func

    def __get__(self, obj, cls=None):
        if obj is None:
            return self
        value = obj.__dict__.get(self.__name__, _missing)
        if value is _missing:
            value = self.func(obj)
            obj.__dict__[self.__name__] = value
        return value


cached_property = CachedProperty


class _Container:
    """
    A container for the Sequence Ontology and related objects.
    This enables lazy-loading of the ontology.
    """

    ontology: "SequenceOntology"

    @cached_property
    def ontology(self) -> "SequenceOntology":
        return SequenceOntology.default()


_C: _Container = _Container()


class Term(str):
    """
    Representation of a Sequence Ontology term.
    Inherits from `str` and adds methods for querying the ontology.
    """

    def get_id(self) -> Id:
        """Get the ID of this term."""
        return _C.ontology.term_to_id[self]

    def parents(self) -> "Consequences":
        """Return the parent terms of this term."""
        # obo file has reversed edges,
        # so networkx descendants are ancestors and vice versa
        return Consequences(
            map(_C.ontology.get_term, _C.ontology.graph.successors(self.get_id()))
        )

    def ancestors(self) -> "Consequences":
        """
        Return the ancestor terms of this term,
        i.e. from the term to the root.
        (Parents, Grandparents, etc.)
        """
        # obo file has reversed edges,
        # so networkx descendants are ancestors and vice versa
        return Consequences(
            map(
                _C.ontology.get_term,
                networkx.descendants(_C.ontology.graph, self.get_id()),
            )
        )

    def children(self) -> "Consequences":
        """Return the child terms of this term."""
        # obo file has reversed edges,
        # so networkx ancestors are descendants and vice versa
        return Consequences(
            map(_C.ontology.get_term, _C.ontology.graph.predecessors(self.get_id()))
        )

    def descendants(self) -> "Consequences":
        """
        Return the descendant terms of this term,
        i.e. from the term to the leaves.
        (Children, Grandchildren, etc.)
        """
        # obo file has reversed edges,
        # so networkx ancestors are descendants and vice versa
        return Consequences(
            map(
                _C.ontology.get_term,
                networkx.ancestors(_C.ontology.graph, self.get_id()),
            )
        )

    def is_ancestor(self, target_term: Self) -> bool:
        """
        Check whether this term is an ancestor of the target term.
        """
        if not target_term:
            return False
        return _C.ontology.path_lengths[target_term].get(self) is not None

    def is_descendant(self, target_term: Self) -> bool:
        """
        Check whether this term is a descendant of the target term.
        """
        if not target_term:
            return False
        return _C.ontology.path_lengths[self].get(target_term) is not None

    def is_a(self, target_term: Self) -> bool:
        """
        Check whether this term is a descendant of the target term.

        This delegates to `is_descendant`, i.e. is just a synonym for convenience.
        """
        return self.is_descendant(target_term)

    def path_length(self, target_term: Self) -> int | None:
        """
        Return the path length between this term and the target term.

        This is irrespective of the direction of the path.
        Returns `None` if there is no path between the terms.
        """
        return _C.ontology.path_lengths[self].get(
            target_term
        ) or _C.ontology.path_lengths[target_term].get(self)


class Consequences(list[Term]):
    def any_is_a(self, target_term: Term) -> bool:
        return any(term.is_a(target_term) for term in self)


class SequenceOntology:
    """
    Representation of the Sequence Ontology.
    Contains the graph, term to id mapping, id to term mapping,
    as well as pre-computed path lengths between all terms.
    """

    graph: networkx.MultiDiGraph
    term_to_id: dict[Term, Id]
    id_to_term: dict[Id, Term]
    path_lengths: dict[Term, dict[Term, int]]

    def __init__(self, graph: networkx.MultiDiGraph):
        self.graph = graph
        self.term_to_id = {
            Term(data["name"]): id_ for id_, data in graph.nodes(data=True)
        }
        self.id_to_term = {v: k for k, v in self.term_to_id.items()}
        self.path_lengths = {
            self.id_to_term[source_id]: {
                self.id_to_term[target_id]: length
                for target_id, length in lengths.items()
            }
            for source_id, lengths in networkx.all_pairs_shortest_path_length(graph)
        }

    @classmethod
    def from_obo(cls, path: str | PathLike | TextIO):
        """Read an Ontology from an OBO file."""
        if isinstance(path, list) and len(path) == 1:
            graph = obonet.read_obo(path[0])
        elif isinstance(path, (str, PathLike, TextIO)):
            graph = obonet.read_obo(path)
        else:
            raise TypeError(
                f"Expected str, PathLike, or TextIO, "
                f"got {type(path)} for value {path}"
            )
        return cls(graph)

    @classmethod
    def default(cls):
        """Load the default Sequence Ontology from bundled resources."""
        import importlib.resources

        so_resource = importlib.resources.files("vembrane.resources").joinpath(
            DEFAULT_SEQUENCE_ONTOLOGY_FILE
        )
        with importlib.resources.as_file(so_resource) as so_file:
            return cls.from_obo(so_file)

    def get_id(self, term: Term) -> Id:
        """Get the ID of a term."""
        return self.term_to_id[term]

    def get_term(self, id_: Id) -> Term:
        """Get the term for an ID."""
        return self.id_to_term[id_]
