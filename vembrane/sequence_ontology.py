from typing import NewType

import networkx
import obonet

DEFAULT_SEQUENCE_ONTOLOGY_PATH = "resources/sequence_ontology.2024-06-06.obo.xz"

Term = NewType("Term", str)
Id = NewType("Id", str)


class SequenceOntology:
    graph: networkx.MultiDiGraph
    term_to_id: dict[Term, Id]
    id_to_term: dict[Id, Term]
    path_lengths: dict[Term, dict[Term, int]]

    def __init__(self, graph: networkx.MultiDiGraph):
        self.graph = graph
        self.term_to_id = {data["name"]: id_ for id_, data in graph.nodes(data=True)}
        self.id_to_term = {v: k for k, v in self.term_to_id.items()}
        self.path_lengths = {
            self.id_to_term[source_id]: {
                self.id_to_term[target_id]: length
                for target_id, length in lengths.items()
            }
            for source_id, lengths in networkx.all_pairs_shortest_path_length(graph)
        }

    @classmethod
    def from_obo(cls, path: str):
        graph = obonet.read_obo(path)
        return cls(graph)

    @classmethod
    def default(cls):
        return cls.from_obo(DEFAULT_SEQUENCE_ONTOLOGY_PATH)

    def get_id(self, term: Term) -> Id:
        return self.term_to_id[term]

    def get_term(self, id_: Id) -> Term:
        return self.id_to_term[id_]

    def get_parents(self, term: Term) -> set[Term]:
        # obo file has reversed edges,
        # so networkx descendants are ancestors and vice versa
        return networkx.descendants(self.graph, self.get_id(term))

    def get_children(self, term: Term) -> set[Term]:
        # obo file has reversed edges,
        # so networkx ancestors are descendants and vice versa
        return networkx.ancestors(self.graph, self.get_id(term))

    def is_descendant(self, source_term: Term, target_term: Term) -> bool:
        if not source_term or not target_term:
            return False
        return self.path_lengths[target_term].get(source_term) is not None

    def is_ancestor(self, source_term: Term, target_term: Term) -> bool:
        if not source_term or not target_term:
            return False
        return self.path_lengths[source_term].get(target_term) is not None

    def is_a(self, term: Term, target_term: Term) -> bool:
        return self.is_descendant(term, target_term)

    def any_is_a(self, terms: set[Term], target_term: Term) -> bool:
        return any(self.is_a(term, target_term) for term in terms)

    def path_length(self, source_term: Term, target_term: Term) -> int | None:
        return self.path_lengths[source_term].get(target_term) or self.path_lengths[
            target_term
        ].get(source_term)

    def shortest_path_length(self, terms: set[Term], target_term: Term) -> int | None:
        return min(
            (self.path_length(term, target_term) for term in terms),
            key=lambda x: x or float("inf"),
        )
