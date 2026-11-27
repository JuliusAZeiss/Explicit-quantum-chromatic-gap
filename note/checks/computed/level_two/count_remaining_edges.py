"""Verify the retained order and edge count after the independent-set deletion."""

from collections import Counter, defaultdict
from itertools import combinations_with_replacement
import json
from pathlib import Path

import count_physical_edges as physical


def main():
    with Path(__file__).with_name("independent_set.json").open() as stream:
        certificate = json.load(stream)
    graph = physical.square.graph_165()
    bottom = frozenset(certificate["bottom"])
    assert bottom <= set(range(graph.n))
    assert all(not (graph.adj[a] & bottom) for a in bottom)
    intersecting = set(map(tuple, certificate["intersecting"]))
    assert all(sum(x * y for x, y in zip(p, q)) > 0
               for p, q in combinations_with_replacement(intersecting, 2))

    partitions = physical.top_partitions()
    adjacency = [sum(1 << j for j in neighbors) for neighbors in graph.adj]
    shapes = {(len(p.spaces), p.needs_scalar) for p in partitions}
    data = {}
    for shape in shapes:
        by_membership = defaultdict(Counter)
        for labels, all_counts, nonscalar_counts in physical.signature_stream(adjacency, *shape):
            membership = sum(1 << i for i, v in enumerate(labels) if v in bottom)
            by_membership[membership][all_counts, nonscalar_counts] += 1
        data[shape] = by_membership

    removed = 0
    incident = 0
    total_vertices = 0
    total_directed_edges = 0
    for p in partitions:
        targets = Counter()
        for q in partitions:
            masks = tuple(sorted(sum(1 << i for i, part in enumerate(p.spaces)
                                     if physical.has_overlap(part, slot)) for slot in q.spaces))
            targets[masks, q.needs_scalar] += 1
        degree_cache = {}
        for membership, signatures in data[len(p.spaces), p.needs_scalar].items():
            top_sum = tuple(sum(part[i][j] for t, part in enumerate(p.spaces)
                                if membership & (1 << t)) for i in range(4) for j in range(4))
            delete = top_sum in intersecting
            for signature, multiplicity in signatures.items():
                if signature not in degree_cache:
                    all_counts, nonscalar_counts = signature
                    degree = 0
                    for (masks, restricted), repeats in targets.items():
                        choices = physical.injective_count(all_counts, masks)
                        if restricted:
                            choices -= physical.injective_count(nonscalar_counts, masks)
                        degree += repeats * choices
                    degree_cache[signature] = degree
                degree = degree_cache[signature]
                total_vertices += multiplicity
                total_directed_edges += degree * multiplicity
                if delete:
                    removed += multiplicity
                    incident += degree * multiplicity
    assert total_vertices == 4330533
    assert total_directed_edges == 2 * 2116320480
    retained = total_vertices - removed
    edges = total_directed_edges // 2 - incident
    print(f"Independent vertices removed: {removed:,}; incident edges: {incident:,}.")
    print(f"Retained graph: {retained:,} vertices, {edges:,} undirected edges.")


if __name__ == "__main__":
    main()
