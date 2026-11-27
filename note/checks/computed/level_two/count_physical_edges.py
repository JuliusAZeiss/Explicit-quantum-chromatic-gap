"""Exact edge count for two levels of the ordinary sparse square construction.

Run with Python 3.8+ in the workspace. The only local dependency is the
original manuscript's exact graph_165() and square bases. No solver is used.
"""

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent))   # verify_small_graphs.py lies next to this file
import verify_small_graphs as square

if hasattr(int, "bit_count"):
    popcount = int.bit_count
else:
    def popcount(x):
        return bin(x).count("1")


@dataclass(frozen=True)
class Partition:
    spaces: tuple
    needs_scalar: bool


def top_partitions():
    identity = tuple(tuple(4 * int(i == j) for j in range(4)) for i in range(4))

    def complement(p):
        return tuple(tuple(x - y for x, y in zip(irow, prow))
                     for irow, prow in zip(identity, p))

    cuts = {}
    triples = []
    for chart, basis in enumerate(square.SQUARE):
        rays = [square.proj4(x) for x in basis]
        for r in range(1, 4):
            plane = square.add(rays[0], rays[r])
            for p in (plane, rays[r]):
                cut = tuple(sorted((p, complement(p))))
                free = p == plane and chart == 5 and r == 3
                cuts[cut] = cuts.get(cut, True) and not free
            parts = (plane, *(rays[i] for i in range(1, 4) if i != r))
            triples.append(Partition(parts, True))
    assert len(cuts) == 27 and sum(not restricted for restricted in cuts.values()) == 1
    assert len({frozenset(t.spaces) for t in triples}) == len(triples) == 18
    return ([Partition((identity,), False)] +
            [Partition(cut, restricted) for cut, restricted in cuts.items()] +
            triples)


@lru_cache(None)
def has_overlap(p, q):
    # Positive semidefinite real symmetric matrices have orthogonal ranges
    # exactly when their trace inner product vanishes.
    return sum(a * b for prow, qrow in zip(p, q) for a, b in zip(prow, qrow)) != 0


def overlap_patterns(partitions):
    """Masks of source slots that constrain each target slot."""
    patterns = Counter()
    for source in partitions:
        for target in partitions:
            masks = tuple(sorted(sum(1 << i for i, p in enumerate(source.spaces)
                                     if has_overlap(p, q)) for q in target.spaces))
            assert all(masks)
            patterns[len(source.spaces), source.needs_scalar, masks, target.needs_scalar] += 1
    return patterns


def allowed_labels(n, arity, needs_scalar):
    # Indices 0,1,2 are the scalar PVMs. Enumerate by position of the first
    # scalar so the restricted ternary case costs O(n²).
    if not needs_scalar:
        yield from permutations(range(n), arity)
        return
    for first in range(arity):
        for scalar in range(3):
            other_labels = tuple(i for i in range(n) if i != scalar)
            for other in permutations(other_labels, arity - 1):
                if all(i >= 3 for i in other[:first]):
                    yield other[:first] + (scalar,) + other[first:]


def signature_stream(adj, arity, needs_scalar):
    """For each source tuple, common-neighbor counts for all subsets of slots."""
    n = len(adj)
    universe = (1 << n) - 1
    for labels in allowed_labels(n, arity, needs_scalar):
        meets = [universe]
        for a in labels:
            neighbors = adj[a]
            meets.extend([x & neighbors for x in meets])
        all_counts = tuple(popcount(x) for x in meets)
        nonscalar_counts = tuple(popcount(x & ~7) for x in meets)
        yield labels, all_counts, nonscalar_counts


def injective_count(c, masks):
    # Inclusion-exclusion on equalities of the target labels. Intersection
    # of admissible sets for target slots with masks i,j has mask i|j.
    if len(masks) == 1:
        return c[masks[0]]
    if len(masks) == 2:
        i, j = masks
        return c[i] * c[j] - c[i | j]
    i, j, k = masks
    return (c[i] * c[j] * c[k] - c[i | j] * c[k] - c[i | k] * c[j]
            - c[j | k] * c[i] + 2 * c[i | j | k])


def graph_counts(partitions, adjacency):
    patterns = overlap_patterns(partitions)
    adj = [sum(1 << j for j in neighbors) for neighbors in adjacency]
    shapes = Counter((len(p.spaces), p.needs_scalar) for p in partitions)
    data = {shape: Counter((all_counts, nonscalar) for _, all_counts, nonscalar
                          in signature_stream(adj, *shape)) for shape in shapes}
    vertices = sum(repeats * sum(data[shape].values()) for shape, repeats in shapes.items())
    directed = 0
    for (arity, restricted, masks, target_restricted), repeats in patterns.items():
        subtotal = 0
        for (all_counts, nonscalar_counts), multiplicity in data[arity, restricted].items():
            choices = injective_count(all_counts, masks)
            if target_restricted:
                choices -= injective_count(nonscalar_counts, masks)
            assert choices >= 0
            subtotal += choices * multiplicity
        directed += subtotal * repeats
    assert directed % 2 == 0
    return vertices, directed // 2


def brute_small(partitions, adjacency):
    # This independent enumeration uses individual bottom assignments rather
    # than counting signatures or inclusion-exclusion.
    vertices = [(p, lab) for p in partitions
                for lab in permutations(range(len(adjacency)), len(p.spaces))
                if not p.needs_scalar or min(lab) < 3]
    pairs = {(p, q): [(i, j) for i, x in enumerate(p.spaces)
                      for j, y in enumerate(q.spaces) if has_overlap(x, y)]
             for p in partitions for q in partitions}
    edges = sum(all(b[j] in adjacency[a[i]] for i, j in pairs[p, q])
                for (p, a), (q, b) in combinations(vertices, 2))
    return len(vertices), edges


def main():
    partitions = top_partitions()
    graph = square.graph_165()
    assert graph.n == 165 and len(graph.edges) == 843
    assert [graph.alias["const", c] for c in square.COLORS] == [0, 1, 2]
    # Full top list with three scalar labels, then a representative selection
    # containing both kinds of binary cut and the genuine ternary condition.
    for n, selected in ((3, partitions), (7, [partitions[0],
                             next(p for p in partitions if len(p.spaces) == 2 and p.needs_scalar),
                             next(p for p in partitions if len(p.spaces) == 2 and not p.needs_scalar),
                             partitions[-1]])):
        adjacency = [neighbors & set(range(n)) for neighbors in graph.adj[:n]]
        result = graph_counts(selected, adjacency)
        assert result == brute_small(selected, adjacency)
        print(f"Direct-enumeration cross-check ({n} labels): {result}", flush=True)
    v, e = graph_counts(partitions, graph.adj)
    assert v == 4330533
    print(f"Full physical graph: {v:,} vertices, {e:,} undirected edges.", flush=True)


if __name__ == "__main__":
    main()
