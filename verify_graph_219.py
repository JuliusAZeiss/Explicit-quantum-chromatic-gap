"""Verify a 219-vertex graph with classical chromatic number 5 and quantum chromatic number 3.

Run with Python 3.8 or later. The adjacent gzip file contains JSON arrays named
vertices, edges, coloring and chains. Vertices are represented by a context
index followed by four labels; the smallest representative is used at a gluing.
An edge is a pair of zero-based vertex indices. A chain is a list of simplices,
each a list of directed edges [tail, head]; all chain operations are mod 2.

The chain certificate lives in the box complex B_edge of the graph, whose
vertices are the arcs (directed edges). It uses a vertex c_0 and chains with
boundary(c_i) = c_{i-1} + T(c_{i-1}), where T reverses every arc. Projecting
to the quotient by T, the connecting homomorphism of the double covering
sends [c_i] to [c_{i-1}]. Three applications give a point, so the cube of
the covering class is nonzero. A 4-coloring would give an equivariant map
to S^2 and force that cube to vanish (the lemma on chain certificates in the
appendix on computer checks).
"""
import gzip
import itertools
import json
import sys
from pathlib import Path

CONTEXTS = (
    (0, 2, 16, 17), (0, 5, 14, 15), (1, 7, 8, 17),
    (1, 9, 10, 16), (2, 5, 12, 13), (3, 6, 9, 15),
    (3, 8, 11, 14), (4, 6, 10, 13), (4, 7, 11, 12),
)
RAYS = (
    (0, 0, 0, 1), (0, 0, 1, 1), (0, 0, 1, 0), (0, 1, 0, -1),
    (0, 1, -1, 0), (0, 1, 0, 0), (1, 1, 1, 1), (1, 1, 1, -1),
    (1, 1, -1, 1), (1, -1, 1, -1), (1, -1, -1, 1), (1, -1, -1, -1),
    (1, 0, 0, 1), (1, 0, 0, -1), (1, 0, 1, 0), (1, 0, -1, 0),
    (1, 1, 0, 0), (1, -1, 0, 0),
)


def check(condition, message):
    if not condition:
        raise ValueError(message)


def vertex_key(context, labels):
    a, b, c, d = labels
    counts = [labels.count(i) for i in range(3)]
    if 4 in counts:
        return ("constant", a)
    if 3 in counts:
        majority = counts.index(3)
        minority = counts.index(1)
        return ("ray", CONTEXTS[context][labels.index(minority)], minority, majority)
    if ((a == b and c == d)
            or labels == (a, a, (a + 1) % 3, (a + 2) % 3)
            or labels == (a, (a + 1) % 3, (a + 2) % 3, (a + 2) % 3)):
        return ("block", context, labels)
    return None


def make_graph():
    """Enumerate the induced graph independently from its nine labeled blocks."""
    classes = {}
    blocks = [[] for _ in CONTEXTS]
    for context, block in enumerate(blocks):
        for labels in itertools.product(range(3), repeat=4):
            key = vertex_key(context, labels)
            if key is not None:
                classes.setdefault(key, []).append((context, *labels))
                block.append((labels, key))
    vertices = sorted(min(reps) for reps in classes.values())
    index = {rep: i for i, rep in enumerate(vertices)}
    ids = {key: index[min(reps)] for key, reps in classes.items()}
    edges = set()
    for block in blocks:
        for (a, key_a), (b, key_b) in itertools.combinations(block, 2):
            if all(x != y for x, y in zip(a, b)):
                edges.add(tuple(sorted((ids[key_a], ids[key_b]))))
    return vertices, edges, {ids[key]: reps for key, reps in classes.items()}


def product(a, b):
    return tuple(tuple(sum(x*y for x, y in zip(row, col)) for col in zip(*b))
                 for row in a)


def check_quantum(representatives, edges):
    # All matrices here represent four times an orthogonal projector.
    projectors = []
    for v in RAYS:
        norm = sum(x*x for x in v)
        check(4 % norm == 0, "nonintegral projector scaling")
        projectors.append(tuple(tuple(x*y*(4//norm) for y in v) for x in v))
    zero = ((0,)*4,)*4
    measurements = {}
    for vertex, reps in representatives.items():
        for context, *labels in reps:
            pvm = tuple(tuple(tuple(sum(projectors[r][i][j]
                                        for r, label in zip(CONTEXTS[context], labels)
                                        if label == k) for j in range(4)) for i in range(4))
                        for k in range(3))
            check(measurements.setdefault(vertex, pvm) == pvm, "inconsistent glued measurement")
        check(all(sum(p[i][j] for p in pvm) == 4*(i == j)
                  for i in range(4) for j in range(4)), "incomplete measurement")
        for p in pvm:
            check(tuple(zip(*p)) == p, "nonsymmetric outcome")
            check(product(p, p) == tuple(tuple(4*x for x in row) for row in p), "nonprojection")
        check(all(product(p, q) == zero for p, q in itertools.combinations(pvm, 2)),
              "nonorthogonal outcomes")
    for a, b in edges:
        left, right = measurements[a], measurements[b]
        check(all(product(p, q) == zero for p, q in zip(left, right)), "monochromatic quantum edge")
        check(all(product(p, q) == product(q, p) for p in left for q in right),
              "incompatible neighboring measurements")


def check_chains(raw, n, edges):
    directed = edges | {(b, a) for a, b in edges}
    check(len(raw) == 4, "wrong number of chains")
    used, previous = set(), set()
    for dimension, raw_chain in enumerate(raw):
        chain = set()
        for raw_simplex in raw_chain:
            simplex = frozenset(tuple(arc) for arc in raw_simplex)
            check(len(simplex) == len(raw_simplex) == dimension + 1, "wrong simplex dimension")
            check(simplex not in chain, "duplicate simplex in chain")
            check(all((a, d) in directed for a, b in simplex for c, d in simplex),
                  "invalid arc-complex simplex")
            chain.add(simplex)
            used.update(vertex for arc in simplex for vertex in arc)
        if dimension:
            boundary = set()
            for simplex in chain:
                for arc in simplex:
                    boundary.symmetric_difference_update([simplex - {arc}])
            reversed_chain = {frozenset((b, a) for a, b in simplex) for simplex in previous}
            check(boundary == previous ^ reversed_chain, "wrong chain boundary")
        else:
            check(len(chain) == 1, "c_0 must be a single arc")
        previous = chain
    check(used == set(range(n)), "chain support differs from graph vertices")
    sizes = tuple(map(len, raw))
    check(sizes == (1, 3, 684, 5832), "unexpected chain sizes")
    return sizes


def main(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    vertices, edges, representatives = make_graph()
    check((len(vertices), len(edges)) == (219, 1353), "wrong graph size")
    check([tuple(v) for v in data["vertices"]] == vertices, "vertex list mismatch")
    check([tuple(e) for e in data["edges"]] == sorted(edges), "edge list mismatch")
    coloring = data["coloring"]
    check(len(coloring) == 219 and set(coloring) == set(range(5)), "invalid color labels")
    check(all(coloring[a] != coloring[b] for a, b in edges), "improper 5-coloring")
    check_quantum(representatives, edges)
    constants = [i for i, rep in enumerate(vertices) if len(set(rep[1:])) == 1]
    check(len(constants) == 3 and all(e in edges for e in itertools.combinations(constants, 2)),
          "constant triangle missing")
    sizes = check_chains(data["chains"], len(vertices), edges)
    print("Verified 219 vertices, 1353 edges, a proper 5-coloring, and a quantum 3-coloring in dimension 4.")
    print("All neighboring measurements commute. The three constant vertices form a triangle.")
    print("Verified all arc simplices and mod-2 boundaries. Chain sizes:", sizes)
    print("Therefore the classical chromatic number is 5 and the quantum chromatic number is 3.")


if __name__ == "__main__":
    default_path = Path(__file__).with_name("graph_219_certificate.json.gz")
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else default_path)
