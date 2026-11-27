"""Construct and check the 207-vertex radius-two quantum coloring graph.

Run with Python 3.8 or later. Only the standard library is needed.
Use --export graph.json to save all vertices, edges, the quantum
measurements, and a proper classical five-coloring.
"""

import argparse
import json
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path

COLORS = range(3)
POSITIONS = range(4)
CELLS = list(product(COLORS, repeat=2))
BASES = (
    ((1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)),
    ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
    ((1, 1, 1, -1), (1, 1, -1, 1), (1, -1, 1, 1), (-1, 1, 1, 1)),
    ((1, 0, 1, 0), (1, 0, -1, 0), (0, 1, 0, 1), (0, 1, 0, -1)),
    ((1, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 1), (0, 0, 1, -1)),
    ((1, 0, 0, 1), (1, 0, 0, -1), (0, 1, 1, 0), (0, 1, -1, 0)),
)
ANCHORS = (0, 0, 3, 3, 3, 3)
SEEDS = ("0001", "0010", "0100", "0012", "0120", "0201")
IDENTITY = tuple(int(i == j) for i, j in product(POSITIONS, repeat=2))
ZERO = (0,) * 16


def grids_and_chart():
    grids = [
        {(a, b): tuple(a if j in (0, i) else b for j in POSITIONS) for a, b in CELLS}
        for i in (1, 2, 3)
    ]
    private = {tuple((int(c) + k) % 3 for c in seed) for seed in SEEDS for k in COLORS}
    boundary = set().union(*(g.values() for g in grids))
    assert len(boundary) == 21 and len(private) == 18 and private.isdisjoint(boundary)
    words = sorted(boundary | private)
    filters = {g[a, (a + 1) % 3]: i for i, g in enumerate(grids, 1) for a in COLORS}

    def edge(v, w):
        if any(a == b for a, b in zip(v, w)):
            return False
        return ((v not in filters or w[0] == w[filters[v]])
                and (w not in filters or v[0] == v[filters[w]]))

    return grids, words, edge, set(filters)


def grid_colorings(colors=()):
    if len(colors) == 9:
        # A row-type coloring has at most two colors in each row.
        row_type = all(len(set(colors[3 * i:3 * i + 3])) < 3 for i in COLORS)
        col_type = all(len(set(colors[i::3])) < 3 for i in COLORS)
        assert row_type != col_type
        yield int(col_type), colors
        return
    r, c = CELLS[len(colors)]
    blocked = {z for (i, j), z in zip(CELLS, colors) if i != r and j != c}
    for z in range(4):
        if z not in blocked:
            yield from grid_colorings(colors + (z,))


def verify_classical_obstruction(grids, words, edge):
    all_maps = list(grid_colorings())
    normalized = [(typ, f) for typ, f in all_maps if (f[0], f[4], f[8]) == (0, 1, 2)]
    assert len(all_maps) == 1056 and len(normalized) == 44
    cross_pairs = [(i, j) for i, (a, b) in enumerate(CELLS)
                   for j, (c, d) in enumerate(CELLS) if a != c and b != d]
    compatible = 0
    for typ, f in normalized:
        for other_type, g in all_maps:
            if all(f[i] != g[j] for i, j in cross_pairs):
                assert typ == other_type
                compatible += 1
    assert compatible == 512

    # A single color permutation normalizes the common scalar triangle
    # in a chart. For each possible grid coloring record the colors it
    # leaves available to each private vertex.
    private = sorted(set(words) - set().union(*(g.values() for g in grids)))
    options = []
    for grid in grids:
        local = []
        for typ, f in normalized:
            masks = []
            for u in private:
                blocked = {c for cell, c in zip(CELLS, f) if edge(u, grid[cell])}
                masks.append(15 ^ sum(1 << c for c in blocked))
            local.append((typ, masks))
        options.append(local)
    peers = [[i for i, v in enumerate(private) if edge(u, v)] for u in private]

    def can_extend(masks):
        if 0 in masks:
            return False
        uncolored = [i for i, m in enumerate(masks) if m != 16]
        if not uncolored:
            return True
        i = min(uncolored, key=lambda j: (bin(masks[j]).count("1"), -len(peers[j])))
        choices = masks[i]
        while choices:
            color = choices & -choices
            choices -= color
            remaining = masks[:]
            remaining[i] = 16  # outside all four color bits
            for j in peers[i]:
                remaining[j] &= ~color
            if can_extend(remaining):
                return True
        return False

    excluded = Counter()
    for (x, a), (y, b), (z, c) in product(*options):
        if (x + y + z) % 2:
            masks = [i & j & k for i, j, k in zip(a, b, c)]
            assert not can_extend(masks)
            excluded[x, y, z] += 1
    assert len(excluded) == 4 and set(excluded.values()) == {10648}
    print(f"Four-color obstruction: {sum(excluded.values())} wrong-parity grid triples excluded.")
    print(f"Bridge: all {compatible} compatible pairs among 44 x 1056 have aligned types.")


def matmul(a, b):
    return tuple(sum(a[4 * i + k] * b[4 * k + j] for k in POSITIONS)
                 for i, j in product(POSITIONS, repeat=2))


def matrix_action(a, v):
    return tuple(sum(a[4 * i + k] * v[k] for k in POSITIONS) for i in POSITIONS)


def projection(ray):
    squared_norm = sum(x * x for x in ray)
    assert 4 % squared_norm == 0
    return tuple((4 // squared_norm) * x * y for x, y in product(ray, repeat=2))


def operator_data(words):
    # All matrices stored here are four times the actual projector.
    result = {}
    for b, basis in enumerate(BASES):
        anchor = ANCHORS[b]
        order = [anchor] + [j for j in POSITIONS if j != anchor]
        projectors = [projection(basis[j]) for j in order]
        for w in words:
            labels = w if b < 3 else tuple((-c) % 3 for c in w)
            result[b, w] = tuple(
                tuple(sum(p[i] for c, p in zip(labels, projectors) if c == color) for i in range(16))
                for color in COLORS
            )
    return result


def classical_color(pvm):
    # Apply the explicit upper coloring for six square bases and n=3.
    # Ordinary colors are 0,1,2. Primed versions of 0 and 1 are 3 and 4.
    def ray_label(ray):
        target = tuple(4 * v for v in ray)
        return next((c for c, p in enumerate(pvm) if matrix_action(p, ray) == target), None)
    for i, basis in enumerate(BASES[:5]):
        c = ray_label(basis[0])
        if c is not None:
            return 3 + c if i == 2 and c < 2 else c
    u = ray_label(BASES[5][0])
    v = ray_label(BASES[5][1])
    assert u is not None and v is not None
    if u < 2:
        return u
    assert v < 2
    return 3 + v


def construct(grids, words, local_edge, filtered):
    pvms = operator_data(words)
    # A filtered word has one identical counterpart across the square.
    # Scalar and unfiltered plane cells retain their original chart tags.
    groups = {}
    ids = {}
    for tag, pvm in pvms.items():
        key = ("shared", pvm) if tag[1] in filtered else ("private", tag)
        if key not in groups:
            groups[key] = []
        groups[key].append(tag)
    copies = list(groups.values())
    for i, tags in enumerate(copies):
        for tag in tags:
            ids[tag] = i
    assert Counter(map(len, copies)) == {1: 180, 2: 27}
    edges = set()

    def join(v, w):
        a, b = ids[v], ids[w]
        assert a != b
        edges.add((min(a, b), max(a, b)))

    for b in range(6):
        for v, w in combinations(words, 2):
            if local_edge(v, w):
                join((b, v), (b, w))
    incidences = []
    for b, d in product(range(3), range(3, 6)):
        matches = []
        for i, g in enumerate(grids):
            for j, h in enumerate(grids):
                # Negation changes (0,2) to the physical pair (0,1)
                # on the column side. Its transpose reverses the named halves.
                m = pvms[b, g[0, 1]]
                if m == pvms[d, h[0, 2]]:
                    matches.append((i, j, False))
                if m == pvms[d, h[2, 0]]:
                    matches.append((i, j, True))
        assert len(matches) == 1
        i, j, flip = matches[0]
        assert flip
        incidences.append((b, d, i, j))
        for (a, c), (x, y) in product(CELLS, repeat=2):
            # Convert column input labels to physical labels, then align the halves.
            if a != (-y) % 3 and c != (-x) % 3:
                join((b, grids[i][a, c]), (d, grids[j][x, y]))
    assert all(len({i for b, d, i, j in incidences if b == k}) == 3 for k in range(3))
    assert all(len({j for b, d, i, j in incidences if d == k}) == 3 for k in range(3, 6))
    return copies, [pvms[tags[0]] for tags in copies], sorted(edges)


def verify_graph(copies, pvms, edges):
    for pvm in set(pvms):
        assert tuple(map(sum, zip(*pvm))) == tuple(4 * n for n in IDENTITY)
        for p in pvm:
            assert p == tuple(p[4 * j + i] for i, j in product(POSITIONS, repeat=2))
            assert matmul(p, p) == tuple(4 * x for x in p)
        assert all(matmul(p, q) == ZERO for p, q in combinations(pvm, 2))
    neighbors = [set() for _ in copies]
    for i, j in edges:
        neighbors[i].add(j)
        neighbors[j].add(i)
        assert all(matmul(p, q) == ZERO for p, q in zip(pvms[i], pvms[j]))

    @lru_cache(None)
    def commute(p, q):
        return all(matmul(a, b) == matmul(b, a) for a, b in product(p, q))

    within_two = set(edges)
    for i in range(len(copies)):
        within_two.update(combinations(sorted(neighbors[i]), 2))
    assert all(commute(pvms[i], pvms[j]) for i, j in within_two)
    one_chart = {w: i for i, tags in enumerate(copies) for b, w in tags if b == 0}
    assert all(one_chart[(a,) * 4] in neighbors[one_chart[(b,) * 4]]
               for a, b in combinations(COLORS, 2))
    coloring = [classical_color(m) for m in pvms]
    assert len(set(coloring)) == 5 and all(coloring[i] != coloring[j] for i, j in edges)
    assert (len(copies), len(edges), len(set(pvms)), len(within_two)) == (207, 1071, 165, 6282)
    print(f"Graph: {len(copies)} vertices, {len(edges)} edges, {len(set(pvms))} distinct measurements.")
    print(f"Radius two: checked all {len(within_two)} unordered pairs at distances one or two.")
    print("A scalar triangle and a proper five-coloring were verified.")
    return coloring


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path, metavar="PATH", help="write the graph and witnesses as JSON")
    parser.add_argument("--check-deletions", type=Path, metavar="GRAPH.json",
                        help="verify vertex-criticality against the certificate rows in GRAPH.json")
    args = parser.parse_args()
    grids, words, edge, filtered = grids_and_chart()
    verify_classical_obstruction(grids, words, edge)
    copies, pvms, edges = construct(grids, words, edge, filtered)
    coloring = verify_graph(copies, pvms, edges)
    if args.check_deletions:
        data = json.loads(args.check_deletions.read_text())
        rows = data["vertex_deletion_colorings"]
        assert len(rows) == len(copies)
        for i, row in enumerate(rows):
            assert isinstance(row, str) and len(row) == len(copies)
            assert row[i] == "x" and row.count("x") == 1 and set(row) <= set("0123x")
            assert all(row[a] != row[b] for a, b in edges)
        print(f"Vertex criticality: checked all {len(rows)} deletion four-colorings.")
    if args.export:
        data = {
            "description": "Vertex and basis IDs start at zero. Each flattened 4x4 matrix is four times the physical projection, in row-major order.",
            "bases_before_reordering": BASES,
            "anchor_indices": ANCHORS,
            "vertices": [{"id": i, "chart_copies": copies[i], "four_times_projectors": pvms[i],
                          "classical_color": coloring[i]} for i in range(len(copies))],
            "edges": edges,
        }
        args.export.write_text(json.dumps(data, separators=(",", ":")) + "\n")
        print(f"Wrote {args.export}")


if __name__ == "__main__":
    main()
