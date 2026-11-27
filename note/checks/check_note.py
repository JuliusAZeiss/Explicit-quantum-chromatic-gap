"""Finite facts of the note (standard library only; python-sat is optional).

Run from this folder:   python3 -B check_note.py

The program imports the exact data of the Mermin--Peres square and helper
functions from the replay programs of the article, which lie in ../../checks.
Every fact checked here has an analytical proof in the note; the program is
an independent replay.  Section numbers refer to the note.

  1. Section 3.1, n = 3: every proper 4-coloring of K_3 x K_3 and of
     K_3 x K_3 x K_3 has exactly one dictator coordinate, identifying
     coordinates carries it along, and two colorings of K_3 x K_3 that
     satisfy the cross condition have the same one (Lemma 3.2).  For n = 4:
     every proper 6-coloring of K_4 x K_4 has exactly one dictator
     coordinate, which is the coordinate it reads, and two such colorings
     that satisfy the cross condition have the same one.  Section 3.2,
     n = 3 (needs python-sat): the graph S'_0(K_3), with fewer edges inside
     the charts, has no proper 4-coloring (Theorem 3.3).
  2. Section 4: the two-edge identity behind Proposition 4.1, on the
     22-vertex graph with its two-qubit coloring.
  3. Section 5.1: the stretched 22-vertex graphs for l = 1, 2: orders and
     sizes, chromatic number four, the edge rule, and commutation up to
     distance 2l + 1 (Proposition 5.1).
  4. Section 5.2: the graphs F_triangle (732 edges) and F_* (678 edges), the
     table of the link, contractibility of all 54 links in twenty orders of
     deletion, and chi(F_*) >= 5 by exhaustive search (Theorem 5.3).
  5. Section 5.3: the graphs F_R for R = 2, 3, 4: orders, the edge rule,
     exact compatibility radius R, and a proper 5-coloring (Theorem 5.4).
  6. Section 5.4: the presentations of weight at most two at level two and
     the count 4,330,533 (Theorem 5.6).
"""
import random
import sys
from itertools import combinations, product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "checks"))

import check_small_graph_22 as g22                                    # noqa: E402
import check_square_graphs as sq                                      # noqa: E402
from check_first_level import CELLS, cells_adjacent, proper_colourings, reads   # noqa: E402
from square_data import BASES, IDENT, INCIDENCE, madd, measurement, orthogonal, plane   # noqa: E402

C3 = range(3)


def neighbors(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def colorable(n, edges, colors, clique):
    """Exact test: SAT solver if python-sat is installed, else backtracking."""
    try:
        from pysat.solvers import Cadical153
    except Exception:                                   # pragma: no cover
        return sq.colourable(n, edges, colors, clique)
    var = lambda v, c: v * colors + c + 1
    solver = Cadical153()
    for v in range(n):
        solver.add_clause([var(v, c) for c in range(colors)])
    for a, b in edges:
        for c in range(colors):
            solver.add_clause([-var(a, c), -var(b, c)])
    for c, v in enumerate(clique):
        solver.add_clause([var(v, c)])
    return solver.solve()


# ------------------------------------------------- 1. dictator coordinates
def dictator_coordinates(f, n, length):
    """Coordinates i with a set S of n - 2 colors and a function gamma such
    that f(x) = gamma(x_i) whenever f(x) is not in S (colors 0, ..., 2n - 3)."""
    out = []
    for i in range(length):
        seen = [set() for _ in range(n)]
        for x, c in f.items():
            seen[x[i]].add(c)
        if any(all(len(seen[a] - set(s)) <= 1 for a in range(n))
               for s in combinations(range(2 * n - 2), n - 2)):
            out.append(i)
    return out


def check_dictators():
    grids = [dict(zip(CELLS, c)) for c in proper_colourings(CELLS, cells_adjacent, 4)]
    coordinate = []
    for f in grids:
        d = dictator_coordinates(f, 3, 2)
        assert len(d) == 1 and (d[0] == 0) == reads(f)[0]
        coordinate.append(d[0])
    pairs = 0
    for f, d in zip(grids, coordinate):
        forbidden = [{f[u] for u in CELLS if cells_adjacent(u, v)} for v in CELLS]
        for c in proper_colourings(CELLS, cells_adjacent, 4, forbidden=forbidden):
            assert dictator_coordinates(dict(zip(CELLS, c)), 3, 2) == [d]
            pairs += 1
    triples = list(product(C3, repeat=3))
    fixed = {(c, c, c): c for c in C3}
    cubes = proper_colourings(triples, lambda u, v: all(a != b for a, b in zip(u, v)), 4, fixed=fixed)
    for c in cubes:
        f = dict(zip(triples, c))
        d = dictator_coordinates(f, 3, 3)
        assert len(d) == 1
        for k in C3:          # identify the two coordinates other than k
            grid = {(x, y): f[tuple(x if m == k else y for m in C3)] for x, y in CELLS}
            assert dictator_coordinates(grid, 3, 2) == [0 if d[0] == k else 1]
    print(f"dictator coordinates (n = 3): unique for all {len(grids)} proper 4-colorings of K_3 x K_3 "
          f"and all {24 * len(cubes)} of K_3 x K_3 x K_3; preserved by identifying coordinates; "
          f"equal for all {pairs} pairs satisfying the cross condition")


def reads_coordinates(f, n):
    """(reads the first coordinate, reads the second coordinate) for a coloring of K_n x K_n."""
    def in_row_class(x, y):
        return any(f[x, yy] == f[x, y] for yy in range(n) if yy != y)

    def in_column_class(x, y):
        return any(f[xx, y] == f[x, y] for xx in range(n) if xx != x)

    first = any(all(in_row_class(x, y) for x in range(n)) for y in range(n))
    second = any(all(in_column_class(x, y) for y in range(n)) for x in range(n))
    return first, second


def check_dictators_k4():
    n, colors = 4, 6
    cells = [(a, a) for a in range(n)] + [(a, b) for a in range(n) for b in range(n) if a != b]
    fixed = {(a, a): a for a in range(n)}          # the diagonal is a clique
    grids = [dict(zip(cells, c)) for c in proper_colourings(cells, cells_adjacent, colors, fixed=fixed)]
    pairs = 0
    for f in grids:
        d = dictator_coordinates(f, n, 2)
        assert len(d) == 1 and reads_coordinates(f, n) == (d[0] == 0, d[0] == 1)
        forbidden = [{f[u] for u in cells if cells_adjacent(u, v)} for v in cells]
        for c in proper_colourings(cells, cells_adjacent, colors, forbidden=forbidden):
            assert dictator_coordinates(dict(zip(cells, c)), n, 2) == d
            pairs += 1
    print(f"dictator coordinates (n = 4): unique, and equal to the coordinate read, for all "
          f"{360 * len(grids)} proper 6-colorings of K_4 x K_4; equal for all {360 * pairs} pairs "
          f"satisfying the cross condition")


def check_fewer_chart_edges():
    """S'_0(K_3): chart edges only between words constant on a common set A_r."""
    try:
        import pysat.solvers  # noqa: F401
    except ImportError:                                 # pragma: no cover
        print("S'_0(K_3): skipped (python-sat is not installed)")
        return
    verts, index, edges, _ = sq.first_level(0)
    kept = set()
    for a, b in edges:
        (i, w), (k, x) = verts[a], verts[b]
        if i == k and not any(w[0] == w[r] and x[0] == x[r] for r in (1, 2, 3)):
            continue
        kept.add((a, b))
    constants = [index[6, (c, c, c, c)] for c in C3]
    assert not colorable(len(verts), kept, 4, constants)
    assert colorable(len(verts), kept, 5, constants)
    print(f"S'_0(K_3): {len(verts)} vertices, {len(kept)} of the {len(edges)} edges of S_0(K_3); "
          "chromatic number 5 by exhaustive search")


# ------------------------------------------------- 2. the two-edge identity
def check_two_edge_identity():
    names, idx, edges = g22.build()
    meas = g22.strategy(names)
    adj = neighbors(len(names), edges)
    walks = 0
    for v in range(len(names)):
        for u in adj[v]:
            for w in adj[v]:
                walks += 1
                for a, c in product(C3, repeat=2):
                    total = g22.Z4
                    for b in C3:
                        if b not in (a, c):
                            term = g22.mul(g22.mul(meas[u][a], meas[v][b]), meas[w][c])
                            total = g22.lin((1, total), (1, term))
                    assert total == g22.mul(meas[u][a], meas[w][c])
    print(f"two-edge identity: holds on all {walks} walks u, v, w of the 22-vertex graph")


# ------------------------------------------------ 3. stretched 22-vertex graph
HEXAGON = [("A", "X"), ("X", "B"), ("B", "Z"), ("Z", "C"), ("C", "Y"), ("Y", "A")]


def stretched(ell):
    """The graph of Proposition 5.1 with chains of ell diamonds."""
    names, idx, edges = g22.build()
    meas = g22.strategy(names)

    def bare(t):
        return t.rstrip("'")

    contexts = []
    for s, t in HEXAGON:
        members = [s, t] + [x + "'" for x in (s, t) if x in "ABC"]
        for r, (left, right) in enumerate(g22.ROWS, 1):
            if {bare(x) for x in left} == {s, t}:
                members.append(f"u{r}")
            if {bare(x) for x in right} == {s, t}:
                members.append(f"v{r}")
        assert len(members) == 5
        if (s, t) in (("A", "X"), ("B", "Z"), ("C", "Y")):
            members.append("o")
        contexts.append(members)
    for r in range(1, 7):
        contexts.append([f"u{r}", f"v{r}"])
    assert len(contexts) == 12
    for members in contexts:                       # each context commutes
        for x, y in combinations(members, 2):
            for p, q in product(meas[idx[x]], meas[idx[y]]):
                assert g22.mul(p, q) == g22.mul(q, p)

    verts, vmeas, color = [], [], []
    base_color = {"o": 3}
    for t in g22.TERMINALS:
        base_color[t] = 2 if t.endswith("'") else 1
    for r in range(1, 7):
        base_color[f"u{r}"], base_color[f"v{r}"] = 3, 4

    def new_vertex(tag, m, c):
        verts.append(tag)
        vmeas.append(m)
        color.append(c)
        return len(verts) - 1

    copy = {}
    for c, members in enumerate(contexts):
        for x in members:
            copy[x, c] = new_vertex((x, c), meas[idx[x]], base_color[x])
    assert len(verts) == 45
    new_edges = set()

    def add(a, b):
        assert a != b
        new_edges.add((min(a, b), max(a, b)))

    for a, b in edges:                             # every edge inside one context
        x, y = names[a], names[b]
        home = [c for c, members in enumerate(contexts) if x in members and y in members]
        assert home
        add(copy[x, home[0]], copy[y, home[0]])
    chains = 0
    for x in names:
        where = [c for c, members in enumerate(contexts) if x in members]
        m = meas[idx[x]]
        inner = sorted(set(range(1, 5)) - {base_color[x]})[:2]
        for c1, c2 in zip(where, where[1:]):
            chains += 1
            tip = copy[x, c1]
            for k in range(ell):
                nxt = copy[x, c2] if k == ell - 1 else new_vertex((x, c1, c2, k, "tip"), m, base_color[x])
                p = new_vertex((x, c1, c2, k, 1), [m[1], m[2], m[0]], inner[0])
                q = new_vertex((x, c1, c2, k, 2), [m[2], m[0], m[1]], inner[1])
                for a, b in ((tip, p), (tip, q), (p, q), (nxt, p), (nxt, q)):
                    add(a, b)
                tip = nxt
    assert chains == 23
    triangle = [copy["o", 0], copy["A", 0], copy["A'", 0]]
    return verts, vmeas, sorted(new_edges), color, triangle


def check_stretched():
    for ell in (1, 2):
        verts, meas, edges, color, triangle = stretched(ell)
        n = len(verts)
        assert n == 22 + 69 * ell and len(edges) == 42 + 115 * ell
        for a, b in edges:
            assert all(g22.mul(meas[a][c], meas[b][c]) == g22.Z4 for c in C3)
            assert color[a] != color[b]
        assert all((min(a, b), max(a, b)) in set(edges) for a, b in combinations(triangle, 2))
        assert not colorable(n, edges, 3, triangle)
        key = [tuple(tuple(tuple(row) for row in p) for p in m) for m in meas]
        cache = {}

        def commute(a, b):
            pair = (key[a], key[b])
            if pair not in cache:
                cache[pair] = all(g22.mul(p, q) == g22.mul(q, p) for p, q in product(meas[a], meas[b]))
            return cache[pair]

        dist = sq.bfs_all(n, edges, n)
        bad = {d for s in range(n) for t, d in dist[s].items() if t > s and not commute(s, t)}
        assert bad and min(bad) > 2 * ell + 1
        print(f"stretched graph, l = {ell}: {n} vertices, {len(edges)} edges, chromatic number 4; "
              f"the coloring is valid and all measurements at distance <= {min(bad) - 1} commute "
              f"(required: {2 * ell + 1})")


# --------------------------------------------- 4. F_triangle and F_* (678 edges)
def build_165():
    simplicial_edges = sq.check_filling()
    n, meas, words_of, chart_vertex, phys, tri = sq.build_165(simplicial_edges)
    per_chart = {}
    for i in BASES:
        chart = set()
        for p in sq.TRIPLES:
            for k in C3:
                a, b = chart_vertex[i, k, p], chart_vertex[i, (k + 1) % 3, p]
                chart.add((min(a, b), max(a, b)))
        for pq in simplicial_edges:
            p, q = tuple(pq)
            if not sq.precedes(p, q):
                p, q = q, p
            for k in C3:
                a, b = chart_vertex[i, k, p], chart_vertex[i, (k + 1) % 3, q]
                chart.add((min(a, b), max(a, b)))
        per_chart[i] = chart
    assert set().union(*per_chart.values()) == tri
    deleted = set()
    for i in BASES:
        for d in sq.DELTA:
            for k in C3:
                a, b = chart_vertex[i, k, d], chart_vertex[i, (k + 1) % 3, d]
                deleted.add((min(a, b), max(a, b)))
    return simplicial_edges, (n, meas, words_of, chart_vertex, phys, tri, per_chart, deleted, tri - deleted)


def check_counts(data):
    n, meas, words_of, chart_vertex, phys, tri, per_chart, deleted, star = data
    assert n == 165 and len(phys) == 843 and tri <= phys
    assert all(len(c) == 147 for c in per_chart.values())
    assert len(tri) == 732 and len(deleted) == 54 and deleted <= tri and len(star) == 678
    for i in BASES:                                   # delta words lie in one basis only
        for d in sq.DELTA:
            for k in C3:
                assert len(words_of[chart_vertex[i, k, d]]) == 1
    for i in BASES:                                   # the end grids are complete in F_*
        for a in sq.ALPHA:
            grid = {}
            for x, y in product(C3, repeat=2):
                p = tuple(0 if a[s] == 0 else (y - x) % 3 for s in C3)
                grid[x, y] = chart_vertex[i, x, p]
            for (x, y), (xx, yy) in combinations(grid, 2):
                if x != xx and y != yy:
                    u, v = grid[x, y], grid[xx, yy]
                    assert (min(u, v), max(u, v)) in star
    print(f"165 measurements: F_phys {len(phys)} edges, F_triangle {len(tri)} edges "
          f"(147 in each chart), F_* {len(star)} edges; all end grids complete")


def link_contractible(adj, u, v):
    """Removal of dominated vertices in the link of the vertex (u, v) of the box complex.

    Left vertices N(v) - u, right vertices N(u) - v; a set of left and right
    vertices is a face when all its left vertices are adjacent to all its
    right vertices.  A vertex whose neighborhood on the other side is contained
    in that of another vertex of its own side is removed.  The link is
    certified contractible when the remaining vertices span one simplex.
    """
    left = set(adj[v]) - {u}
    right = set(adj[u]) - {v}
    changed = True
    while changed:
        changed = False
        for side, other in ((left, right), (right, left)):
            for a in sorted(side):
                na = adj[a] & other
                if any(b != a and na <= (adj[b] & other) for b in side):
                    side.discard(a)
                    changed = True
                    break
            if changed:
                break
    if not left and not right:
        return False
    return all(b in adj[a] for a in left for b in right)


def check_links(data, orders=20, seed=7):
    n, meas, words_of, chart_vertex, phys, tri, per_chart, deleted, star = data
    adj = neighbors(n, tri)
    cv = lambda p, k: chart_vertex[1, k % 3, p]
    u, v = cv(sq.DELTA[0], 0), cv(sq.DELTA[0], 1)
    rows = [cv(sq.ALPHA[0], 2), cv(sq.BETA[2], 0), cv(sq.DELTA[0], 2), cv(sq.O, 2), cv(sq.UPS[1], 2)]
    cols = [cv(sq.ALPHA[0], 1), cv(sq.BETA[2], 2), cv(sq.DELTA[0], 2), cv(sq.O, 1), cv(sq.UPS[1], 1)]
    assert set(rows) == adj[v] - {u} and set(cols) == adj[u] - {v}
    table = [[int(c in adj[r]) for c in cols] for r in rows]
    assert table == [[1, 0, 0, 1, 1], [0, 1, 0, 1, 1], [0, 0, 0, 0, 0],
                     [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]
    rng = random.Random(seed)
    todo = sorted(deleted)
    for trial in range(orders):
        order = todo[:]
        if trial:
            rng.shuffle(order)
        current = set(tri)
        for a, b in order:
            adj = neighbors(n, current)
            assert link_contractible(adj, a, b) and link_contractible(adj, b, a)
            current.discard((a, b))
        assert current == star
    print(f"links: the table of Section 5.2 is confirmed; all 54 deletions have contractible links "
          f"in {orders} orders of deletion")


def check_chromatic(data):
    n, meas, words_of, chart_vertex, phys, tri, per_chart, deleted, star = data
    scalars = [chart_vertex[1, k, sq.O] for k in C3]
    assert not colorable(n, star, 4, scalars)
    assert colorable(n, phys, 5, scalars)
    print("chi(F_*) >= 5 and chi(F_phys) <= 5 by exhaustive search; hence chi = 5 for F_*, F_triangle, F_phys")


# ------------------------------------------------------------ 5. the graphs F_R
def build_FR(R, simplicial_edges):
    """Six separate charts of 39 words and nine wires with R - 2 interior layers."""
    lam = R - 2
    verts, meas, chart = [], [], {}
    for i in BASES:
        for k in C3:
            for p in sq.TRIPLES:
                chart[i, k, p] = len(verts)
                verts.append((i, k, p))
                meas.append(measurement(i, sq.word(k, p)))
    edges = set()

    def add(a, b):
        assert a != b
        edges.add((min(a, b), max(a, b)))

    for i in BASES:
        for p in sq.TRIPLES:
            for k in C3:
                add(chart[i, k, p], chart[i, (k + 1) % 3, p])
        for pq in simplicial_edges:
            p, q = tuple(pq)
            if not sq.precedes(p, q):
                p, q = q, p
            for k in C3:
                add(chart[i, k, p], chart[i, (k + 1) % 3, q])
    chart_edges = len(edges)

    def end_cell(basis, positions, x, y):
        w = tuple(x if s in positions else y for s in range(4))
        p = tuple((w[s] - w[0]) % 3 for s in (1, 2, 3))
        assert p in sq.TRIPLES
        return chart[basis, w[0], p]

    for (i, k), (pos_i, pos_k) in INCIDENCE.items():
        layers = [{(x, y): end_cell(i, pos_i, x, y) for x in C3 for y in C3}]
        for ell in range(1, lam + 1):
            missing = (ell + 1) % 2                  # omit row 0 in odd layers, row 1 in even layers
            cells = {}
            for x in C3:
                if x == missing:
                    continue
                for y in C3:
                    cells[x, y] = len(verts)
                    verts.append(((i, k), ell, x, y))
                    meas.append(meas[layers[0][x, y]])
            layers.append(cells)
        layers.append({(x, y): end_cell(k, pos_k, x, y) for x in C3 for y in C3})
        for a in range(len(layers)):
            for b in (a, a + 1):
                if b >= len(layers):
                    continue
                for (x, y), u in layers[a].items():
                    for (xx, yy), v in layers[b].items():
                        if x != xx and y != yy:
                            add(u, v)
    return verts, meas, sorted(edges), chart_edges


def check_FR(simplicial_edges, data):
    n165, meas165, words_of, chart_vertex, phys, tri, per_chart, deleted, star = data
    ident = {m: v for v, m in enumerate(meas165)}
    scalars = [chart_vertex[1, k, sq.O] for k in C3]
    five = None
    try:                                              # a proper 5-coloring of F_triangle
        from pysat.solvers import Cadical153
        var = lambda v, c: v * 5 + c + 1
        solver = Cadical153()
        for v in range(n165):
            solver.add_clause([var(v, c) for c in range(5)])
        for a, b in tri:
            for c in range(5):
                solver.add_clause([-var(a, c), -var(b, c)])
        assert solver.solve()
        model = set(x for x in solver.get_model() if x > 0)
        five = [next(c for c in range(5) if var(v, c) in model) for v in range(n165)]
    except ImportError:                               # pragma: no cover
        assert sq.colourable(n165, tri, 5, scalars)
    for R in (2, 3, 4):
        verts, meas, edges, chart_edges = build_FR(R, simplicial_edges)
        n = len(verts)
        assert n == 234 + 54 * (R - 2) and chart_edges == 6 * 147
        assert len(set(meas)) == 165
        assert all(orthogonal(meas[a], meas[b]) for a, b in edges)
        dist = sq.bfs_all(n, edges, R + 1)
        assert sq.noncommuting_distances(meas, dist) == {R + 1}
        fold = [ident[m] for m in meas]               # every vertex to its measurement
        assert all((min(fold[a], fold[b]), max(fold[a], fold[b])) in tri for a, b in edges)
        if five is not None:
            assert all(five[fold[a]] != five[fold[b]] for a, b in edges)
        print(f"F_{R}: {n} vertices, {len(edges)} edges; the coloring is valid, all measurements at "
              f"distance <= {R} commute and some at distance {R + 1} do not; the graph maps "
              f"homomorphically to F_triangle and is 5-colorable")


# --------------------------------------------- 6. weight at most two at level two
def check_level_two_count():
    templates = []
    for i in BASES:
        for r in (1, 2, 3):
            s, t = [x for x in (1, 2, 3) if x != r]
            templates.append((i, ((0, r), (s,), (t,))))
    templates.append((6, ((0, 3), (1, 2))))
    assert len(templates) == 19

    def count(n):
        """Distinct measurements of weight <= 2 over n first-level labels, 0, 1, 2 scalar."""
        seen = {frozenset({(c, IDENT)}) for c in range(3)}          # padded scalars
        for i, parts in templates:
            projections = [plane(i, p) for p in parts]
            for labels in product(range(n), repeat=len(parts)):
                internal = sum(x >= 3 for x in labels)
                if internal + (1 if internal <= 1 else 0) > 2:
                    continue
                function = {}
                for x, pr in zip(labels, projections):
                    function[x] = madd(function[x], pr) if x in function else pr
                seen.add(frozenset(function.items()))
        return len(seen)

    def formula(n):
        b = n - 3
        return (n + 27 * (n * (n - 1) - b * (b - 1)) + b * (b - 1)
                + 18 * (n * (n - 1) * (n - 2) - b * (b - 1) * (b - 2)))

    for n in (5, 6, 7, 8):
        assert count(n) == formula(n)
    assert formula(165) == 4330533
    print("level two: enumeration with 5, 6, 7 and 8 first-level labels agrees with the count of "
          "Section 5.4; for 165 labels the count is 4,330,533")


if __name__ == "__main__":
    check_dictators()
    check_dictators_k4()
    check_fewer_chart_edges()
    check_two_edge_identity()
    check_stretched()
    simplicial_edges, data = build_165()
    check_counts(data)
    check_links(data)
    check_chromatic(data)
    check_FR(simplicial_edges, data)
    check_level_two_count()
    print("check_note: all checks passed")
