"""Verify the small graphs of the paper and the finite statements about them.

Run with Python 3.8 or later, from the directory of this file:

    python3 verify_small_graphs.py

Only the standard library is used, and all linear algebra is exact integer
arithmetic. The adjacent file small_graphs_data.json.gz contains the
five-coloring and all vertex-deletion witnesses for the 165-vertex graph.
The program reconstructs the graphs, checks the tables of the analytical proofs and
coloring rules, and compares the 219-vertex graph with its own certificate.

Conventions. The colors of the quantum coloring are 1, 2, 3, read modulo 3,
and 0 is a fourth color. A measurement is a triple (E_1, E_2, E_3) of pairwise
orthogonal projections of C^4 with sum 1; it is stored as the three integer
matrices 4*E_x. For an orthogonal basis b = (b_0, b_1, b_2, b_3) and a word l
in {1,2,3}^4, the measurement b[l] has the projection onto b_i at outcome l_i.
Two measurements are adjacent when E_x F_x = 0 for x = 1, 2, 3.
"""
import gzip
import itertools
import json
import sys
from pathlib import Path

COLORS = (1, 2, 3)
PAIRS = [(x, y) for x in COLORS for y in COLORS if x != y]


def mod3(x):
    return (x - 1) % 3 + 1


def check(condition, message):
    if not condition:
        raise ValueError(message)


# ----------------------------------------------------------------------------
# exact linear algebra

def proj4(v):
    """four times the orthogonal projection onto the line spanned by v"""
    n = sum(a * a for a in v)
    check(4 % n == 0, "unexpected norm")
    return tuple(tuple(4 // n * a * b for b in v) for a in v)


ZERO = tuple(tuple(0 for _ in range(4)) for _ in range(4))


def add(A, B):
    return tuple(tuple(a + b for a, b in zip(r, s)) for r, s in zip(A, B))


def is_zero_product(A, B):
    columns = list(zip(*B))
    return all(sum(x * y for x, y in zip(row, col)) == 0 for row in A for col in columns)


def measurement(basis, word):
    for u, v in itertools.combinations(basis, 2):
        check(sum(a * b for a, b in zip(u, v)) == 0, "basis is not orthogonal")
    parts = []
    for x in COLORS:
        M = ZERO
        for b, label in zip(basis, word):
            if label == x:
                M = add(M, proj4(b))
        parts.append(M)
    return tuple(parts)


def orthogonal(E, F):
    return all(is_zero_product(E[x], F[x]) for x in range(3))


class Graph:
    """orthogonality graph of a list of measurements, each given as a labeling of a basis"""

    def __init__(self):
        self.names, self.meas, self.index, self.alias = [], [], {}, {}

    def add(self, name, basis, word):
        E = measurement(basis, word)
        if E not in self.index:
            self.index[E] = len(self.meas)
            self.meas.append(E)
            self.names.append(name)
        self.alias[name] = self.index[E]
        return self.index[E]

    def finish(self):
        n = len(self.meas)
        self.n = n
        self.adj = [set() for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if orthogonal(self.meas[i], self.meas[j]):
                    self.adj[i].add(j)
                    self.adj[j].add(i)
        self.edges = sorted((i, j) for i in range(n) for j in self.adj[i] if i < j)
        return self

    def induced(self, names):
        idx = [self.alias[name] for name in names]
        check(len(set(idx)) == len(idx), "repeated vertex")
        pos = {v: i for i, v in enumerate(idx)}
        return [set(pos[u] for u in self.adj[v] if u in pos) for v in idx]


# ----------------------------------------------------------------------------
# colorings

def colorings(adj, lists):
    """all proper colorings with the given color lists, vertices in index order"""
    n = len(adj)
    col = [None] * n

    def rec(v):
        if v == n:
            yield list(col)
            return
        for c in lists[v]:
            if all(col[u] != c for u in adj[v] if u < v):
                col[v] = c
                yield from rec(v + 1)
        col[v] = None

    yield from rec(0)


def forced_steps_to_conflict(adj, allowed):
    """Return the number of forced assignments before a vertex runs out of colors."""
    if any(not options for options in allowed.values()):
        return 0
    steps = 0
    while allowed:
        v = next((v for v, options in allowed.items() if len(options) == 1), None)
        check(v is not None, "forced assignments stalled on an even grid pattern")
        color, = allowed.pop(v)
        steps += 1
        for u in adj[v]:
            if u in allowed:
                allowed[u].discard(color)
                if not allowed[u]:
                    return steps
    raise ValueError("an even grid pattern extends to a proper coloring")


def find_coloring(adj, k, removed=()):
    """a proper coloring with colors 0..k-1 of the graph minus the removed vertices, or None"""
    n = len(adj)
    live = [v for v in range(n) if v not in removed]
    col = {}

    def options(v):
        return [c for c in range(k) if all(col.get(u) != c for u in adj[v])]

    def rec():
        rest = [v for v in live if v not in col]
        if not rest:
            return True
        v = min(rest, key=lambda u: (len(options(u)), -len(adj[u])))
        for c in options(v):
            if c > 1 + max(col.values(), default=-1):
                break                      # colors are interchangeable
            col[v] = c
            if rec():
                return True
            del col[v]
        return False

    return dict(col) if rec() else None


def is_proper(adj, col, removed=()):
    return all(col[u] != col[v] for u in range(len(adj)) if u not in removed
               for v in adj[u] if v not in removed)


# ----------------------------------------------------------------------------
# the configurations

# joint eigenbases of the rows and columns of the Mermin-Peres square
SQUARE = (
    ((1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)),
    ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
    ((1, 1, 1, -1), (1, 1, -1, 1), (1, -1, 1, 1), (-1, 1, 1, 1)),
    ((1, 0, 1, 0), (1, 0, -1, 0), (0, 1, 0, 1), (0, 1, 0, -1)),
    ((1, 1, 0, 0), (1, -1, 0, 0), (0, 0, 1, 1), (0, 0, 1, -1)),
    ((1, 0, 0, 1), (1, 0, 0, -1), (0, 1, 1, 0), (0, 1, -1, 0)),
)
# the relations between the planes P^i_j = span(b^i_0, b^i_j): (i, j, i', j', same plane?)
RELATIONS = (
    (1, 1, 4, 2, True), (1, 2, 5, 2, True), (1, 3, 6, 2, True),
    (2, 1, 5, 1, True), (2, 2, 4, 1, True), (2, 3, 6, 1, True),
    (3, 1, 5, 3, True), (3, 2, 4, 3, True), (3, 3, 6, 3, False),
)
# the eighteen vectors and nine contexts of Section 2
RAYS = (
    (0, 0, 0, 1), (0, 0, 1, 1), (0, 0, 1, 0), (0, 1, 0, -1),
    (0, 1, -1, 0), (0, 1, 0, 0), (1, 1, 1, 1), (1, 1, 1, -1),
    (1, 1, -1, 1), (1, -1, 1, -1), (1, -1, -1, 1), (1, -1, -1, -1),
    (1, 0, 0, 1), (1, 0, 0, -1), (1, 0, 1, 0), (1, 0, -1, 0),
    (1, 1, 0, 0), (1, -1, 0, 0),
)
CONTEXTS = (
    (0, 2, 16, 17), (0, 5, 14, 15), (1, 7, 8, 17),
    (1, 9, 10, 16), (2, 5, 12, 13), (3, 6, 9, 15),
    (3, 8, 11, 14), (4, 6, 10, 13), (4, 7, 11, 12),
)
STANDARD = SQUARE[1]


def plane_word(j, x, y):
    """the cell <x,y> of P_j = span(b_0, b_j) as a word"""
    return tuple(x if i in (0, j) else y for i in range(4))


def ray_word(j, x, y):
    """the cell <x,y> of the projection onto b_j as a word"""
    return tuple(x if i == j else y for i in range(4))


def add_constants(G, which=COLORS):
    for x in which:
        G.add(("const", x), STANDARD, (x,) * 4)


def gadget_words():
    w = [(a, k, mod3(k + 1), mod3(k + 2)) for k in COLORS for a in COLORS]
    r = [ray_word(j, k, mod3(k + 1)) for j in COLORS for k in COLORS]
    return w, r


def graph_165():
    G = Graph()
    add_constants(G)
    for i in (1, 2, 3):
        for j in COLORS:
            for x, y in PAIRS:
                G.add(("cell", i, j, x, y), SQUARE[i - 1], plane_word(j, x, y))
    w, r = gadget_words()
    for i in range(1, 7):
        for word in w:
            G.add(("w", i) + word, SQUARE[i - 1], word)
        for word in r:
            G.add(("r", i) + word, SQUARE[i - 1], word)
    return G.finish()


def graph_32():
    G = Graph()
    add_constants(G, (1, 3))
    for i in (1, 2, 3):
        a, b, c = (1, 2, 3) if i < 3 else (3, 2, 1)
        for j, (x, y) in ((1, (a, b)), (2, (a, b)), (3, (b, c))):
            G.add(("cell", i, j, x, y), SQUARE[i - 1], plane_word(j, x, y))
            G.add(("cell", i, j, y, x), SQUARE[i - 1], plane_word(j, y, x))
    for i in range(1, 7):
        a, b, c = (3, 2, 1) if i in (3, 6) else (1, 2, 3)
        G.add(("e", i), SQUARE[i - 1], (a, c, c, b))
        G.add(("e'", i), SQUARE[i - 1], (c, a, b, c))
    return G.finish()


# a vertex-critical graph with chromatic number 4 found by a computer search:
# labelings (i, word) of the bases of the square
GRAPH_28 = (
    (1, 1111), (1, 1212), (1, 1332), (1, 2121), (1, 2211), (1, 2332), (1, 3123),
    (1, 3213), (1, 3333), (2, 1212), (2, 1332), (2, 2121), (2, 2211), (2, 2332),
    (2, 3123), (2, 3223), (3, 1321), (3, 2233), (3, 2323), (3, 3112), (3, 3232),
    (4, 1332), (4, 3213), (5, 1332), (5, 3123), (5, 3213), (6, 1231), (6, 3112),
)


def graph_28():
    G = Graph()
    for i, word in GRAPH_28:
        G.add((i, word), SQUARE[i - 1], tuple(int(ch) for ch in str(word)))
    return G.finish()


def system_219():
    """constants, cells and triple measurements of the eighteen-ray system"""
    G = Graph()
    add_constants(G)
    for c, ctx in enumerate(CONTEXTS):
        basis = tuple(RAYS[r] for r in ctx)
        for i in range(4):
            for x, y in PAIRS:
                G.add(("ray", ctx[i], x, y), basis, ray_word(i, x, y))
    for c, ctx in enumerate(CONTEXTS):
        basis = tuple(RAYS[r] for r in ctx)
        for x, y in PAIRS:
            G.add(("plane", c, x, y), basis, (x, x, y, y))
        for a in COLORS:
            G.add(("m", c, 1, a), basis, (a, a, mod3(a + 1), mod3(a + 2)))
            G.add(("m", c, 2, a), basis, (a, mod3(a + 1), mod3(a + 2), mod3(a + 2)))
    return G.finish()


# ----------------------------------------------------------------------------
# the lemmas

def row_and_column_sets(cells):
    """cells[(x,y)] = color of <x,y>; rows with a row cell, columns with a column cell"""
    I = {x for x in COLORS if any(cells[(x, y)] == x for y in COLORS if y != x)}
    J = {y for y in COLORS if any(cells[(x, y)] == y for x in COLORS if x != y)}
    return I, J


def is_true(cells):
    I, J = row_and_column_sets(cells)
    check(not (I & J) and (len(I) >= 2) != (len(J) >= 2), "grid lemma fails")
    return len(I) >= 2


def check_lemmas():
    G = Graph()
    add_constants(G)
    b = STANDARD
    for x, y in PAIRS:
        G.add(("Q", x, y), b, (x, x, y, y))        # Q = span(b_0, b_1)
        G.add(("b2", x, y), b, ray_word(2, x, y))
        G.add(("b3", x, y), b, ray_word(3, x, y))
    for k in COLORS:
        G.add(("m", k), b, (k, k, mod3(k + 1), mod3(k + 2)))
    G.finish()
    const = [("const", x) for x in COLORS]
    pinned = [[x] for x in COLORS]
    free = [0, 1, 2, 3]

    # the grid lemma
    names = const + [("Q", x, y) for x, y in PAIRS]
    count = {True: 0, False: 0}
    for col in colorings(G.induced(names), pinned + [free] * 6):
        cells = {p: col[3 + i] for i, p in enumerate(PAIRS)}
        check(all(cells[p] in (p[0], p[1], 0) for p in PAIRS), "cell with a wrong color")
        I, J = row_and_column_sets(cells)
        t = is_true(cells)
        count[t] += 1
        # comparison with the row-type colorings of K_3 x K_3 (classes with two cells in a row)
        full = dict(cells)
        full.update({(x, x): x for x in COLORS})
        classes = {}
        for cell, color in full.items():
            classes.setdefault(color, []).append(cell)
        row_classes = {c for c, members in classes.items()
                       if any(u[0] == v[0] for u, v in itertools.combinations(members, 2))}
        row_type = any(all(full[(x, y)] in row_classes for x in COLORS) for y in COLORS)
        check(row_type == t, "true differs from row-type")
        rows, cols = (I, J) if t else (J, I)
        col_cells = [p for p in PAIRS if cells[p] == (p[1] if t else p[0])]
        null = sorted(p for p in PAIRS if cells[p] == 0)
        line = 0 if t else 1               # coordinate that indexes the special line
        if len(rows) == 3:
            check(not col_cells, "grid lemma (ii)")
        else:
            s = (set(COLORS) - rows).pop()
            check(null == sorted(p for p in PAIRS if p[line] == s), "grid lemma (ii)")
            check(all(p[1 - line] == s for p in col_cells), "grid lemma (ii)")
    check(count == {True: 22, False: 22}, "wrong number of grid colorings")
    print("grid lemma: 44 colorings of a grid, 22 true and 22 false, structure as stated; "
          "true coincides with row-type")

    # two orthogonal projections (the lines of b_2 and b_3)
    names = const + [("b2", x, y) for x, y in PAIRS] + [("b3", x, y) for x, y in PAIRS]
    total = both = 0
    for col in colorings(G.induced(names), pinned + [free] * 12):
        total += 1
        A = {p: col[3 + i] for i, p in enumerate(PAIRS)}
        B = {p: col[9 + i] for i, p in enumerate(PAIRS)}
        if is_true(A) and is_true(B):
            both += 1
            found = False
            for P, Q in ((A, B), (B, A)):
                for s in COLORS:
                    x, t = [v for v in COLORS if v != s]
                    found |= (P[(s, x)] == P[(s, t)] == 0 and P[(x, s)] == P[(t, s)] == s
                              and P[(x, t)] == x and P[(t, x)] == t
                              and Q[(x, s)] == Q[(t, s)] == 0 and Q[(s, x)] == Q[(s, t)] == s
                              and Q[(x, t)] == x and Q[(t, x)] == t)
            check(found, "two true orthogonal projections in an unexpected configuration")
    check((total, both) == (720, 6), "unexpected counts for two grids")
    print("two orthogonal noncomplementary projections: 720 colorings, both true in 6")

    # triples: (Q, b_2, b_3) with the measurements m_k
    names = (const + [("Q", x, y) for x, y in PAIRS] + [("b2", x, y) for x, y in PAIRS]
             + [("b3", x, y) for x, y in PAIRS] + [("m", k) for k in COLORS])
    patterns = {}
    for col in colorings(G.induced(names), pinned + [free] * 21):
        key = tuple(is_true({p: col[3 + 6 * g + i] for i, p in enumerate(PAIRS)}) for g in range(3))
        patterns[key] = patterns.get(key, 0) + 1
    check(all(sum(key) == 1 for key in patterns) and sum(patterns.values()) == 7986,
          "triple lemma fails")
    print("triple lemma: 7986 colorings of the 24-vertex gadget, exactly one true part in each")


def check_gadget(G165):
    """the 39-vertex gadget of one basis forces an odd number of true planes"""
    w, r = gadget_words()
    names = ([("const", x) for x in COLORS] + [("cell", 1, j, x, y) for j in COLORS for x, y in PAIRS]
             + [("w", 1) + word for word in w] + [("r", 1) + word for word in r])
    adj = G165.induced(names)
    check(len(adj) == 39, "gadget size")
    grid_adj = [set(u for u in adj[v] if u < 21) for v in range(21)]
    check(all(not any(21 > u >= 3 and (u - 3) // 6 != (v - 3) // 6 for u in adj[v])
              for v in range(3, 21)), "cells of different planes are adjacent")

    def truth_pattern(col):
        return tuple(is_true(dict(zip(PAIRS, col[start:start + 6]))) for start in (3, 9, 15))

    seen = {}
    immediate = propagated = 0
    for col in colorings(grid_adj, [[x] for x in COLORS] + [[0, 1, 2, 3]] * 18):
        key = truth_pattern(col)
        seen[key] = seen.get(key, 0) + 1
        if sum(key) % 2:
            continue
        allowed = {v: set(range(4)) - {col[u] for u in adj[v] if u < 21} for v in range(21, 39)}
        if forced_steps_to_conflict(adj, allowed) == 0:
            immediate += 1
        else:
            propagated += 1
    check(len(seen) == 8 and set(seen.values()) == {44 ** 3 // 8}, "grid coloring counts")
    check((immediate, propagated) == (18804, 23788), "unexpected counts of forced conflicts")

    # Evaluating all 39 words in any one of the four basis coordinates gives a proper coloring.
    words = ([(x,) * 4 for x in COLORS] + [plane_word(j, x, y) for j in COLORS for x, y in PAIRS]
             + w + r)
    odd_patterns = set()
    for pos in range(4):
        col = [word[pos] for word in words]
        check(is_proper(adj, col), "basis-coordinate witness fails")
        odd_patterns.add(truth_pattern(col))
    check(odd_patterns == {key for key in seen if sum(key) % 2}, "missing odd pattern")
    print("39-vertex gadget: 42592 even grid colorings excluded (18804 immediate conflicts, "
          "23788 by forced assignments); the four basis coordinates witness all odd patterns")


# ----------------------------------------------------------------------------
# the square

def plane4(basis, j):
    return add(proj4(basis[0]), proj4(basis[j]))


def check_square():
    four = tuple(tuple(4 * int(i == j) for j in range(4)) for i in range(4))
    for i, j, i2, j2, same in RELATIONS:
        P, Q = plane4(SQUARE[i - 1], j), plane4(SQUARE[i2 - 1], j2)
        check((P == Q) if same else (add(P, Q) == four), "relation between planes fails")
    solutions = 0
    for choice in itertools.product(range(4), repeat=6):
        t = {(i, j): choice[i - 1] in (0, j) for i in range(1, 7) for j in COLORS}
        solutions += all((t[(i, j)] == t[(i2, j2)]) == same for i, j, i2, j2, same in RELATIONS)
    check(solutions == 0, "the square has a consistent choice")
    print("square: the nine relations between the planes hold; no consistent choice of basis vectors")


# ----------------------------------------------------------------------------
# the graphs

def check_graph_165(G, data):
    check((G.n, len(G.edges)) == (165, 843), "wrong size of the 165-vertex graph")
    check(data["names_165"] == [repr(name) for name in G.names], "vertex order differs from the data file")
    col = data["coloring_165"]
    check(len(col) == 165 and set(col) <= set(range(5)) and is_proper(G.adj, col), "bad 5-coloring")
    rows = data["critical_165"]
    check(isinstance(rows, list) and len(rows) == G.n, "expected one coloring per deleted vertex")
    for v, line in enumerate(rows):
        check(isinstance(line, str) and len(line) == G.n and line[v] == "." and line.count(".") == 1,
              f"invalid deletion marker or row length for vertex {v}")
        check(set(line) <= set("0123."), f"unexpected color in deletion witness for vertex {v}")
        check(is_proper(G.adj, line, removed=(v,)), "bad 4-coloring of a vertex-deleted subgraph")
    print("165-vertex graph: 165 vertices, 843 edges; the stored 5-coloring is proper; "
          "every vertex-deleted subgraph has a proper 4-coloring")


def outcome_containing(E, vector):
    """Return the outcome whose range contains vector, if any. Each matrix is scaled by four."""
    for x, matrix in zip(COLORS, E):
        if all(sum(a * b for a, b in zip(row, vector)) == 4 * v for row, v in zip(matrix, vector)):
            return x
    return None


def check_vector_coloring(G, palette_size):
    """Check the explicit coloring using the five pairwise nonorthogonal vectors b^i_0."""
    colors = []
    residual = 0
    for E in G.meas:
        for i, basis in enumerate(SQUARE[:5], start=1):
            x = outcome_containing(E, basis[0])
            if x is not None:
                color = 4 if (palette_size == 5 and i == 3 and x == 3) else x
                break
        else:
            residual += 1
            if palette_size == 4:
                color = 4
            else:
                x = outcome_containing(E, SQUARE[5][0])
                check(x is not None, "uncovered measurement has no label on the first basis-6 vector")
                color = {1: 4, 2: 5, 3: 3}[x]
        colors.append(color)
    if palette_size == 4:
        check(residual == 1, "expected one measurement needing the fourth color")
    check(is_proper(G.adj, colors), "improper vector coloring")
    print(f"{G.n}-vertex graph: verified the explicit {palette_size}-coloring by vectors")


def check_three_colors(G, expected, name):
    check((G.n, len(G.edges)) == expected, "wrong size of " + name)
    check(find_coloring(G.adj, 3) is None, name + " is 3-colorable")
    col = find_coloring(G.adj, 4)
    check(col is not None and is_proper(G.adj, col), name + " is not 4-colorable")
    critical = all(find_coloring(G.adj, 3, removed=(v,)) is not None for v in range(G.n))
    print(f"{name}: {G.n} vertices, {len(G.edges)} edges, not 3-colorable, 4-colorable, "
          f"vertex-critical: {critical}")
    return critical


def check_local_tables(G, removed=()):
    """Reconstruct the local tables of the analytical proof for the 28-vertex graph and its 29-vertex variant."""
    live = {E: v for E, v in G.index.items() if v not in removed}
    size = len(live)
    check(size in (28, 29), "unexpected size for a parity-table check")
    # Each tuple is a column of the local table, in lexicographic bit order.
    # A string gives the colors of the words in the key. An empty string excludes the bit pattern.
    columns = {
        (1332, 3123, 3213): ("", "233", "321", "", "312", "", "", "133"),
        (1332, 3123): ("", "23", "32", "", "31", "", "", "13"),
        (1332, 3213): ("", "23", "31", "", "32", "", "", "13"),
        (1321, 3112): ("12", "21", "31", "13"),
        (1231, 3112): ("12", "31", "21", "13"),
    }
    if size == 28:
        word_lists = [(1332, 3123, 3213), (1332, 3123), (1321, 3112),
                      (1332, 3213), (1332, 3123, 3213), (1231, 3112)]
    else:
        word_lists = [(1332, 3123), (1332, 3123), (1321, 3112),
                      (1332, 3123), (1332, 3123), (1321, 3112)]

    plane_cells = {}
    words_by_basis = {}
    for i, (basis, words) in enumerate(zip(SQUARE, word_lists), start=1):
        cells = []
        for j in COLORS:
            members = []
            for x, y in PAIRS:
                v = live.get(measurement(basis, plane_word(j, x, y)))
                if v is None:
                    continue
                z = 6 - x - y
                const = live.get(measurement(basis, (z,) * 4))
                check(const in G.adj[v], "cell's omitted color is not excluded by a constant")
                members.append(v)
                cells.append((v, j, x, y))
            if members:
                check(len(members) == 1 or (len(members) == 2 and members[1] in G.adj[members[0]]),
                      "cell orientations do not define a common bit")
                plane_cells[i, j] = set(members)
        plane_indices = [j for j in COLORS if (i, j) in plane_cells]
        check(plane_indices == ([1, 2] if i in (3, 6) else [1, 2, 3]), "missing or unexpected plane")
        vertices = [live[measurement(basis, tuple(int(x) for x in str(word)))] for word in words]
        words_by_basis[i] = vertices
        local_vertices = [v for v, j, x, y in cells] + vertices
        adj = G.induced([G.names[v] for v in local_vertices])
        patterns = itertools.product((0, 1), repeat=len(plane_indices))
        check(len(columns[words]) == 2 ** len(plane_indices), "incomplete parity table")
        for bits, entry in zip(patterns, columns[words]):
            pinned = [[x if bits[j - 1] else y] for v, j, x, y in cells]
            extensions = colorings(adj, pinned + [COLORS] * len(words))
            outputs = ["".join(map(str, col[len(cells):])) for col in extensions]
            check(outputs == ([entry] if entry else []), f"table mismatch at basis {i}, bits {bits}")

    for i, j, other_i, other_j, same in RELATIONS:
        if same:
            check(plane_cells[i, j] == plane_cells[other_i, other_j], "plane occurrences use different cells")
    for v, w in zip(words_by_basis[3], words_by_basis[6]):
        check(w in G.adj[v], "missing edge between corresponding special words of bases 3 and 6")
    print(f"{size}-vertex graph: verified the local parity tables, eight shared planes and two cross edges")


def check_graph_32(G):
    check_vector_coloring(G, 4)
    check(not check_three_colors(G, (32, 72), "32-vertex graph"), "unexpected criticality")
    removed = [G.alias[name] for name in (("cell", 2, 1, 1, 2), ("cell", 3, 3, 2, 1), ("cell", 3, 3, 1, 2))]
    check_local_tables(G, removed)
    keep = [v for v in range(G.n) if v not in removed]
    pos = {v: i for i, v in enumerate(keep)}
    H = Graph()
    H.n, H.adj = len(keep), [set(pos[u] for u in G.adj[v] if u in pos) for v in keep]
    H.edges = [(u, v) for u in range(H.n) for v in H.adj[u] if u < v]
    check(check_three_colors(H, (29, 63), "its 29-vertex subgraph"), "the subgraph is not critical")


def check_F(G):
    """compare the eighteen-ray system with the certificate of the 219-vertex graph"""
    path = Path(__file__).with_name("graph_219_certificate.json.gz")
    if not path.exists():
        print("certificate of the 219-vertex graph not found: comparison skipped")
        return
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    idx = []
    for c, *labels in data["vertices"]:
        basis = tuple(RAYS[r] for r in CONTEXTS[c])
        idx.append(G.index.get(measurement(basis, tuple(x + 1 for x in labels))))
    check(G.n == 219 and None not in idx and len(set(idx)) == 219, "vertex sets differ")
    edges = {(min(idx[u], idx[v]), max(idx[u], idx[v])) for u, v in data["edges"]}
    bases = [tuple(RAYS[r] for r in ctx) for ctx in CONTEXTS]

    def is_labeling(E, basis):
        # every basis vector lies in the range of exactly one of the three projections
        return all(sum(not is_zero_product(proj4(b), E[y]) for y in range(3)) == 1 for b in basis)

    contexts = [{c for c, basis in enumerate(bases) if is_labeling(E, basis)} for E in G.meas]
    inside = {(u, v) for u, v in G.edges if contexts[u] & contexts[v]}
    check(edges == inside, "edges of F are not the orthogonal pairs inside a context")
    check((len(inside), len(G.edges) - len(inside)) == (1353, 111), "orthogonality counts for F")
    print("eighteen rays: the 219 measurements of the one-in-three system are the vertices of F, "
          "and the 1353 edges of F are the orthogonal pairs inside a common context")


def check_level_one():
    """all labelings of the nine contexts: measurements and orthogonal pairs"""
    G = Graph()
    contexts = {}
    for c, ctx in enumerate(CONTEXTS):
        basis = tuple(RAYS[r] for r in ctx)
        for word in itertools.product(COLORS, repeat=4):
            contexts.setdefault(G.add((c, word), basis, word), set()).add(c)
    G.finish()
    inside = sum(1 for u, v in G.edges if contexts[u] & contexts[v])
    check((G.n, len(G.edges), inside) == (597, 6024, 5538), "unexpected counts at the first level")
    print("first level: the labelings of the nine contexts are 597 measurements with 6024 "
          "orthogonal pairs, 5538 of them inside a common context")


def main():
    path = Path(__file__).with_name("small_graphs_data.json.gz")
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    check_square()
    check_lemmas()
    G165 = graph_165()
    check_gadget(G165)
    check_graph_165(G165, data)
    check_vector_coloring(G165, 5)
    check_graph_32(graph_32())
    G28 = graph_28()
    degree_sums = [sum(len(G28.adj[v]) for v, (basis, word) in enumerate(G28.names) if basis == i)
                   for i in range(1, 7)]
    check(degree_sums == [45, 29, 21, 8, 11, 8], "row degree sums for the 28-vertex graph")
    triangle = [G28.alias[1, word] for word in (3333, 1212, 2121)]
    check(all(v in G28.adj[u] for u, v in itertools.combinations(triangle, 2)), "missing triangle")
    check_local_tables(G28)
    check_vector_coloring(G28, 4)
    check(check_three_colors(G28, (28, 61), "28-vertex graph"), "the 28-vertex graph is not critical")
    check_F(system_219())
    check_level_one()
    print("verified: all statements checked")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
