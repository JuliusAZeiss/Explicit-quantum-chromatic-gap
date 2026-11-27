"""The 22-vertex graph of Section 6.1 (standard library only).

Checks:
  1. 22 vertices, 42 edges, chromatic number four (exhaustive), and
     vertex- and edge-criticality (an additional fact, not used in the paper);
  2. the two-qubit strategy: a quantum 3-coloring in dimension four;
  3. its only 4-cycles are the twelve cycles (o, T, p, T'), and under the
     displayed 4-coloring no complete bipartite pair of shores has color sets
     {1,3} and {2,4}: this is the hypothesis of the analytical proof that the
     box complex maps equivariantly to a circle (Proposition 6.4);
  4. an independent cellular computation of w^2 = 0 on the quotient of the box
     complex (existence of a Z/4 lift of the cocycle of the double cover), with
     K_4 (w^2 != 0) and C_5 (w^2 = 0) as controls for the routine.
"""
from itertools import combinations, product

TERMINALS = ["A", "A'", "B", "B'", "C", "C'", "X", "Y", "Z"]      # ' = bar
ROWS = [(("A", "X"), ("Z", "C'")), (("A'", "X"), ("C", "Z")),
        (("X", "B"), ("Y", "C'")), (("X", "B'"), ("Y", "C")),
        (("A", "Y"), ("B", "Z")), (("A'", "Y"), ("B'", "Z"))]


def build():
    names = ["o"] + TERMINALS + [f"u{r}" for r in range(1, 7)] + [f"v{r}" for r in range(1, 7)]
    idx = {x: n for n, x in enumerate(names)}
    edges = set()

    def add(a, b):
        a, b = idx[a], idx[b]
        edges.add((min(a, b), max(a, b)))
    for t in TERMINALS:
        add("o", t)
    for s in "ABC":
        add(s, s + "'")
    for r, (left, right) in enumerate(ROWS, 1):
        add(f"u{r}", f"v{r}")
        for t in left:
            add(f"u{r}", t)
        for t in right:
            add(f"v{r}", t)
    return names, idx, sorted(edges)


def proper_colourings_exist(n, edges, k):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    colour = [-1] * n

    def go(v):
        if v == n:
            return True
        for c in range(k):
            if all(colour[w] != c for w in adj[v]):
                colour[v] = c
                if go(v + 1):
                    return True
                colour[v] = -1
        return False
    return go(0)


# ------------------------------------------------------------------ strategy
def kron(a, b):
    return [[a[i][j] * b[k][l] for j in range(2) for l in range(2)]
            for i in range(2) for k in range(2)]


def mul(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def lin(*terms):
    n = len(terms[0][1])
    return [[sum(c * m[i][j] for c, m in terms) for j in range(n)] for i in range(n)]


I2 = [[1, 0], [0, 1]]
SX = [[0, 1], [1, 0]]
SY = [[0, -1j], [1j, 0]]
SZ = [[1, 0], [0, -1]]
I4 = kron(I2, I2)
Z4 = [[0] * 4 for _ in range(4)]


def strategy(names):
    obs = {"A": kron(I2, SZ), "B": kron(SX, SX), "C": kron(SZ, SX),
           "X": kron(SX, I2), "Y": kron(SZ, I2), "Z": kron(SY, SY)}
    for s in "ABC":
        obs[s + "'"] = lin((-1, obs[s]))
    assert mul(obs["A"], obs["X"]) == mul(obs["C"], obs["Z"])
    assert mul(obs["X"], obs["B"]) == mul(obs["Y"], obs["C"])
    assert mul(obs["A"], obs["Y"]) == lin((-1, mul(obs["B"], obs["Z"])))
    meas = {"o": [I4, Z4, Z4]}
    for t in TERMINALS:
        u = obs[t]
        meas[t] = [Z4, lin((0.5, I4), (0.5, u)), lin((0.5, I4), (-0.5, u))]
    for r, (left, right) in enumerate(ROWS, 1):
        for name, (s, t) in ((f"u{r}", left), (f"v{r}", right)):
            u, v = obs[s], obs[t]
            assert mul(u, v) == mul(v, u)
            uv = mul(u, v)
            meas[name] = [lin((0.5, I4), (-0.5, uv)),
                          lin((0.25, I4), (-0.25, u), (-0.25, v), (0.25, uv)),
                          lin((0.25, I4), (0.25, u), (0.25, v), (0.25, uv))]
    for m in meas.values():
        assert lin((1, m[0]), (1, m[1]), (1, m[2])) == I4
        for a in range(3):
            assert mul(m[a], m[a]) == m[a]
    return [meas[x] for x in names]


# ---------------------------------------------------- box complex: w^2 = 0 ?
def w_squared_vanishes(n, edges):
    """Does the cocycle of the double cover B(G) -> B(G)/Z_2 lift to Z/4?

    Quotient cells: vertices = graph edges; 1-cells = pairs of edges sharing an
    endpoint; 2-cells = (vertex, three neighbors) and 4-cycles.  With the 0/1
    lift wt of w (antisymmetric under reversing a 1-cell), w^2 = Sq^1 w is
    represented by (delta wt)/2.  It vanishes iff this 2-cochain is a mod 2
    coboundary, which is a linear system over F_2.
    """
    edges = [tuple(e) for e in edges]
    eid = {e: i for i, e in enumerate(edges)}
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    one = {}
    for i, e in enumerate(edges):
        for j in range(i):
            f = edges[j]
            shared = set(e) & set(f)
            if shared:
                (v,) = shared
                one[j, i] = (len(one), int((v != e[0]) ^ (v != f[0])))

    def value(a, b):                       # oriented 1-cell a -> b
        if a < b:
            return one[a, b][1], one[a, b][0]
        return -one[b, a][1], one[b, a][0]

    rows = []
    def face(cycle):
        total, cols = 0, 0
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            val, col = value(a, b)
            total += val
            cols ^= 1 << col
        assert total % 2 == 0
        rows.append((cols, (total // 2) % 2))
    E = lambda a, b: eid[(min(a, b), max(a, b))]
    for v in range(n):
        for b1, b2, b3 in combinations(sorted(adj[v]), 3):
            face([E(v, b1), E(v, b2), E(v, b3)])
    for a, b in combinations(range(n), 2):
        common = sorted(adj[a] & adj[b])
        for c, d in combinations(common, 2):
            if a < c:
                face([E(a, c), E(c, b), E(b, d), E(d, a)])
    # Gaussian elimination over F_2 on [cols | rhs]
    pivots = {}
    for cols, rhs in rows:
        while cols:
            p = cols.bit_length() - 1
            if p in pivots:
                pc, pr = pivots[p]
                cols ^= pc
                rhs ^= pr
            else:
                pivots[p] = (cols, rhs)
                break
        else:
            if rhs:
                return False
    return True


if __name__ == "__main__":
    names, idx, edges = build()
    n = len(names)
    assert n == 22 and len(edges) == 42
    assert not proper_colourings_exist(n, edges, 3) and proper_colourings_exist(n, edges, 4)
    print("22 vertices, 42 edges, chromatic number 4")
    for e in edges:
        assert proper_colourings_exist(n, [f for f in edges if f != e], 3)
    for v in range(n):
        rest = [(a - (a > v), b - (b > v)) for a, b in edges if v not in (a, b)]
        assert proper_colourings_exist(n - 1, rest, 3)
    print("deleting any edge or any vertex leaves a 3-colorable graph")

    meas = strategy(names)
    for a, b in edges:
        for c in range(3):
            assert mul(meas[a][c], meas[b][c]) == Z4
    assert all((idx["o"], idx[t]) in edges for t in ("A", "A'")) and (idx["A"], idx["A'"]) in edges
    print("two-qubit strategy: a quantum 3-coloring in dimension four; o, A, A' is a triangle")

    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    four_cycles = set()
    for a, b in combinations(range(n), 2):
        for c, d in combinations(sorted(adj[a] & adj[b]), 2):
            four_cycles.add((frozenset((a, b)), frozenset((c, d))))
            four_cycles.add((frozenset((c, d)), frozenset((a, b))))
    unordered = {frozenset(p) for p in four_cycles}
    assert len(unordered) == 12
    private = {idx[f"{s}{r}"] for s in "uv" for r in range(1, 7)}
    for pair in unordered:
        shores = sorted(pair, key=lambda s: idx["o"] not in s)
        assert idx["o"] in shores[0] and len(shores[0] & private) == 1
        assert all(names[t] in TERMINALS for t in shores[1])
    # the displayed 4-coloring: terminals 1 or 2 (bars complementary), o and u_r 3, v_r 4
    colour = {idx["o"]: 3}
    for t in TERMINALS:
        colour[idx[t]] = 2 if t.endswith("'") else 1
    for r in range(1, 7):
        colour[idx[f"u{r}"]] = 3
        colour[idx[f"v{r}"]] = 4
    assert all(colour[a] != colour[b] for a, b in edges)
    for s1, s2 in four_cycles:
        assert not ({colour[v] for v in s1} == {1, 3} and {colour[v] for v in s2} == {2, 4})
    print("the twelve 4-cycles all pass through o; no shores colored {1,3} against {2,4}")

    k4 = [(a, b) for a, b in combinations(range(4), 2)]
    c5 = [(i, (i + 1) % 5) if i < 4 else (0, 4) for i in range(5)]
    assert not w_squared_vanishes(4, k4) and w_squared_vanishes(5, c5)
    assert w_squared_vanishes(n, edges)
    print("cellular computation: w^2 = 0 on the box complex (controls: K_4 nonzero, C_5 zero)")
    print("check_small_graph_22: all checks passed")
