"""Finite facts about the graphs built on the square at the first level.

Checks, with exact arithmetic and the standard library only:
  1. the six bases and the nine plane identities of Section 2.1;
  2. the first level S_lambda(K_3) for lambda = 0, 1, 2: 342 + 81 lambda
     vertices (273 distinct measurements), the worked example of Section 2.2,
     validity of the strategy and its exact compatibility radius lambda + 2;
  3. the 165 measurements of Section 6.2 and the graph F_phys with 843 edges;
  4. the 21 triangles of the appendix on the tower over these measurements:
     they are simplices of the triangulation of the 3-torus by ordered chains,
     their mod 2 boundary is the three loops, and the graph edges used by the
     sampling map of that appendix are edges of F_phys;
  5. chi(F_phys) = 5 by exhaustive search (independent of the proof), and
     chi(S_0(K_3)) = 5 when the optional SAT solver is installed.
"""
import sys
from itertools import combinations, product

from square_data import BASES, INCIDENCE, check_square, commute, measurement, orthogonal

C3 = range(3)


def tuple_mod(v):
    return tuple(x % 3 for x in v)


# ---------------------------------------------------------------- level one
def first_level(lam=0):
    """S_lambda(K_3): six basis charts and nine wires with lambda interior layers."""
    verts = []
    for i in BASES:
        for w in product(C3, repeat=4):
            if w[0] in w[1:]:
                verts.append((i, w))
    index = {v: n for n, v in enumerate(verts)}
    edges = set()
    for (i, w), (k, x) in combinations(verts, 2):
        if i == k:
            if all(a != b for a, b in zip(w, x)):
                edges.add((index[i, w], index[k, x]))

    def endpoint_word(positions, x, y):
        return tuple(x if s in positions else y for s in range(4))
    plane_word = {}
    for (i, k), (pos_i, pos_k) in INCIDENCE.items():
        layer = {}
        for x, y in product(C3, repeat=2):
            layer[0, x, y] = index[i, endpoint_word(pos_i, x, y)]
            layer[lam + 1, x, y] = index[k, endpoint_word(pos_k, x, y)]
            for ell in range(1, lam + 1):
                layer[ell, x, y] = len(verts)
                index[(i, k), ell, x, y] = len(verts)
                verts.append(((i, k), ell, x, y))
                plane_word[len(verts) - 1] = (i, endpoint_word(pos_i, x, y))
        for ell in range(lam + 2):
            for x, y, xx, yy in product(C3, repeat=4):
                if x != xx and y != yy:
                    a = layer[ell, x, y]
                    edges.add((min(a, layer[ell, xx, yy]), max(a, layer[ell, xx, yy])))
                    if ell <= lam:
                        b = layer[ell + 1, xx, yy]
                        edges.add((min(a, b), max(a, b)))
    return verts, index, edges, plane_word


def bfs_all(n, edges, cutoff):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    within = {}
    for s in range(n):
        dist = {s: 0}
        frontier = [s]
        for d in range(1, cutoff + 1):
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = d
                        nxt.append(v)
            frontier = nxt
        within[s] = dist
    return within


def noncommuting_distances(meas, dist):
    """Distances at which some pair of vertices has noncommuting measurements.

    Commutation is tested once for every pair of distinct measurements.
    """
    ids, key = {}, []
    for m in meas:
        key.append(ids.setdefault(m, len(ids)))
    rep = sorted(ids, key=ids.get)
    cache = {}
    out = set()
    for s in range(len(meas)):
        for t, d in dist[s].items():
            if t <= s or key[s] == key[t]:
                continue
            pair = (min(key[s], key[t]), max(key[s], key[t]))
            if pair not in cache:
                cache[pair] = commute(rep[pair[0]], rep[pair[1]])
            if not cache[pair]:
                out.add(d)
    return out


def check_first_level():
    for lam in (1, 2):
        verts, index, edges, plane_word = first_level(lam)
        assert len(verts) == 342 + 81 * lam
        meas = [measurement(*(plane_word[n] if n in plane_word else v)) for n, v in enumerate(verts)]
        assert all(orthogonal(meas[a], meas[b]) for a, b in edges)
        dist = bfs_all(len(verts), edges, lam + 3)
        assert noncommuting_distances(meas, dist) == {lam + 3}
        print(f"S_{lam}(K_3): {len(verts)} vertices; strategy valid and compatible exactly up to distance {lam + 2}")
    verts, index, edges, plane_word = first_level(0)
    assert not plane_word
    assert len(verts) == 342 and all(sum(1 for v in verts if v[0] == i) == 57 for i in BASES)
    meas = [measurement(i, w) for i, w in verts]
    assert len(set(meas)) == 273
    for a, b in edges:
        assert orthogonal(meas[a], meas[b])
    # worked example: row 3 = basis 3, column 3 = basis 6; r, g, b = 0, 1, 2
    r, g = 0, 1
    u = index[3, (r, g, g, r)]
    v = index[6, (r, g, g, r)]
    twin = index[6, (g, r, r, g)]
    assert (min(u, v), max(u, v)) in edges                      # equal words, adjacent
    assert meas[u] == meas[twin] and u != twin                  # different words, equal measurement
    assert (min(u, twin), max(u, twin)) not in edges
    # compatibility: exactly radius 2 for the displayed strategy
    dist = bfs_all(len(verts), edges, 3)
    assert noncommuting_distances(meas, dist) == {3}
    print(f"S_0(K_3): 342 vertices, 273 distinct measurements, {len(edges)} edges; strategy valid, "
          "2-compatible and not 3-compatible; worked example confirmed")


# ---------------------------------------------------------- 165 measurements
O = (0, 0, 0)
ALPHA = [tuple(0 if s == i else 1 for s in C3) for i in C3]
BETA = [tuple_mod(-x for x in a) for a in ALPHA]
DELTA = [tuple(2 if s == i else 0 for s in C3) for i in C3]
UPS = [tuple((s - i) % 3 for s in C3) for i in C3]
TRIPLES = [O] + ALPHA + BETA + DELTA + UPS
assert len(set(TRIPLES)) == 13


def word(k, p):
    return (k % 3,) + tuple((k + x) % 3 for x in p)


def triangles_21():
    out = []
    for i in C3:
        j, l = (i + 1) % 3, (i + 2) % 3
        out += [(O, DELTA[i], ALPHA[i]), (DELTA[i], ALPHA[i], UPS[j]),
                (ALPHA[i], UPS[j], UPS[i]), (ALPHA[i], UPS[i], BETA[i]),
                (BETA[i], UPS[i], UPS[l]), (BETA[i], UPS[l], DELTA[j]),
                (BETA[i], DELTA[j], O)]
    return out


def precedes(p, q):
    """q - p is a nonzero 0/1 vector modulo three."""
    d = tuple((b - a) % 3 for a, b in zip(p, q))
    return d != (0, 0, 0) and all(x in (0, 1) for x in d)


def check_filling():
    tris = triangles_21()
    assert len(tris) == 21 and len({frozenset(t) for t in tris}) == 21
    edge_count = {}
    for t in tris:
        assert all(v in TRIPLES for v in t)
        # an ordered chain: the three vertices can be ordered p < q < r with
        # all three differences 0/1 vectors
        orders = [(p, q, r) for p in t for q in t for r in t
                  if len({p, q, r}) == 3 and precedes(p, q) and precedes(q, r) and precedes(p, r)]
        assert len(orders) == 1, t
        for e in combinations(t, 2):
            edge_count[frozenset(e)] = edge_count.get(frozenset(e), 0) + 1
    boundary = {e for e, c in edge_count.items() if c % 2}
    loops = set()
    for i in C3:
        loops |= {frozenset((O, ALPHA[i])), frozenset((ALPHA[i], BETA[i])), frozenset((BETA[i], O))}
        # the loop u -> u * alpha_i passes O, alpha_i, beta_i by unit steps
        assert precedes(O, ALPHA[i]) and precedes(ALPHA[i], BETA[i]) and precedes(BETA[i], O)
    assert boundary == loops and len(edge_count) == 36
    assert all(c in (1, 2) for c in edge_count.values())
    print("filling: 21 triangles, 36 edges, mod 2 boundary = the three loops O, alpha_i, beta_i")
    return list(edge_count)


def build_165(simplicial_edges):
    ids, meas, words_of = {}, [], []
    chart_vertex = {}
    for i in BASES:
        for k in C3:
            for p in TRIPLES:
                m = measurement(i, word(k, p))
                if m not in ids:
                    ids[m] = len(meas)
                    meas.append(m)
                    words_of.append([])
                words_of[ids[m]].append((i, k, p))
                chart_vertex[i, k, p] = ids[m]
    n = len(meas)
    phys = {(a, b) for a, b in combinations(range(n), 2) if orthogonal(meas[a], meas[b])}

    def e(a, b):
        return (min(a, b), max(a, b))

    # edges used by the sampling map: equal triples with different shifts, and a
    # simplicial edge p < q lifted from shift k to shift k + 1
    sampled = set()
    for i in BASES:
        for p in TRIPLES:
            for k in C3:
                sampled.add(e(chart_vertex[i, k, p], chart_vertex[i, (k + 1) % 3, p]))
        for pq in simplicial_edges:
            p, q = tuple(pq)
            if not precedes(p, q):
                p, q = q, p
            assert precedes(p, q) and not precedes(q, p)
            for k in C3:
                sampled.add(e(chart_vertex[i, k, p], chart_vertex[i, (k + 1) % 3, q]))
    return n, meas, words_of, chart_vertex, phys, sampled


def check_counts(data):
    n, meas, words_of, chart_vertex, phys, sampled = data
    assert n == 165
    multiplicity = sorted(len({w[0] for w in ws}) for ws in words_of)
    assert multiplicity.count(6) == 3 and multiplicity.count(2) == 54 and multiplicity.count(1) == 108
    assert len(phys) == 843 and sampled <= phys
    print(f"165 measurements (3 scalar, 54 on planes, 108 with a rank-one outcome): "
          f"F_phys has {len(phys)} edges and contains the {len(sampled)} edges used by the sampling map")


# -------------------------------------------------------- chromatic numbers
def neighbours(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def colourable(n, edges, colours, clique):
    """Exact backtracking search with the colors of a clique fixed."""
    adj = neighbours(n, edges)
    colour = [-1] * n
    for c, v in enumerate(clique):
        colour[v] = c
    sat = [set() for _ in range(n)]
    for v in clique:
        for w in adj[v]:
            sat[w].add(colour[v])

    sys.setrecursionlimit(10000)

    def solve():
        best, best_key = -1, None
        for v in range(n):
            if colour[v] < 0:
                key = (len(sat[v]), len(adj[v]))
                if best_key is None or key > best_key:
                    best, best_key = v, key
        if best < 0:
            return True
        v = best
        for c in range(colours):
            if c in sat[v]:
                continue
            colour[v] = c
            touched = [w for w in adj[v] if c not in sat[w]]
            for w in touched:
                sat[w].add(c)
            if all(len(sat[w]) < colours or colour[w] >= 0 for w in adj[v]) and solve():
                return True
            for w in touched:
                sat[w].discard(c)
            colour[v] = -1
        return False

    ok = solve()
    if ok:
        assert all(colour[a] != colour[b] for a, b in edges)
    return ok


def check_chromatic(data):
    n, meas, words_of, chart_vertex, phys, sampled = data
    scalars = [chart_vertex[1, k, O] for k in C3]
    try:
        from pysat.solvers import Cadical153
    except Exception:                                   # pragma: no cover
        Cadical153 = None

    def sat_colourable(order, edges, colours, clique):
        if Cadical153 is None:
            return colourable(order, edges, colours, clique)
        var = lambda v, c: v * colours + c + 1
        s = Cadical153()
        for v in range(order):
            s.add_clause([var(v, c) for c in range(colours)])
        for a, b in edges:
            for c in range(colours):
                s.add_clause([-var(a, c), -var(b, c)])
        for c, v in enumerate(clique):
            s.add_clause([var(v, c)])
        return s.solve()

    assert not sat_colourable(n, phys, 4, scalars)
    assert sat_colourable(n, phys, 5, scalars)
    print("chi(F_phys) = 5 by exhaustive search (165 vertices)")
    if Cadical153 is not None:
        verts, index, edges, _ = first_level(0)
        constants = [index[6, (c, c, c, c)] for c in C3]
        assert not sat_colourable(len(verts), edges, 4, constants)
        assert sat_colourable(len(verts), edges, 5, constants)
        print("chi(S_0(K_3)) = 5 by exhaustive search (342 vertices)")


if __name__ == "__main__":
    assert check_square()
    print("square: six orthogonal bases and the nine plane identities hold")
    check_first_level()
    simplicial_edges = check_filling()
    data = build_165(simplicial_edges)
    check_counts(data)
    check_chromatic(data)
    print("check_square_graphs: all checks passed")
