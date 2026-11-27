"""Dimension three: the five-color ceiling and the graphs Q_R (Section 5.3).

Exact arithmetic in Q(sqrt 2), standard library only.  Checks:
  1. Peres' 33 rays: 57 lines and 40 orthogonal triads after completing every
     orthogonal pair, 63 wires, and no choice of exactly one line per triad;
  2. the graph Q_R for R = 1, ..., 4, with wires of max(R - 1, 1) edges:
     |V| = 1080 + 567 max(R - 2, 0), every edge is orthogonal at each outcome,
     the measurements commute up to distance max(R, 2), and some pair at
     distance max(R, 2) + 1 does not commute;
  3. the 5-coloring produced by the proof of the ceiling (Theorem 5.5(a)) is
     proper on Q_R, with every case of the proof exercised on an auxiliary graph.

The lower bound chi(Q_R) >= 5 follows from item 1 and the lemma on proper
4-colorings, which check_first_level.py enumerates.
"""
import sys
from fractions import Fraction as Fr
from itertools import combinations, permutations, product


# ------------------------------------------------------------- Q(sqrt 2)
class S2:
    """a + b sqrt(2) with rational a, b."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a, self.b = Fr(a), Fr(b)

    def __add__(self, o):
        o = lift(o)
        return S2(self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __neg__(self):
        return S2(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-lift(o))

    def __rsub__(self, o):
        return lift(o) - self

    def __mul__(self, o):
        o = lift(o)
        return S2(self.a * o.a + 2 * self.b * o.b, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__

    def inv(self):
        n = self.a * self.a - 2 * self.b * self.b
        return S2(self.a / n, -self.b / n)

    def __truediv__(self, o):
        return self * lift(o).inv()

    def __eq__(self, o):
        o = lift(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def sign(self):
        a, b = self.a, self.b
        if a >= 0 and b >= 0:
            return int(a > 0 or b > 0)
        if a <= 0 and b <= 0:
            return -1
        # opposite signs: compare a^2 with 2 b^2
        big = a * a - 2 * b * b
        if a > 0:
            return 1 if big > 0 else -1
        return -1 if big > 0 else 1

    def __repr__(self):
        return f"({self.a}+{self.b}r2)"


def lift(x):
    return x if isinstance(x, S2) else S2(x)


ZERO, ONE, R2 = S2(0), S2(1), S2(0, 1)


def dot(u, v):
    return sum((a * b for a, b in zip(u, v)), ZERO)


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def normal_form(v):
    """Projective normal form: first nonzero coordinate equal to one."""
    lead = next(c for c in v if c != ZERO)
    return tuple(c / lead for c in v)


def proj(v):
    n = dot(v, v)
    return tuple(tuple(v[i] * v[j] / n for j in range(3)) for i in range(3))


def madd(x, y):
    return tuple(tuple(x[i][j] + y[i][j] for j in range(3)) for i in range(3))


def mmul(x, y):
    return tuple(tuple(sum((x[i][k] * y[k][j] for k in range(3)), ZERO) for j in range(3))
                 for i in range(3))


Z3 = tuple(tuple(ZERO for _ in range(3)) for _ in range(3))
I3 = tuple(tuple(ONE if i == j else ZERO for j in range(3)) for i in range(3))


def is_zero(m):
    return all(c == ZERO for row in m for c in row)


def commutes(m, n):
    return all(mmul(m[a], n[b]) == mmul(n[b], m[a]) for a in range(3) for b in range(3))


def trace(m):
    return m[0][0] + m[1][1] + m[2][2]


# ------------------------------------------------------------ Peres' rays
def peres():
    seeds = [(ZERO, ZERO, ONE), (ZERO, ONE, ONE), (ZERO, ONE, R2), (ONE, ONE, R2)]
    rays = set()
    for s in seeds:
        for p in permutations(s):
            for signs in product((1, -1), repeat=3):
                rays.add(normal_form(tuple(c * e for c, e in zip(p, signs))))
    assert len(rays) == 33
    ordered = sorted(rays, key=repr)
    pairs = [(u, v) for u, v in combinations(ordered, 2) if dot(u, v) == ZERO]
    ray_triads = [t for t in combinations(ordered, 3)
                  if all(dot(a, b) == ZERO for a, b in combinations(t, 2))]
    assert len(pairs) == 72 and len(ray_triads) == 16      # 24 pairs lie in no triad of rays
    lines = set(rays)
    triads = set()
    for u, v in combinations(sorted(rays, key=repr), 2):
        if dot(u, v) == ZERO:
            w = normal_form(cross(u, v))
            lines.add(w)
            triads.add(frozenset((u, v, w)))
    lines = sorted(lines, key=repr)
    index = {l: i for i, l in enumerate(lines)}
    triads = sorted(tuple(sorted(index[l] for l in t)) for t in triads)
    assert len(lines) == 57 and len(triads) == 40
    for t in triads:
        assert all(dot(lines[a], lines[b]) == ZERO for a, b in combinations(t, 2))
    return lines, triads


def no_valuation(nlines, triads):
    """No set of lines meets every triad in exactly one line."""
    by_line = [[] for _ in range(nlines)]
    for t in triads:
        for l in t:
            by_line[l].append(t)
    state = [None] * nlines                 # None / True (chosen) / False

    def assign(l, val, trail):
        if state[l] is not None:
            return state[l] == val
        state[l] = val
        trail.append(l)
        if val:
            for t in by_line[l]:
                for m in t:
                    if m != l and not assign(m, False, trail):
                        return False
        else:
            for t in by_line[l]:
                rest = [m for m in t if state[m] is not False]
                if not rest:
                    return False
                if len(rest) == 1 and state[rest[0]] is None:
                    if not assign(rest[0], True, trail):
                        return False
        return True

    def search():
        open_triads = [t for t in triads if not any(state[l] for l in t)]
        if not open_triads:
            return True
        t = min(open_triads, key=lambda t: sum(state[l] is None for l in t))
        for l in t:
            if state[l] is None:
                trail = []
                if assign(l, True, trail) and search():
                    return True
                for m in trail:
                    state[m] = None
        return False
    return not search()


# -------------------------------------------------------------- the graph
def build_q(R, lines, triads):
    """Wires with layers 0..R, hence R edges per wire (R >= 1)."""
    verts, meas = [], []
    chart = {}
    Qp = [proj(l) for l in lines]
    for c, t in enumerate(triads):
        for w in product(range(3), repeat=3):
            chart[c, w] = len(verts)
            verts.append(("chart", c, w))
            m = []
            for a in range(3):
                acc = Z3
                for pos, l in enumerate(t):
                    if w[pos] == a:
                        acc = madd(acc, Qp[l])
                m.append(acc)
            meas.append(tuple(m))
    edges = set()

    def add(a, b):
        if a != b:
            edges.add((min(a, b), max(a, b)))
    for c in range(len(triads)):
        for w, x in combinations(list(product(range(3), repeat=3)), 2):
            if all(a != b for a, b in zip(w, x)):
                add(chart[c, w], chart[c, x])
    links = 0
    for l in range(len(lines)):
        charts = [c for c, t in enumerate(triads) if l in t]
        assert charts
        root = charts[0]
        for other in charts[1:]:
            links += 1
            layer = {}
            for end, c in ((0, root), (R, other)):
                pos = triads[c].index(l)
                for x, y in product(range(3), repeat=2):
                    w = tuple(x if p == pos else y for p in range(3))
                    layer[end, x, y] = chart[c, w]
            comp = tuple(tuple((ONE if i == j else ZERO) - Qp[l][i][j] for j in range(3))
                         for i in range(3))
            for ell in range(1, R):
                for x, y in product(range(3), repeat=2):
                    layer[ell, x, y] = len(verts)
                    verts.append(("link", l, root, other, ell, x, y))
                    m = []
                    for a in range(3):
                        acc = Z3
                        if a == x:
                            acc = madd(acc, Qp[l])
                        if a == y:
                            acc = madd(acc, comp)
                        m.append(acc)
                    meas.append(tuple(m))
            for ell in range(R + 1):
                for (x, y), (xx, yy) in product(product(range(3), repeat=2), repeat=2):
                    if x != xx and y != yy:
                        add(layer[ell, x, y], layer[ell, xx, yy])
                        if ell < R:
                            add(layer[ell, x, y], layer[ell + 1, xx, yy])
    return verts, meas, sorted(edges), links


def check_strategy(R, verts, meas, edges):
    """Commutation at distance <= R; returns the distances <= R + 1 of noncommuting pairs."""
    n = len(verts)
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
        assert all(is_zero(mmul(meas[a][c], meas[b][c])) for c in range(3))
    # identical measurements commute trivially; memoise on measurement ids
    ids, key = {}, []
    for m in meas:
        key.append(ids.setdefault(m, len(ids)))
    rep = {k: m for m, k in ids.items()}
    memo = {}
    bad_at = set()
    for s in range(n):
        dist = {s: 0}
        frontier = [s]
        for d in range(1, R + 2):
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = d
                        nxt.append(v)
            frontier = nxt
        for t, d in dist.items():
            if t <= s:
                continue
            pair = (min(key[s], key[t]), max(key[s], key[t]))
            if pair not in memo:
                memo[pair] = commutes(rep[pair[0]], rep[pair[1]])
            if not memo[pair]:
                assert d == R + 1, (R, d)
                bad_at.add(d)
    return len(ids), bad_at


# ------------------------------------------------- the ceiling, constructively
def five_colouring(meas, edges, zs=((1, 2, 4), (1, 3, 7), (2, 3, 5))):
    """The coloring from the proof of the five-color ceiling.

    Colors 4 and 5 go to the vertices whose first, respectively second, outcome
    has rank at least two.  The others are colored 1, 2, 3 using the marks of a
    generic vector z on lines.
    """
    n = len(meas)
    rank = [[trace(p) for p in m] for m in meas]
    for r in rank:
        assert all(x.b == 0 and x.a.denominator == 1 for x in r) and sum(x.a for x in r) == 3
    rank = [[int(x.a) for x in r] for r in rank]

    def line_of(p):                     # rank-one projector -> normal form of its range
        col = next(tuple(p[i][j] for i in range(3)) for j in range(3)
                   if any(p[i][j] != ZERO for i in range(3)))
        return normal_form(col)

    def kernel_line(p):                 # rank-two projector -> its kernel line
        comp = tuple(tuple((ONE if i == j else ZERO) - p[i][j] for j in range(3)) for i in range(3))
        return line_of(comp)

    for z in zs:
        zv = tuple(S2(c) for c in z)
        zz = dot(zv, zv)

        def mark(line):
            w = dot(line, zv)
            val = 3 * w * w - zz * dot(line, line)       # sign of weight - 1/3
            s = val.sign()
            return None if s == 0 else (1 if s < 0 else 2)

        colour = [None] * n
        ok = True
        triangles = {}
        for v in range(n):
            r = rank[v]
            if r[0] >= 2:
                colour[v] = 4
            elif r[1] >= 2:
                colour[v] = 5
            elif r[2] == 3:
                colour[v] = 3
            elif r[2] == 2:
                colour[v] = mark(kernel_line(meas[v][2]))
                ok &= colour[v] is not None
            else:
                ls = tuple(line_of(p) for p in meas[v])
                basis = frozenset(ls)
                order = sorted(basis, key=repr)
                perm = tuple(order.index(l) for l in ls)
                even = sum(1 for i in range(3) for j in range(i) if perm[j] > perm[i]) % 2 == 0
                triangles.setdefault((basis, even), {}).setdefault(perm, []).append(v)
        if not ok:
            continue
        for (basis, even), members in triangles.items():
            marks = {l: mark(l) for l in basis}
            if None in marks.values():
                ok = False
                break
            assert set(marks.values()) == {1, 2}       # both marks occur in a basis
            order = sorted(basis, key=repr)
            coset = [p for p in permutations(range(3))
                     if (sum(1 for i in range(3) for j in range(i) if p[j] > p[i]) % 2 == 0) == even]
            lists = {p: {1, 2, 3} - {marks[order[p[2]]]} for p in coset}
            choice = None
            for cols in permutations((1, 2, 3)):
                if all(cols[i] in lists[p] for i, p in enumerate(coset)):
                    choice = dict(zip(coset, cols))
                    break
            assert choice is not None                    # Hall's condition
            for p, vs in members.items():
                for v in vs:
                    colour[v] = choice[p]
        if ok:
            assert all(c in (1, 2, 3, 4, 5) for c in colour)
            assert all(colour[a] != colour[b] for a, b in edges)
            return colour
    raise AssertionError("no generic vector among the candidates")


def auxiliary_graph(lines, triads):
    """All three-outcome coarse-grainings of the 40 triads, all rank patterns."""
    Qp = [proj(l) for l in lines]
    seen = {}
    for t in triads:
        for w in product(range(3), repeat=3):
            m = []
            for a in range(3):
                acc = Z3
                for pos, l in enumerate(t):
                    if w[pos] == a:
                        acc = madd(acc, Qp[l])
                m.append(acc)
            seen.setdefault(tuple(m), len(seen))
    meas = sorted(seen, key=seen.get)
    edges = [(a, b) for a, b in combinations(range(len(meas)), 2)
             if all(is_zero(mmul(meas[a][c], meas[b][c])) for c in range(3))]
    return meas, edges


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    lines, triads = peres()
    assert no_valuation(len(lines), triads)
    degree = [sum(l in t for t in triads) for l in range(len(lines))]
    assert min(degree) >= 1 and sum(d - 1 for d in degree) == 63
    print("Peres: 33 rays, 57 lines, 40 triads, 63 wires, no one-in-three selection")

    aux_meas, aux_edges = auxiliary_graph(lines, triads)
    patterns = {tuple(int(trace(p).a) for p in m) for m in aux_meas}
    assert patterns == {r for r in product(range(4), repeat=3) if sum(r) == 3}
    five_colouring(aux_meas, aux_edges)
    print(f"ceiling: the proof's 5-coloring is proper on all {len(aux_meas)} coarse-grainings "
          f"of the 40 triads ({len(aux_edges)} orthogonal pairs, every rank pattern)")

    for R in (1, 2, 3, 4):
        length = max(R - 1, 1)                      # edges per wire
        verts, meas, edges, links = build_q(length, lines, triads)
        assert links == 63 and len(verts) == 1080 + 567 * max(R - 2, 0)
        # commutation up to distance length + 1 = max(R, 2), failure at length + 2
        distinct, bad_at = check_strategy(length + 1, verts, meas, edges)
        assert bad_at == {length + 2}
        colour = five_colouring(meas, edges)
        assert set(colour) <= {1, 2, 3, 4, 5}
        print(f"Q_{R}: {len(verts)} vertices, {len(edges)} edges, {distinct} distinct measurements; "
              f"compatible exactly up to distance {length + 1}; the proof's 5-coloring is proper")

    print("check_dimension_three: all checks passed")
