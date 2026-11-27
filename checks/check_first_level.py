"""Finite facts behind Section 2 and Appendix A of the paper.

Exact arithmetic, standard library only.

  1. The description of the square by two-qubit Pauli observables: every basis
     is a common eigenbasis of its row or column, every plane P_e is the
     eigenspace for +1 of the observable at its edge of K_3,3, and the products of
     the rows and columns are I, I, I and I, I, -I.
  2. The parity lemma of Section 2.1, by enumeration of all 4^6 choices of one
     position in each basis.
  3. The lemma on proper 4-colorings of powers of the triangle (Section 2.4,
     proved in Appendix A.1), by enumeration:
       (a) every proper 4-coloring of K_3 x K_3 reads exactly one coordinate;
       (b) two such colorings f, g with f(u) != g(v) whenever u, v differ in
           both coordinates read the same coordinate;
       (c) for a proper 4-coloring of K_3 x K_3 x K_3 exactly one of the three
           restrictions f_i reads the first coordinate.
  4. The coloring with 2n - 1 colors from Appendix A.2 of all measurements that
     label the rays of the six bases with n labels, for n = 3 (417 measurements,
     among them the 273 measurements of S_0(K_3)) and n = 4.
  5. The graph W_lambda(D) of Section 2.5 for the eighteen triads of the square
     and lambda = 0, 1: 27 wires, one of them between complementary parts, the
     identification of the wire ends, the edge rule, compatibility up to distance
     lambda + 2, and the absence of a one-in-three selection.
"""
from fractions import Fraction as Fr
from itertools import combinations, product

from square_data import (BASES, IDENT, INCIDENCE, Q, ZERO, check_square, madd, mmul, msub,
                         orthogonal, plane)

C3 = range(3)


# ------------------------------------------------------------ 1. Pauli square
def kron(a, b):
    return tuple(tuple(a[i][k] * b[j][l] for k in range(2) for l in range(2))
                 for i in range(2) for j in range(2))


def check_pauli_square():
    one = ((1, 0), (0, 1))
    x = ((0, 1), (1, 0))
    z = ((1, 0), (0, -1))
    iy = ((0, 1), (-1, 0))                       # i times the Pauli matrix Y
    yy = tuple(tuple(-c for c in row) for row in kron(iy, iy))   # Y (x) Y, real
    square = [[kron(x, one), kron(one, x), kron(x, x)],
              [kron(one, z), kron(z, one), kron(z, z)],
              [kron(x, z), kron(z, x), yy]]

    def apply(m, v):
        return tuple(sum(m[a][b] * v[b] for b in range(4)) for a in range(4))

    def eigenvalue(m, v):
        w = apply(m, v)
        for sign in (1, -1):
            if w == tuple(sign * c for c in v):
                return sign
        raise AssertionError("not an eigenvector")

    lines = {i: square[i - 1] for i in (1, 2, 3)}
    lines.update({k: [square[r][k - 4] for r in range(3)] for k in (4, 5, 6)})
    for i, observables in lines.items():
        patterns = [tuple(eigenvalue(m, v) for m in observables) for v in BASES[i]]
        assert len(set(patterns)) == 4, i          # a common eigenbasis, no degeneracy
    for (i, k), (pos_i, _) in INCIDENCE.items():
        obs = square[i - 1][k - 4]
        proj = tuple(tuple(Fr(int(a == b) + obs[a][b], 2) for b in range(4)) for a in range(4))
        assert plane(i, pos_i) == proj, (i, k)     # P_e = (I + O_e) / 2

    def prod3(ms):
        out = IDENT
        for m in ms:
            out = mmul(out, tuple(tuple(Fr(c) for c in row) for row in m))
        return out

    minus = tuple(tuple(-c for c in row) for row in IDENT)
    assert [prod3(lines[i]) for i in (1, 2, 3, 4, 5)] == [IDENT] * 5
    assert prod3(lines[6]) == minus
    print("Pauli square: six common eigenbases, P_e = (I + O_e)/2 at all nine edges of K_3,3, "
          "products I, I, I and I, I, -I")


# ------------------------------------------------------------ 2. parity lemma
def check_parity():
    edges = sorted(INCIDENCE)
    consistent = 0
    for choice in product(range(4), repeat=6):
        s = dict(zip(range(1, 7), choice))
        pairs = 0
        ok = True
        for (i, k) in edges:
            pos_i, pos_k = INCIDENCE[i, k]
            a, b = s[i] in pos_i, s[k] in pos_k
            pairs += a + b
            ok = ok and (a == b)
        assert pairs % 2 == 1                       # the count in the proof is always odd
        consistent += ok
    assert consistent == 0
    print("parity lemma: none of the 4^6 = 4096 choices of positions is consistent on all nine wires")


# ---------------------------------------------------- 3. proper 4-colorings
def proper_colourings(vertices, adjacent, colours, fixed=None, forbidden=None):
    """All proper colorings, as tuples indexed like `vertices` (backtracking)."""
    n = len(vertices)
    nbrs = [[b for b in range(n) if adjacent(vertices[a], vertices[b])] for a in range(n)]
    colour = [-1] * n
    out = []

    def rec(a):
        if a == n:
            out.append(tuple(colour))
            return
        v = vertices[a]
        if fixed and v in fixed:
            options = (fixed[v],)
        else:
            options = range(colours)
        for c in options:
            if forbidden and c in forbidden[a]:
                continue
            if any(colour[b] == c for b in nbrs[a] if b < a):
                continue
            colour[a] = c
            rec(a + 1)
        colour[a] = -1

    rec(0)
    return out


CELLS = list(product(C3, repeat=2))


def cells_adjacent(u, v):
    return u[0] != v[0] and u[1] != v[1]


def reads(f):
    """(reads the first coordinate, reads the second coordinate) for f: cell -> color."""
    def in_row_class(x, y):
        return any(f[x, yy] == f[x, y] for yy in C3 if yy != y)

    def in_column_class(x, y):
        return any(f[xx, y] == f[x, y] for xx in C3 if xx != x)

    first = any(all(in_row_class(x, y) for x in C3) for y in C3)
    second = any(all(in_column_class(x, y) for y in C3) for x in C3)
    return first, second


def check_four_colourings():
    grids = [dict(zip(CELLS, c)) for c in proper_colourings(CELLS, cells_adjacent, 4)]
    kinds = []
    for f in grids:
        first, second = reads(f)
        assert first != second                                           # part (a)
        kinds.append(first)
    print(f"lemma (a): all {len(grids)} proper 4-colorings of K_3 x K_3 read exactly one "
          f"coordinate ({sum(kinds)} the first, {len(kinds) - sum(kinds)} the second)")

    pairs = 0
    for f, kind in zip(grids, kinds):
        forbidden = [{f[u] for u in CELLS if cells_adjacent(u, v)} for v in CELLS]
        for c in proper_colourings(CELLS, cells_adjacent, 4, forbidden=forbidden):
            g = dict(zip(CELLS, c))
            assert reads(g)[0] == kind                                   # part (b)
            pairs += 1
    print(f"lemma (b): all {pairs} pairs of colorings on consecutive layers read the same coordinate")

    triples = list(product(C3, repeat=3))

    def triples_adjacent(u, v):
        return all(a != b for a, b in zip(u, v))

    # the three constant triples are pairwise adjacent; fixing their colors
    # divides the number of colorings by 4 * 3 * 2 and does not change which
    # coordinates are read
    fixed = {(c, c, c): c for c in C3}
    cubes = proper_colourings(triples, triples_adjacent, 4, fixed=fixed)
    for c in cubes:
        f = dict(zip(triples, c))
        true_parts = 0
        for i in C3:
            grid = {(x, y): f[tuple(x if k == i else y for k in C3)] for x, y in CELLS}
            first, second = reads(grid)
            assert first != second
            true_parts += first
        assert true_parts == 1                                           # part (c)
    print(f"lemma (c): in all {24 * len(cubes)} proper 4-colorings of K_3 x K_3 x K_3 exactly one "
          "of the three restrictions reads the first coordinate")


# ------------------------------------- 4. the coloring with 2n - 1 colors
def check_colouring(n):
    rays = {(i, s): BASES[i][s] for i in BASES for s in range(4)}

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def contains(i, w, a, vec):
        """vec lies in the outcome space of label a of the labeling (i, w)."""
        return any(w[s] == a for s in range(4)) and all(
            dot(vec, rays[i, s]) == 0 for s in range(4) if w[s] != a)

    witnesses = [BASES[i][0] for i in (1, 2, 3, 4, 5)]            # p_1, ..., p_5
    u, v, w6, _ = BASES[6]
    # the facts about these vectors used in the proof
    assert all(dot(p, q) != 0 for p, q in combinations(witnesses, 2))
    assert tuple(a + b for a, b in zip(u, w6)) == witnesses[0]
    assert tuple(a + b for a, b in zip(u, v)) == tuple(2 * c for c in witnesses[1])
    assert tuple(a + b for a, b in zip(v, w6)) == witnesses[2]
    assert all(dot(u, witnesses[i]) != 0 for i in (0, 1, 3, 4)) and dot(v, witnesses[2]) != 0

    last = n - 1                                                   # the label called n in the paper

    def colour_of(i, w):
        for idx, p in enumerate(witnesses, start=1):
            for a in range(n):
                if contains(i, w, a, p):
                    return ("primed", a) if idx == 3 and a != last else ("ordinary", a)
        assert i == 6 and len({w[0], w[1], w[2]}) == 3
        if w[0] != last:
            return ("ordinary", w[0])
        assert w[1] != last
        return ("primed", w[1])

    def key(i, w):
        return tuple(plane(i, [s for s in range(4) if w[s] == a]) for a in range(n))

    measurements = {}
    for i in BASES:
        for w in product(range(n), repeat=4):
            k = key(i, w)
            c = colour_of(i, w)
            if k in measurements:
                assert measurements[k][1] == c      # the color depends on the measurement only
            else:
                measurements[k] = ((i, w), c)
    reps = [rep for rep, _ in measurements.values()]
    cols = [c for _, c in measurements.values()]
    assert len(set(cols)) <= 2 * n - 1

    orth = {(a, b): dot(rays[a], rays[b]) == 0 for a in rays for b in rays}
    edges = 0
    for (x, (i, w)), (y, (k, z)) in combinations(enumerate(reps), 2):
        if all(orth[(i, s), (k, t)] for s in range(4) for t in range(4) if w[s] == z[t]):
            assert cols[x] != cols[y]
            edges += 1
    print(f"coloring: n = {n}: {len(reps)} measurements, {edges} orthogonal pairs, "
          f"properly colored with {len(set(cols))} <= {2 * n - 1} colors")
    return len(reps)


# ------------------------------ 5. the one-in-three graph of the square triads
def check_one_in_three_square(lam):
    from check_square_graphs import bfs_all, noncommuting_distances

    triads = []                                   # (basis, r, parts)
    for i in BASES:
        for r in (1, 2, 3):
            s, t = [x for x in (1, 2, 3) if x != r]
            triads.append((i, r, (plane(i, (0, r)), Q[i, s], Q[i, t])))
    assert len(triads) == 18

    def word_measurement(parts, w):
        out = []
        for a in C3:
            acc = ZERO
            for p in range(3):
                if w[p] == a:
                    acc = madd(acc, parts[p])
            out.append(acc)
        return tuple(out)

    verts, meas, chart = [], [], {}
    for c, (_, _, parts) in enumerate(triads):
        for w in product(C3, repeat=3):
            chart[c, w] = len(verts)
            verts.append(("chart", c, w))
            meas.append(word_measurement(parts, w))
    edges = set()

    def add(a, b):
        assert a != b
        edges.add((min(a, b), max(a, b)))
    for c in range(18):
        for w, x in combinations(list(product(C3, repeat=3)), 2):
            if all(a != b for a, b in zip(w, x)):
                add(chart[c, w], chart[c, x])

    classes = {}                                  # {P, I - P} -> [(triad, part index)]
    for c, (_, _, parts) in enumerate(triads):
        for p, proj in enumerate(parts):
            classes.setdefault(frozenset((proj, msub(IDENT, proj))), []).append((c, p))
    assert all(len({c for c, _ in occ}) == len(occ) for occ in classes.values())
    links = complementary = 0
    constraints = []                              # (triad, part, triad', part', equal?)
    for occ in classes.values():
        (c0, p0) = occ[0]
        first = triads[c0][2][p0]
        for (c1, p1) in occ[1:]:
            links += 1
            equal = triads[c1][2][p1] == first
            complementary += not equal
            constraints.append((c0, p0, c1, p1, equal))
            layer = {}
            for x, y in product(C3, repeat=2):
                layer[0, x, y] = chart[c0, tuple(x if q == p0 else y for q in range(3))]
                if equal:
                    end = tuple(x if q == p1 else y for q in range(3))
                else:
                    end = tuple(y if q == p1 else x for q in range(3))
                layer[lam + 1, x, y] = chart[c1, end]
                wanted = tuple(madd(first if a == x else ZERO,
                                    msub(IDENT, first) if a == y else ZERO) for a in C3)
                # both ends carry the measurement "x on P_1, y on its complement"
                assert meas[layer[0, x, y]] == wanted and meas[layer[lam + 1, x, y]] == wanted
                for ell in range(1, lam + 1):
                    layer[ell, x, y] = len(verts)
                    verts.append(("link", c0, c1, ell, x, y))
                    meas.append(wanted)
            for ell in range(lam + 2):
                for (x, y), (xx, yy) in product(product(C3, repeat=2), repeat=2):
                    if x != xx and y != yy:
                        add(layer[ell, x, y], layer[ell, xx, yy])
                        if ell <= lam:
                            add(layer[ell, x, y], layer[ell + 1, xx, yy])
    assert links == 27 and complementary == 1
    assert len(verts) == 18 * 27 + 27 * 9 * lam
    assert all(orthogonal(meas[a], meas[b]) for a, b in edges)
    dist = bfs_all(len(verts), edges, lam + 3)
    bad = noncommuting_distances(meas, dist)
    assert bad and min(bad) == lam + 3            # compatible exactly up to distance lam + 2

    # one-in-three selections: one part per triad, equal parts alike, complementary parts opposite
    by_basis = {i: [c for c, (b, _, _) in enumerate(triads) if b == i] for i in BASES}
    local = {}
    for i, cs in by_basis.items():
        inner = [k for k in constraints if k[0] in cs and k[2] in cs]
        assert len(inner) == 3 and all(k[4] for k in inner)
        good = []
        for choice in product(range(3), repeat=3):
            sel = dict(zip(cs, choice))
            if all((sel[c0] == p0) == (sel[c1] == p1) for c0, p0, c1, p1, _ in inner):
                good.append(sel)
        assert len(good) == 4                     # the four positions of the basis
        local[i] = good
    outer = [k for k in constraints if triads[k[0]][0] != triads[k[2]][0]]
    assert len(outer) == 9
    selections = 0
    for combo in product(*(local[i] for i in sorted(local))):
        sel = {}
        for part in combo:
            sel.update(part)
        if all(((sel[c0] == p0) == (sel[c1] == p1)) == equal for c0, p0, c1, p1, equal in outer):
            selections += 1
    assert selections == 0
    print(f"W_{lam} for the 18 triads of the square: {len(verts)} vertices, {len(edges)} edges, 27 wires "
          f"(one between complementary parts); strategy valid and compatible exactly up to distance "
          f"{lam + 2}; no one-in-three selection")


if __name__ == "__main__":
    assert check_square()
    check_pauli_square()
    check_parity()
    check_four_colourings()
    assert check_colouring(3) == 417
    check_colouring(4)
    for lam in (0, 1):
        check_one_in_three_square(lam)
    print("check_first_level: all checks passed")
