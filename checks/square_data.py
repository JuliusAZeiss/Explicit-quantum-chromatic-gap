"""Exact data for the Mermin--Peres square used in the paper.

Standard library only.  Everything is exact (fractions.Fraction).

Conventions follow Section 2 of the paper:
  * six bases i = 1..6 of C^4, positions s = 0..3,
  * A_r = {0, r} and P^i_r = Q^i_0 + Q^i_r for r = 1, 2, 3,
  * the nine incidences e = (i, k), i <= 3 < k, with the position sets
    A_{i,e} and A_{k,e} spanning the common plane P_e.
"""
from fractions import Fraction as Fr
from itertools import product

TABLE = {
    1: ["++++", "+-+-", "++--", "+--+"],
    2: ["+000", "0+00", "00+0", "000+"],
    3: ["+++-", "++-+", "+-++", "-+++"],
    4: ["+0+0", "+0-0", "0+0+", "0+0-"],
    5: ["++00", "+-00", "00++", "00+-"],
    6: ["+00+", "+00-", "0++0", "0+-0"],
}
SIGN = {"+": 1, "-": -1, "0": 0}
BASES = {i: [tuple(SIGN[c] for c in w) for w in rows] for i, rows in TABLE.items()}

ZERO = tuple(tuple(Fr(0) for _ in range(4)) for _ in range(4))
IDENT = tuple(tuple(Fr(int(a == b)) for b in range(4)) for a in range(4))


def madd(x, y):
    return tuple(tuple(x[a][b] + y[a][b] for b in range(4)) for a in range(4))


def msub(x, y):
    return tuple(tuple(x[a][b] - y[a][b] for b in range(4)) for a in range(4))


def mmul(x, y):
    return tuple(tuple(sum(x[a][c] * y[c][b] for c in range(4)) for b in range(4))
                 for a in range(4))


def ray_projector(v):
    n = sum(c * c for c in v)
    return tuple(tuple(Fr(v[a] * v[b], n) for b in range(4)) for a in range(4))


Q = {(i, s): ray_projector(BASES[i][s]) for i in BASES for s in range(4)}


def plane(i, positions):
    out = ZERO
    for s in positions:
        out = madd(out, Q[i, s])
    return out


A = {1: (0, 1), 2: (0, 2), 3: (0, 3)}
COMPLEMENT_A3 = (1, 2)

# The nine identities (1) of the paper: (row basis, r) = (column basis, r'),
# with the complemented incidence (3,3) = I - (6,3).
IDENTITIES = [
    ((1, 1), (4, 2)), ((1, 2), (5, 2)), ((1, 3), (6, 2)),
    ((2, 1), (5, 1)), ((2, 2), (4, 1)), ((2, 3), (6, 1)),
    ((3, 1), (5, 3)), ((3, 2), (4, 3)), ((3, 3), (6, 3)),
]
# incidence e = (i, k) -> (positions at i, positions at k)
INCIDENCE = {}
for (i, r), (k, rp) in IDENTITIES:
    if (i, k) == (3, 6):
        INCIDENCE[(i, k)] = (A[r], COMPLEMENT_A3)
    else:
        INCIDENCE[(i, k)] = (A[r], A[rp])


def check_square():
    """Orthonormality of the six bases and the nine plane identities."""
    for i in BASES:
        for s, t in product(range(4), repeat=2):
            dot = sum(a * b for a, b in zip(BASES[i][s], BASES[i][t]))
            assert (dot == 0) == (s != t)
        assert plane(i, range(4)) == IDENT
    for (i, k), (pos_i, pos_k) in INCIDENCE.items():
        assert plane(i, pos_i) == plane(k, pos_k), (i, k)
    # the complemented incidence is the only one whose column positions miss 0
    assert [e for e, (_, pk) in INCIDENCE.items() if 0 not in pk] == [(3, 6)]
    assert plane(3, A[3]) == msub(IDENT, plane(6, A[3]))
    # every row basis meets every column basis in exactly one plane, and the
    # three planes at a basis use A_1, A_2, A_3 (A_3^c at basis 6 on (3,6))
    for i in (1, 2, 3):
        assert sorted(INCIDENCE[i, k][0] for k in (4, 5, 6)) == [A[1], A[2], A[3]]
    for k in (4, 5):
        assert sorted(INCIDENCE[i, k][1] for i in (1, 2, 3)) == [A[1], A[2], A[3]]
    assert sorted(INCIDENCE[i, 6][1] for i in (1, 2, 3)) == [A[1], A[2], COMPLEMENT_A3]
    return True


def measurement(i, word):
    """Exact projective measurement of a basis word with colors 0, 1, 2."""
    out = []
    for a in range(3):
        out.append(plane(i, [s for s in range(4) if word[s] == a]))
    return tuple(out)


def orthogonal(m, n):
    return all(mmul(m[a], n[a]) == ZERO for a in range(3))


def commute(m, n):
    return all(mmul(m[a], n[b]) == mmul(n[b], m[a]) for a in range(3) for b in range(3))
