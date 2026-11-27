"""Finite core of Appendix F (random descents; standard library only).

Checks at the first level, with the triangle K_3 as the substituted graph:
  1. no choice of positions is consistent on all nine edges of K_3,3 (the
     parity lemma of Section 2.1), and for each edge exactly 16 choices have
     that edge as their only defective edge; under the uniform distribution
     on these 16, each end of each edge points to its first entry in exactly
     half of them (Lemma F.1);
  2. for lambda = 1, 2: the descent map of every cut is a proper partial
     3-coloring of S_lambda(K_3) outside the cut layer; every interior vertex
     is outside the domain in a fraction 1/(9 lambda) of the cuts and is
     otherwise sent to its first and to its second entry equally often
     (Lemma F.2);
  3. for lambda = 0: every word with exactly two values, each at two
     positions, is two-valued along exactly one edge at its basis; the descent
     map of every cut is a proper partial 3-coloring of S_0(K_3) on its
     domain; each such word is outside the domain in a fraction 1/18 of the
     cuts and is otherwise sent to each of its two values equally often; all
     other words are always in the domain (Lemma F.5);
  4. the two identities behind the bounds: delta + (1 - delta) p/2 = p for
     delta = 1/(9 lambda), p = 2/(9 lambda + 1), and (17/18) 2^beta = 1 for
     beta = log2(18/17).
"""
from fractions import Fraction
from itertools import product
from math import log2

from check_square_graphs import first_level
from square_data import INCIDENCE

EDGES = sorted(INCIDENCE)          # the nine edges (i, k) of K_3,3, i <= 3 < k


def points_first(sigma, edge, end):
    """Does sigma point to the first entry at this end of this edge?"""
    i, k = edge
    pos_i, pos_k = INCIDENCE[edge]
    return sigma[i] in pos_i if end == i else sigma[k] in pos_k


def defective(sigma):
    return [e for e in EDGES if points_first(sigma, e, e[0]) != points_first(sigma, e, e[1])]


def choices_by_defect():
    """Choices of positions whose only defective edge is a given edge."""
    out = {e: [] for e in EDGES}
    for values in product(range(4), repeat=6):
        sigma = dict(zip(range(1, 7), values))
        bad = defective(sigma)
        assert bad, "a choice of positions consistent on all nine edges"
        if len(bad) == 1:
            out[bad[0]].append(sigma)
    return out


def check_choices(by_defect):
    for e0, sigmas in by_defect.items():
        assert len(sigmas) == 16
        for e in EDGES:
            for end in e:
                assert sum(points_first(s, e, end) for s in sigmas) == 8
    print("choices of positions: every choice has a defective edge; 16 choices for each "
          "single defective edge; each end of each edge points to its first entry in 8 of them")


def check_interior(lam, by_defect):
    verts, index, edges, _ = first_level(lam)
    interior = [n for n, v in enumerate(verts) if len(v) == 4]
    outside = {n: 0 for n in interior}
    first = {n: 0 for n in interior}
    second = {n: 0 for n in interior}
    cuts = 0
    for e0, sigmas in by_defect.items():
        for sigma in sigmas:
            for l0 in range(1, lam + 1):
                cuts += 1
                image = {}
                for n, v in enumerate(verts):
                    if len(v) == 2:                       # a word of a basis chart
                        b, w = v
                        image[n] = w[sigma[b]]
                        continue
                    e, ell, x, y = v                      # an interior vertex of a wire
                    if e == e0 and ell == l0:
                        outside[n] += 1
                        continue
                    if e == e0:
                        end = e[0] if ell < l0 else e[1]
                    else:
                        end = e[0]
                        assert points_first(sigma, e, e[0]) == points_first(sigma, e, e[1])
                    if points_first(sigma, e, end):
                        image[n] = x
                        first[n] += 1
                    else:
                        image[n] = y
                        second[n] += 1
                for a, b in edges:
                    if a in image and b in image:
                        assert image[a] != image[b], "the descent map is not a homomorphism"
    assert cuts == 9 * 16 * lam
    for n in interior:
        assert Fraction(outside[n], cuts) == Fraction(1, 9 * lam)
        assert first[n] == second[n]
    print(f"lambda = {lam}: {cuts} cuts; every descent map is a proper partial 3-coloring of "
          f"S_{lam}(K_3) ({len(verts)} vertices, {len(edges)} edges); each interior vertex is "
          f"cut with frequency 1/{9 * lam} and otherwise sent to either entry equally often")


def check_zero(by_defect):
    verts, index, edges, _ = first_level(0)
    assert all(len(v) == 2 for v in verts)

    def two_valued_along(v, edge):
        b, w = v
        if b not in edge:
            return False
        pos = INCIDENCE[edge][0] if b == edge[0] else INCIDENCE[edge][1]
        inside = {w[s] for s in pos}
        rest = {w[s] for s in range(4) if s not in pos}
        return len(inside) == 1 and len(rest) == 1 and inside != rest

    balanced = [n for n, (b, w) in enumerate(verts)
                if sorted(w.count(c) for c in set(w)) == [2, 2]]
    assert len(balanced) == 6 * 18
    for n in balanced:
        assert sum(two_valued_along(verts[n], e) for e in EDGES) == 1
    for n, v in enumerate(verts):
        if n not in set(balanced):
            assert not any(two_valued_along(v, e) for e in EDGES)
    outside = {n: 0 for n in balanced}
    value_count = {n: {} for n in balanced}
    cuts = 0
    for e0, sigmas in by_defect.items():
        for sigma in sigmas:
            for b0 in e0:
                cuts += 1
                image = {}
                for n, (b, w) in enumerate(verts):
                    if b == b0 and two_valued_along((b, w), e0):
                        outside[n] += 1
                        continue
                    image[n] = w[sigma[b]]
                    if n in value_count:
                        value_count[n][image[n]] = value_count[n].get(image[n], 0) + 1
                assert all(n in image for n in range(len(verts)) if n not in outside)
                for a, b in edges:
                    if a in image and b in image:
                        assert image[a] != image[b], "the descent map is not a homomorphism"
    assert cuts == 9 * 16 * 2
    for n in balanced:
        assert Fraction(outside[n], cuts) == Fraction(1, 18)
        counts = sorted(value_count[n].values())
        assert len(counts) == 2 and counts[0] == counts[1]
    print(f"lambda = 0: {cuts} cuts; every descent map is a proper partial 3-coloring of S_0(K_3) "
          f"on its domain; the {len(balanced)} words with two values at two positions each are "
          "removed with frequency 1/18 and otherwise sent to either value equally often")


def check_identities():
    for lam in range(1, 51):
        delta = Fraction(1, 9 * lam)
        p = Fraction(2, 9 * lam + 1)
        assert delta + (1 - delta) * p / 2 == p
    beta = log2(18 / 17)
    assert abs(17 / 18 * 2 ** beta - 1) < 1e-12 and beta < 0.083
    print(f"identities: delta + (1 - delta) p/2 = p for lambda = 1, ..., 50; beta = {beta:.5f}")


def main():
    by_defect = choices_by_defect()
    check_choices(by_defect)
    for lam in (1, 2):
        check_interior(lam, by_defect)
    check_zero(by_defect)
    check_identities()
    print("check_descents: all checks passed")


if __name__ == "__main__":
    main()
