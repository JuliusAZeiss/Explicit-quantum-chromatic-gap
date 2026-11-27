"""Numerical companion to Appendix D (rounding to large outcomes).

For 8 <= d <= 400 the probability
    p = Pr[ Beta(s, d - s) > 1/2 + 1/d ],   s = ceil(d/3),
is computed exactly as a rational number through the identity
    Pr[Beta(a, b) > x] = Pr[Bin(a + b - 1, x) <= a - 1],
and compared with the lower bound theta^d / (768 d), theta = 3 / 2^(5/3), that is
proved in the appendix.  The script also reports p^(1/d), which tends to theta.
"""
from fractions import Fraction as Fr
from math import comb, log, exp

THETA = 3 / 2 ** (5 / 3)


def tail(d):
    s = -(-d // 3)
    x = Fr(1, 2) + Fr(1, d)
    n = d - 1
    return sum(comb(n, k) * x ** k * (1 - x) ** (n - k) for k in range(s))


if __name__ == "__main__":
    worst = None
    for d in range(8, 401):
        p = tail(d)
        # exact comparison through logarithms of exact rationals is unnecessary:
        # the margin is several orders of magnitude
        lp = log(p.numerator) - log(p.denominator)
        bound = d * log(THETA) - log(768 * d)
        assert lp > bound + 1, d
        ratio = lp - bound
        worst = ratio if worst is None else min(worst, ratio)
        if d in (8, 50, 100, 200, 400):
            print(f"d = {d:3d}: p^(1/d) = {exp(lp / d):.5f}   log(p / bound) = {ratio:.2f}")
    print(f"theta = {THETA:.5f}; p >= theta^d/(768 d) for 8 <= d <= 400 "
          f"(smallest log-margin {worst:.2f})")
    print(f"1/log2(2^(5/3)/3) = {1 / (5 / 3 - log(3) / log(2)):.4f}")
    print("check_rounding: all checks passed")
