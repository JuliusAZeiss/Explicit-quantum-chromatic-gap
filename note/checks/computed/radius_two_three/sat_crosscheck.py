"""Independent test that the graphs with 207 and 279 vertices have no proper 4-coloring.

Needs python-sat.  The two programs next to this file prove the lower bound
chi >= 5 by a parity argument with an exhaustive local check.  This program
ignores that argument: it hands the two graphs to a SAT solver, with the
colors of one triangle fixed.  It takes about twenty minutes (one minute for
the first graph) and is therefore not part of ../../run_all.py.

    python3 -B sat_crosscheck.py
"""
import contextlib
import io
import json
import time
from pathlib import Path

from pysat.solvers import Cadical153

import verify_square_radius_three_replacements as radius_three


def colorable(n, edges, colors, clique):
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


def triangle(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    for a, b in edges:
        common = adj[a] & adj[b]
        if common:
            return [a, b, min(common)]
    raise AssertionError("no triangle")


def main():
    data = json.loads(Path(__file__).with_name("square_radius_two_207.json").read_text())
    graphs = [(len(data["vertices"]), [tuple(e) for e in data["edges"]])]
    with contextlib.redirect_stdout(io.StringIO()):
        result = radius_three.construct(279)
    graphs.append((len(result["vertices"]), [tuple(e) for e in result["edges"]]))
    for n, edges in graphs:
        start = time.time()
        fixed = triangle(n, edges)
        assert colorable(n, edges, 5, fixed)
        assert not colorable(n, edges, 4, fixed)
        print(f"{n} vertices, {len(edges)} edges: a proper 5-coloring exists and no proper "
              f"4-coloring ({time.time() - start:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
