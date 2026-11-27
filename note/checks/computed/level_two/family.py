"""Symbolic generator and class counts for the two-level square graphs.

A vertex is (template index, tuple of distinct bottom labels). The optional
threshold argument selects the six- or seven-color lower bound. The main program
counts vertices by their nine classical colors without expanding the graphs.
"""
from collections import Counter
from functools import lru_cache
import gzip
import itertools
import json
from pathlib import Path
import sys

SUPPLEMENT = Path(__file__).resolve().parent   # verify_small_graphs.py and its data lie next to this file
sys.path.insert(0, str(SUPPLEMENT))
import verify_small_graphs as square


def templates():
    # (square basis, slot at each of the four positions, unrestricted labels)
    yield (0, (0, 0, 0, 0), True)
    for basis in range(6):
        for ray in (1, 2, 3):
            yield (basis, tuple(int(i != ray) for i in range(4)), False)
    for basis in range(3):
        for ray in (1, 2, 3):
            # The complemented incidence is basis 3 to basis 6, axis 3.
            yield (basis, tuple(int(i not in (0, ray)) for i in range(4)), (basis, ray) == (2, 3))
    for basis in range(6):
        for ray in (1, 2, 3):
            other = [i for i in range(4) if i not in (0, ray)]
            yield (basis, tuple(0 if i in (0, ray) else 1 + other.index(i) for i in range(4)), False)


TEMPLATES = tuple(templates())


def dot(x, y):
    return sum(a*b for a,b in zip(x,y))


def supports(source, target):
    return [i for i,v in enumerate(square.SQUARE[source]) if dot(v,target)]


@lru_cache(None)
def upper_color(template, colors):
    basis, slots, _ = TEMPLATES[template]
    word = tuple(colors[x] for x in slots)
    for i, witness_basis in enumerate(square.SQUARE[:5]):
        palette = {word[pos] for pos in supports(basis,witness_basis[0])}
        if len(palette) == 1:
            a, = palette
            return a+5 if i == 2 and a < 4 else a
    assert basis == 5
    a,b = word[:2]
    return a if a < 4 else b+5


@lru_cache(None)
def bottom():
    graph = square.graph_165()
    with gzip.open(SUPPLEMENT / "small_graphs_data.json.gz", "rt") as stream:
        coloring = tuple(json.load(stream)["coloring_165"])
    assert len(graph.meas) == 165 and square.is_proper(graph.adj, coloring)
    assert list(range(3)) == [graph.alias[("const",c)] for c in square.COLORS]
    return graph, coloring


def class_counts():
    graph, color = bottom()
    groups = Counter((color[v], v < 3) for v in range(graph.n))
    counts = Counter()
    for template, (_, slots, unrestricted) in enumerate(TEMPLATES):
        for assignments in itertools.product(groups, repeat=1+max(slots)):
            if not unrestricted and not any(scalar for _,scalar in assignments):
                continue
            ways = 1
            for i,group in enumerate(assignments):
                ways *= max(0, groups[group] - assignments[:i].count(group))
            counts[upper_color(template,tuple(c for c,_ in assignments))] += ways
    return counts


@lru_cache(None)
def independent_certificate():
    data = json.loads(Path(__file__).with_name("independent_set.json").read_text())
    labels = frozenset(data["bottom"])
    projections = frozenset(tuple(p) for p in data["intersecting"])
    graph,_ = bottom()
    assert all(not graph.adj[v].intersection(labels) for v in labels)
    assert all(dot(p,q)>0 for p in projections for q in projections)
    return labels, projections


@lru_cache(None)
def removed_pattern(template, inside):
    basis,slots,_ = TEMPLATES[template]
    selected = [square.proj4(square.SQUARE[basis][i]) for i,s in enumerate(slots) if inside[s]]
    projection = tuple(sum(p[i][j] for p in selected) for i in range(4) for j in range(4))
    return projection in independent_certificate()[1]


def count_removed():
    graph,_ = bottom()
    labels,_ = independent_certificate()
    groups = Counter((v in labels,v<3) for v in range(graph.n))
    size = 0
    for t,(_,slots,unrestricted) in enumerate(TEMPLATES):
        for flags in itertools.product(groups,repeat=1+max(slots)):
            if not unrestricted and not any(s for _,s in flags):
                continue
            if removed_pattern(t,tuple(inside for inside,_ in flags)):
                ways = 1
                for i,group in enumerate(flags):
                    ways *= max(0,groups[group]-flags[:i].count(group))
                size += ways
    return size


def vertices(threshold=7):
    if threshold not in (6,7):
        raise ValueError("supported thresholds are 6 and 7")
    graph,_ = bottom()
    labels,_ = independent_certificate()
    for t, (_, slots, unrestricted) in enumerate(TEMPLATES):
        k = 1 + max(slots)
        for mask in range(1 << k):
            if not unrestricted and not mask:
                continue
            pools = [range(3) if mask & (1 << i) else range(3,graph.n) for i in range(k)]
            for v in itertools.product(*pools):
                if len(set(v)) < k:
                    continue
                if threshold == 7 or not removed_pattern(t,tuple(x in labels for x in v)):
                    yield t,v


def adjacent(left,right):
    graph,_ = bottom()
    lt,lv = left
    rt,rv = right
    lb,ls,_ = TEMPLATES[lt]
    rb,rs,_ = TEMPLATES[rt]
    return all(not dot(u,v) or rv[rs[j]] in graph.adj[lv[ls[i]]]
               for i,u in enumerate(square.SQUARE[lb]) for j,v in enumerate(square.SQUARE[rb]))


if __name__ == "__main__":
    counts = class_counts()
    assert sum(counts.values()) == 4330533
    removed = count_removed()
    assert removed == 1407626
    print("Class sizes for colors 0 to 8:", [counts[c] for c in range(9)])
    print("Order at threshold 7:", sum(counts.values()))
    print("Order at threshold 6:", sum(counts.values())-removed)
