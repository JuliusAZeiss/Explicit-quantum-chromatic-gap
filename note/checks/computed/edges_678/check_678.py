"""Exact finite check for the six-chart, 165-vertex, 678-edge square.
Plain Python 3. No external packages or SAT solver are required.
The local test enumerates all proper colorings of K_3 x K_3 using a
four-color palette (which includes the three-color possibilities).
"""
from itertools import combinations, product

C = range(3)
O = (0, 0, 0)
A = [tuple(int(t != i) for t in C) for i in C]
B = [tuple((-x) % 3 for x in a) for a in A]
D = [tuple((-int(t == i)) % 3 for t in C) for i in C]
U = [tuple((t - i) % 3 for t in C) for i in C]
difference = [O] + A + B + D + U
phase_edges = set()
for i in C:
    a, b = A[i], B[i]
    d, e = D[i], D[(i + 1) % 3]
    u, v, w = U[i], U[(i + 1) % 3], U[(i + 2) % 3]
    triangles = [(O,d,a), (d,a,v), (a,v,u), (a,u,b),
                 (b,u,w), (b,w,e), (b,e,O)]
    for triangle in triangles:
        phase_edges.update(frozenset(pair) for pair in combinations(triangle, 2))

def word(p, k):
    return (k,) + tuple((k + x) % 3 for x in p)

vertices = [word(p, k) for p in difference for k in C]
def phase(v):
    return tuple((x - v[0]) % 3 for x in v[1:])

def edge(v, w):
    p, q = phase(v), phase(w)
    return all(a != b for a, b in zip(v, w)) and (
        (p == q and p not in D) or frozenset((p, q)) in phase_edges)

cells = [(i,j) for i in C for j in C]
grids = [[(r,) + tuple(r if t == i else c for t in C) for r,c in cells] for i in C]
W = [(a, k, (k+1)%3, (k+2)%3) for a,k in cells]
R = [word(p, k) for p in D for k in C]
assert sum(edge(v,w) for v,w in combinations(vertices, 2)) == 138
assert not any(edge(v,w) for v,w in combinations(R, 2))
for block in [W] + grids:
    assert all(edge(block[i],block[j]) == all(x != y for x,y in zip(cells[i],cells[j]))
               for i,j in combinations(range(9),2))
def four_colorings():
    """Enumerate all proper colorings on the nine positions, in row order."""
    colors = []
    def extend():
        if len(colors) == 9:
            yield tuple(colors)
            return
        i,j = cells[len(colors)]
        forbidden = {c for c,(a,b) in zip(colors,cells) if a != i and b != j}
        for c in range(4):
            if c not in forbidden:
                colors.append(c)
                yield from extend()
                colors.pop()
    yield from extend()

graph_maps = list(four_colorings())
# The diagonal triangle in each plane has the three fixed colors.
plane_maps = [f for f in graph_maps if (f[0],f[4],f[8]) == (0,1,2)]
def row_type(f):
    row = all(len(set(f[3*r:3*r+3])) < 3 for r in C)
    col = all(len(set(f[c::3])) < 3 for c in C)
    assert row != col
    return row

def touches(left, right):
    return [(i,j) for i,v in enumerate(left) for j,w in enumerate(right) if edge(v,w)]
grid_W = [touches(grid,W) for grid in grids]
def ray_domains(block, f):
    # A bit (4*i+c) records whether color c at ray word i is allowed.
    bits=(1 << 36)-1
    for i,j in touches(R,block):
        bits &= ~(1 << (4*i + f[j]))
    return bits
plane_records = [[(f,row_type(f),ray_domains(grid,f)) for f in plane_maps] for grid in grids]
outcomes = set()
for g in graph_maps:
    wr = ray_domains(W,g)
    options = [[(t,m) for f,t,m in rec if all(f[i] != g[j] for i,j in adj)]
               for rec,adj in zip(plane_records,grid_W)]
    for choice in product(*options):
        allowed=wr
        for t,m in choice:
            allowed &= m
        if all((allowed >> (4*i)) & 15 for i in range(9)):
            outcomes.add(tuple(int(t) for t,m in choice))
assert len(graph_maps) == 1056 and len(plane_maps) == 44
assert outcomes == {(1,0,0),(0,1,0),(0,0,1),(1,1,1)}
print('Chart edges: 138. W-grid colorings: 1056. Normalized plane colorings: 44.')
print('Possible plane-type triples:', sorted(outcomes))
print('Identified graph edges:', 138 + 6*(138-48))

# Restore the exact vertex order of the edge-deletion witnesses.
import gzip
from pathlib import Path
bases = [
    [(1,1,1,1),(1,-1,1,-1),(1,1,-1,-1),(1,-1,-1,1)],
    [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)],
    [(1,1,1,-1),(1,1,-1,1),(1,-1,1,1),(-1,1,1,1)],
    [(1,0,1,0),(1,0,-1,0),(0,1,0,1),(0,1,0,-1)],
    [(1,1,0,0),(1,-1,0,0),(0,0,1,1),(0,0,1,-1)],
    [(1,0,0,1),(1,0,0,-1),(0,1,1,0),(0,1,-1,0)]]

def measurement(basis, v):
    # All matrix entries of all projectors become integral when scaled by four.
    return tuple(tuple(sum((4 // sum(z*z for z in ray)) * ray[a] * ray[b]
                           for color,ray in zip(v,basis) if color == c)
                       for a,b in product(range(4),repeat=2)) for c in C)

measurements = sorted({measurement(b,v) for b in bases for v in vertices})
ids = {m:i for i,m in enumerate(measurements)}
graph_edges = set()
for b in bases:
    local = {v:ids[measurement(b,v)] for v in vertices}
    for v,w in combinations(vertices,2):
        if edge(v,w):
            graph_edges.add(tuple(sorted((local[v],local[w]))))
assert (len(measurements),len(graph_edges)) == (165,678)
triangular_edges = set(graph_edges)
for b in bases:
    for p in D:
        three = [ids[measurement(b,word(p,k))] for k in C]
        triangular_edges.update(tuple(sorted(pair)) for pair in combinations(three,2))
assert len(triangular_edges) == 732
forced = set()
with gzip.open(Path(__file__).with_name('edge_forced_732.txt.gz'),'rt') as certificate:
    for line in certificate:
        aa,bb,f = line.split()
        e = (int(aa),int(bb))
        assert len(f) == 165 and set(f) <= set('0123')
        assert {(a,b) for a,b in triangular_edges if f[a] == f[b]} == {e}
        forced.add(e)
assert forced == graph_edges
print('All 678 edges are necessary even in the 732-edge triangular host.')

# Two measurements are adjacent precisely when at every outcome their
# projections are orthogonal. By positive semidefiniteness, zero trace
# of the product is equivalent to orthogonality. The flattened matrices
# are symmetric, so the trace is their ordinary dot product, up to scale.
orth_edges = {(i,j) for i,j in combinations(range(len(measurements)),2)
              if all(sum(a*b for a,b in zip(p,q)) == 0
                     for p,q in zip(measurements[i],measurements[j]))}
assert len(orth_edges) == 843 and triangular_edges <= orth_edges
taken = set()
sizes = []
with gzip.open(Path(__file__).with_name('edge_packing_843.txt.gz'),'rt') as certificate:
    for line in certificate:
        f = line.strip()
        assert len(f) == 165 and set(f) <= set('0123')
        bad = {(i,j) for i,j in orth_edges if f[i] == f[j]}
        assert bad and not (taken & bad)
        assert len(graph_edges & bad) == 1
        sizes.append(len(bad))
        taken.update(bad)
assert sorted(sizes) == [1]*621 + [2]*57
print('843-edge host: 678 disjoint nonempty monochromatic edge sets, occupying 735 edges.')
print('Every five-chromatic subgraph of the full host has at least 678 edges.')
