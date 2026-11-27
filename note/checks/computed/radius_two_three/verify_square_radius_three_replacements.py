"""Explicit radius-three graphs in dimension four. Requires the companion module
verify_square_radius_two_207.py and Python 3.8+. No third-party packages.

construct(279) and construct(287) return the full labeled graphs, measurements
and five-colorings. The local parity theorem and type-preserving bridges
rule out four classical colors: each of six charts has an odd number of
row-type ends, whereas exactly one end of each of nine wires is row-type.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path

import verify_square_radius_two_207 as square


CELLS = square.CELLS
OFF_DIAGONAL = {xy for xy in CELLS if xy[0] != xy[1]}
DIAGONAL = set(CELLS) - OFF_DIAGONAL
STRIPE = {xy for xy in CELLS if xy[0] != 0}
GRID_EDGES = [(x,y) for x in CELLS for y in CELLS if x[0] != y[0] and x[1] != y[1]]


def port_at_each_end(operators, grids):
    result = {}
    for b, d in product(range(3), range(3, 6)):
        matches = [(i,j) for i,g in enumerate(grids) for j,h in enumerate(grids)
                   if operators[b, g[0,1]] == operators[d, h[2,0]]]
        assert len(matches) == 1
        result[b,d] = matches[0]
    # The nine identifications use all three ports at each chart exactly once.
    ends = [(b,i) for (b,d),(i,j) in result.items()] + [(d,j) for (b,d),(i,j) in result.items()]
    assert len(set(ends)) == 18
    return result


def rename_tail(word, port):
    # A cyclic permutation of the last three coordinates sends {0,1} to {0,port+1}.
    return (word[0], *(word[1 + (j - port) % 3] for j in range(3)))


def chart_template(port, mode):
    _, old, _, _ = square.grids_and_chart()
    if mode == 'plain':
        return sorted(old), set()
    if mode == 'all':
        remove = ['0100','1211','2022']
        add = ['0020','1101','2212']
        filtered = OFF_DIAGONAL
    elif mode == 'two':
        remove = ['1211','2022']
        add = ['0020','1101']
        filtered = {(0,1),(0,2)}
    else:
        assert mode == 'old'
        remove, add = [], []
        filtered = {(a, (a+1)%3) for a in range(3)}
    to_word = lambda s: rename_tail(tuple(map(int,s)), port)
    words = set(old) - set(map(to_word, remove)) | set(map(to_word, add))
    assert len(words) == 39
    return sorted(words), filtered


def local_rule(marked, half):
    def adjacent(u, v):
        return (all(a != b for a,b in zip(u,v)) and
                (u not in marked or v[0] == v[half]) and
                (v not in marked or u[0] == u[half]))
    return adjacent


def check_maps(sequence, edges):
    def neighbor(u,v):
        return tuple(sorted((u,v))) in edges
    for m in sequence:
        assert set(m) == set(CELLS)
        assert all(neighbor(m[x],m[y]) for x,y in GRID_EDGES)
    for left,right in zip(sequence,sequence[1:]):
        assert all(neighbor(left[x],right[y]) for x,y in GRID_EDGES)


def construct(order=279):
    if order not in (279,287):
        raise ValueError('Choose 279 or 287.')
    grids, old, _, _ = square.grids_and_chart()
    ops = square.operator_data(list(product(range(3),repeat=4)))
    ports = port_at_each_end(ops, grids)
    graph_keys, physical, edges = [], [], set()
    end, filters = [], []
    checked_charts = set()
    def add_vertex(key,operator):
        v=len(physical)
        graph_keys.append(key)
        physical.append(operator)
        return v
    def add(u,v):
        assert u!=v
        edges.add(tuple(sorted((u,v))))

    for b in range(6):
        matching = b%3
        port = ports[matching, matching+3][int(b>=3)]
        mode = 'all' if order==279 else ('two' if b==0 else 'old' if b==3 else 'plain')
        words, active = chart_template(port, mode)
        filters.append(active)
        chart = {w:add_vertex(['chart',b,list(w)],ops[b,w]) for w in words}
        end.append(chart)
        marked = {grids[port][xy] for xy in active}
        edge = local_rule(marked,port+1)
        for u,v in combinations(words,2):
            if edge(u,v):
                add(chart[u],chart[v])
        # Every distinct chart edge rule is audited from scratch.
        key = (tuple(words),tuple(sorted(marked)),port)
        if key not in checked_charts:
            square.verify_classical_obstruction(grids, words, edge)
            checked_charts.add(key)

    for (b,d),(i,j) in ports.items():
        left = {xy:end[b][grids[i][xy]] for xy in CELLS}
        right = {(x,y):end[d][grids[j][-y%3,-x%3]] for x,y in CELLS}
        assert all(physical[left[x]] == physical[right[x]] for x in CELLS)
        matched = (d==b+3 if order==279 else (b,d)==(0,3))
        if matched:
            active_left = filters[b]
            active_right = {(x,y) for x,y in CELLS if (-y%3,-x%3) in filters[d]}
            middle = DIAGONAL if order==279 else {(1,0),(1,1),(1,2),(2,1),(2,2)}
        else:
            active_left=active_right=set()
            middle=STRIPE
        relay={xy:add_vertex(['relay',b,d,list(xy)],physical[left[xy]]) for xy in sorted(middle)}
        for x,y in GRID_EDGES:
            if x in active_left and y in active_right:
                add(left[x],right[y])
            if y in relay:
                add(left[x],relay[y])
                add(right[x],relay[y])
                if x in relay:
                    add(relay[x],relay[y])
        # Physical coordinates at the right end exchange the two local axes.
        # Equal physical types across all nine wires imply nine row-type ends.
        left_with_relay = {**left, **relay}
        right_with_relay = {**right, **relay}
        sequence = [left, left_with_relay]
        if matched and order == 287:
            # Among the four unbuffered cells, only (2,0) is adjacent to
            # (0,1) or (0,2). Both of the latter have direct edges from
            # their left copies to the right copy of (2,0). Switch it first.
            sequence.append({**left_with_relay, (2, 0): right[2, 0]})
        sequence.extend([right_with_relay, right])
        check_maps(sequence, edges)
    assert len(physical)==order
    colors=[square.classical_color(p) for p in physical]
    assert len(set(colors))==5
    phys=sorted(set(physical))
    lookup={p:i for i,p in enumerate(phys)}
    which=[lookup[p] for p in physical]
    prior=set(square.operator_data(old).values())
    assert len(phys)==165 and len(set(phys)-prior)==(18 if order==279 else 2)
    mm=lru_cache(None)(square.matmul)
    for m in phys:
        for p in m:
            assert mm(p,p)==tuple(4*t for t in p)
            assert all(p[4*i+j]==p[4*j+i] for i,j in product(range(4),repeat=2))
        assert tuple(map(sum,zip(*m))) == tuple(4*t for t in square.IDENTITY)
        assert all(mm(p,q)==square.ZERO for p,q in combinations(m,2))
    @lru_cache(None)
    def commute(a,b):
        return all(mm(p,q)==mm(q,p) for p,q in product(phys[a],phys[b]))
    neighbors=[set() for m in physical]
    for u,v in edges:
        assert colors[u]!=colors[v]
        assert all(mm(p,q)==square.ZERO for p,q in zip(physical[u],physical[v]))
        neighbors[u].add(v)
        neighbors[v].add(u)
    dist=Counter()
    for src in range(len(physical)):
        seen={src}
        layer=seen
        for d in (1,2,3):
            layer={w for v in layer for w in neighbors[v]}-seen
            seen=seen|layer
            for w in layer:
                if src<w:
                    assert commute(*sorted((which[src],which[w])))
                    dist[d]+=1
    tri=[end[0][(a,)*4] for a in range(3)]
    assert all(tuple(sorted((u,v))) in edges for u,v in combinations(tri,2))
    result=dict(vertices=graph_keys,edges=sorted(edges),measurements_times_four=phys,measurement_index=which,
                classical_five_coloring=colors,pairs_by_distance=dict(sorted(dist.items())))
    print(f'{order} vertices, {len(edges)} edges, 165 physical measurements, {len(set(phys)-prior)} replacements.')
    print('Commuting pairs at distances 1, 2, 3:', [dist[i] for i in (1,2,3)])
    return result


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--order',type=int,default=279,choices=[279,287])
    ap.add_argument('--export',type=Path)
    args=ap.parse_args()
    graph=construct(args.order)
    if args.export:args.export.write_text(json.dumps(graph,separators=(',',':'))+'\n')
