# Programs of the note

Two kinds of programs accompany the note (`../main.pdf`).

- `check_note.py` independently recomputes finite facts that are **proved** in Sections 3–5 of the note. No statement depends on it.
- The programs in `computed/` establish the facts of Section 6, which are **known only by computation**.

All programs use exact arithmetic and the Python standard library (Python 3.10 or newer). `check_note.py` imports the data of the square and helper functions from the replay programs of the article in `../../checks`, and uses `python-sat` for its colorability tests when it is installed (otherwise a slower exact search, and the test of $S'_0(K_3)$ is skipped).

To run everything in sequence (about five minutes):

```sh
python3 -B run_all.py
```

`logs/run_all_2026-11-27.log` records a complete passing run (CPython 3.12.14 on macOS, with `python-sat` installed). `logs/sat_crosscheck_2026-11-27.log` records a run of `computed/radius_two_three/sat_crosscheck.py`.

## `check_note.py`: replay of proved facts

| What it recomputes | Statement of the note |
|---|---|
| Every proper 4-coloring of $K_3^{\times2}$ (1,056) and of $K_3^{\times3}$ (373,104) has exactly one dictator coordinate; identifying coordinates carries it along; the 12,288 pairs of colorings that satisfy the cross condition have the same one | Lemma 3.2, for $n=3$ and $L\le3$ |
| Every proper 6-coloring of $K_4^{\times2}$ (2,754,000) has exactly one dictator coordinate, which is the coordinate it reads; all pairs of such colorings that satisfy the cross condition have the same one | Lemma 3.2(a) and (c), for $n=4$ and $L=2$ |
| The graph $S'_0(K_3)$, in which two words of a chart are adjacent only if both are constant on a common set $A_r$, has no proper 4-coloring and has a proper 5-coloring (needs `python-sat`) | Theorem 3.3, for $n=3$ |
| The two-edge identity (4.2) on all 360 walks $u,v,w$ of the 22-vertex graph with its two-qubit coloring | Proposition 4.1 |
| The stretched graphs for $l=1,2$: 91 and 160 vertices, 157 and 272 edges, chromatic number four, the edge rule, commutation up to distance $2l+1$ (and not beyond) | Proposition 5.1 |
| $F_\triangle$ (732 edges, 147 in each chart) and $F_*$ (678 edges); complete end grids; the table of the complex $\Lambda$; contractibility of all 54 complexes in twenty orders of deletion; $\chi(F_*)\ge5$ by exhaustive search | Theorem 5.3 |
| $F_R$ for $R=2,3,4$: 234, 288, 342 vertices, the edge rule, compatibility radius exactly $R$, the homomorphism to $F_\triangle$, a proper 5-coloring | Theorem 5.4 |
| The measurements of weight at most two at level two, enumerated for 5 to 8 first-level labels and compared with the counting formula; the value 4,330,533 for 165 labels | Theorem 5.6(e) |

## `computed/`: facts of Section 6

| Folder and program | Fact |
|---|---|
| `edges_678/check_678.py` with `edge_forced_732.txt.gz`, `edge_packing_843.txt.gz` | (C2): every subgraph of $F_\triangle$ with chromatic number five contains the 678 edges of $F_*$; every subgraph of $F_{\rm phys}$ with chromatic number five has at least 678 edges |
| `radius_two_three/verify_square_radius_two_207.py` with `square_radius_two_207.json` | (C3): the graph with 207 vertices and 1,071 edges for $R=2$ |
| `radius_two_three/verify_square_radius_three_replacements.py`, `verify_replacement_deletions.py` with `square_radius_three_replacement_deletions.json.gz` | (C3): the graph with 279 vertices and 1,449 edges for $R=3$ (and a variant with 287 vertices) |
| `radius_two_three/sat_crosscheck.py` (new, needs `python-sat`, about twenty minutes, not part of `run_all.py`) | (C3): an independent confirmation by a SAT solver that the graphs with 207 and 279 vertices have no proper 4-coloring |
| `level_two/family.py`, `count_physical_edges.py`, `count_remaining_edges.py` with `independent_set.json` | (C4): the order and the 2,116,320,480 edges of $H_2$; the independent set of 1,407,626 vertices; the remaining graph with 2,922,907 vertices and 732,610,036 edges |

Fact (C1), the criticality of $F_{22}$, is checked by `../../checks/check_small_graph_22.py`.

Apart from `sat_crosscheck.py`, the programs in `computed/` were written earlier in the project and are included unchanged, with one exception: in `level_two/family.py` and `level_two/count_physical_edges.py` the line that locates `verify_small_graphs.py` now points to the folder of the program itself. That file and `small_graphs_data.json.gz` belong to an earlier supplement of the project and provide the graph $F_{\rm phys}$ for the counts. The comments and messages of these programs use the terms of the time they were written (for example "host" for $F_{\rm phys}$ and "threshold 7" and "threshold 6" for $H_2$ and for the graph after deletion of the independent set).
