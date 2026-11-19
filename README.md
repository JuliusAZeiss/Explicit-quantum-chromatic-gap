# Verification programs for "Explicit graphs with quantum chromatic number three and large chromatic number"

Author of the paper: Julius A. Zeiss.

Every theorem, proposition and lemma of the paper is proved by hand. The
programs here are independent checks of the finite statements. They use only
the Python standard library (Python 3.8 or later) and exact integer
arithmetic. No numerical solver is involved.

## Run

    python3 verify_graph_219.py
    python3 verify_small_graphs.py

Each program prints what it has verified and exits with status 0. Together
they take a few seconds.

## Files

| File | Content |
|---|---|
| `verify_graph_219.py` | Checks the graph F on 219 vertices: its 1353 edges, a proper 5-coloring, the quantum 3-coloring in dimension four with commuting measurements on edges, and the mod 2 chain certificate that excludes a 4-coloring. |
| `graph_219_certificate.json.gz` | Vertices, edges, 5-coloring and chains for F. |
| `verify_small_graphs.py` | Rebuilds the graphs on 28, 29, 32 and 165 vertices from the Mermin–Peres square, checks their vertex and edge counts, chromatic numbers and vertex-criticality, the parity rule for the 39-vertex gadget, the lemmas on four-colorings of grids and triples, and the counts for the first-level graph on 597 vertices. |
| `small_graphs_data.json.gz` | The 5-coloring of the 165-vertex graph and a 4-coloring of each of its vertex-deleted subgraphs. |

The programs and their output are described in the appendix "A topological
certificate and computer checks" of the paper.
