# Replay programs

These programs independently recompute finite combinatorial and algebraic facts from the accompanying manuscript (`../paper/main.pdf`), together with descriptive counts reported in the text. Every theorem in the paper has an analytical proof, so no statement depends on these scripts.

The scripts use exact arithmetic and only the Python standard library (Python 3.10 or newer). When `python-sat` is installed, `check_square_graphs.py` uses its SAT solver to check both $\chi(F_{\mathrm{phys}}) = 5$ and $\chi(S_0(K_3)) = 5$; without `python-sat`, $\chi(F_{\mathrm{phys}}) = 5$ falls back to a slower backtracking search and the check of $\chi(S_0(K_3)) = 5$ is skipped.

To run all six verification programs in sequence:

```sh
python3 -B run_all.py
```

| Program | What it recomputes | Location in `paper/main.pdf` | Typical runtime |
|---|---|---|---|
| `check_first_level.py` | Pauli description (2) of the six bases and nine planes; Lemma 2.1 over all $4^6 = 4{,}096$ position choices; Lemma 2.3 over all $1{,}056$ proper $4$-colorings of $K_3^{\times2}$, all $12{,}288$ consecutive-layer pairs, and all $373{,}104$ proper $4$-colorings of $K_3^{\times3}$; Proposition A.1 for $n=3$ ($417$ measurements) and $n=4$ ($1{,}408$ measurements); $W_\lambda$ on the $18$ triads of the square for $\lambda=0,1$ ($27$ wires, compatibility radius $\lambda+2$, no one-in-three selection). | §2.1 (eq. (2), Lem. 2.1); §2.4 (Lem. 2.3); App. A (Lem. 2.3, Prop. A.1); §2.5 (Prop. 2.5) | ~30 s |
| `check_square_graphs.py` | Nine plane identities (1); $S_\lambda(K_3)$ for $\lambda=0,1,2$ ($342+81\lambda$ vertices, $273$ distinct measurements at $\lambda=0$, worked example, compatibility radius $\lambda+2$); the $165$ measurements (12) and $843$ edges of $F_{\mathrm{phys}}$; the $21$ triangles of the tower over $F_{\mathrm{phys}}$ and their mod-$2$ boundary; $\chi(F_{\mathrm{phys}})=5$; $\chi(S_0(K_3))=5$ (requires `python-sat`). | §2.1 (eq. (1)); §2.2; §2.3 (Lem. 2.2); §2.4 (Thm. 2.4); §6.2 (eq. (12), Thm. 6.2); App. E | ~1 min |
| `check_descents.py` | The $16$ choices of positions for each single defective edge of $K_{3,3}$ and their uniform marginals; for $\lambda=1,2$, all $144\lambda$ cuts: every descent map is a proper partial $3$-coloring of $S_\lambda(K_3)$, with each interior vertex cut with frequency $1/(9\lambda)$; for $\lambda=0$, all $288$ cuts and the $108$ words removed with frequency $1/18$; the two identities behind the bounds. | §3.4 (Cor. 3.7); App. F (Lems. F.1, F.2, F.5) | < 1 s |
| `check_small_graph_22.py` | The $22$-vertex, $42$-edge graph $F_{22}$: $\chi=4$, two-qubit strategy, its twelve $4$-cycles and the shore condition in Proposition 6.4, and an independent cellular computation confirming $w^2=0$ (with $K_4$ and $C_5$ as controls). Also verifies vertex- and edge-criticality (not used in the paper). | §6.1 (Prop. 6.1); §6.3 (Prop. 6.4) | < 1 s |
| `check_dimension_three.py` | Peres' $33$ rays ($57$ lines, $40$ triads, $63$ wires, no one-in-three selection); $Q_R$ for $R \le 4$ ($1080 + 567\max\{R-2,0\}$ vertices, compatibility radius $\max\{R,2\}$); the $5$-coloring of Theorem 5.5(a) on $Q_R$ and on all $585$ coarse-grainings of the $40$ triads. | §5.3 (Thm. 5.5) | ~4–5 min |
| `check_rounding.py` | Exact verification of the beta-tail bound $p \ge \theta^d/(768\,d)$ in the proof of Theorem 5.7 for $8 \le d \le 400$. | §5.4 (Thm. 5.7); App. D | a few s |
| `square_data.py` | Exact rational coordinates of the six bases of $\mathbb{C}^4$ (imported by the scripts above). | §2.1 | — |
| `run_all.py` | Executes the six verification scripts in order and stops on the first failure. | — | ~6 min |

`logs/run_all_2026-11-26.log` records a complete passing run of `python3 -B run_all.py` on 26 November 2026 (CPython 3.12.14 on macOS, with `python-sat` installed).

The supplementary note in `../note/` has its own programs in `../note/checks/`, which import `square_data.py`, `check_first_level.py`, `check_square_graphs.py` and `check_small_graph_22.py` from this folder.
