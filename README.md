# Explicit graphs with quantum chromatic number three and large chromatic number

**Author:** Julius A. Zeiss  
**Affiliation:** Institute for Quantum Information, RWTH Aachen University, Otto-Blumenthal-Straße 20, 52074 Aachen, Germany

This repository accompanies the manuscript *"Explicit graphs with quantum chromatic number three and large chromatic number"* and contains the compiled PDF, the complete LaTeX sources, and exact-arithmetic Python scripts that independently replay the finite combinatorial and algebraic facts proved analytically in the paper. The folder `note/` holds a supplementary note with further results of the same project (see [Supplementary note](#supplementary-note) below).

## Abstract

Ciardo and Mottet proved that graphs with quantum chromatic number three can have arbitrarily large classical chromatic number even under $R$-compatible strategies, where measurements within distance $R$ commute. Their construction relies on seven hardness-of-approximation reductions, leading them to ask for explicit families and a direct proof. We provide both by iterating a single graph operation built from the Mermin–Peres square on the triangle. For every level $j \ge 1$ and compatibility parameter $R \ge 1$, both the $j$th iterate and an explicit sparse induced subgraph have quantum chromatic number three and a box-complex-certified classical chromatic number at least $2j+3$. At fixed $j$ the sparse graph's size grows polynomially in $R$, and at fixed $R$ it grows singly exponentially in $j$. These graphs admit an $R$-compatible strategy of local dimension $4^j$ using adaptive Pauli measurements on $2j$ EPR pairs, and both resources are rigid: *every* quantum $3$-coloring of these graphs has local dimension divisible by $4^j$ and compatibility radius at most $\max\{R,2\}$. Over all graphs with an $R$-compatible quantum $3$-coloring, our sparse family achieves the optimal polynomial scaling in $R$ at fixed classical chromatic number; as that number grows, any box-complex lower bound forces exponentially many edges. We further prove that any graph admitting a quantum $3$-coloring in local dimension three is classically $5$-colorable—attained explicitly at every compatibility radius from Peres' rays—while quantum $3$-colorings in local dimension four force classical $6$-colorability, and those in higher local dimensions imply a tighter exponential upper bound on the classical chromatic number. Finally, we exhibit the smallest known witnesses on $22$ and $165$ vertices for four and five classical colors, lift Lalonde's computer-assisted $15$-vertex lower bound to $21$ vertices for five colors, and extend the construction to every nonbipartite target graph.

## Repository contents

```text
.
├── LICENSE.md                   # Which license applies to which files
├── LICENSES/
│   └── CC0-1.0.txt              # Legal code of CC0 1.0 Universal
├── paper/
│   ├── main.pdf                 # Compiled manuscript (71 pages)
│   ├── main.tex                 # Root LaTeX source
│   ├── quantumarticle.cls       # Quantum journal LaTeX class
│   ├── references.bib           # BibTeX bibliography
│   ├── sections/                # Sections 1–8
│   └── appendices/              # Appendices A–F
├── checks/
│   ├── README.md                # Detailed index of replay programs
│   ├── run_all.py               # Runs all six verification scripts in sequence
│   ├── square_data.py           # Exact rational coordinates of the six bases of C^4
│   ├── check_first_level.py     # Pauli square, Lemma 2.1, Lemma 2.3, Prop. A.1, W_lambda
│   ├── check_square_graphs.py   # S_lambda(K_3), F_phys (165 vertices, 843 edges), 21-triangle tower
│   ├── check_descents.py        # Cuts and descent maps of Appendix F at the first level
│   ├── check_small_graph_22.py  # F_22 (22 vertices, 42 edges), strategy, and w^2 = 0 check
│   ├── check_dimension_three.py # Peres' 33 rays (57 lines, 40 triads), Q_R, and 5-coloring
│   ├── check_rounding.py        # Beta-tail inequality in Appendix D for 8 <= d <= 400
│   └── logs/
│       └── run_all_2026-11-26.log
└── note/
    ├── README.md                # Short description of the supplementary note
    ├── main.pdf                 # Compiled note
    ├── main.tex, references.bib, sections/
    └── checks/
        ├── README.md            # Index of the programs of the note
        ├── run_all.py           # Runs all programs of the note in sequence
        ├── check_note.py        # Replay of finite facts proved in Sections 3–5 of the note
        ├── computed/            # Programs and certificates for the facts of Section 6
        └── logs/
            └── run_all_2026-11-27.log
```

## Building the paper

The paper compiles cleanly with `latexmk` (TeX Live / pdfLaTeX + BibTeX):

```sh
cd paper
latexmk -pdf main.tex
```

## Running the verification scripts

All results in the paper have analytical proofs; no theorem depends on a computation by the author (Lalonde's published computer-assisted classification is the sole computer-assisted external import, used only in Proposition 6.5). The Python scripts in `checks/` use exact arithmetic and the Python 3.10+ standard library (with `python-sat` optional for SAT-based chromatic-number checks in `check_square_graphs.py`):

```sh
cd checks
python3 -B run_all.py
```

See [`checks/README.md`](checks/README.md) for a table mapping each script to the corresponding sections, equations, and statements in [`paper/main.pdf`](paper/main.pdf).

## Supplementary note

[`note/main.pdf`](note/main.pdf), *Further results on the square substitution: complete graphs, maps from spheres, and small graphs*, collects results of the same project that the paper does not contain:

- one square substitution over the complete graph $K_n$ has quantum chromatic number $n$ and chromatic number exactly $2n-1$;
- if the box complex of a graph receives an equivariant map from a sphere of dimension $2h+1$, then $j$ substitutions give chromatic number at least $4h+2j+3$; for complete graphs this is $2n-3+2j$ for odd $n$ and $2n-4+2j$ for even $n$;
- the box complex of a graph with a quantum $3$-coloring receives no equivariant map from the two-sphere;
- further small graphs with quantum chromatic number three: chromatic number four at a prescribed compatibility parameter, the $165$ measurements with only $678$ edges, chromatic number five at a prescribed compatibility radius on $234+54(R-2)$ vertices, and a family of graphs of measurements whose first level has $165$ vertices;
- a list of facts known only by computation, with the programs and certificates that verify them.

Every numbered statement of the note has an analytical proof; it imports one lemma of Brakensiek and Guruswami and one theorem of Simons, Tardif and Wehlau, both stated in full. The facts known only by computation are set apart in its last section. See [`note/checks/README.md`](note/checks/README.md) for the programs.

```sh
cd note
latexmk -pdf main.tex
cd checks
python3 -B run_all.py
```

## Earlier version of this repository

Until 27 November 2026 this repository held the verification programs for the version of 19 November 2026 of the paper (the graph on 219 vertices and the small graphs on 28, 29, 32 and 165 vertices). The present paper is organized around the square substitution and has its own programs. The earlier files are kept under the tag [`programs-2026-11-19`](https://github.com/JuliusAZeiss/explicit-quantum-chromatic-gap/tree/programs-2026-11-19).

## License

- The paper and the supplementary note (PDFs, LaTeX sources, bibliography files) and this README are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- The programs, data files and logs in `checks/` and `note/checks/` are dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
- `paper/quantumarticle.cls` is the LaTeX class of the journal Quantum and is distributed under the LaTeX Project Public License, as its header states.

See [`LICENSE.md`](LICENSE.md) for the details.
