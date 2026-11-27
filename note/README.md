# Further results on the square substitution

Supplementary note to the article in `../paper/`. It collects results of the same project that the article does not contain: one substitution over a complete graph (chromatic number exactly 2n-1), maps from spheres, the two-sphere, further small graphs, and facts known only by computation.

- `main.pdf` (28 pages), `main.tex`, `sections/`, `references.bib`. Build with `latexmk -pdf main.tex`.
- `checks/check_note.py` replays finite facts that are proved in Sections 3-5 of the note. It imports helpers from `../../checks`.
- `checks/computed/` holds the programs and certificates for Section 6 (facts with no analytical proof): `edges_678/`, `radius_two_three/`, `level_two/`. They are copies of earlier programs of the project; only the two lines that locate `verify_small_graphs.py` in `level_two/` were changed, and `radius_two_three/sat_crosscheck.py` is new.
- `checks/run_all.py` runs everything (about five minutes); `checks/logs/run_all_2026-11-27.log` is a complete passing run (CPython 3.12.14, with `python-sat`).

Status on 27 November 2026: complete draft. The proofs were checked line by line and then read independently a second time. The second read found no error; the one small gap and the minor points it raised are fixed.
