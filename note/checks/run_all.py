"""Run the programs of the note in sequence and stop at the first failure.

First the replay program check_note.py (facts with analytical proofs), then
the programs in computed/ (facts of Section 6, known only by computation).
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = [
    (".", ["check_note.py"]),
    ("computed/edges_678", ["check_678.py"]),
    ("computed/radius_two_three", ["verify_square_radius_two_207.py",
                                   "--check-deletions", "square_radius_two_207.json"]),
    ("computed/radius_two_three", ["verify_square_radius_three_replacements.py", "--order", "279"]),
    ("computed/radius_two_three", ["verify_replacement_deletions.py"]),
    ("computed/level_two", ["family.py"]),
    ("computed/level_two", ["count_physical_edges.py"]),
    ("computed/level_two", ["count_remaining_edges.py"]),
]
for folder, command in RUNS:
    label = command[0] if folder == "." else f"{folder}/{command[0]}"
    print(f"== {label} {' '.join(command[1:])}".rstrip(), flush=True)
    subprocess.run([sys.executable, "-B", *command], cwd=HERE / folder, check=True)
print("all programs of the note passed")
