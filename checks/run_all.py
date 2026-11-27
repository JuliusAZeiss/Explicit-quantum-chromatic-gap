"""Run the replay programs of the paper and stop at the first failure."""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for name in ("check_first_level.py", "check_square_graphs.py", "check_descents.py",
             "check_small_graph_22.py", "check_dimension_three.py", "check_rounding.py"):
    print(f"== {name}", flush=True)
    subprocess.run([sys.executable, "-B", str(HERE / name)], cwd=HERE, check=True)
print("all replay programs passed")
