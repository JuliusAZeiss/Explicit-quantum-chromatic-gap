"""Verify all single-vertex deletions of the 279- and 287-vertex examples."""
import gzip
import json
from pathlib import Path

from verify_square_radius_three_replacements import construct


def verify_deletions(order, rows):
    graph = construct(order)
    if len(rows) != order:
        raise ValueError("Expected one coloring for each vertex deletion")
    for removed, colors in enumerate(rows):
        if len(colors) != order or colors[removed] != "x" or colors.count("x") != 1:
            raise ValueError(f"Invalid omission at vertex {removed}")
        if not set(colors) <= set("0123x"):
            raise ValueError("Unexpected color")
        for u, v in graph["edges"]:
            if colors[u] == colors[v]:
                raise ValueError(f"Monochromatic edge {u, v} after deleting {removed}")
    print(f"All {order} deletions have proper classical four-colorings.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    data = Path(__file__).with_name("square_radius_three_replacement_deletions.json.gz")
    with gzip.open(data, "rt") as handle:
        rows_by_order = json.load(handle)
    for order in (279, 287):
        verify_deletions(order, rows_by_order[str(order)])
