#!/usr/bin/env python3
"""
patch_pyrocko.py — Apply a one-line fix to Snuffler's marker click detection.

Snuffler's marker_under_cursor() uses exact tuple matching, which means
phase markers with wildcard channels (e.g. YI.QCB01..*) are never
selectable by click. This script patches it to use match_nslc() instead,
consistent with how the hoovering() function already works.

Run once after creating the conda environment:

    conda activate qcf_eq
    python patch_pyrocko.py
"""

import sys
from pathlib import Path

OLD = (
    "                    for nslc_id in marker_nslc_ids:\n"
    "                        if nslc_id in relevant_nslc_ids:\n"
    "                            return marker\n"
)

NEW = (
    "                    for nslc in relevant_nslc_ids:\n"
    "                        if marker.match_nslc(nslc):\n"
    "                            return marker\n"
)

def find_pile_viewer():
    import pyrocko
    base = Path(pyrocko.__file__).parent
    path = base / "gui" / "snuffler" / "pile_viewer.py"
    if not path.exists():
        sys.exit(f"Could not find pile_viewer.py at {path}")
    return path


def main():
    path = find_pile_viewer()
    src = path.read_text()

    if NEW in src:
        print(f"Already patched: {path}")
        return

    if OLD not in src:
        print(f"WARNING: Expected code not found in {path}")
        print("Pyrocko may have been updated — check pile_viewer.py manually.")
        sys.exit(1)

    patched = src.replace(OLD, NEW)
    path.write_text(patched)
    print(f"Patched: {path}")


if __name__ == "__main__":
    main()
