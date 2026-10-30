#!/usr/bin/env python3
"""
patch_pyrocko.py — Apply known one-line fixes to Snuffler.

Two independent bugs patched here:

1. Marker click detection. marker_under_cursor() uses exact tuple
   matching, which means phase markers with wildcard channels (e.g.
   YI.QCB01..*) are never selectable by click. Patched to use
   match_nslc() instead, consistent with how the hoovering() function
   already works.

2. "Sort by Names" is a no-op. In pile_viewer.py's sort-menu definition,
   every other entry (Distance, Azimuth, Backazimuth) looks up a real
   per-trace sort key via self.station_attrib(...), but "Sort by Names"
   is wired to `lambda tr: ()` -- the same empty tuple for every trace,
   so it can never reorder anything. Confirmed directly in the installed
   pyrocko 2025.12.04 source, not something in this project's own
   scripts. Patched to sort by the trace's own embedded station code
   (`tr.station`), which -- unlike Distance/Azimuth -- doesn't need an
   external stations file or an active event origin at all.

Run once after creating the conda environment:

    conda activate qcf_eq
    python patch_pyrocko.py
"""

import sys
from pathlib import Path

PATCHES = [
    (
        '"Sort by Distance/Azimuth" never sets origin (dist_m always None)',
        (
            "        def s_sortingmode_change(self, ignore=None):\n"
            "            for menuitem, valfunc in self.menuitems_ssorting:\n"
            "                if menuitem.isChecked():\n"
            "                    self._ssort = valfunc\n"
            "\n"
            "            self.sortingmode_change()\n"
        ),
        (
            "        def s_sortingmode_change(self, ignore=None):\n"
            "            for menuitem, valfunc in self.menuitems_ssorting:\n"
            "                if menuitem.isChecked():\n"
            "                    self._ssort = valfunc\n"
            "\n"
            "            if self.get_active_event() is None:\n"
            "                for marker in self.markers:\n"
            "                    if isinstance(marker, EventMarker):\n"
            "                        self.set_active_event_marker(marker)\n"
            "                        break\n"
            "\n"
            "            self.sortingmode_change()\n"
        ),
    ),
    (
        "marker click detection (wildcard-channel markers unselectable)",
        (
            "                    for nslc_id in marker_nslc_ids:\n"
            "                        if nslc_id in relevant_nslc_ids:\n"
            "                            return marker\n"
        ),
        (
            "                    for nslc in relevant_nslc_ids:\n"
            "                        if marker.match_nslc(nslc):\n"
            "                            return marker\n"
        ),
    ),
    (
        '"Sort by Names" is a no-op (always returns the same empty key)',
        (
            "                ('Sort by Names',\n"
            "                    lambda tr: (),\n"
            "                    qg.QKeySequence(qc.Qt.Key_S, qc.Qt.Key_N)),\n"
        ),
        (
            "                ('Sort by Names',\n"
            "                    lambda tr: (tr.station,),\n"
            "                    qg.QKeySequence(qc.Qt.Key_S, qc.Qt.Key_N)),\n"
        ),
    ),
]


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
    changed = False

    for label, old, new in PATCHES:
        if new in src:
            print(f"Already patched: {label}")
            continue
        if old not in src:
            print(f"WARNING: expected code for '{label}' not found in {path}")
            print("Pyrocko may have been updated -- check pile_viewer.py manually.")
            continue
        src = src.replace(old, new)
        changed = True
        print(f"Patched: {label}")

    if changed:
        path.write_text(src)
        print(f"\nWrote {path}")
    else:
        print("\nNo changes needed.")


if __name__ == "__main__":
    main()
