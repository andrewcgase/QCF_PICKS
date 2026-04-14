#!/usr/bin/env python3
"""
03_snuffler_prep.py — Convert the associated catalog to Pyrocko marker format
                       for manual QC in Snuffler.

Reads:  config.EVENTS_ASSOC    (events_associated.csv)
        config.PICKS_ASSOC     (picks_associated.csv)
Writes: config.MARKERS_FOR_QC  (events_to_review.markers)

After running this script, open Snuffler to review picks:

  snuffler \\
    /path/to/DATA/corrected_clock_final/QCBxx/HHZ/*.msd \\
    /path/to/DATA/raw_data/QCBxx/Data/*/*/*.msd \\
    /path/to/DATA/onshore/NET/STA/Data/*/*/*.msd \\
    --markers=markers/events_to_review.markers

In Snuffler:
  - Use 'e' key to jump between events
  - Right-click picks to edit phase name / delete
  - 'p'/'s' keys to add P/S picks
  - File → Save Markers to save your changes as markers/events_reviewed.markers

Then run 04_export.py to convert the reviewed markers to output catalogs.

Usage
-----
  conda activate pyrocko      # environment with pyrocko installed
  python 03_snuffler_prep.py
"""

import sys
import numpy as np
import pandas as pd
from obspy import UTCDateTime

try:
    from pyrocko.gui.snuffler.marker import EventMarker, PhaseMarker, save_markers
    from pyrocko import model as pmodel
except ImportError:
    sys.exit(
        "Pyrocko not found.  Install it with:\n"
        "  conda install -c pyrocko pyrocko\n"
        "or activate the correct conda environment."
    )

import config


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def utc_to_pyrocko(time_str):
    """Convert ISO time string to Pyrocko float time (seconds since epoch)."""
    from pyrocko.util import str_to_time
    # Normalise format
    t = str(time_str).replace("T", " ").rstrip("Z")
    return str_to_time(t, format="%Y-%m-%d %H:%M:%S.OPTFRAC")


def station_id_to_nslc(station_id):
    """
    Convert GaMMA station_id 'NET.STA.LOC.' to (net, sta, loc, '*') tuple.
    Snuffler will match against all channels.
    """
    parts = station_id.rstrip(".").split(".")
    net = parts[0] if len(parts) > 0 else ""
    sta = parts[1] if len(parts) > 1 else ""
    loc = parts[2] if len(parts) > 2 else ""
    return (net, sta, loc, "*")


# Phase kind codes used by Snuffler (0 = default colour)
PHASE_KIND = {"P": 0, "S": 1}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    config.MARKERS_DIR.mkdir(parents=True, exist_ok=True)

    for f in [config.EVENTS_ASSOC, config.PICKS_ASSOC]:
        if not f.exists():
            sys.exit(f"Run 02_associate.py first — {f} not found")

    events_df = pd.read_csv(config.EVENTS_ASSOC)
    picks_df  = pd.read_csv(config.PICKS_ASSOC)

    # Only keep picks associated to events
    picks_df = picks_df[picks_df["event_index"] >= 0].copy()
    print(f"Events: {len(events_df)},  associated picks: {len(picks_df)}")

    markers = []

    for _, ev in events_df.iterrows():
        ev_time = utc_to_pyrocko(ev["time"])
        ev_depth_m = float(ev["depth_km"]) * 1000.0

        pyrocko_event = pmodel.Event(
            lat=float(ev["latitude"]),
            lon=float(ev["longitude"]),
            depth=ev_depth_m,
            time=ev_time,
            name=f"ev{int(ev['event_index'])}",
            magnitude=None if ev.get("magnitude", 999) == 999 else float(ev["magnitude"]),
        )
        event_marker = EventMarker(pyrocko_event)
        markers.append(event_marker)

        # Phase markers linked to this event
        ev_picks = picks_df[picks_df["event_index"] == ev["event_index"]]
        for _, pick in ev_picks.iterrows():
            phase = str(pick["phase_type"]).upper()
            nslc = station_id_to_nslc(pick["station_id"])
            pick_time = utc_to_pyrocko(pick["phase_time"])

            phase_marker = PhaseMarker(
                nslc_ids=[nslc],
                tmin=pick_time,
                tmax=pick_time,
                kind=PHASE_KIND.get(phase, 0),
                phasename=phase,
                event=pyrocko_event,
                event_time=ev_time,
            )
            markers.append(phase_marker)

    save_markers(markers, str(config.MARKERS_FOR_QC))
    print(f"\nMarker file written: {config.MARKERS_FOR_QC}")
    print(f"  {len(events_df)} event markers")
    print(f"  {len(picks_df)} phase markers")

    # Print Snuffler launch command
    print("\n--- Launch Snuffler for manual QC ---")
    print("Add your data directories to the command below, then run it:\n")
    print(f"snuffler \\")
    print(f"  {config.OBS_CORRECTED}/QCBxx/HHZ/*.msd \\")
    print(f"  {config.OBS_RAW}/QCBxx/Data/*/*/*.msd \\")
    print(f"  {config.ONSHORE}/NET/STA/Data/*/*/*.msd \\")
    print(f"  --markers={config.MARKERS_FOR_QC}\n")
    print("When done reviewing, save markers to:")
    print(f"  {config.MARKERS_REVIEWED}")
    print("Then run: python 04_export.py")


if __name__ == "__main__":
    main()
