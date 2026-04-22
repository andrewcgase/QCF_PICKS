#!/usr/bin/env python3
"""
02_picks_to_markers.py — Convert the existing SeisBench/GaMMA pick catalog
                          into a Pyrocko markers file for QC in Snuffler.

Reads:  config.INPUT_EVENTS   (gamma_events.csv)
        config.INPUT_PICKS    (gamma_picks_refined.csv)
Writes: config.MARKERS_FOR_QC (markers/events_to_review.markers)

After running, open Snuffler with the packaged waveforms:

  snuffler waveforms/*.mseed --markers=markers/events_to_review.markers

Snuffler workflow
-----------------
  e / E       — jump to next / previous event
  f           — toggle bandpass filter (2–20 Hz recommended)
  p / s       — add P or S pick at cursor
  Right-click pick → Delete, or rename phase
  File → Save Markers  — save your edits

Save reviewed markers to markers/events_reviewed.markers (or any name),
then run 03_export.py to produce the final catalog tables.

Usage
-----
  conda activate qcf_eq
  python 02_picks_to_markers.py                    # all events with picks files
  python 02_picks_to_markers.py --min-picks 5      # skip poorly-constrained events
  python 02_picks_to_markers.py --only-extracted   # only events with ev*.mseed files
"""

import argparse
from pathlib import Path

import pandas as pd
from pyrocko.gui import marker as pm
from pyrocko.model import Event
from pyrocko import util as putil

import config


def gamma_score_to_kind(score):
    """
    Map a GaMMA score (0–1) to a Pyrocko marker kind index (0–7).
    Higher score → warmer colour.
    """
    if score is None or pd.isna(score):
        return 3   # yellow
    if score >= 0.8:
        return 1   # green
    if score >= 0.5:
        return 3   # yellow
    return 0       # red


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-picks", type=int, default=5,
                        help="Skip events with fewer than N associated picks (default: 5).")
    parser.add_argument("--min-score", type=float, default=0.0,
                        help="Skip events with gamma_score below this value (default: 0).")
    parser.add_argument("--only-extracted", action="store_true",
                        help="Only include events whose ev*.mseed file exists in waveforms/.")
    parser.add_argument("--month", default=None,
                        help="Filter to a single month YYYY-MM (e.g. 2022-04). "
                             "Sets default --out to markers/YYYY-MM.markers.")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Limit to the first N events (for working in batches).")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output markers file (default: markers/events_to_review.markers).")
    args = parser.parse_args()

    config.MARKERS_DIR.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load catalogs
    # -----------------------------------------------------------------------
    events_df = pd.read_csv(config.INPUT_EVENTS)
    picks_df  = pd.read_csv(config.INPUT_PICKS)

    # Filter events
    if "num_picks" in events_df.columns:
        events_df = events_df[events_df["num_picks"] >= args.min_picks]
    if "gamma_score" in events_df.columns and args.min_score > 0:
        events_df = events_df[events_df["gamma_score"] >= args.min_score]
    # Add _month column for filtering
    events_df["_dt"]    = pd.to_datetime(events_df["time"])
    events_df["_month"] = events_df["_dt"].dt.strftime("%Y-%m")

    if args.month:
        events_df = events_df[events_df["_month"] == args.month]
        if args.out is None:
            args.out = config.MARKERS_DIR / f"{args.month}.markers"

    if args.only_extracted:
        # Scan all monthly subdirectories under waveforms/
        extracted = set()
        for sub in config.WAVEFORMS_DIR.iterdir():
            if sub.is_dir():
                for p in sub.glob("ev*.mseed"):
                    extracted.add(int(p.stem.lstrip("ev")))
        # Also catch any flat ev*.mseed still in the root (old-style)
        for p in config.WAVEFORMS_DIR.glob("ev*.mseed"):
            extracted.add(int(p.stem.lstrip("ev")))
        events_df = events_df[events_df["event_index"].isin(extracted)]

    if args.max_events:
        events_df = events_df.head(args.max_events)

    print(f"Events: {len(events_df)}")

    # Keep only picks that belong to filtered events
    valid_idxs = set(events_df["event_index"].astype(int))
    picks_df = picks_df[
        (picks_df["event_index"] >= 0) &
        (picks_df["event_index"].astype(int).isin(valid_idxs))
    ]
    print(f"Picks:  {len(picks_df)}")

    # -----------------------------------------------------------------------
    # Build station_id → SEED network.station.location.channel string
    # gamma station_id format: "YI.QCB01."  (network.station.location.)
    # -----------------------------------------------------------------------
    def parse_nslc(station_id):
        """Return (net, sta, loc, chan) best guess from gamma station_id."""
        parts = str(station_id).rstrip(".").split(".")
        if len(parts) >= 2:
            net, sta = parts[0], parts[1]
            loc = parts[2] if len(parts) > 2 else ""
        else:
            net, sta, loc = "", parts[0], ""
        return net, sta, loc

    # -----------------------------------------------------------------------
    # Build Pyrocko markers
    # -----------------------------------------------------------------------
    markers = []

    for _, ev in events_df.iterrows():
        ev_idx  = int(ev["event_index"])
        ot      = putil.str_to_time(str(ev["time"]).replace("T", " "))

        lat  = float(ev["latitude"])  if "latitude"  in ev else 0.0
        lon  = float(ev["longitude"]) if "longitude" in ev else 0.0
        dep  = float(ev["depth_km"])  * 1000.0 if "depth_km" in ev else 0.0
        mag  = float(ev["magnitude"]) if "magnitude" in ev and not pd.isna(ev.get("magnitude")) else None

        score = ev.get("gamma_score", None)
        kind  = gamma_score_to_kind(score)

        event_obj = Event(
            time=ot,
            lat=lat,
            lon=lon,
            depth=dep,
            magnitude=mag,
            name=f"ev{ev_idx}",
        )

        ev_marker = pm.EventMarker(event=event_obj)
        markers.append(ev_marker)

        # Phase markers for this event's picks
        ev_picks = picks_df[picks_df["event_index"].astype(int) == ev_idx]
        for _, pk in ev_picks.iterrows():
            try:
                net, sta, loc = parse_nslc(pk["station_id"])
                phase = str(pk.get("phase_type", "P")).upper()
                pt    = putil.str_to_time(str(pk["phase_time"]).replace("T", " "))
                nslc  = [(net, sta, loc, "*")]
                pm_marker = pm.PhaseMarker(
                    nslc_ids=nslc,
                    tmin=pt,
                    tmax=pt,
                    phasename=phase,
                    event=event_obj,
                    kind=kind,
                )
                markers.append(pm_marker)
            except Exception:
                continue

    # -----------------------------------------------------------------------
    # Write markers file
    # -----------------------------------------------------------------------
    out_path = args.out if args.out else config.MARKERS_FOR_QC
    pm.save_markers(markers, str(out_path))
    n_ev  = sum(1 for m in markers if isinstance(m, pm.EventMarker))
    n_ph  = sum(1 for m in markers if isinstance(m, pm.PhaseMarker))
    print(f"\nWrote {n_ev} event markers and {n_ph} phase markers to:")
    print(f"  {out_path}")

    print(f"\nLaunch Snuffler:")
    waveform_path = (config.WAVEFORMS_DIR / args.month
                     if args.month else config.WAVEFORMS_DIR)
    print(f"  snuffler {waveform_path}/ "
          f"--markers={out_path}")


if __name__ == "__main__":
    main()
