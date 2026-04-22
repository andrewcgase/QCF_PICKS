#!/usr/bin/env python3
"""
04_export.py — Export the QC'd catalog to multiple formats.

Reads:  config.MARKERS_REVIEWED (events_reviewed.markers — saved from Snuffler)
Writes:
  config.CATALOG_CSV    — flat CSV (one row per event)
  config.CATALOG_QML    — QuakeML (ObsPy Catalog)
  config.HYPODD_PHASE   — HypoDD phase.dat (input to ph2dt)
  config.TOMODD_PHASE   — TomoDD absolute travel-time file

Usage
-----
  conda activate obspy        # environment with obspy + pyrocko
  python 04_export.py

  # Export from a specific markers file (e.g. after a second QC pass):
  python 04_export.py --markers markers/my_revision.markers
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.core.event import (
    Catalog, Event, Origin, Pick, Arrival, WaveformStreamID,
    QuantityError, ResourceIdentifier
)
from obspy.geodetics import gps2dist_azimuth

try:
    from pyrocko.gui.snuffler.marker import load_markers, EventMarker, PhaseMarker
    from pyrocko.util import time_to_str
except ImportError:
    sys.exit("Pyrocko not found — conda install -c pyrocko pyrocko")

import config


# ---------------------------------------------------------------------------
# Load reviewed markers
# ---------------------------------------------------------------------------

def load_reviewed_catalog(markers_path):
    """
    Parse Pyrocko marker file → (events_df, picks_df).

    events_df columns: event_id, time, latitude, longitude, depth_km, magnitude
    picks_df  columns: event_id, station_id, network, station, phase, pick_time, score
    """
    markers = load_markers(str(markers_path))

    # Link phase markers to their events via hash (not done automatically on load)
    hash_to_event = {}
    for m in markers:
        if isinstance(m, EventMarker):
            hash_to_event[m.get_event_hash()] = m.get_event()
    for m in markers:
        if isinstance(m, PhaseMarker) and m._event_hash in hash_to_event:
            m.set_event(hash_to_event[m._event_hash])

    ev_rows, ph_rows = [], []
    event_map = {}   # pyrocko event name → row index

    for m in markers:
        if isinstance(m, EventMarker):
            e = m.get_event()
            ev_rows.append({
                "event_id":  e.name,
                "time":      time_to_str(e.time),
                "latitude":  e.lat,
                "longitude": e.lon,
                "depth_km":  (e.depth or 0.0) / 1000.0,
                "magnitude": e.magnitude if e.magnitude != 999.0 else None,
            })
            event_map[e.name] = e

        elif isinstance(m, PhaseMarker):
            ev = m.get_event()
            if ev is None:
                continue
            nslc = list(m.nslc_ids)[0] if m.nslc_ids else ("", "", "", "")
            ph_rows.append({
                "event_id":   ev.name,
                "network":    nslc[0],
                "station":    nslc[1],
                "location":   nslc[2],
                "channel":    nslc[3],
                "station_id": f"{nslc[0]}.{nslc[1]}.{nslc[2]}.",
                "phase":      m._phasename or "?",
                "pick_time":  time_to_str(m.tmin),
                "score":      None,
            })

    events_df = pd.DataFrame(ev_rows)
    picks_df  = pd.DataFrame(ph_rows)
    return events_df, picks_df


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def export_csv(events_df, picks_df):
    events_df.to_csv(config.CATALOG_CSV, index=False,
                     float_format="%.4f")
    print(f"CSV written: {config.CATALOG_CSV}  ({len(events_df)} events)")


# ---------------------------------------------------------------------------
# QuakeML export
# ---------------------------------------------------------------------------

def export_quakeml(events_df, picks_df):
    catalog = Catalog()

    for _, ev_row in events_df.iterrows():
        event_id = ev_row["event_id"]
        ot = UTCDateTime(ev_row["time"])

        origin = Origin(
            time=ot,
            latitude=ev_row["latitude"],
            longitude=ev_row["longitude"],
            depth=ev_row["depth_km"] * 1000.0,
            resource_id=ResourceIdentifier(f"quakeml:qcf/{event_id}/origin"),
        )

        ev_picks_df = picks_df[picks_df["event_id"] == event_id]
        obspy_picks, arrivals = [], []

        for i, ph_row in ev_picks_df.iterrows():
            wfid = WaveformStreamID(
                network_code=ph_row["network"],
                station_code=ph_row["station"],
                location_code=ph_row["location"],
                channel_code=ph_row.get("channel", "*"),
            )
            pick = Pick(
                resource_id=ResourceIdentifier(f"quakeml:qcf/{event_id}/pick/{i}"),
                time=UTCDateTime(ph_row["pick_time"]),
                waveform_id=wfid,
                phase_hint=ph_row["phase"],
            )
            arrival = Arrival(
                resource_id=ResourceIdentifier(f"quakeml:qcf/{event_id}/arrival/{i}"),
                pick_id=pick.resource_id,
                phase=ph_row["phase"],
            )
            obspy_picks.append(pick)
            arrivals.append(arrival)

        origin.arrivals = arrivals

        obspy_event = Event(
            resource_id=ResourceIdentifier(f"quakeml:qcf/{event_id}"),
            origins=[origin],
            picks=obspy_picks,
            preferred_origin_id=origin.resource_id,
        )
        catalog.append(obspy_event)

    catalog.write(str(config.CATALOG_QML), format="QUAKEML")
    print(f"QuakeML written: {config.CATALOG_QML}  ({len(catalog)} events)")


# ---------------------------------------------------------------------------
# HypoDD phase.dat export
# ---------------------------------------------------------------------------

def export_hypodd(events_df, picks_df):
    """
    Write HypoDD phase.dat format (input to ph2dt).

    Format:
      # YR MO DY HR MN SC  LAT       LON        DEP   MAG  EH  EZ  RMS  ID
      # STA  TT  W  PHA
    """
    lines = []

    for _, ev_row in events_df.iterrows():
        event_id = ev_row["event_id"]
        ev_id_int = _event_id_int(event_id)
        t = UTCDateTime(ev_row["time"])
        mag = ev_row["magnitude"]
        mag_str = f"{mag:.1f}" if pd.notna(mag) else "0.0"

        lines.append(
            f"# {t.year:4d} {t.month:02d} {t.day:02d} "
            f"{t.hour:02d} {t.minute:02d} {t.second + t.microsecond/1e6:05.2f} "
            f"{ev_row['latitude']:8.4f} {ev_row['longitude']:9.4f} "
            f"{ev_row['depth_km']:6.3f} {mag_str:>5} "
            f"0.00 0.00 0.00 {ev_id_int:>8d}"
        )

        ev_picks = picks_df[picks_df["event_id"] == event_id]
        ot = UTCDateTime(ev_row["time"])

        for _, ph_row in ev_picks.iterrows():
            tt = UTCDateTime(ph_row["pick_time"]) - ot
            if tt < 0 or tt > 1800:
                continue   # sanity check
            weight = _score_to_weight(ph_row.get("score"))
            phase  = ph_row["phase"].upper()
            lines.append(f"{ph_row['station']:<6s} {tt:7.3f} {weight:.1f} {phase}")

    with open(config.HYPODD_PHASE, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    n_events = len(events_df)
    print(f"HypoDD phase.dat written: {config.HYPODD_PHASE}  ({n_events} events)")


# ---------------------------------------------------------------------------
# TomoDD absolute travel-time export
# ---------------------------------------------------------------------------

def export_tomodd(events_df, picks_df):
    """
    Write TomoDD absolute travel-time file.

    Format (same layout as HypoDD phase.dat but TomoDD reads it directly):
      # EV_ID  OT_ISO  LAT  LON  DEP_KM
      STA  PHASE  TT_ABS  WEIGHT
    """
    lines = []

    for _, ev_row in events_df.iterrows():
        event_id = ev_row["event_id"]
        ev_id_int = _event_id_int(event_id)
        t = UTCDateTime(ev_row["time"])
        ot_str = t.strftime("%Y-%m-%dT%H:%M:%S.%f")

        lines.append(
            f"# {ev_id_int:>8d}  {ot_str}  "
            f"{ev_row['latitude']:9.4f}  {ev_row['longitude']:10.4f}  "
            f"{ev_row['depth_km']:7.3f}"
        )

        ev_picks = picks_df[picks_df["event_id"] == event_id]
        ot = UTCDateTime(ev_row["time"])

        for _, ph_row in ev_picks.iterrows():
            tt = UTCDateTime(ph_row["pick_time"]) - ot
            if tt < 0 or tt > 1800:
                continue
            weight = _score_to_weight(ph_row.get("score"))
            phase  = ph_row["phase"].upper()
            lines.append(f"{ph_row['station']:<6s}  {phase}  {tt:9.4f}  {weight:.3f}")

    with open(config.TOMODD_PHASE, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    n_events = len(events_df)
    print(f"TomoDD phase file written: {config.TOMODD_PHASE}  ({n_events} events)")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _event_id_int(event_id_str):
    """Extract integer from 'evXXXX' or return hash."""
    try:
        return int(str(event_id_str).lstrip("ev"))
    except ValueError:
        return abs(hash(event_id_str)) % 10**8


def _score_to_weight(score):
    """Map PhaseNet confidence score to pick weight [0–1]."""
    if score is None or (isinstance(score, float) and np.isnan(score)):
        return 1.0
    s = float(score)
    if s >= 0.9:
        return 1.0
    elif s >= 0.7:
        return 0.8
    elif s >= 0.5:
        return 0.6
    elif s >= 0.3:
        return 0.4
    return 0.2


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--markers", type=Path, default=config.MARKERS_REVIEWED,
        help=f"Reviewed Pyrocko markers file (default: {config.MARKERS_REVIEWED})"
    )
    args = parser.parse_args()

    if not args.markers.exists():
        sys.exit(
            f"{args.markers} not found.\n"
            "Complete the Snuffler QC step and save markers before exporting."
        )

    config.CATALOGS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading reviewed markers from {args.markers}")
    events_df, picks_df = load_reviewed_catalog(args.markers)
    print(f"  {len(events_df)} events,  {len(picks_df)} picks")

    if events_df.empty:
        sys.exit("No events found in marker file.")

    export_csv(events_df, picks_df)
    export_quakeml(events_df, picks_df)
    export_hypodd(events_df, picks_df)
    export_tomodd(events_df, picks_df)

    print("\nAll exports complete.")
    print(f"  CSV     : {config.CATALOG_CSV}")
    print(f"  QuakeML : {config.CATALOG_QML}")
    print(f"  HypoDD  : {config.HYPODD_PHASE}")
    print(f"  TomoDD  : {config.TOMODD_PHASE}")


if __name__ == "__main__":
    main()
