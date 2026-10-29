#!/usr/bin/env python3
"""
02_nor_to_markers.py — Convert the TOQUES SEISAN Nordic (.nor) catalog into
                        a Pyrocko markers file for QC in Snuffler.

Station codes in the Nordic file are CB01–CB28; they are mapped to QCB01–QCB28
so phase markers align with the waveform traces in Snuffler.

Reads:  config.TOQUES_NOR            (QCF_Pickblue_NLL_liu.nor)
Writes: config.MARKERS_TOQUES_FOR_QC (markers/toques_to_review.markers)

Snuffler workflow
-----------------
  e / E       — jump to next / previous event
  f           — toggle bandpass filter (2–20 Hz recommended)
  p / s       — add P or S pick at cursor
  Right-click pick → Delete, or rename phase
  File → Save Markers  — save your edits

Save reviewed markers to markers/toques_reviewed.markers (or any name),
then run 03_export.py --toques to produce the final catalog tables.

Usage
-----
  conda activate qcf_eq
  python 02_nor_to_markers.py
  python 02_nor_to_markers.py --only-extracted   # only events with toq*.mseed
  python 02_nor_to_markers.py --month 2021-10
  python 02_nor_to_markers.py --max-events 20    # for testing
"""

import argparse
from pathlib import Path

from pyrocko.gui import marker as pm
from pyrocko.model import Event
from pyrocko import util as putil

import config
import nordic


def nordic_weight_to_kind(weight_str):
    """Map Nordic pick weight (0–4, 9) to Pyrocko marker kind (colour)."""
    w = weight_str.strip() if weight_str else ""
    if w in ("", "0", "1"):
        return 1   # green — reliable pick
    if w in ("2", "3"):
        return 3   # yellow — use with caution
    return 0       # red — weight 4/9, do not use


def obs_station(cb_code):
    """Convert Nordic CB01 station code to QCB01 to match raw data and waveforms."""
    if cb_code.startswith("CB"):
        return "Q" + cb_code
    return cb_code


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--month", default=None,
                        help="Filter to YYYY-MM (e.g. 2021-10). "
                             "Sets default output to markers/YYYY-MM_toques.markers.")
    parser.add_argument("--only-extracted", action="store_true",
                        help="Only include events whose toq*.mseed exists in waveforms_toques/.")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Limit to first N events (for testing).")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output markers file path.")
    args = parser.parse_args()

    config.MARKERS_DIR.mkdir(parents=True, exist_ok=True)

    events_all = nordic.read_nor_file(str(config.TOQUES_NOR))
    print(f"Loaded {len(events_all)} events from {config.TOQUES_NOR.name}")

    # Tag each event with its original 1-based catalog index *before* any
    # filtering. That index is what 01_extract_toques.py baked into the
    # waveform filenames (toq{idx:04d}.mseed, numbered across the full
    # 574-event catalog, not per month). Re-enumerating from 1 after
    # --month already shrank the list (the old bug here) silently produces
    # local indices that only coincidentally match the real filenames for
    # the catalog's first chronological chunk (Sep 2021) -- every other
    # month comes up empty, which is exactly what Anne hit with
    # `--only-extracted --month 2021-12` (real files toq0126-toq0169, but
    # the re-enumerated local indices were 1-44).
    events_raw = list(enumerate(events_all, start=1))  # [(idx, event), ...]

    # Filter by month
    if args.month:
        yr_s, mo_s = args.month.split("-")
        events_raw = [(i, e) for i, e in events_raw
                      if e["header"]["year"] == int(yr_s)
                      and e["header"]["month"] == int(mo_s)]
        if args.out is None:
            args.out = config.MARKERS_DIR / f"{args.month}_toques.markers"

    # Filter to events with extracted waveforms
    if args.only_extracted:
        extracted = set()
        wdir = config.WAVEFORMS_TOQUES_DIR
        if wdir.exists():
            for sub in wdir.iterdir():
                if sub.is_dir():
                    for p in sub.glob("toq*.mseed"):
                        extracted.add(int(p.stem.lstrip("toq")))
            for p in wdir.glob("toq*.mseed"):
                extracted.add(int(p.stem.lstrip("toq")))
        events_raw = [(i, e) for i, e in events_raw if i in extracted]

    if args.max_events:
        events_raw = events_raw[:args.max_events]

    print(f"Events after filtering: {len(events_raw)}")

    markers = []

    for ev_idx, ev in events_raw:
        h = ev["header"]
        try:
            dt = nordic.event_datetime(h)
            ot = putil.str_to_time(dt.strftime("%Y-%m-%d %H:%M:%S.%f"))
        except Exception:
            continue

        lat = h.get("lat") or 0.0
        lon = h.get("lon") or 0.0
        dep = (h.get("depth") or 0.0) * 1000.0
        mag = h.get("mag1") or None

        event_obj = Event(
            time=ot,
            lat=lat,
            lon=lon,
            depth=dep,
            magnitude=mag,
            name=f"toq{ev_idx:04d}",
        )
        markers.append(pm.EventMarker(event=event_obj))

        for pk in ev["picks"]:
            try:
                sta_cb = pk.get("station", "")
                if not sta_cb:
                    continue
                phase = (pk.get("phase") or "").upper()[:1]
                if phase not in ("P", "S"):
                    continue

                tt = nordic.pick_travel_time(h, pk)
                if tt is None or tt < 0 or tt > 1800:
                    continue

                pt = ot + tt
                sta = obs_station(sta_cb)
                kind = nordic_weight_to_kind(pk.get("weight", ""))

                markers.append(pm.PhaseMarker(
                    nslc_ids=[(config.OBS_NETWORK, sta, "", "*")],
                    tmin=pt,
                    tmax=pt,
                    phasename=phase,
                    event=event_obj,
                    kind=kind,
                ))
            except Exception:
                continue

    out_path = args.out if args.out else config.MARKERS_TOQUES_FOR_QC
    pm.save_markers(markers, str(out_path))

    n_ev = sum(1 for m in markers if isinstance(m, pm.EventMarker))
    n_ph = sum(1 for m in markers if isinstance(m, pm.PhaseMarker))
    print(f"\nWrote {n_ev} event markers and {n_ph} phase markers to:")
    print(f"  {out_path}")

    waveform_path = config.WAVEFORMS_TOQUES_DIR
    if args.month:
        waveform_path = waveform_path / args.month
    stations_path = Path(__file__).parent / "stations.pf"
    stations_flag = f" --stations={stations_path}" if stations_path.exists() else ""
    if not stations_flag:
        print(f"\nNOTE: {stations_path} not found -- run make_stations_file.py first, "
              f"or Snuffler's 'sort by distance' will silently do nothing (it needs "
              f"station coordinates to compute distance).")
    print(f"\nLaunch Snuffler:")
    print(f"  snuffler {waveform_path}/ --markers={out_path}{stations_flag}")


if __name__ == "__main__":
    main()
