#!/usr/bin/env python3
"""
01_repick.py — Re-pick P and S phases for events in the detection catalog
               using SeisBench PhaseNet with PickBlue weights.

Reads:  config.INPUT_EVENTS   (gamma_events.csv — 10,471 detected events)
Writes: config.PICKS_REPICKED (picks_repicked.csv)

Strategy
--------
Events are grouped by calendar day.  Each station's daily miniSEED file is
loaded once per day, then all event windows falling within that day are
sliced from it and passed to PhaseNet in a single batch.  This keeps disk
I/O proportional to (stations × days) rather than (stations × events).

Data hierarchy
--------------
  OBS raw    : raw_data/QCBxx/Data/YYYY/DDD/YI_QCBxx_CHAN__YYYY_DDD.msd
  OBS clk-corr: corrected_clock_final/QCBxx/CHAN/YI_QCBxx_CHAN__YYYY_DDD.c1.msd
  Onshore    : onshore/NET/STA/Data/YYYY/DDD/NET_STA_CHAN__YYYY_DDD.msd

Usage
-----
  conda activate obspy        # environment with obspy + seisbench
  python 01_repick.py
  python 01_repick.py --max-events 100   # test run on first 100 events
"""

import argparse
import glob
import sys
from collections import defaultdict
from math import radians, sin, cos, sqrt, atan2
from pathlib import Path

import numpy as np
import pandas as pd
import obspy
from obspy import UTCDateTime, read, read_inventory
import seisbench.models as sbm
from tqdm import tqdm

import config


# ---------------------------------------------------------------------------
# Haversine distance
# ---------------------------------------------------------------------------

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = radians(lat1), radians(lat2)
    a = sin(radians(lat2 - lat1) / 2)**2 + \
        cos(p1) * cos(p2) * sin(radians(lon2 - lon1) / 2)**2
    return R * 2 * atan2(sqrt(a), sqrt(1 - a))


# ---------------------------------------------------------------------------
# Station inventory
# ---------------------------------------------------------------------------

def load_stations():
    """
    Return a DataFrame with columns: network, station, latitude, longitude,
    elevation_m, data_path_root.

    Merges OBS metadata (EarthScope CSV) with onshore metadata (StationXML).
    """
    rows = []

    # OBS stations
    obs_meta = pd.read_csv(config.OBS_META_CSV)
    for _, row in obs_meta.iterrows():
        sta = row["Station"]
        if sta in config.SKIP_STATIONS:
            continue
        if sta in config.CLOCK_CORRECTED_STATIONS:
            data_root = config.OBS_CORRECTED / sta
        else:
            data_root = config.OBS_RAW / sta / "Data"
        rows.append({
            "network": config.OBS_NETWORK,
            "station": sta,
            "latitude": row["Latitude"],
            "longitude": row["Longitude"],
            "elevation_m": row["Depth (m)"],   # negative = below sea level
            "data_root": data_root,
            "corrected": sta in config.CLOCK_CORRECTED_STATIONS,
        })

    # Onshore stations (from StationXML)
    if config.ONSHORE_XML.exists():
        inv = read_inventory(str(config.ONSHORE_XML))
        for network in inv:
            for station in network:
                sta_dir = config.ONSHORE / network.code / station.code / "Data"
                if not sta_dir.exists():
                    continue
                rows.append({
                    "network": network.code,
                    "station": station.code,
                    "latitude": station.latitude,
                    "longitude": station.longitude,
                    "elevation_m": station.elevation,
                    "data_root": sta_dir,
                    "corrected": False,
                })
    else:
        print(f"WARNING: {config.ONSHORE_XML} not found — onshore stations skipped")

    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} stations  "
          f"({(df.network == config.OBS_NETWORK).sum()} OBS, "
          f"{(df.network != config.OBS_NETWORK).sum()} onshore)")
    return df


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_daily_files(data_root, network, station, year, jday, corrected=False):
    """
    Return a list of miniSEED file paths for the given station/day.
    Tries HHZ/HH1/HH2 (OBS) and BHZ/BHN/BHE/HHZ/HHN/HHE (onshore).
    """
    files = []
    tag = f"{year}_{jday:03d}"

    if corrected:
        # corrected_clock_final/QCBxx/CHAN/YI_QCBxx_CHAN__YYYY_DDD*.c1.msd
        for chan in ["HHZ", "HH1", "HH2", "EDH", "HDH"]:
            chan_dir = data_root / chan
            if chan_dir.exists():
                files += glob.glob(str(chan_dir / f"{network}_{station}_{chan}__{tag}*.msd"))
    elif network == config.OBS_NETWORK:
        # raw_data/QCBxx/Data/YYYY/DDD/YI_QCBxx_CHAN__YYYY_DDD.msd
        day_dir = data_root / str(year) / f"{jday:03d}"
        if day_dir.exists():
            for chan in ["HHZ", "HH1", "HH2", "HDH"]:
                files += glob.glob(str(day_dir / f"{network}_{station}_{chan}__{tag}*.msd"))
    else:
        # onshore: onshore/NET/STA/Data/YYYY/DDD/NET_STA_CHAN__YYYY_DDD.msd
        day_dir = data_root / str(year) / f"{jday:03d}"
        if day_dir.exists():
            for chan in ["HHZ", "HHN", "HHE", "HH1", "HH2",
                         "BHZ", "BHN", "BHE", "BH1", "BH2"]:
                files += glob.glob(str(day_dir / f"{network}_{station}_{chan}__{tag}*.msd"))

    return files


# ---------------------------------------------------------------------------
# Core picking
# ---------------------------------------------------------------------------

def pick_windows(model, stream, event_rows):
    """
    For each event in event_rows, slice a window from `stream`,
    run PhaseNet, and return a list of pick dicts.
    """
    picks_out = []

    for _, ev in event_rows.iterrows():
        t0 = UTCDateTime(ev["time"]) - config.WINDOW_BEFORE_S
        t1 = UTCDateTime(ev["time"]) + config.WINDOW_AFTER_S

        win = stream.slice(starttime=t0, endtime=t1)
        if not win or max(tr.stats.npts for tr in win) < 10:
            continue

        try:
            win.detrend("demean")
            result = model.classify(
                win,
                batch_size=config.SEISBENCH_BATCH,
                P_threshold=config.PICK_THRESHOLD_P,
                S_threshold=config.PICK_THRESHOLD_S,
            )
            for pick in result.picks:
                net  = win[0].stats.network
                sta  = win[0].stats.station
                loc  = win[0].stats.location
                picks_out.append({
                    "station_id":  f"{net}.{sta}.{loc}.",
                    "phase_time":  str(pick.peak_time),
                    "phase_score": round(float(pick.peak_value), 6),
                    "phase_type":  pick.phase.lower(),
                    "event_index": int(ev["event_index"]),
                })
        except Exception as exc:
            pass   # skip windows that fail (gaps, short traces, etc.)

    return picks_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-events", type=int, default=None,
                        help="Process only the first N events (testing).")
    args = parser.parse_args()

    config.CATALOGS_DIR.mkdir(parents=True, exist_ok=True)

    # Load model
    print(f"Loading PhaseNet ({config.SEISBENCH_WEIGHTS})...")
    model = sbm.PhaseNet.from_pretrained(config.SEISBENCH_WEIGHTS)
    model.eval()

    # Load catalog
    events = pd.read_csv(config.INPUT_EVENTS)
    events["time"] = pd.to_datetime(events["time"])
    if args.max_events:
        events = events.head(args.max_events)
        print(f"TEST MODE: processing first {args.max_events} events")
    print(f"Events to process: {len(events)}")

    # Load station inventory
    stations = load_stations()

    # Group events by (year, julian_day)
    events["_year"] = events["time"].dt.year
    events["_jday"] = events["time"].dt.day_of_year
    day_groups = events.groupby(["_year", "_jday"])
    print(f"Calendar days spanned: {len(day_groups)}")

    all_picks = []
    skipped_stations = 0

    for (year, jday), day_events in tqdm(day_groups, desc="Days"):

        for _, sta_row in stations.iterrows():
            net     = sta_row["network"]
            sta     = sta_row["station"]

            # Distance filter: skip stations far from all events on this day
            min_dist = min(
                haversine_km(sta_row["latitude"], sta_row["longitude"],
                             ev["latitude"], ev["longitude"])
                for _, ev in day_events.iterrows()
            )
            if min_dist > config.MAX_STATION_DIST_KM:
                continue

            # Find files for this station/day
            files = find_daily_files(
                sta_row["data_root"], net, sta,
                year, jday, sta_row["corrected"]
            )
            if not files:
                skipped_stations += 1
                continue

            # Load and merge all channels for the day
            try:
                stream = obspy.Stream()
                for f in files:
                    stream += read(f, dtype="float64")
                stream.merge(fill_value="interpolate")
            except Exception:
                skipped_stations += 1
                continue

            # Per-event windowing + picking
            picks = pick_windows(model, stream, day_events)
            all_picks.extend(picks)

    print(f"\nTotal picks: {len(all_picks)}  "
          f"(station-days with no data: {skipped_stations})")

    df_picks = pd.DataFrame(all_picks)
    if df_picks.empty:
        print("No picks produced — check data paths and thresholds in config.py")
        sys.exit(1)

    df_picks.to_csv(config.PICKS_REPICKED, index=False,
                    date_format="%Y-%m-%dT%H:%M:%S.%f")
    print(f"Picks written to {config.PICKS_REPICKED}")


if __name__ == "__main__":
    main()
