#!/usr/bin/env python3
"""
Process Anne's revised .tbl pick files into master CSVs.

For each .tbl file in the tables/ directory the script:
  1. Parses the Reveal JSON pick table (FFID, Time pairs)
  2. Splits picks into segments wherever consecutive FFIDs are separated by
     more than GAP_THRESHOLD shots (avoids interpolating through water-wave
     shadow zones or phase boundaries)
  3. Within each segment: linearly interpolates over missing FFIDs, then
     applies a boxcar (moving-average) smoothing filter
  4. Undoes the 6 km/s reduction velocity applied during picking:
     T_actual = T_picked + offset / 6, where offset is the great-circle
     source–receiver distance in km
  5. Assigns pick uncertainties (in seconds) by FFID range from
     OBS uncertainties.csv — rows outside any defined range are dropped
  6. Merges FFID source positions from the cruise shotlog
  7. Attaches OBS receiver metadata (lat, lon, depth) from EarthScope
  8. Writes one CSV per station to tables/processed/

Usage
-----
  python process_anne_picks.py                      # default gap threshold = 20 FFIDs
  python process_anne_picks.py --gap-threshold 50   # stricter
  python process_anne_picks.py --boxcar-window 5    # wider smoothing
  python process_anne_picks.py --station PO01A_12S01  # single file for testing

Output columns (same schema as Josh's master files)
----------------------------------------------------
  FFID, Time (Seconds), Time Boxcar (Seconds), Uncertainty (Seconds),
  FFID Latitude, FFID Longitude, FFID Depth (km),
  OBS Number, OBS Line, OBS Depth (km), OBS Latitude, OBS Longitude
"""

import argparse
import json
import glob
import os
import sys
from math import floor, log10, radians, sin, cos, sqrt, atan2
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Paths (relative to this script's location)
# ---------------------------------------------------------------------------
BASE = Path(__file__).parent
TABLES_DIR  = BASE / "tables"
OUTPUT_DIR  = BASE / "tables" / "processed"
UNC_FILE    = BASE / "OBS uncertainties.csv"
ES_FILE     = BASE / "EarthScope_OBS_Deployments.csv"
ALLSHOTS    = BASE / "allshots_formatted.csv"

# Uncertainty rank → seconds
RANK_TO_SEC = {1: 0.024, 2: 0.050, 3: 0.100}

FFID_DEPTH_KM = 0.012       # fixed source depth per project convention
REDUCTION_VELOCITY = 6.0    # km/s — reduction velocity used during picking


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two (lat, lon) points."""
    R = 6371.0
    phi1, phi2 = radians(lat1), radians(lat2)
    dphi = radians(lat2 - lat1)
    dlam = radians(lon2 - lon1)
    a = sin(dphi / 2) ** 2 + cos(phi1) * cos(phi2) * sin(dlam / 2) ** 2
    return R * 2 * atan2(sqrt(a), sqrt(1 - a))


def round_sig(x, sig=7):
    """Round to sig significant figures; pass through zero."""
    if not isinstance(x, (int, float)) or x == 0:
        return x
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def parse_tbl(path):
    """Return sorted list of [FFID, Time] pairs from a Reveal .tbl file."""
    with open(path) as f:
        data = json.load(f)
    picks = data.get("f_of_x_picks", [])
    return sorted(picks, key=lambda p: p[0])


def split_segments(picks, gap_threshold):
    """
    Partition picks into contiguous segments.
    A new segment begins whenever two consecutive FFIDs differ by more than
    gap_threshold.  Each segment is returned as a list of [FFID, Time] pairs.
    """
    if not picks:
        return []
    segments, current = [], [picks[0]]
    for i in range(1, len(picks)):
        if picks[i][0] - picks[i - 1][0] > gap_threshold:
            segments.append(current)
            current = []
        current.append(picks[i])
    segments.append(current)
    return segments


def interpolate_segment(picks, boxcar_window=3):
    """
    Linearly interpolate Time over every integer FFID within the segment's
    range, then apply a rolling-mean (boxcar) filter.
    Returns a DataFrame with columns [FFID, Time, Time Boxcar].
    """
    df = pd.DataFrame(picks, columns=["FFID", "Time"])
    df["FFID"] = pd.to_numeric(df["FFID"])
    df["Time"] = pd.to_numeric(df["Time"]) / 1000.0   # Reveal picks are in ms → convert to s
    df = df.set_index("FFID")

    full_range = range(int(df.index.min()), int(df.index.max()) + 1)
    df = df.reindex(full_range).interpolate()

    df["Time"]        = df["Time"].apply(lambda x: round_sig(x, 7))
    df["Time Boxcar"] = (
        df["Time"]
        .rolling(window=boxcar_window, min_periods=1)
        .mean()
        .apply(lambda x: round_sig(x, 7))
    )
    return df.reset_index().rename(columns={"index": "FFID"})


def load_uncertainty_table(path):
    """
    Load OBS uncertainties.csv and return a clean DataFrame with columns:
      cruise, station, ffid_start, ffid_end, rank

    The file has inconsistent column naming (trailing spaces, duplicate names)
    so we select by position: cols 0-4 are cruise, station, ffid_start,
    ffid_end, rank.
    """
    df = pd.read_csv(path, header=0, usecols=[0, 1, 2, 3, 4])
    df.columns = ["cruise", "station", "ffid_start", "ffid_end", "rank"]
    df = df.dropna(subset=["cruise", "station", "ffid_start", "ffid_end", "rank"])
    df["cruise"]     = df["cruise"].astype(str).str.strip()
    df["station"]    = df["station"].astype(str).str.strip()
    df["ffid_start"] = df["ffid_start"].astype(int)
    df["ffid_end"]   = df["ffid_end"].astype(int)
    df["rank"]       = df["rank"].astype(int)
    return df


def load_shotlog():
    """
    Load allshots_formatted.csv — the unified shot table for the entire survey.
    Returns a DataFrame indexed by shot number with sourceLat, sourceLon columns.
    """
    df = pd.read_csv(ALLSHOTS, low_memory=False,
                     usecols=["shot", "sourceLat", "sourceLon"])
    df = df.dropna(subset=["shot", "sourceLat", "sourceLon"])
    df["shot"] = df["shot"].astype(int)
    return df.set_index("shot")


def get_obs_line(station):
    """Extract OBS line from station name, e.g. '12S01' -> '12S', '1AS03' -> '1AS'."""
    # Station number is always the last 2 characters
    return station[:-2]


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------

def process_file(tbl_path, unc_df, shotlog, es_df, gap_threshold, boxcar_window):
    name   = Path(tbl_path).stem                     # e.g. "PO01A_12S01"
    cruise = name.split("_")[0]                      # e.g. "PO01A"
    station = "_".join(name.split("_")[1:])          # e.g. "12S01"
    obs_line = get_obs_line(station)                 # e.g. "12S"

    # --- EarthScope OBS metadata ---
    obs_row = es_df[es_df["Station"] == station]
    if obs_row.empty:
        print(f"  SKIP {name}: no EarthScope record for '{station}'")
        return None
    obs_lat   = obs_row["Latitude"].iloc[0]
    obs_lon   = obs_row["Longitude"].iloc[0]
    obs_depth = obs_row["Depth (km)"].iloc[0]

    # --- Uncertainty ranges for this cruise + station ---
    unc_rows = unc_df[
        (unc_df["cruise"] == cruise) & (unc_df["station"] == station)
    ]
    if unc_rows.empty:
        print(f"  SKIP {name}: no uncertainty data for {cruise} / {station}")
        return None
    unc_ranges = list(zip(
        unc_rows["ffid_start"], unc_rows["ffid_end"], unc_rows["rank"]
    ))

    # --- Parse picks ---
    picks = parse_tbl(tbl_path)
    if not picks:
        print(f"  SKIP {name}: empty pick table")
        return None

    # --- Split and interpolate ---
    segments = split_segments(picks, gap_threshold)
    seg_dfs = [interpolate_segment(seg, boxcar_window) for seg in segments]
    df = pd.concat(seg_dfs, ignore_index=True)

    # --- Assign uncertainties; drop rows outside defined ranges ---
    def lookup_unc(ffid):
        for start, end, rank in unc_ranges:
            if start <= ffid <= end:
                return RANK_TO_SEC.get(rank)
        return np.nan

    df["Uncertainty"] = df["FFID"].apply(lookup_unc)
    n_dropped = df["Uncertainty"].isna().sum()
    df = df.dropna(subset=["Uncertainty"])
    if n_dropped:
        print(f"  {name}: dropped {n_dropped} interpolated rows outside uncertainty ranges")

    # --- Merge FFID positions from unified shot table ---
    df = df.join(shotlog[["sourceLat", "sourceLon"]], on="FFID", how="left")

    n_missing = df["sourceLat"].isna().sum()
    if n_missing:
        print(f"  WARNING {name}: {n_missing} FFIDs not found in shotlog")

    # --- Undo reduction velocity (picks were made at V_red = 6 km/s) ---
    # T_picked = T_actual - offset / V_red  →  T_actual = T_picked + offset / V_red
    def _offset(row):
        if pd.isna(row["sourceLat"]):
            return np.nan
        return haversine_km(row["sourceLat"], row["sourceLon"], obs_lat, obs_lon)

    df["offset_km"] = df.apply(_offset, axis=1)
    correction = df["offset_km"] / REDUCTION_VELOCITY
    df["Time"]        = (df["Time"]        + correction).apply(lambda x: round_sig(x, 7))
    df["Time Boxcar"] = (df["Time Boxcar"] + correction).apply(lambda x: round_sig(x, 7))

    # --- Build final output ---
    final = pd.DataFrame({
        "FFID":                   df["FFID"],
        "Time (Seconds)":         df["Time"],
        "Time Boxcar (Seconds)":  df["Time Boxcar"],
        "Uncertainty (Seconds)":  df["Uncertainty"],
        "FFID Latitude":          df["sourceLat"],
        "FFID Longitude":         df["sourceLon"],
        "FFID Depth (km)":        FFID_DEPTH_KM,
        "OBS Number":             station,
        "OBS Line":               obs_line,
        "OBS Depth (km)":         obs_depth,
        "OBS Latitude":           obs_lat,
        "OBS Longitude":          obs_lon,
    })

    n_segs = len(segments)
    seg_info = ", ".join(
        f"{len(s)} picks → {len(d)} rows"
        for s, d in zip(segments, seg_dfs)
    )
    print(f"  {name}: {n_segs} segment(s) [{seg_info}] → {len(final)} output rows")

    return final


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--gap-threshold", type=int, default=20, metavar="N",
        help="Max consecutive-FFID gap to interpolate across (default: 20). "
             "Gaps larger than this mark a segment boundary."
    )
    parser.add_argument(
        "--boxcar-window", type=int, default=3, metavar="W",
        help="Rolling-mean window size for boxcar smoothing (default: 3)."
    )
    parser.add_argument(
        "--station", type=str, default=None, metavar="NAME",
        help="Process a single station only, e.g. PO01A_12S01 (for testing)."
    )
    args = parser.parse_args()

    print(f"Settings: gap_threshold={args.gap_threshold} FFIDs, "
          f"boxcar_window={args.boxcar_window}")
    print()

    # Load shared resources
    unc_df  = load_uncertainty_table(UNC_FILE)
    es_df   = pd.read_csv(ES_FILE)
    shotlog = load_shotlog()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Find .tbl files to process
    if args.station:
        tbl_files = [TABLES_DIR / f"{args.station}.tbl"]
        missing = [f for f in tbl_files if not f.exists()]
        if missing:
            sys.exit(f"File not found: {missing[0]}")
    else:
        tbl_files = sorted(
            f for f in TABLES_DIR.glob("*.tbl")
            if f.stem not in ("test", "test1")
        )

    print(f"Processing {len(tbl_files)} file(s)...\n")

    success, skipped = 0, 0
    for tbl_path in tbl_files:
        result = process_file(
            tbl_path, unc_df, shotlog, es_df,
            args.gap_threshold, args.boxcar_window
        )
        if result is not None:
            out_path = OUTPUT_DIR / f"{Path(tbl_path).stem}.csv"
            result.to_csv(out_path, index=False)
            success += 1
        else:
            skipped += 1

    print(f"\nDone: {success} files written to {OUTPUT_DIR}, {skipped} skipped.")


if __name__ == "__main__":
    main()
