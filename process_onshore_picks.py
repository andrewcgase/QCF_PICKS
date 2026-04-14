#!/usr/bin/env python3
"""
Process Anne's onshore/offshore .tbl pick files into master CSVs.

For each .tbl file in the onshore_tables/ directory the script:
  1. Parses the Reveal JSON pick table (FFID, Time pairs)
  2. Splits picks into segments wherever consecutive FFIDs are separated by
     more than GAP_THRESHOLD shots (avoids interpolating through shadow zones
     or phase boundaries)
  3. Within each segment: linearly interpolates over missing FFIDs, then
     applies a boxcar (moving-average) smoothing filter
  4. Assigns pick uncertainties (in seconds) by FFID range from
     onshore-offshore uncertainties.csv — rows outside any defined range
     are dropped
  5. Merges FFID source positions from the cruise shotlog
  6. Attaches station receiver metadata (lat, lon, elevation) from
     onshore_station_locations.csv
  7. Writes one CSV per station×line to onshore_tables/processed/

File naming convention
----------------------
  Each .tbl file must be named {STATION}_{LINE}.tbl, where STATION is the
  seismic station code (e.g. BNAB) and LINE is the seismic line name
  (e.g. PM01B2).  The two parts are split on the first underscore only.

  Examples:
    BNAB_PM01B2.tbl → station=BNAB, line=PM01B2
    SIT_T12A14B.tbl → station=SIT, line=T12A14B

Uncertainty ranks (onshore)
---------------------------
  rank 1 → 0.100 s   (good)
  rank 2 → 0.250 s   (okay)
  rank 3 → 0.500 s   (bad, but picked)

Required input files
--------------------
  onshore_tables/*.tbl              Pick files from Reveal (Anne provides)
  onshore-offshore uncertainties.csv  Uncertainty table (Anne provides)
  onshore_station_locations.csv     Station lat/lon/elevation (Anne provides)
  allshots_formatted.csv            Unified cruise shot table

Usage
-----
  python process_onshore_picks.py                      # all files, gap threshold = 20 FFIDs
  python process_onshore_picks.py --gap-threshold 50   # stricter segmenting
  python process_onshore_picks.py --boxcar-window 5    # wider smoothing
  python process_onshore_picks.py --station BNAB_PM01B2  # single file

Output columns
--------------
  FFID, Time (Seconds), Time Boxcar (Seconds), Uncertainty (Seconds),
  FFID Latitude, FFID Longitude, FFID Depth (km),
  Station, Latitude, Longitude, Elevation (km)
"""

import argparse
import json
import sys
from math import floor, log10
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Paths (relative to this script's location)
# ---------------------------------------------------------------------------
BASE          = Path(__file__).parent
ONSHORE_DIR   = BASE / "onshore_tables"
OUTPUT_DIR    = BASE / "onshore_tables" / "processed"
UNC_FILE      = BASE / "onshore-offshore uncertainties.csv"
STATIONS_FILE = BASE / "onshore_station_locations.csv"
ALLSHOTS      = BASE / "allshots_formatted.csv"

# Uncertainty rank → seconds (onshore convention)
RANK_TO_SEC = {1: 0.100, 2: 0.250, 3: 0.500}

FFID_DEPTH_KM = 0.012   # fixed source depth per project convention


# ---------------------------------------------------------------------------
# Helpers (identical to process_anne_picks.py)
# ---------------------------------------------------------------------------

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
    Load onshore-offshore uncertainties.csv and return a clean DataFrame with
    columns: station, line, ffid_start, ffid_end, rank.

    The file has trailing spaces in column headers and extra blank columns, so
    we select by position (cols 0-4) and rename manually.
    """
    df = pd.read_csv(path, header=0, usecols=[0, 1, 2, 3, 4])
    df.columns = ["station", "line", "ffid_start", "ffid_end", "rank"]
    df = df.dropna(subset=["station", "line", "ffid_start", "ffid_end", "rank"])
    df["station"]    = df["station"].astype(str).str.strip()
    df["line"]       = df["line"].astype(str).str.strip()
    df["ffid_start"] = pd.to_numeric(df["ffid_start"], errors="coerce").astype("Int64")
    df["ffid_end"]   = pd.to_numeric(df["ffid_end"],   errors="coerce").astype("Int64")
    df["rank"]       = pd.to_numeric(df["rank"],        errors="coerce").astype("Int64")
    df = df.dropna()

    # Warn about inverted ranges (data entry errors)
    bad = df[df["ffid_start"] > df["ffid_end"]]
    if not bad.empty:
        print("WARNING: inverted FFID ranges in uncertainty table (ffid_start > ffid_end):")
        for _, row in bad.iterrows():
            print(f"  {row['station']}, {row['line']}: "
                  f"{row['ffid_start']} > {row['ffid_end']}")
        print("  These rows will be skipped during uncertainty lookup.")

    return df


def load_station_locations(path):
    """
    Load onshore_station_locations.csv.
    Returns a DataFrame indexed by Station with Latitude, Longitude, Elevation (km).
    """
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df["Station"] = df["Station"].astype(str).str.strip()
    missing = df[df[["Latitude", "Longitude", "Elevation (km)"]].isna().any(axis=1)]
    if not missing.empty:
        print("WARNING: missing coordinates in onshore_station_locations.csv:")
        for _, row in missing.iterrows():
            print(f"  {row['Station']}")
    return df.set_index("Station")


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


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------

def process_file(tbl_path, unc_df, shotlog, stations_df, gap_threshold, boxcar_window):
    stem    = Path(tbl_path).stem                   # e.g. "BNAB_PM01B2"
    station, line = stem.split("_", 1)              # e.g. "BNAB", "PM01B2"

    # --- Station metadata ---
    if station not in stations_df.index:
        print(f"  SKIP {stem}: '{station}' not found in onshore_station_locations.csv")
        return None
    sta_row = stations_df.loc[station]
    if pd.isna(sta_row["Latitude"]) or pd.isna(sta_row["Longitude"]):
        print(f"  SKIP {stem}: coordinates for '{station}' are missing — fill in "
              f"onshore_station_locations.csv first")
        return None
    sta_lat  = float(sta_row["Latitude"])
    sta_lon  = float(sta_row["Longitude"])
    sta_elev = float(sta_row["Elevation (km)"])

    # --- Uncertainty ranges for this station + line ---
    unc_rows = unc_df[
        (unc_df["station"] == station) & (unc_df["line"] == line)
    ]
    if unc_rows.empty:
        print(f"  SKIP {stem}: no uncertainty data for station={station}, line={line}")
        return None

    # Exclude inverted ranges
    valid_unc = unc_rows[unc_rows["ffid_start"] <= unc_rows["ffid_end"]]
    unc_ranges = list(zip(
        valid_unc["ffid_start"], valid_unc["ffid_end"], valid_unc["rank"]
    ))
    if not unc_ranges:
        print(f"  SKIP {stem}: all uncertainty ranges are inverted for "
              f"station={station}, line={line}")
        return None

    # --- Parse picks ---
    picks = parse_tbl(tbl_path)
    if not picks:
        print(f"  SKIP {stem}: empty pick table")
        return None

    # --- Split and interpolate ---
    segments = split_segments(picks, gap_threshold)
    seg_dfs  = [interpolate_segment(seg, boxcar_window) for seg in segments]
    df = pd.concat(seg_dfs, ignore_index=True)

    # --- Assign uncertainties; drop rows outside defined ranges ---
    def lookup_unc(ffid):
        for start, end, rank in unc_ranges:
            if start <= ffid <= end:
                return RANK_TO_SEC.get(int(rank))
        return np.nan

    df["Uncertainty"] = df["FFID"].apply(lookup_unc)
    n_dropped = df["Uncertainty"].isna().sum()
    df = df.dropna(subset=["Uncertainty"])
    if n_dropped:
        print(f"  {stem}: dropped {n_dropped} interpolated rows outside "
              f"uncertainty ranges")

    if df.empty:
        print(f"  SKIP {stem}: no rows remain after uncertainty filtering")
        return None

    # --- Merge FFID positions from unified shot table ---
    df = df.join(shotlog[["sourceLat", "sourceLon"]], on="FFID", how="left")

    n_missing = df["sourceLat"].isna().sum()
    if n_missing:
        print(f"  WARNING {stem}: {n_missing} FFIDs not found in shotlog")

    # --- Build final output ---
    final = pd.DataFrame({
        "FFID":                   df["FFID"],
        "Time (Seconds)":         df["Time"],
        "Time Boxcar (Seconds)":  df["Time Boxcar"],
        "Uncertainty (Seconds)":  df["Uncertainty"],
        "FFID Latitude":          df["sourceLat"],
        "FFID Longitude":         df["sourceLon"],
        "FFID Depth (km)":        FFID_DEPTH_KM,
        "Station":                station,
        "Latitude":               sta_lat,
        "Longitude":              sta_lon,
        "Elevation (km)":         sta_elev,
    })

    n_segs   = len(segments)
    seg_info = ", ".join(
        f"{len(s)} picks → {len(d)} rows"
        for s, d in zip(segments, seg_dfs)
    )
    print(f"  {stem}: {n_segs} segment(s) [{seg_info}] → {len(final)} output rows")

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
             "Gaps larger than this start a new segment."
    )
    parser.add_argument(
        "--boxcar-window", type=int, default=3, metavar="W",
        help="Rolling-mean window size for boxcar smoothing (default: 3)."
    )
    parser.add_argument(
        "--station", type=str, default=None, metavar="NAME",
        help="Process a single file only, e.g. BNAB_PM01B2 (for testing)."
    )
    args = parser.parse_args()

    print(f"Settings: gap_threshold={args.gap_threshold} FFIDs, "
          f"boxcar_window={args.boxcar_window}")
    print()

    # Load shared resources
    unc_df      = load_uncertainty_table(UNC_FILE)
    stations_df = load_station_locations(STATIONS_FILE)
    shotlog     = load_shotlog()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Find .tbl files to process
    if args.station:
        tbl_files = [ONSHORE_DIR / f"{args.station}.tbl"]
        missing = [f for f in tbl_files if not f.exists()]
        if missing:
            sys.exit(f"File not found: {missing[0]}")
    else:
        tbl_files = sorted(ONSHORE_DIR.glob("*.tbl"))

    if not tbl_files:
        print(f"No .tbl files found in {ONSHORE_DIR}")
        print("Place Anne's pick files there, named {{STATION}}_{{LINE}}.tbl")
        print("  e.g.  BNAB_PM01B2.tbl,  SIT_T12A14B.tbl")
        return

    print(f"Processing {len(tbl_files)} file(s)...\n")

    success, skipped = 0, 0
    for tbl_path in tbl_files:
        result = process_file(
            tbl_path, unc_df, shotlog, stations_df,
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
