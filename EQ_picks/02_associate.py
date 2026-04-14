#!/usr/bin/env python3
"""
02_associate.py — Re-associate picks from 01_repick.py using GaMMA.

Reads:  config.PICKS_REPICKED  (picks_repicked.csv)
Writes: config.EVENTS_ASSOC    (events_associated.csv)
        config.PICKS_ASSOC     (picks_associated.csv)

The GaMMA configuration mirrors the one used in the original detection run
(seisbench/apply/apply_SBM.ipynb): BGMM method, eikonal 1-D velocity model,
Vp/Vs = 1.73.

Usage
-----
  conda activate obspy        # environment with obspy + gamma
  python 02_associate.py
"""

import sys
import numpy as np
import pandas as pd
from pyproj import Proj
from obspy import read_inventory

from gamma.utils import association, estimate_eps

import config


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_stations():
    """
    Build the GaMMA-format stations DataFrame:
    columns: id, longitude, latitude, elevation(m), x(km), y(km), z(km)
    """
    rows = []

    # OBS
    obs_meta = pd.read_csv(config.OBS_META_CSV)
    for _, row in obs_meta.iterrows():
        sta = row["Station"]
        if sta in config.SKIP_STATIONS:
            continue
        rows.append({
            "id":           f"{config.OBS_NETWORK}.{sta}..",
            "longitude":    row["Longitude"],
            "latitude":     row["Latitude"],
            "elevation(m)": float(row["Depth (m)"]),
        })

    # Onshore
    if config.ONSHORE_XML.exists():
        inv = read_inventory(str(config.ONSHORE_XML))
        for network in inv:
            for station in network:
                rows.append({
                    "id":           f"{network.code}.{station.code}..",
                    "longitude":    station.longitude,
                    "latitude":     station.latitude,
                    "elevation(m)": station.elevation,
                })

    stations = pd.DataFrame(rows)

    # Stereographic projection centred on the network
    x0 = stations["longitude"].median()
    y0 = stations["latitude"].median()
    proj = Proj(f"+proj=sterea +lon_0={x0} +lat_0={y0} +units=km")

    stations[["x(km)", "y(km)"]] = stations.apply(
        lambda r: pd.Series(proj(longitude=r["longitude"], latitude=r["latitude"])),
        axis=1,
    )
    stations["z(km)"] = stations["elevation(m)"].apply(lambda e: -e / 1e3)

    return stations, proj, x0, y0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    config.CATALOGS_DIR.mkdir(parents=True, exist_ok=True)

    if not config.PICKS_REPICKED.exists():
        sys.exit(f"Run 01_repick.py first — {config.PICKS_REPICKED} not found")

    picks_df = pd.read_csv(config.PICKS_REPICKED)
    print(f"Loaded {len(picks_df)} picks from {config.PICKS_REPICKED.name}")

    # Keep only associated-quality picks
    picks_df = picks_df.rename(columns={
        "station_id": "id",
        "phase_time":  "timestamp",
        "phase_score": "prob",
        "phase_type":  "type",
    })
    picks_df["timestamp"] = pd.to_datetime(picks_df["timestamp"])
    picks_df["type"] = picks_df["type"].str.upper()

    stations, proj, x0, y0 = load_stations()
    print(f"Stations: {len(stations)}")

    # Build GaMMA config
    cfg = dict(config.GAMMA_CONFIG)
    cfg["center"] = (x0, y0)

    x_lims = (stations["x(km)"].min() - 50, stations["x(km)"].max() + 50)
    y_lims = (stations["y(km)"].min() - 50, stations["y(km)"].max() + 50)
    cfg["x(km)"] = x_lims
    cfg["y(km)"] = y_lims
    cfg["xlim_degree"] = (
        proj(longitude=x_lims[0], latitude=y0, inverse=True)[0],
        proj(longitude=x_lims[1], latitude=y0, inverse=True)[0],
    )
    cfg["ylim_degree"] = (
        proj(longitude=x0, latitude=y_lims[0], inverse=True)[1],
        proj(longitude=x0, latitude=y_lims[1], inverse=True)[1],
    )
    cfg["bfgs_bounds"] = (
        (x_lims[0] - 1, x_lims[1] + 1),
        (y_lims[0] - 1, y_lims[1] + 1),
        (0, cfg["z(km)"][1] + 1),
        (None, None),
    )
    cfg["dbscan_eps"] = estimate_eps(stations, cfg["vel"]["p"])

    # Eikonal travel-time table
    v = config.EIKONAL_VEL
    cfg["eikonal"] = {
        "vel": v,
        "h": config.EIKONAL_DZ,
        "xlim": x_lims,
        "ylim": y_lims,
        "zlim": cfg["z(km)"],
    }

    print("Running GaMMA association...")
    print(f"  method={cfg['method']}, min_picks={cfg['min_picks_per_eq']}, "
          f"dbscan_eps={cfg['dbscan_eps']:.1f} km")

    events_list, assignments = association(picks_df, stations, cfg, 0, cfg["method"])
    print(f"  → {len(events_list)} events associated")

    # Build output DataFrames
    events = pd.DataFrame(events_list)
    if events.empty:
        print("No events associated — try relaxing config.GAMMA_CONFIG thresholds")
        sys.exit(1)

    events[["longitude", "latitude"]] = events.apply(
        lambda r: pd.Series(proj(longitude=r["x(km)"], latitude=r["y(km)"], inverse=True)),
        axis=1,
    )
    events["depth_km"] = events["z(km)"]
    events.to_csv(config.EVENTS_ASSOC, index=False,
                  float_format="%.3f", date_format="%Y-%m-%dT%H:%M:%S.%f")
    print(f"Events written to {config.EVENTS_ASSOC}")

    assignments_df = pd.DataFrame(assignments,
                                   columns=["pick_index", "event_index", "gamma_score"])
    picks_out = picks_df.join(
        assignments_df.set_index("pick_index")
    ).fillna(-1).astype({"event_index": int})
    picks_out = picks_out.rename(columns={
        "id": "station_id", "timestamp": "phase_time",
        "prob": "phase_score", "type": "phase_type",
    })
    picks_out.to_csv(config.PICKS_ASSOC, index=False,
                     date_format="%Y-%m-%dT%H:%M:%S.%f")
    print(f"Picks written to {config.PICKS_ASSOC}")

    # Summary
    assoc = picks_out[picks_out["event_index"] >= 0]
    print(f"\nAssociated picks: {len(assoc)} / {len(picks_out)}  "
          f"({100*len(assoc)/len(picks_out):.1f}%)")
    print(f"P picks: {(assoc['phase_type']=='P').sum()},  "
          f"S picks: {(assoc['phase_type']=='S').sum()}")


if __name__ == "__main__":
    main()
