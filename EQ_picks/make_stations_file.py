#!/usr/bin/env python3
"""
make_stations_file.py -- Write a Pyrocko stations file covering every OBS +
onshore station, so Snuffler can compute per-trace distance (needed for its
"sort by distance" menu option, which silently does nothing without this --
Snuffler has no station coordinates loaded at all currently, since neither
01_extract_windows.py/01_extract_toques.py nor 02_picks_to_markers.py/
02_nor_to_markers.py ever write or load a --stations file).

Station locations are static (don't change per event/month), so this only
needs to be (re)run once, or whenever the station lists in
EarthScope_OBS_Deployments.csv / onshore.xml change.

Usage:
    python3 make_stations_file.py
    snuffler waveforms_toques/2021-09/ --markers=markers/2021-09_toques.markers \\
        --stations=stations.pf
"""
import csv
from pathlib import Path

from pyrocko import model

import config

OUT_PATH = Path(__file__).parent / "stations.pf"


def main():
    stations = []

    with open(config.OBS_META_CSV) as f:
        for row in csv.DictReader(f):
            stations.append(model.Station(
                network=config.OBS_NETWORK, station=row["Station"], location="",
                lat=float(row["Latitude"]), lon=float(row["Longitude"]),
                elevation=float(row["Depth (m)"]),  # already negative-down == positive-up convention
            ))

    if config.ONSHORE_XML.exists():
        import obspy
        inv = obspy.read_inventory(str(config.ONSHORE_XML))
        for net in inv:
            for sta in net:
                stations.append(model.Station(
                    network=net.code, station=sta.code, location="",
                    lat=sta.latitude, lon=sta.longitude,
                    elevation=sta.elevation if sta.elevation is not None else 0.0,
                ))
    else:
        print(f"WARNING: {config.ONSHORE_XML} not found -- onshore stations omitted")

    model.dump_stations(stations, str(OUT_PATH))
    print(f"Wrote {len(stations)} stations ({sum(1 for s in stations if s.network == config.OBS_NETWORK)} "
          f"OBS, {sum(1 for s in stations if s.network != config.OBS_NETWORK)} onshore) to {OUT_PATH}")
    print(f"\nPass this to Snuffler: --stations={OUT_PATH}")


if __name__ == "__main__":
    main()
