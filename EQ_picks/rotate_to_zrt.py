#!/usr/bin/env python3
"""
rotate_to_zrt.py -- Rotate extracted TOQUES waveform windows (Z + two
horizontals) into Z/R/T using each event's preliminary Nordic-catalog
location, for picking in Snuffler.

Two different rotations are needed depending on station type:

  Onshore (BHN/BHE, already true north/east):
    N/E -> R/T directly via ObsPy's rotate_ne_rt(back_azimuth) -- standard,
    unambiguous, no station-specific correction needed.

  OBS (HH1/HH2, NOT aligned to true north -- these free-fall to the seafloor
  and settle at an arbitrary heading):
    H1/H2 -> N/E first, using each station's real measured orientation angle
    (see data/obs_orientation/README.md for provenance: OBSIC/WHOI, Doran &
    Laske 2017 Rayleigh-wave method), THEN N/E -> R/T as above.
    H1/H2 -> N/E uses the convention confirmed against DLOPy's own source
    (data/obs_orientation/README.md): theta = clockwise angle from true
    North to H1 (the "Orientation, degrees" column), and
        N = H1*cos(theta) - H2*sin(theta)
        E = H1*sin(theta) + H2*cos(theta)
    OBS stations with no measured orientation, or whose reported
    uncertainty exceeds --max-orientation-error (default 30 deg -- e.g.
    QCB10 at 345+/-174 deg is essentially unconstrained), are left
    unrotated (Z/1/2 passed through as-is) rather than rotated with a
    near-meaningless angle -- printed as a warning, not silently dropped.

Z is passed through unchanged for every station. Any other channel present
(e.g. HDH/EDH hydrophone) is passed through unchanged too.

Output channel naming: the component letter is replaced the usual way
(...N->...R, ...E->...T for onshore; ...1->...R, ...2->...T for OBS),
keeping the same band+instrument code, e.g. BHN->BHR, HH1->HHR.

Usage:
    python3 rotate_to_zrt.py                      # all months, all events
    python3 rotate_to_zrt.py --month 2021-09
    python3 rotate_to_zrt.py --max-orientation-error 20
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt

import config
import nordic

OBS_ORIENTATION_CSV = Path(__file__).parent / "data" / "obs_orientation" / "YI_experiment_metrics.csv"
WAVEFORMS_ZRT_DIR = config.WAVEFORMS_TOQUES_DIR.parent / "waveforms_toques_zrt"


def load_obs_orientations():
    """station -> (theta_deg, error_deg). Only stations with a real
    (non '--') measured orientation are included."""
    import csv
    out = {}
    with open(OBS_ORIENTATION_CSV) as f:
        for row in csv.DictReader(f):
            theta, err = row["Orientation, degrees"], row["Error, degrees"]
            if theta in ("--", "") or err in ("--", ""):
                continue
            out[row["Station"]] = (float(theta), float(err))
    return out


def load_obs_locations():
    """OBS station -> (lat, lon), from EarthScope_OBS_Deployments.csv."""
    import csv
    out = {}
    with open(config.OBS_META_CSV) as f:
        for row in csv.DictReader(f):
            out[row["Station"]] = (float(row["Latitude"]), float(row["Longitude"]))
    return out


def load_onshore_locations():
    """onshore station -> (lat, lon), from the station XML inventory."""
    out = {}
    if not config.ONSHORE_XML.exists():
        return out
    inv = obspy.read_inventory(str(config.ONSHORE_XML))
    for net in inv:
        for sta in net:
            out[sta.code] = (sta.latitude, sta.longitude)
    return out


def h1h2_to_ne(h1, h2, theta_deg):
    """DLOPy convention -- see module docstring / data/obs_orientation/README.md."""
    th = np.radians(theta_deg)
    n = h1 * np.cos(th) - h2 * np.sin(th)
    e = h1 * np.sin(th) + h2 * np.cos(th)
    return n, e


def rotate_stream(st, ev_lat, ev_lon, sta_lat, sta_lon, obs_orientation, max_orient_err):
    """Rotate one event's stream (all channels for one station) to Z/R/T in
    place where possible. Returns (rotated_stream, status) where status is
    'onshore', 'obs', or 'skipped:<reason>'."""
    _, _, back_azimuth = gps2dist_azimuth(ev_lat, ev_lon, sta_lat, sta_lon)

    by_comp = {tr.stats.channel[-1]: tr for tr in st}
    others = [tr for tr in st if tr.stats.channel[-1] not in ("1", "2", "N", "E")]

    if "N" in by_comp and "E" in by_comp:
        n_tr, e_tr = by_comp["N"], by_comp["E"]
        r_data, t_data = rotate_ne_rt(n_tr.data.astype(float), e_tr.data.astype(float), back_azimuth)
        status = "onshore"
    elif "1" in by_comp and "2" in by_comp:
        h1_tr, h2_tr = by_comp["1"], by_comp["2"]
        station = h1_tr.stats.station
        if station not in obs_orientation:
            return st, "skipped:no orientation measured"
        theta, err = obs_orientation[station]
        if err > max_orient_err:
            return st, f"skipped:orientation error {err:g} deg > max {max_orient_err:g} deg"
        n_data, e_data = h1h2_to_ne(h1_tr.data.astype(float), h2_tr.data.astype(float), theta)
        r_data, t_data = rotate_ne_rt(n_data, e_data, back_azimuth)
        n_tr, e_tr = h1_tr, h2_tr  # reuse for stats template below
        status = "obs"
    else:
        return st, "skipped:no horizontal pair found"

    r_tr, t_tr = n_tr.copy(), e_tr.copy()
    r_tr.data, t_tr.data = r_data.astype(n_tr.data.dtype), t_data.astype(e_tr.data.dtype)
    r_tr.stats.channel = n_tr.stats.channel[:-1] + "R"
    t_tr.stats.channel = e_tr.stats.channel[:-1] + "T"

    out = obspy.Stream(others)
    if "Z" in by_comp:
        out += by_comp["Z"]
    out += r_tr
    out += t_tr
    return out, status


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--month", default=None, help="only process this YYYY-MM subfolder")
    p.add_argument("--max-orientation-error", type=float, default=30.0,
                   help="OBS stations with reported orientation error (deg) above this "
                        "are left unrotated (default: 30)")
    args = p.parse_args()

    events = nordic.read_nor_file(str(config.TOQUES_NOR))
    obs_orientation = load_obs_orientations()
    obs_loc = load_obs_locations()
    onshore_loc = load_onshore_locations()
    print(f"{len(obs_orientation)} OBS stations have a usable measured orientation "
          f"(error <= {args.max_orientation_error:g} deg checked per-station below)")

    month_dirs = ([config.WAVEFORMS_TOQUES_DIR / args.month] if args.month
                  else sorted(p for p in config.WAVEFORMS_TOQUES_DIR.iterdir() if p.is_dir()))

    counts = {}
    for month_dir in month_dirs:
        if not month_dir.exists():
            sys.exit(f"{month_dir} not found")
        out_dir = WAVEFORMS_ZRT_DIR / month_dir.name
        out_dir.mkdir(parents=True, exist_ok=True)
        for mseed_path in sorted(month_dir.glob("toq*.mseed")):
            ev_idx = int(mseed_path.stem.lstrip("toq"))
            if ev_idx < 1 or ev_idx > len(events):
                continue
            h = events[ev_idx - 1]["header"]
            ev_lat, ev_lon = h.get("lat"), h.get("lon")
            if ev_lat is None or ev_lon is None:
                counts["skipped:no event location"] = counts.get("skipped:no event location", 0) + 1
                continue

            st = obspy.read(str(mseed_path))
            out_st = obspy.Stream()
            for station in sorted({tr.stats.station for tr in st}):
                sta_st = st.select(station=station)
                if station in obs_loc:
                    sta_lat, sta_lon = obs_loc[station]
                elif station in onshore_loc:
                    sta_lat, sta_lon = onshore_loc[station]
                else:
                    out_st += sta_st
                    counts["skipped:no station location"] = counts.get("skipped:no station location", 0) + 1
                    continue
                rotated, status = rotate_stream(sta_st, ev_lat, ev_lon, sta_lat, sta_lon,
                                                 obs_orientation, args.max_orientation_error)
                out_st += rotated
                counts[status] = counts.get(status, 0) + 1
            out_st.write(str(out_dir / mseed_path.name), format="MSEED", encoding="STEIM2")

    print("\nPer-station rotation outcome counts:")
    for k, v in sorted(counts.items()):
        print(f"  {k}: {v}")
    print(f"\nWrote rotated waveforms to {WAVEFORMS_ZRT_DIR}")
    print(f"Point 02_nor_to_markers.py / Snuffler at {WAVEFORMS_ZRT_DIR} instead of "
          f"{config.WAVEFORMS_TOQUES_DIR} to pick on Z/R/T.")


if __name__ == "__main__":
    main()
