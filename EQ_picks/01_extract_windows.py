#!/usr/bin/env python3
"""
01_extract_windows.py — Extract waveform windows around each detected event
                         and save as per-event miniSEED files for sharing.

Only stations that contributed picks to each event are extracted.
Files are saved as Steim2-compressed miniSEED, one file per event.

Reads:  config.INPUT_EVENTS   (gamma_events.csv)
        config.INPUT_PICKS    (gamma_picks_refined.csv)
Writes: config.WAVEFORMS_DIR/ ev{event_index}.mseed  (one per event)

Estimated output size
---------------------
  All 10,471 events, 3 channels, 120 s windows: ~1–2 GB compressed
  With --min-picks 8 (~5,000 events):            ~500 MB compressed
  With --channels Z, --min-picks 8:              ~150 MB compressed

Usage
-----
  conda activate qcf_eq
  python 01_extract_windows.py                          # all events, HHZ/HH1/HH2
  python 01_extract_windows.py --min-picks 8            # quality filter
  python 01_extract_windows.py --channels Z             # vertical only
  python 01_extract_windows.py --max-events 50          # test run
  python 01_extract_windows.py --dry-run                # estimate size, no writing
  python 01_extract_windows.py --workers 8              # parallel workers (default: 4)
"""

import argparse
import glob
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import obspy
from obspy import UTCDateTime, read
import pandas as pd
from tqdm import tqdm

import config


# ---------------------------------------------------------------------------
# Data access helpers
# ---------------------------------------------------------------------------

def find_daily_files(data_root, network, station, year, jday,
                     corrected=False, channels=None):
    """
    Return a list of miniSEED file paths for station/day.
    channels: list of channel codes to search for, e.g. ['HHZ', 'HH1', 'HH2']
    """
    if channels is None:
        channels = ["HHZ", "HH1", "HH2", "EDH", "HDH",
                    "BHZ", "BHN", "BHE", "BH1", "BH2"]
    tag = f"{year}_{jday:03d}"
    files = []

    if corrected:
        # corrected_clock_final/QCBxx/CHAN/YI_QCBxx_CHAN__YYYY_DDD*.c1.msd
        for chan in channels:
            chan_dir = data_root / chan
            if chan_dir.exists():
                files += glob.glob(str(chan_dir / f"{network}_{station}_{chan}__{tag}*.msd"))
    elif network == config.OBS_NETWORK:
        # raw_data/QCBxx/Data/YYYY/DDD/YI_QCBxx_CHAN__YYYY_DDD.msd
        day_dir = data_root / str(year) / f"{jday:03d}"
        if day_dir.exists():
            for chan in channels:
                files += glob.glob(str(day_dir / f"{network}_{station}_{chan}__{tag}*.msd"))
    else:
        # onshore/NET/STA/Data/YYYY/DDD/NET_STA_CHAN__YYYY_DDD.msd
        day_dir = data_root / str(year) / f"{jday:03d}"
        if day_dir.exists():
            for chan in channels:
                files += glob.glob(str(day_dir / f"{network}_{station}_{chan}__{tag}*.msd"))

    return files


def build_station_index():
    """
    Return dict: station_code → {network, data_root, corrected}

    OBS stations are discovered by scanning the data directories on disk
    (raw_data/ and corrected_clock_final/) rather than relying on the
    metadata CSV, whose station-name column may differ from directory names.
    """
    index = {}

    # OBS: clock-corrected stations (QCB06, QCB28)
    if config.OBS_CORRECTED.exists():
        for sta_dir in sorted(config.OBS_CORRECTED.iterdir()):
            sta = sta_dir.name
            if not sta_dir.is_dir() or sta in config.SKIP_STATIONS:
                continue
            index[sta] = {
                "network":   config.OBS_NETWORK,
                "data_root": sta_dir,
                "corrected": True,
            }

    # OBS: raw stations (everything else under raw_data/)
    if config.OBS_RAW.exists():
        for sta_dir in sorted(config.OBS_RAW.iterdir()):
            sta = sta_dir.name
            if not sta_dir.is_dir() or sta in config.SKIP_STATIONS:
                continue
            if sta in config.CLOCK_CORRECTED_STATIONS:
                continue   # already added from corrected_clock_final
            data_root = sta_dir / "Data"
            if data_root.exists():
                index[sta] = {
                    "network":   config.OBS_NETWORK,
                    "data_root": data_root,
                    "corrected": False,
                }

    # Onshore stations from StationXML inventory
    if config.ONSHORE_XML.exists():
        from obspy import read_inventory
        inv = read_inventory(str(config.ONSHORE_XML))
        for network in inv:
            for station in network:
                sta_dir = config.ONSHORE / network.code / station.code / "Data"
                if sta_dir.exists():
                    index[station.code] = {
                        "network":   network.code,
                        "data_root": sta_dir,
                        "corrected": False,
                    }

    print(f"Station index: {len(index)} stations "
          f"({', '.join(sorted(index)[:8])}{'...' if len(index) > 8 else ''})")
    return index


# ---------------------------------------------------------------------------
# Per-day worker (runs in a subprocess)
# ---------------------------------------------------------------------------

def process_day(args_tuple):
    """
    Process all events that fall on a single calendar day.
    Returns (written, skipped, already_done) counts.
    """
    (day_events,        # list of dicts (one per event)
     station_index,     # dict: sta_code → {network, data_root, corrected}
     channels,          # list of channel codes
     before_s,          # float
     after_s,           # float
     out_dir,           # Path
     skip_existing,     # bool
     ) = args_tuple

    stream_cache = {}
    written = skipped = already_done = 0

    for ev in day_events:
        ev_idx   = int(ev["event_index"])
        ot       = UTCDateTime(ev["time"])
        # Monthly subdirectory based on event origin time (window may cross
        # month boundary but the file is always placed by origin month)
        ev_month = f"{ot.year:04d}-{ot.month:02d}"
        out_path = Path(out_dir) / ev_month / f"ev{ev_idx}.mseed"

        if skip_existing and out_path.exists():
            already_done += 1
            continue

        t0     = ot - before_s
        t1     = ot + after_s
        year   = int(ev["_year"])
        jday   = int(ev["_jday"])

        days_needed = {(year, jday)}
        if t0.julday != jday or t0.year != year:
            days_needed.add((t0.year, t0.julday))
        if t1.julday != jday or t1.year != year:
            days_needed.add((t1.year, t1.julday))

        event_stream = obspy.Stream()

        for sta, info in station_index.items():

            for (yr, jd) in days_needed:
                cache_key = (sta, yr, jd)
                if cache_key not in stream_cache:
                    files = find_daily_files(
                        info["data_root"], info["network"],
                        sta, yr, jd,
                        corrected=info["corrected"],
                        channels=channels,
                    )
                    if files:
                        try:
                            st = obspy.Stream()
                            for f in files:
                                st += read(f, dtype="float64")
                            st.merge(fill_value="interpolate")
                            stream_cache[cache_key] = st
                        except Exception:
                            stream_cache[cache_key] = obspy.Stream()
                    else:
                        stream_cache[cache_key] = obspy.Stream()

                event_stream += stream_cache[cache_key]

        if not event_stream:
            skipped += 1
            continue

        try:
            win = event_stream.slice(starttime=t0, endtime=t1).copy()
            win.detrend("demean")
            if not win:
                skipped += 1
                continue
            for tr in win:
                tr.data = tr.data.astype(np.int32)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            win.write(str(out_path), format="MSEED", encoding="STEIM2")
            written += 1
        except Exception:
            skipped += 1

    return written, skipped, already_done


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-picks", type=int, default=5,
                        help="Skip events with fewer than N picks (default: 5).")
    parser.add_argument("--min-score", type=float, default=0.0,
                        help="Skip events with gamma_score below this value (default: 0).")
    parser.add_argument("--channels", default="all",
                        choices=["all", "Z"],
                        help="'all' → HHZ/HH1/HH2/BH*; 'Z' → vertical only (default: all).")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Process only the first N qualifying events (testing).")
    parser.add_argument("--before", type=float, default=config.WINDOW_BEFORE_S,
                        help=f"Seconds before origin time (default: {config.WINDOW_BEFORE_S}).")
    parser.add_argument("--after", type=float, default=config.WINDOW_AFTER_S,
                        help=f"Seconds after origin time (default: {config.WINDOW_AFTER_S}).")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel worker processes (default: 4).")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip events whose output file already exists (resume a partial run).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Estimate output size and exit without writing files.")
    args = parser.parse_args()

    config.WAVEFORMS_DIR.mkdir(parents=True, exist_ok=True)

    if args.channels == "Z":
        channels = ["HHZ", "BHZ"]
    else:
        channels = ["HHZ", "HH1", "HH2", "HDH", "EDH",
                    "BHZ", "BHN", "BHE", "BH1", "BH2"]

    # Load catalog and picks
    events_df = pd.read_csv(config.INPUT_EVENTS)
    picks_df  = pd.read_csv(config.INPUT_PICKS)

    # Filter events
    if "num_picks" in events_df.columns:
        events_df = events_df[events_df["num_picks"] >= args.min_picks]
    if "gamma_score" in events_df.columns and args.min_score > 0:
        events_df = events_df[events_df["gamma_score"] >= args.min_score]
    if args.max_events:
        events_df = events_df.head(args.max_events)

    print(f"Events to extract: {len(events_df)}"
          f"  (min_picks≥{args.min_picks}, min_score≥{args.min_score})")

    def sta_code(station_id):
        parts = str(station_id).rstrip(".").split(".")
        return parts[1] if len(parts) > 1 else parts[0]

    event_stations = defaultdict(set)
    for _, row in picks_df[picks_df["event_index"] >= 0].iterrows():
        event_stations[int(row["event_index"])].add(sta_code(row["station_id"]))

    station_index = build_station_index()

    # Dry run: estimate size (all stations × all events)
    win_len = args.before + args.after
    n_stations = len(station_index)
    n_chan = 1 if args.channels == "Z" else 3
    est_bytes = len(events_df) * n_stations * n_chan * win_len * 100 * 4
    est_compressed = est_bytes * 0.25
    print(f"Estimated output: {est_bytes/1e9:.2f} GB uncompressed, "
          f"~{est_compressed/1e9:.2f} GB compressed (Steim2)")
    if args.dry_run:
        return

    # Group events by calendar day for cache-efficient I/O
    events_df = events_df.copy()
    events_df["_dt"]   = pd.to_datetime(events_df["time"])
    events_df["_year"] = events_df["_dt"].dt.year
    events_df["_jday"] = events_df["_dt"].dt.day_of_year

    day_groups = {}
    for _, ev in events_df.iterrows():
        key = (int(ev["_year"]), int(ev["_jday"]))
        day_groups.setdefault(key, []).append(ev.to_dict())

    print(f"Calendar days spanned: {len(day_groups)}  |  workers: {args.workers}")

    work_items = [
        (evs, station_index, channels, args.before, args.after,
         config.WAVEFORMS_DIR, args.skip_existing)
        for evs in day_groups.values()
    ]

    total_written = total_skipped = total_already_done = 0

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_day, item): item for item in work_items}
        with tqdm(total=len(events_df), desc="Events") as pbar:
            for fut in as_completed(futures):
                try:
                    w, s, d = fut.result()
                    total_written      += w
                    total_skipped      += s
                    total_already_done += d
                    pbar.update(w + s + d)
                except Exception as exc:
                    day_evs = futures[fut][0]
                    day_key = (day_evs[0].get("_year"), day_evs[0].get("_jday"))
                    print(f"\nWorker crashed on day {day_key} ({len(day_evs)} events): {exc}",
                          flush=True)
                    pbar.update(len(day_evs))

    print(f"\nDone: {total_written} events written to {config.WAVEFORMS_DIR}, "
          f"{total_skipped} skipped (no waveform data), "
          f"{total_already_done} skipped (already existed)")
    print(f"Next step — build markers:")
    print(f"  python 02_picks_to_markers.py --only-extracted")
    print(f"  snuffler {config.WAVEFORMS_DIR}/*.mseed "
          f"--markers={config.MARKERS_FOR_QC}")


if __name__ == "__main__":
    main()
