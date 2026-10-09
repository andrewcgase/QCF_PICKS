#!/usr/bin/env python3
"""
01_extract_toques.py — Extract waveform windows for TOQUES OBS catalog events.

Instructor-only step. Reads the TOQUES SEISAN Nordic catalog and packages
per-event 3-minute waveform windows from the raw OBS + onshore data,
saving one miniSEED file per event for transfer to the student.

Station codes in the Nordic file are CB01–CB28; these are mapped to QCB01–QCB28
to match the raw data directory names.

Reads:  config.TOQUES_NOR        (QCF_Pickblue_NLL_liu.nor)
Writes: config.WAVEFORMS_TOQUES_DIR/ toq{N:04d}.mseed  (one per event)

Usage
-----
  conda activate qcf_eq
  python 01_extract_toques.py                  # all 574 events
  python 01_extract_toques.py --max-events 20  # test run
  python 01_extract_toques.py --skip-existing  # resume a partial run
  python 01_extract_toques.py --dry-run        # estimate size, no writing
  python 01_extract_toques.py --workers 8      # parallel workers (default: 4)
"""

import argparse
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import obspy
from obspy import UTCDateTime, read
from tqdm import tqdm

import config
import nordic


# ---------------------------------------------------------------------------
# Data access helpers  (same logic as 01_extract_windows.py)
# ---------------------------------------------------------------------------

def find_daily_files(data_root, network, station, year, jday,
                     corrected=False, channels=None):
    if channels is None:
        channels = ["HHZ", "HH1", "HH2", "EDH", "HDH",
                    "BHZ", "BHN", "BHE", "BH1", "BH2"]
    tag = f"{year}_{jday:03d}"
    files = []
    if corrected:
        for chan in channels:
            chan_dir = data_root / chan
            if chan_dir.exists():
                files += glob.glob(str(chan_dir / f"{network}_{station}_{chan}__{tag}*.msd"))
    else:
        day_dir = data_root / str(year) / f"{jday:03d}"
        if day_dir.exists():
            for chan in channels:
                files += glob.glob(str(day_dir / f"{network}_{station}_{chan}__{tag}*.msd"))
    return files


def build_station_index():
    index = {}
    if config.OBS_CORRECTED.exists():
        for sta_dir in sorted(config.OBS_CORRECTED.iterdir()):
            sta = sta_dir.name
            if not sta_dir.is_dir() or sta in config.SKIP_STATIONS:
                continue
            index[sta] = {"network": config.OBS_NETWORK,
                          "data_root": sta_dir, "corrected": True}
    if config.OBS_RAW.exists():
        for sta_dir in sorted(config.OBS_RAW.iterdir()):
            sta = sta_dir.name
            if not sta_dir.is_dir() or sta in config.SKIP_STATIONS:
                continue
            if sta in config.CLOCK_CORRECTED_STATIONS:
                continue
            data_root = sta_dir / "Data"
            if data_root.exists():
                index[sta] = {"network": config.OBS_NETWORK,
                              "data_root": data_root, "corrected": False}
    if config.ONSHORE_XML.exists():
        from obspy import read_inventory
        inv = read_inventory(str(config.ONSHORE_XML))
        for network in inv:
            for station in network:
                sta_dir = config.ONSHORE / network.code / station.code / "Data"
                if sta_dir.exists():
                    index[station.code] = {"network": network.code,
                                           "data_root": sta_dir, "corrected": False}
    print(f"Station index: {len(index)} stations "
          f"({', '.join(sorted(index)[:8])}{'...' if len(index) > 8 else ''})")
    return index


# ---------------------------------------------------------------------------
# Per-day worker (runs in a subprocess)
# ---------------------------------------------------------------------------

def process_day_toques(args_tuple):
    """Process all TOQUES events that fall on a single calendar day."""
    (day_events, station_index, channels, before_s, after_s,
     out_dir, skip_existing) = args_tuple

    stream_cache = {}
    written = skipped = already_done = 0

    for ev in day_events:
        ev_label = ev["ev_label"]
        ot = UTCDateTime(ev["time"])
        ev_month = f"{ot.year:04d}-{ot.month:02d}"
        out_path = Path(out_dir) / ev_month / f"{ev_label}.mseed"

        if skip_existing and out_path.exists():
            already_done += 1
            continue

        t0 = ot - before_s
        t1 = ot + after_s

        days_needed = {(ot.year, ot.julday)}
        if (t0.year, t0.julday) != (ot.year, ot.julday):
            days_needed.add((t0.year, t0.julday))
        if (t1.year, t1.julday) != (ot.year, ot.julday):
            days_needed.add((t1.year, t1.julday))

        event_stream = obspy.Stream()

        for sta, info in station_index.items():
            for (yr, jd) in days_needed:
                cache_key = (sta, yr, jd)
                if cache_key not in stream_cache:
                    files = find_daily_files(
                        info["data_root"], info["network"], sta, yr, jd,
                        corrected=info["corrected"], channels=channels,
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
    parser.add_argument("--channels", default="all", choices=["all", "Z"],
                        help="'all' → HHZ/HH1/HH2/BH*; 'Z' → vertical only (default: all).")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Process only the first N events (for testing).")
    parser.add_argument("--before", type=float, default=config.WINDOW_BEFORE_S,
                        help=f"Seconds before origin time (default: {config.WINDOW_BEFORE_S}).")
    parser.add_argument("--after", type=float, default=config.WINDOW_AFTER_S,
                        help=f"Seconds after origin time (default: {config.WINDOW_AFTER_S}).")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel worker processes (default: 4).")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip events whose output file already exists (resume).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Estimate output size and exit without writing files.")
    args = parser.parse_args()

    config.WAVEFORMS_TOQUES_DIR.mkdir(parents=True, exist_ok=True)

    channels = (["HHZ", "BHZ"] if args.channels == "Z"
                else ["HHZ", "HH1", "HH2", "HDH", "EDH",
                      "BHZ", "BHN", "BHE", "BH1", "BH2"])

    events_raw = nordic.read_nor_file(str(config.TOQUES_NOR))
    print(f"TOQUES catalog: {len(events_raw)} events in {config.TOQUES_NOR.name}")

    if args.max_events:
        events_raw = events_raw[:args.max_events]

    events = []
    for i, ev in enumerate(events_raw, start=1):
        try:
            dt = nordic.event_datetime(ev["header"])
            events.append({
                "ev_label": f"toq{i:04d}",
                "time": dt.strftime("%Y-%m-%dT%H:%M:%S.%f"),
            })
        except Exception:
            continue

    print(f"Events to extract: {len(events)}")

    station_index = build_station_index()

    if args.dry_run:
        win_len = args.before + args.after
        n_chan = 1 if args.channels == "Z" else 3
        est = len(events) * len(station_index) * n_chan * win_len * 100 * 4
        print(f"Estimated output: {est/1e9:.2f} GB uncompressed, "
              f"~{est*0.25/1e9:.2f} GB compressed (Steim2)")
        return

    # Group events by calendar day for cache-efficient I/O
    day_groups = {}
    for ev in events:
        ot = UTCDateTime(ev["time"])
        key = (ot.year, ot.julday)
        day_groups.setdefault(key, []).append(ev)

    print(f"Calendar days spanned: {len(day_groups)}  |  workers: {args.workers}")

    work_items = [
        (evs, station_index, channels, args.before, args.after,
         config.WAVEFORMS_TOQUES_DIR, args.skip_existing)
        for evs in day_groups.values()
    ]

    total_written = total_skipped = total_already_done = 0

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_day_toques, item): item for item in work_items}
        with tqdm(total=len(events), desc="Events") as pbar:
            for fut in as_completed(futures):
                try:
                    w, s, d = fut.result()
                    total_written      += w
                    total_skipped      += s
                    total_already_done += d
                    pbar.update(w + s + d)
                except Exception as exc:
                    day_evs = futures[fut][0]
                    print(f"\nWorker crashed ({len(day_evs)} events): {exc}", flush=True)
                    pbar.update(len(day_evs))

    print(f"\nDone: {total_written} written to {config.WAVEFORMS_TOQUES_DIR}, "
          f"{total_skipped} skipped (no data), "
          f"{total_already_done} skipped (already existed)")
    print(f"\nNext step:")
    print(f"  python 02_nor_to_markers.py --only-extracted")
    print(f"  snuffler {config.WAVEFORMS_TOQUES_DIR}/ "
          f"--markers={config.MARKERS_TOQUES_FOR_QC}")


if __name__ == "__main__":
    main()
