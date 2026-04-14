# QCF First-Arrival Picks

Processing pipeline for seismic first-arrival picks from the Queen Charlotte Fault (QCF) survey.
Picks are exported from Reveal as `.tbl` files and processed into master CSVs for use in
seismic tomography.

---

## Directory structure

```
QCF_PICKS/
├── tables/                          # OBS pick files from Anne (128 × .tbl)
│   └── processed/                   # Output CSVs, one per OBS station
├── onshore_tables/                  # Onshore/offshore pick files from Anne (.tbl)
│   └── processed/                   # Output CSVs, one per station × line
├── picks/                           # Josh's original processed data (reference)
├── OBS uncertainties.csv            # Uncertainty table for OBS stations
├── onshore-offshore uncertainties.csv  # Uncertainty table for onshore stations
├── EarthScope_OBS_Deployments.csv   # OBS receiver metadata (lat, lon, depth)
├── onshore_station_locations.csv    # Onshore receiver metadata (Anne fills in)
├── allshots_formatted.csv           # Unified shot table for all survey lines
├── process_anne_picks.py            # Processes OBS picks → tables/processed/
├── process_onshore_picks.py         # Processes onshore picks → onshore_tables/processed/
└── Interpolation_Example.ipynb      # Notebook demonstrating the processing pipeline
```

---

## Output column schema

Both scripts produce CSVs with the following columns:

| Column | Description |
|---|---|
| `FFID` | Field File Identifier (shot number) |
| `Time (Seconds)` | Interpolated first-arrival travel time |
| `Time Boxcar (Seconds)` | Boxcar-smoothed travel time (rolling mean) |
| `Uncertainty (Seconds)` | Pick uncertainty from uncertainty table |
| `FFID Latitude` | Shot latitude (from allshots_formatted.csv) |
| `FFID Longitude` | Shot longitude |
| `FFID Depth (km)` | Source depth (fixed at 0.012 km per project convention) |
| `OBS Number` *(OBS only)* | Station identifier, e.g. `12S01` |
| `OBS Line` *(OBS only)* | Line identifier, e.g. `12S` |
| `OBS Depth (km)` *(OBS only)* | Receiver depth (negative, from EarthScope) |
| `OBS Latitude` *(OBS only)* | Receiver latitude |
| `OBS Longitude` *(OBS only)* | Receiver longitude |
| `Station` *(onshore only)* | Station code, e.g. `BNAB` |
| `Latitude` *(onshore only)* | Receiver latitude |
| `Longitude` *(onshore only)* | Receiver longitude |
| `Elevation (km)` *(onshore only)* | Receiver elevation |

---

## Processing OBS picks (`process_anne_picks.py`)

### What Anne needs to provide
- `.tbl` files in `tables/`, named `{CRUISE}_{STATION}.tbl`
  (e.g. `PO01A_12S01.tbl`, `PO02A_2S21.tbl`)
- Updated rows in `OBS uncertainties.csv` if any uncertainty ranges change

### Uncertainty ranks (OBS)
| Rank | Uncertainty |
|---|---|
| 1 | 24 ms |
| 2 | 50 ms |
| 3 | 100 ms |

### Running
```bash
python process_anne_picks.py                       # all 128 stations
python process_anne_picks.py --station PO01A_12S01 # single station (testing)
python process_anne_picks.py --gap-threshold 50    # stricter gap splitting
python process_anne_picks.py --boxcar-window 5     # wider smoothing
```

Output goes to `tables/processed/`.

---

## Processing onshore picks (`process_onshore_picks.py`)

### What Anne needs to provide
1. **Station coordinates** — fill in `onshore_station_locations.csv` with the
   latitude, longitude, and elevation (in km) for each land station:
   BNAB, CRAG, NDB, PCLB, RUBB, S32K, SIT, U33K, V35K, WRAK

2. **Pick files** — place `.tbl` files in `onshore_tables/`, named
   `{STATION}_{LINE}.tbl` where STATION is the seismic station code and
   LINE is the seismic line name.
   Examples:
   ```
   BNAB_PM01B2.tbl
   SIT_T12A14B.tbl
   U33K_PO01A.tbl
   ```

### Uncertainty ranks (onshore)
| Rank | Uncertainty |
|---|---|
| 1 | 100 ms |
| 2 | 250 ms |
| 3 | 500 ms |

### Running
```bash
python process_onshore_picks.py                        # all files in onshore_tables/
python process_onshore_picks.py --station BNAB_PM01B2  # single file (testing)
python process_onshore_picks.py --gap-threshold 50     # stricter gap splitting
python process_onshore_picks.py --boxcar-window 5      # wider smoothing
```

Output goes to `onshore_tables/processed/`.

---

## Processing pipeline overview

For each `.tbl` file, both scripts apply the same steps:

1. **Parse** — load `[FFID, Time]` pairs from Reveal's JSON format
2. **Split** — divide picks into segments at FFID gaps larger than `--gap-threshold`
   (default 100). This avoids interpolating across shadow zones or phase boundaries
   (e.g. the transition from direct water wave to refracted first arrival)
3. **Interpolate** — linearly fill in every integer FFID within each segment
4. **Smooth** — apply a boxcar (rolling-mean) filter of width `--boxcar-window`
   (default 3) to produce `Time Boxcar`
5. **Undo reduction velocity** *(OBS only)* — OBS picks were made on data displayed
   with a 6 km/s reduction velocity. The stored times are reduced times:
   `T_reduced = T_actual - offset / 6`. The script recovers actual travel times:
   `T_actual = T_reduced + offset / 6`, where offset is the great-circle
   source–receiver distance computed from the merged FFID and OBS positions.
6. **Uncertainty** — assign uncertainty in seconds based on FFID range from the
   appropriate uncertainty CSV; rows outside any defined range are dropped
7. **Shot positions** — merge FFID latitude/longitude from `allshots_formatted.csv`
8. **Receiver metadata** — attach station position and depth/elevation

See `Interpolation_Example.ipynb` for a visual walkthrough of steps 1–4.

---

## Known data quality issues in `onshore-offshore uncertainties.csv`

The following entries appear to contain typos. **Confirm with Anne before processing:**

| Station | Line | Issue |
|---|---|---|
| BNAB | T0607 | `ffid_end = 4900` — should probably be `49000` or `49300` (ffid_start is 48731) |
| PCLB | PM04A | `ffid_start = 40375 > ffid_end = 40122` — range is inverted |
| CRAG | PO01A | `ffid_end = 222650` — likely `22650` (extra leading `2`) |
| U33K | PO01A | `ffid_end = 222639` — likely `22639` (extra leading `2`) |

The script will warn about inverted ranges and skip them automatically.
The `222xxx` values are not inverted so they will not trigger a warning but will
produce an incorrect (oversized) uncertainty range.

---

## Reference files

| File | Description |
|---|---|
| `EarthScope_OBS_Deployments.csv` | Authoritative OBS deployment positions and depths from EarthScope DMC. 148 stations across lines 12S, 1AS, 1BS, 2S, 3S, 8S, QCB. Depths are in km (negative). |
| `allshots_formatted.csv` | 74,069 shots from all survey cruises (PO01A, PO02A, PO12A, PO12, PM*, T*, RM* lines). Columns: shot, date, time, sourceLat, sourceLon, shipLat, shipLon, waterDepth, line. |
| `OBS uncertainties.csv` | Per-cruise, per-station FFID ranges with uncertainty rank for OBS picks. |
| `onshore-offshore uncertainties.csv` | Per-station, per-line FFID ranges with uncertainty rank for onshore station picks. |
