# QCF Earthquake Pick Revision Pipeline

QC first-arrival phases for the QCF passive-source earthquake catalog, then
export the reviewed catalog for seismic tomography.

---

## Overview

The pipeline takes an existing detection catalog (SeisBench + GaMMA picks) and
produces a carefully QC'd catalog suitable for double-difference tomography
(HypoDD / TomoDD).  Manual review is done in Snuffler with the existing
PhaseNet/PickBlue picks pre-loaded as markers.

```
[detection catalog + SeisBench picks]
        │
        ▼
01_extract_windows.py   Package per-event waveforms as miniSEED for sharing
        │
        ▼
02_picks_to_markers.py  Convert SeisBench picks → Pyrocko markers file
        │
    [Snuffler]          Review waveforms, delete bad picks, adjust times
        │
        ▼
03_export.py            Export reviewed catalog to CSV / QuakeML / HypoDD / TomoDD
```

---

## Setup

### 1. Clone and create environment

```bash
git clone https://github.com/andrewcgase/QCF_PICKS.git
cd QCF_PICKS/EQ_picks
conda env create -f environment.yml
conda activate qcf_eq
```

### 2. Edit `config.py`

Open `config.py` and set the paths for your local copy of the data:

```python
DATA_ROOT     = Path("/your/path/to/DATA")      # contains raw_data/, corrected_clock_final/, onshore/
INPUT_EVENTS  = Path("/your/path/to/gamma_events.csv")   # detection catalog
INPUT_PICKS   = Path("/your/path/to/gamma_picks_refined.csv")
ONSHORE_XML   = Path("/your/path/to/onshore.xml")
OBS_META_CSV  = Path("/your/path/to/EarthScope_OBS_Deployments.csv")
```

---

## Running the pipeline

### Step 1 — Extract waveform windows

```bash
conda activate qcf_eq
python 01_extract_windows.py                    # all events (≥5 picks)
python 01_extract_windows.py --min-picks 8      # stricter quality filter
python 01_extract_windows.py --channels Z       # vertical only (smaller)
python 01_extract_windows.py --max-events 50    # test run
python 01_extract_windows.py --dry-run          # estimate size, no writing
```

Reads `INPUT_EVENTS` and `INPUT_PICKS`, extracts a waveform window around each
event's origin time from every station that contributed picks, and saves one
Steim2-compressed miniSEED file per event: `waveforms/ev{event_index}.mseed`.

**Estimated output sizes:**
| Filter | Approx. compressed size |
|---|---|
| All 10,471 events, 3 channels, 120 s | ~1–2 GB |
| `--min-picks 8` (~5,000 events) | ~500 MB |
| `--min-picks 8 --channels Z` | ~150 MB |

**Key parameters in `config.py`:**
| Parameter | Default | Description |
|---|---|---|
| `WINDOW_BEFORE_S` | `30` | Seconds before origin time |
| `WINDOW_AFTER_S` | `90` | Seconds after origin time |

### Step 2 — Build Snuffler markers

```bash
python 02_picks_to_markers.py                   # all filtered events
python 02_picks_to_markers.py --only-extracted  # only events with mseed files
python 02_picks_to_markers.py --min-picks 8
```

Converts `gamma_events.csv` + `gamma_picks_refined.csv` to a Pyrocko
`.markers` file at `markers/events_to_review.markers`.

Then open Snuffler:

```bash
snuffler waveforms/*.mseed --markers=markers/events_to_review.markers
```

**In Snuffler:**
- `e` / `E` — jump to next / previous event
- `p` / `s` — add a P or S pick at cursor position
- Right-click a pick → Delete, or change phase name
- `f` — toggle filter (use bandpass 2–20 Hz for most phases)
- `File → Save Markers` — save your changes

Save reviewed markers to `markers/events_reviewed.markers`.

### Step 3 — Export

```bash
python 03_export.py
# or specify a different markers file:
python 03_export.py --markers markers/my_revised.markers
```

Exports the QC'd catalog to four formats:

| File | Format | Use for |
|---|---|---|
| `catalogs/catalog_final.csv` | CSV | General analysis, plotting |
| `catalogs/catalog_final.qml` | QuakeML | ObsPy / SeisComp workflows |
| `catalogs/hypodd_phase.dat` | HypoDD `phase.dat` | Input to `ph2dt` → HypoDD |
| `catalogs/tomodd_phase.dat` | TomoDD | Input to TomoDD inversion |

---

## Data layout expected

```
DATA/
├── corrected_clock_final/
│   ├── QCB06/
│   │   ├── HHZ/   YI_QCB06_HHZ__YYYY_DDD.c1.msd
│   │   ├── HH1/
│   │   └── HH2/
│   └── QCB28/  ...
├── raw_data/
│   ├── QCB01/
│   │   └── Data/
│   │       └── YYYY/
│   │           └── DDD/  YI_QCB01_HHZ__YYYY_DDD.msd
│   └── QCBxx/ ...
└── onshore/
    ├── AK/
    │   └── S32K/
    │       └── Data/YYYY/DDD/  AK_S32K_BHZ__YYYY_DDD.msd
    └── NET/STA/Data/...
```

**OBS stations excluded from raw_data** (use clock-corrected versions):
QCB06, QCB28

**OBS stations skipped** (bad clock / data):
QCB10, QCB14, QCB20

---

## Output column descriptions

### `catalog_final.csv`
| Column | Description |
|---|---|
| `event_id` | Internal ID from Snuffler marker |
| `time` | Origin time (UTC ISO) |
| `latitude` | Decimal degrees |
| `longitude` | Decimal degrees |
| `depth_km` | Hypocentral depth (km) |
| `magnitude` | Magnitude (if assigned) |

### `hypodd_phase.dat`
Standard HypoDD `phase.dat` format (see HypoDD manual §3.1).
Feed directly to `ph2dt`.  Pick weights are derived from PhaseNet scores.

### `tomodd_phase.dat`
Absolute travel times in the same layout as HypoDD `phase.dat`.
Each event block begins with `#` followed by the event header line.

---

## Acknowledgements

- [SeisBench](https://github.com/seisbench/seisbench) — ML phase picker framework
- [PickBlue](https://doi.org/10.1029/2023JB026518) — OBS-trained PhaseNet weights
- [GaMMA](https://github.com/AI4EPS/GaMMA) — Gaussian Mixture Model Association
- [Pyrocko](https://pyrocko.org) — Seismology toolkit and Snuffler waveform viewer
