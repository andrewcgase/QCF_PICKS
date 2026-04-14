# QCF Earthquake Pick Revision Pipeline

Re-pick and QC first-arrival phases for the QCF passive-source earthquake catalog,
then export the reviewed catalog for seismic tomography.

---

## Overview

The pipeline takes an existing detection catalog (from an initial SeisBench + GaMMA
run) and produces a carefully QC'd catalog suitable for double-difference tomography
(HypoDD / TomoDD).

```
[detection catalog]
        │
        ▼
01_repick.py        Re-pick P/S on event windows using PhaseNet + PickBlue
        │
        ▼
02_associate.py     Re-associate picks → events using GaMMA (BGMM + eikonal)
        │
        ▼
03_snuffler_prep.py Convert to Pyrocko markers for manual QC in Snuffler
        │
    [Snuffler]      Review waveforms, delete bad picks, adjust times
        │
        ▼
04_export.py        Export reviewed catalog to CSV / QuakeML / HypoDD / TomoDD
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

### Step 1 — Re-pick

```bash
conda activate qcf_eq
python 01_repick.py                    # process all events
python 01_repick.py --max-events 50   # test run on first 50 events
```

Reads `INPUT_EVENTS`, extracts waveform windows around each event from all nearby
stations, and runs PhaseNet + PickBlue.  Output: `catalogs/picks_repicked.csv`.

**Key parameters in `config.py`:**
| Parameter | Default | Description |
|---|---|---|
| `SEISBENCH_WEIGHTS` | `"PickBlue"` | SeisBench model weights |
| `PICK_THRESHOLD_P` | `0.3` | Minimum P probability |
| `PICK_THRESHOLD_S` | `0.3` | Minimum S probability |
| `WINDOW_BEFORE_S` | `30` | Seconds before origin time |
| `WINDOW_AFTER_S` | `90` | Seconds after origin time |
| `MAX_STATION_DIST_KM` | `400` | Distance cutoff for station inclusion |

### Step 2 — Re-associate

```bash
python 02_associate.py
```

Runs GaMMA (BGMM method, eikonal 1-D velocity model) on the new picks.
Output: `catalogs/events_associated.csv`, `catalogs/picks_associated.csv`.

**Key parameters in `config.py`:**
| Parameter | Default | Description |
|---|---|---|
| `GAMMA_CONFIG["min_picks_per_eq"]` | `5` | Minimum picks for association |
| `GAMMA_CONFIG["max_sigma11"]` | `3.0` | Max origin-time uncertainty (s) |
| `EIKONAL_VEL` | `vp=[5.5,5.5,6.7,7.8]` | 1-D P velocity model (km/s) |

### Step 3 — Snuffler QC

```bash
python 03_snuffler_prep.py
```

Creates `markers/events_to_review.markers`.  The script prints the Snuffler
launch command for your data.  Typical usage:

```bash
snuffler \
  /path/to/DATA/corrected_clock_final/QCBxx/HHZ/*.msd \
  /path/to/DATA/raw_data/QCBxx/Data/*/*/*.msd \
  /path/to/DATA/onshore/NET/STA/Data/*/*/*.msd \
  --markers=markers/events_to_review.markers
```

**In Snuffler:**
- `e` / `E` — jump to next / previous event
- `p` / `s` — add a P or S pick at cursor position
- Right-click a pick → Delete, or change phase name
- `f` — toggle filter (use bandpass 2–20 Hz for most phases)
- `File → Save Markers` — save your changes

Save reviewed markers to `markers/events_reviewed.markers`.

### Step 4 — Export

```bash
python 04_export.py
# or specify a different markers file:
python 04_export.py --markers markers/my_revised.markers
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
