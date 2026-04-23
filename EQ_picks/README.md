# QCF Earthquake Pick Revision Pipeline

QC first-arrival phases for the QCF passive-source earthquake catalog, then
export the reviewed catalog for seismic tomography.

---

## Overview

The pipeline takes an existing detection catalog (SeisBench + GaMMA picks) and
produces a QC'd catalog suitable for double-difference tomography
(HypoDD / TomoDD). Manual review is done in Snuffler with the PhaseNet/PickBlue
picks pre-loaded as markers.

```
[GaMMA detection catalog — gamma_events.csv, gamma_picks_refined.csv]
        │
        ▼
01_extract_windows.py     (instructor only) Package per-event waveforms
        │
        ▼
     [transfer waveforms/ to student via rsync or scp]
        │
        ▼
02_picks_to_markers.py    Convert GaMMA picks → Pyrocko markers file
        │
    [Snuffler]            Review waveforms, delete bad picks, adjust times
        │
        ▼
03_export.py              Export reviewed catalog to CSV / QuakeML / HypoDD / TomoDD
```

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/andrewcgase/QCF_PICKS.git
cd QCF_PICKS/EQ_picks
```

### 2. Create and activate the conda environment

```bash
conda env create -f environment.yml
conda activate qcf_eq
```

### 3. Apply the Snuffler patch

Snuffler has a bug in pick-click detection that prevents loaded phase markers
from being selectable. Run this once after creating the environment:

```bash
conda activate qcf_eq
python patch_pyrocko.py
```

The script is idempotent — running it a second time does nothing.

### 4. Edit `config.py`

Open `config.py` and update the two path variables for your machine.
The student workflow only needs `INPUT_EVENTS` and `INPUT_PICKS`
(the GaMMA catalog CSVs); you do **not** need access to the raw seismograph
data (that was only required for Step 1, which the instructor already ran).

```python
SEISBENCH_DIR = Path("/your/path/to/seisbench/apply")   # contains gamma_events.csv etc.
```

Everything else (`WAVEFORMS_DIR`, `MARKERS_DIR`, `CATALOGS_DIR`) is relative to
the repo directory and requires no changes.

---

## Receiving waveforms from the instructor

The waveform package is too large for git. The instructor will transfer
`waveforms/` to you separately.

**Option A — rsync (recommended, resumable):**

```bash
rsync -avz --progress \
    instructor@hostname:/path/to/QCF_PICKS/EQ_picks/waveforms/ \
    waveforms/
```

**Option B — scp:**

```bash
scp -r instructor@hostname:/path/to/QCF_PICKS/EQ_picks/waveforms/ \
    waveforms/
```

After transfer, the directory should look like:

```
waveforms/
├── 2021-10/
│   ├── ev1234.mseed
│   └── ev1235.mseed
├── 2021-11/
│   └── ...
└── 2022-04/
    └── ...
```

Each file is a 3-minute window (60 s before + 120 s after origin time) for all
available stations, saved as Steim2-compressed miniSEED.

---

## Step 2 — Build Snuffler markers

Convert the GaMMA pick catalog to a Pyrocko markers file.

```bash
conda activate qcf_eq

# All events that have an extracted waveform file:
python 02_picks_to_markers.py --only-extracted

# To work one month at a time (recommended — faster Snuffler load):
python 02_picks_to_markers.py --only-extracted --month 2022-04
```

Options:

| Flag | Default | Description |
|---|---|---|
| `--only-extracted` | off | Only include events with a waveform file in `waveforms/` |
| `--month YYYY-MM` | all | Restrict to a single calendar month |
| `--min-picks N` | 5 | Skip events with fewer than N associated picks |
| `--max-events N` | all | Cap at N events (useful for testing) |
| `--out PATH` | auto | Output markers file path |

Output: `markers/events_to_review.markers` (or `markers/YYYY-MM.markers` when
`--month` is used).

---

## Step 3 — QC in Snuffler

### Launch

```bash
# All waveforms at once:
snuffler waveforms/ --markers=markers/events_to_review.markers

# One month only (faster):
snuffler waveforms/2022-04/ --markers=markers/2022-04.markers
```

### Navigation

| Key | Action |
|---|---|
| `e` | Jump to next event |
| `E` | Jump to previous event |
| `f` | Toggle bandpass filter (2–20 Hz recommended) |
| `+` / `-` | Zoom time axis in / out |
| `s` (window) | Scale traces |
| `a` / `A` | Scroll earlier / later |

### Reviewing picks

Loaded phase markers appear as coloured vertical lines on the traces.
**Green** = high PhaseNet confidence (≥0.8), **yellow** = medium, **red** = low.

**To delete a bad pick:**

1. Press `Escape` to exit picking mode (if active).
2. Left-click the pick marker to select it (it will highlight).
3. Press `Backspace` to delete.

**To add a new pick:**

1. Press `p` (P-wave) or `s` (S-wave) to enter picking mode.
2. Left-click on the waveform at the desired arrival time.
3. Press `Escape` to exit picking mode.

**To move a pick:**

Delete the existing pick and add a new one at the correct time.

### Saving

When finished reviewing a session:

```
File → Save Markers
```

Save to `markers/events_reviewed.markers` (or append a date/suffix if working
in batches, e.g. `markers/2022-04_reviewed.markers`).

**Your edits are not auto-saved.** Save before closing Snuffler.

---

## Step 4 — Export

```bash
conda activate qcf_eq

# Export using the default reviewed markers file:
python 03_export.py

# Or specify a custom file:
python 03_export.py --markers markers/2022-04_reviewed.markers
```

Output files written to `catalogs/`:

| File | Format | Use for |
|---|---|---|
| `catalog_final.csv` | CSV | General analysis, plotting |
| `catalog_final.qml` | QuakeML | ObsPy / SeisComp workflows |
| `hypodd_phase.dat` | HypoDD `phase.dat` | Input to `ph2dt` → HypoDD |
| `tomodd_phase.dat` | TomoDD | Input to TomoDD inversion |

Pick weights in HypoDD/TomoDD output are derived from PhaseNet confidence
scores (high confidence → weight 1.0; low confidence → weight 0.2).

---

## File structure

```
EQ_picks/
├── config.py                   Edit paths for your machine
├── environment.yml             Conda environment definition
├── patch_pyrocko.py            Run once after conda env setup
│
├── 01_extract_windows.py       (instructor) Extract waveform windows
├── 02_picks_to_markers.py      Build Pyrocko markers from GaMMA picks
├── 03_export.py                Export reviewed catalog
│
├── waveforms/                  NOT in git — transfer from instructor
│   ├── 2021-10/
│   │   └── ev{N}.mseed
│   └── YYYY-MM/
│
├── markers/                    Created locally
│   ├── events_to_review.markers   Generated by 02_picks_to_markers.py
│   └── events_reviewed.markers    Saved from Snuffler after QC
│
└── catalogs/                   Created locally
    ├── catalog_final.csv
    ├── catalog_final.qml
    ├── hypodd_phase.dat
    └── tomodd_phase.dat
```

---

## Troubleshooting

**Snuffler loads but I see no waveforms for an event**
- The event may not have had a waveform extracted (low pick count, data gap).
  Use `--only-extracted` when building markers.

**Phase markers are visible but clicking them does nothing**
- The pyrocko patch was not applied. Run `python patch_pyrocko.py`.

**`03_export.py` reports 0 picks**
- Make sure you saved markers from Snuffler before running export.
- Check that the markers file path matches what `--markers` expects.

**Snuffler is very slow to load**
- Work one month at a time: `--month YYYY-MM` when building markers, then
  point Snuffler at `waveforms/YYYY-MM/` instead of `waveforms/`.

---

## Acknowledgements

- [SeisBench](https://github.com/seisbench/seisbench) — ML phase picker framework
- [PickBlue](https://doi.org/10.1029/2023JB026518) — OBS-trained PhaseNet weights
- [GaMMA](https://github.com/AI4EPS/GaMMA) — Gaussian Mixture Model Association
- [Pyrocko](https://pyrocko.org) — Seismology toolkit and Snuffler waveform viewer
