# EQ_picks — QCF Earthquake Pick Revision Pipeline

## Project goal
QC first-arrival phases from earthquake catalogs for the Queen Charlotte Fault (QCF) passive-source OBS dataset, then export reviewed catalogs for double-difference seismic tomography (HypoDD / TomoDD).

Manual review is done in **Snuffler** (Pyrocko waveform viewer) using pre-loaded phase markers.

## Student context
Student **Anne** is performing the Snuffler picking review for both catalogs. The instructor (andrewgase) runs Step 1 to extract waveform windows and transfers them to Anne via rsync/scp. Anne works Steps 2–4.

---

## Pipeline A — GaMMA catalog (SeisBench + PhaseNet/PickBlue ML picks)

```
[gamma_events.csv + gamma_picks_refined.csv]
        │
        ▼
01_extract_windows.py     (instructor) Extract per-event waveforms → waveforms/
        │
        ▼
02_picks_to_markers.py    Convert GaMMA picks → Pyrocko markers
        │
    [Snuffler]            Anne reviews, saves markers/events_reviewed.markers
        │
        ▼
03_export.py              → catalogs/catalog_final.{csv,qml,hypodd,tomodd}
```

## Pipeline B — TOQUES catalog (SEISAN Nordic format, 574 events, from Liu)

```
[QCF_Pickblue_NLL_liu.nor]
        │
        ▼
01_extract_toques.py      (instructor) Extract per-event waveforms → waveforms_toques/
        │
        ▼
02_nor_to_markers.py      Convert Nordic picks → Pyrocko markers
        │
    [Snuffler]            Anne reviews, saves markers/toques_reviewed.markers
        │
        ▼
03_export.py --toques     → catalogs/toques_catalog_final.{csv,qml,hypodd,tomodd}
```

The TOQUES catalog uses the same raw OBS recordings as Pipeline A but is a different set of
earthquakes. Station codes in the Nordic file (`CB01`–`CB28`) are automatically mapped to
`QCB01`–`QCB28` to match the raw data and waveform trace IDs in Snuffler.

---

## Key files

| File | Purpose |
|---|---|
| `config.py` | All paths and parameters — edit here only |
| `environment.yml` | Conda env (`qcf_eq`) |
| `patch_pyrocko.py` | One-time fix for Snuffler pick-click bug |
| `nordic.py` | SEISAN Nordic format parser (copied from tomoDD_student project) |
| `01_extract_windows.py` | (A) Instructor: extract waveforms from raw OBS + onshore data |
| `01_extract_toques.py` | (B) Instructor: extract waveforms for TOQUES events |
| `02_picks_to_markers.py` | (A) GaMMA CSV → Pyrocko markers |
| `02_nor_to_markers.py` | (B) Nordic .nor → Pyrocko markers |
| `03_export.py` | Export reviewed markers → CSV / QuakeML / HypoDD / TomoDD |
| `STUDENT_GUIDE.md` | Step-by-step guide written for Anne |
| `README.md` | Full technical documentation |

---

## Data paths (from config.py)

- Raw OBS data: `/home/andrewgase/research/QCF/DATA/raw_data/`
- Clock-corrected OBS (QCB06, QCB28): `/home/andrewgase/research/QCF/DATA/corrected_clock_final/`
- GaMMA catalog: `/home/andrewgase/research/QCF/seisbench/apply/`
- TOQUES Nordic catalog: `/home/andrewgase/research/QCF/tomoDD_student/raw_data/toques_catalog/QCF_Pickblue_NLL_liu.nor`
- GaMMA waveforms: `waveforms/` (not in git; transferred separately)
- TOQUES waveforms: `waveforms_toques/` (not in git; transferred separately)
- Markers: `markers/`
- Exports: `catalogs/`

---

## Conda environment

```bash
conda activate qcf_eq
```

## Known issues / patches

- **Snuffler pick-click bug**: loaded phase markers are not selectable without `patch_pyrocko.py`. Run once after env creation.
- Clock corrections applied only to QCB06 and QCB28; QCB10, QCB14, QCB20 are skipped entirely (bad clock/data).

## Seismic network

- OBS network: `YI` (Queen Charlotte Fault deployment)
- Onshore networks: AK, AT, CN, US
- GaMMA pipeline: PhaseNet with PickBlue weights (OBS-trained), GaMMA associator
- TOQUES pipeline: NLLoc-located picks from Liu (Garza et al. 2024 experiment)
