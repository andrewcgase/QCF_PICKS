# Student Guide: Reviewing Earthquake Phase Picks

This guide walks you through the science and the mechanics of manually
reviewing automated earthquake picks for the Queen Charlotte Fault (QCF)
catalog. The README has the full setup instructions; this guide tells you
**what to look for** once you're inside Snuffler.

---

## Why are we doing this?

The catalog was built by a machine-learning picker (PhaseNet/PickBlue) followed
by an association algorithm (GaMMA). That pipeline is fast and consistent, but
it makes mistakes:

- Picks placed on noise, instrument glitches, or secondary arrivals
- Good arrivals that were missed entirely
- Systematically early or late picks on a particular station

Your job is to catch those errors. The reviewed picks feed directly into HypoDD
and TomoDD inversions, where pick timing errors translate directly into location
errors. A 0.1 s timing error on a P-wave at 50 km is ~0.5 km of mislocated
hypocentre.

---

## Before you start: the setup

Follow the steps in `README.md` first. The short version:

```bash
conda activate qcf_eq
python patch_pyrocko.py          # only needed once
python 02_picks_to_markers.py --only-extracted --month YYYY-MM
snuffler waveforms/YYYY-MM/ --markers=markers/YYYY-MM.markers
```

Work **one month at a time** — it keeps Snuffler fast and lets you save
progress in batches.

---

## Snuffler orientation

When Snuffler opens you will see a waterfall of seismic traces, one per
station-channel, arranged by time. Coloured vertical lines are the
**pre-loaded phase markers**:

| Colour | Meaning |
|--------|---------|
| Green  | High-confidence pick (PhaseNet score ≥ 0.8) — probably correct |
| Yellow | Medium confidence (0.5–0.8) — check carefully |
| Red    | Two cases: (1) pre-loaded pick with low ML confidence (< 0.5) — assume wrong until proven right; (2) event markers and any pick you add manually always appear red because they have no ML score — red on a manual pick does **not** mean it is bad |

Press `e` to jump to the first event. Each jump centres the view on the
event origin time and loads the associated waveform window.

**Recommended filter:** press `f` to toggle the 2–20 Hz bandpass. Leave it on
for the entire session unless an event is very large (then try 1–10 Hz).

---

## What does a good pick look like?

### P-wave (compressional)

- Sudden onset on **vertical component** (channel ending in `Z` or `HHZ`)
- The first motion is a sharp upward or downward kick — not a gradual increase
- Arrives **before** the S-wave
- At near-station distances (< 50 km), the gap between P and S is just a few
  seconds; at far stations (> 200 km) it can be 30+ seconds

### S-wave (shear)

- Larger amplitude arrival on **horizontal components** (`N` / `E`, or `HHN` /
  `HHE`)
- Often looks like the waveform "wakes up" — amplitude grows suddenly
- Can be hard to pick if the P coda is still ringing

### Both phases

A correct pick sits at the **leading edge of the arrival** — the very first
moment the ground moves, not the peak. If the marker is on the peak or in the
coda, it is late; drag it (or delete and re-pick) to the true onset.

---

## Decision tree for each pick

Work through each coloured marker in the event window:

```
Is there a clear arrival in the waveform near this marker?
│
├── YES → Is the marker on the onset (not the peak)?
│         ├── YES → Leave it. Move on.
│         └── NO  → Delete the marker. Add a new one at the true onset.
│
└── NO  → Is this a quiet station (low SNR) or a dead channel?
          ├── YES, quiet → Delete the pick. Low-SNR picks hurt more than they help.
          └── YES, dead  → Delete the pick. Note the station for the instructor.
```

**When in doubt, delete.** One bad pick in HypoDD pulls the hypocenter
toward the wrong station. A missing pick just means less constraint.

---

## Adding a new pick

1. Press `Escape` to make sure you are not already in picking mode.
2. Press `p` (P-wave) or `s` (S-wave) to enter picking mode — the cursor changes.
3. **Click on the waveform at the onset.** Be precise: zoom in with `+` if needed.
4. Press `Escape` to exit picking mode.

The new marker will be white (uncoloured) and is associated with the current
event. It will be exported just like the original picks.

---

## Efficient workflow tips

- **Keyboard-first.** `e` (next event), `E` (prev event), `f` (filter), `+/-`
  (zoom), `Backspace` (delete selected marker). Mouse only for picking and
  selecting.
- **Set a rhythm.** Give each event ~30–60 seconds. If an event is ambiguous
  after that, mark it for later (or just leave the green picks and move on).
- **Work by station.** For each event, scan top to bottom. Stations are sorted
  by distance — near stations arrive earliest and have the sharpest onsets.
- **Watch for the pattern.** Good events show a moveout: picks step later as
  you go to farther stations. A pick that is wildly early or late compared
  to its neighbours is suspect.
- **Save often.** `File → Save Markers` every 15–20 events. Snuffler does not
  auto-save. Save to `markers/YYYY-MM_reviewed.markers`.

---

## Common mistakes to avoid

| Mistake | Consequence | Fix |
|---------|-------------|-----|
| Picking the peak instead of the onset | Systematically late picks → biased locations | Zoom in; pick the first deviation from the noise floor |
| Leaving a red pick on a flat trace | Pulls hypocenter toward that station | Delete it |
| Adding an S-pick on the vertical component | Phase label wrong → bad Vp/Vs | Only pick S on horizontal channels |
| Forgetting to save before closing | Hours of work lost | Save every 20 events |
| Working on all months at once | Snuffler is slow, easy to lose track | One month per session |

---

## After QC: exporting the catalog

Once you have reviewed a month and saved your markers file:

```bash
conda activate qcf_eq
python 03_export.py --markers markers/YYYY-MM_reviewed.markers
```

Output files appear in `catalogs/`:

| File | Use |
|------|-----|
| `catalog_final.csv` | Quick look, plotting |
| `catalog_final.qml` | ObsPy / SeisComp |
| `hypodd_phase.dat` | Input to `ph2dt` → HypoDD |
| `tomodd_phase.dat` | Input to TomoDD |

---

## Getting help

- **Snuffler crashes / no waveforms:** check that you ran `--only-extracted`
  when building markers.
- **Can't click on markers:** the pyrocko patch was not applied —
  run `python patch_pyrocko.py`.
- **Pick count looks too low after export:** make sure you saved markers in
  Snuffler *before* running `03_export.py`.
- **Anything else:** ask the instructor or open an issue on the repo.
