# OBS horizontal component orientations (network YI)

Downloaded 2026-10-29 from https://obsic-metrics.whoi.edu/YI (OBSIC metrics
hub, WHOI) -- the `YI_experiment_metrics.csv` link on that page, plus its
linked methodology PDF.

Covers all 25 long-term QCB stations (the ones feeding the TOQUES catalog --
see `nordic.py`'s CB01-CB28 -> QCB01-QCB28 mapping). The other ~120 rows
(short-term 12S*/etc. deployments) have no orientation measured (`--`).

## Method and convention (from the PDF)

Rayleigh-wave arrival-angle method (DLOPy; Doran & Laske 2017, BSSA
107(2):691-708, doi:10.1785/0120160165).

**Orientation, degrees**: the clockwise angle *theta* from true North to the
H1 component, 0-360 deg. **Error, degrees**: its uncertainty.

Confirmed directly against DLOPy's own source (`pysave/locfuns.py`'s
`rot2d(N, E, alpha)`, used internally to test trial rotation angles) that
applying this reported theta to raw H1/H2 recovers true geographic N/E as:

    N = H1*cos(theta) - H2*sin(theta)
    E = H1*sin(theta) + H2*cos(theta)

`scripts/rotate_to_zrt.py` uses this formula for OBS stations, then ObsPy's
standard `rotate_ne_rt` (unambiguous, no OBS-specific convention needed) for
the N/E -> R/T step using each event's preliminary location.

## Data quality

Several stations have very large uncertainty and should probably not be
rotated with confidence -- e.g. QCB10 (345 +/- 174 deg, i.e. essentially
unconstrained), QCB13 (108 +/- 116), QCB14 (248 +/- 281). See
`scripts/rotate_to_zrt.py --max-orientation-error` (default 30 deg) -- any
OBS station whose Error exceeds this is left as Z/1/2 (unrotated), not
silently rotated with a near-meaningless angle.

## Regenerating

```
curl -sL https://obsic-metrics.whoi.edu/YI/csv -o YI_experiment_metrics.csv
curl -sL https://obsic-metrics.whoi.edu/static/images/AACSE_orientations_Readme.pdf -o AACSE_orientations_Readme.pdf
```
