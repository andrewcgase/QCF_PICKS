"""
config.py — Central configuration for the EQ_picks pipeline.

All scripts import from here.  Edit paths and parameters in one place.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Data directories
# ---------------------------------------------------------------------------
DATA_ROOT     = Path("/Users/andrew/research/QCF/DATA")

# OBS miniSEED: raw_data/QCBxx/Data/YYYY/DDD/YI_QCBxx_CHAN__YYYY_DDD.msd
OBS_RAW       = DATA_ROOT / "raw_data"

# Clock-corrected OBS (QCB06, QCB28 only):
# corrected_clock_final/QCBxx/CHAN/YI_QCBxx_CHAN__YYYY_DDD.c1.msd
OBS_CORRECTED = DATA_ROOT / "corrected_clock_final"

# Onshore: onshore/NET/STA/Data/YYYY/DDD/NET_STA_CHAN__YYYY_DDD.msd
ONSHORE       = DATA_ROOT / "onshore"

# Stations to pull from corrected_clock_final instead of raw_data
CLOCK_CORRECTED_STATIONS = {"QCB06", "QCB28"}

# Stations to skip entirely (bad clock / bad data)
SKIP_STATIONS = {"QCB10", "QCB14", "QCB20"}

# OBS network code
OBS_NETWORK = "YI"

# Onshore networks present in DATA/onshore/
ONSHORE_NETWORKS = ["AK", "AT", "CN", "US"]

# ---------------------------------------------------------------------------
# Station metadata
# ---------------------------------------------------------------------------
SEISBENCH_DIR = Path("/Users/andrew/research/QCF/seisbench/apply")
ONSHORE_XML   = SEISBENCH_DIR / "xml" / "onshore.xml"
OBS_META_CSV  = Path(__file__).parent.parent / "WA_picks" / "EarthScope_OBS_Deployments.csv"

# ---------------------------------------------------------------------------
# Input detection catalog (from initial SeisBench + GaMMA run)
# ---------------------------------------------------------------------------
INPUT_EVENTS = SEISBENCH_DIR / "gamma_events.csv"
INPUT_PICKS  = SEISBENCH_DIR / "gamma_picks_refined.csv"

# ---------------------------------------------------------------------------
# Output paths (all under EQ_picks/)
# ---------------------------------------------------------------------------
BASE         = Path(__file__).parent
CATALOGS_DIR = BASE / "catalogs"
MARKERS_DIR  = BASE / "markers"

# Step 1 output
PICKS_REPICKED   = CATALOGS_DIR / "picks_repicked.csv"

# Step 2 output
EVENTS_ASSOC     = CATALOGS_DIR / "events_associated.csv"
PICKS_ASSOC      = CATALOGS_DIR / "picks_associated.csv"

# Step 1 waveform packages
WAVEFORMS_DIR    = BASE / "waveforms"

# Step 2: Pyrocko marker files
MARKERS_FOR_QC   = MARKERS_DIR  / "events_to_review.markers"
MARKERS_REVIEWED = MARKERS_DIR  / "events_reviewed.markers"   # saved by Snuffler after QC

# Step 4 final exports
CATALOG_CSV      = CATALOGS_DIR / "catalog_final.csv"
CATALOG_QML      = CATALOGS_DIR / "catalog_final.qml"
HYPODD_PHASE     = CATALOGS_DIR / "hypodd_phase.dat"
TOMODD_PHASE     = CATALOGS_DIR / "tomodd_phase.dat"

# ---------------------------------------------------------------------------
# SeisBench parameters
# ---------------------------------------------------------------------------
SEISBENCH_WEIGHTS  = "PickBlue"    # sbm.PhaseNet.from_pretrained("PickBlue")
PICK_THRESHOLD_P   = 0.3
PICK_THRESHOLD_S   = 0.3
SEISBENCH_BATCH    = 256
WINDOW_BEFORE_S    = 30.0          # seconds before origin time
WINDOW_AFTER_S     = 90.0          # seconds after origin time
MAX_STATION_DIST_KM = 400          # skip stations beyond this distance

# ---------------------------------------------------------------------------
# GaMMA association parameters
# Reproduced from seisbench/apply/apply_SBM.ipynb
# ---------------------------------------------------------------------------

# 1-D velocity model for eikonal travel-time table
EIKONAL_VEL = {
    "z":  [0.0, 5.5, 16.0, 32.0],   # km
    "p":  [5.5, 5.5,  6.7,  7.8],   # km/s
    "s":  [v / 1.73 for v in [5.5, 5.5, 6.7, 7.8]],
}
EIKONAL_DZ = 1.0   # km grid spacing for eikonal table

GAMMA_CONFIG = {
    "method": "BGMM",
    "oversample_factor": 5,          # BGMM oversampling
    "use_amplitude": False,
    "vel": {"p": 6.0, "s": 6.0 / 1.75},
    "dims": ["x(km)", "y(km)", "z(km)"],
    "use_dbscan": True,
    "dbscan_min_samples": 3,
    # dbscan_eps is estimated at runtime from station geometry via estimate_eps()
    "min_picks_per_eq": 5,
    "min_p_picks_per_eq": 0,
    "min_s_picks_per_eq": 0,
    "max_sigma11": 3.0,   # seconds
    "max_sigma22": 1.0,
    "max_sigma12": 1.0,
    "ncpu": 4,            # reduce if memory is limited
    "z(km)": (0, 50),
}

# ---------------------------------------------------------------------------
# HypoDD / TomoDD export
# ---------------------------------------------------------------------------
# Weight assigned to picks by uncertainty rank (used in phase.dat)
PHASE_WEIGHTS = {
    "P": {0.3: 1.0, 0.5: 1.0, 0.7: 1.0, 0.9: 1.0},  # score → weight (linear fallback)
    "S": {0.3: 0.5, 0.5: 0.8, 0.7: 1.0, 0.9: 1.0},
}
