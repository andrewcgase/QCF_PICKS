"""SEISAN Nordic-format parsers, for the earthquake pick/catalog files Liu
supplied in raw_data/{toques_catalog,craig_aftershock_catalog}/.

Two related formats appear in this project:

1. Full Nordic bulletin (.nor) -- one file, many events, each event a type-1
   hypocenter header line, some type-3/E/I comment/error lines (skipped),
   a type-6 waveform-file line (skipped), a type-7 phase-column-header line
   (skipped, but see NOTE below), then one phase-pick line per station/phase,
   ending at a blank line. Only raw_data/toques_catalog/*.nor has this --
   raw_data/craig_aftershock_catalog/ only has the summary format below (no
   picks yet).

2. Hypocenter summary (*_headers_*.txt) -- one line per event, a simplified
   whitespace-delimited restatement of the same type-1 header (verified by
   cross-checking several events' fields against their .nor type-1 line --
   they match field-for-field). Both raw_data/ catalogs have this.

NOTE on column positions: fixed-width Nordic files are notorious for the
printed column-header LABEL not landing exactly on the data field's actual
start column (the label text is placed for human readability, not as a
strict ruler). Learned this the hard way here -- an initial header-label-
based guess put the SECON and TRES fields one column early, silently
dropping the seconds' leading digit (e.g. "31.080" truncated to "1.080").
Caught and fixed by cross-validating against physical sanity (pick time
compared to the event's own origin time) before trusting any of this. The
column positions below are verified this way, not just read off the label
line -- see git history / session notes if this needs re-deriving for a
differently-formatted Nordic export.
"""

from __future__ import annotations

import datetime


def _to_int(s, default=None):
    s = s.strip()
    return int(s) if s else default


def _to_float(s, default=None):
    s = s.strip()
    return float(s) if s else default


def parse_type1_line(line):
    """Parse a Nordic type-1 (hypocenter header) line -- verified column-by-
    column against real data (see module docstring)."""
    return dict(
        year=_to_int(line[1:5]),
        month=_to_int(line[6:8]),
        day=_to_int(line[8:10]),
        hour=_to_int(line[11:13]),
        minute=_to_int(line[13:15]),
        second=_to_float(line[16:20]),
        lat=_to_float(line[23:30]),
        lon=_to_float(line[30:38]),
        depth=_to_float(line[38:43]),
        agency=line[45:48].strip(),
        n_stations=_to_int(line[48:51], 0),
        rms=_to_float(line[51:55]),
        mag1=_to_float(line[55:59]),
        mag1_type=line[59].strip(),
        mag1_agency=line[60:63].strip(),
        mag2=_to_float(line[63:67]),
        mag2_type=line[67].strip(),
        mag2_agency=line[68:71].strip(),
    )


def parse_pick_line(line):
    """Parse a Nordic phase-pick line -- verified column-by-column against
    real data (see module docstring; SECON/TRES columns corrected from the
    header-label positions)."""
    station = line[1:6].strip()
    component = line[6:8].strip()
    onset = line[9].strip()
    phase = line[10:14].strip()
    weight = line[14].strip()
    first_motion = line[16].strip()
    hour = _to_int(line[18:20])
    minute = _to_int(line[20:22])
    second = _to_float(line[22:28])
    residual = _to_float(line[64:69])
    distance_km = _to_float(line[72:75])
    azimuth = _to_float(line[76:79])
    return dict(station=station, component=component, onset=onset, phase=phase,
                weight=weight, first_motion=first_motion, hour=hour, minute=minute,
                second=second, residual=residual, distance_km=distance_km, azimuth=azimuth)


def read_nor_file(path, encoding="latin-1"):
    """Read a full Nordic bulletin (.nor) file: returns a list of dicts, one
    per event, each {'header': {...from parse_type1_line...}, 'picks': [...]}."""
    events = []
    current = None
    with open(path, encoding=encoding) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip():
                if current is not None:
                    events.append(current)
                    current = None
                continue
            if len(line) >= 80 and line[79] == "1" and line[1:5].strip().isdigit():
                current = {"header": parse_type1_line(line), "picks": []}
                continue
            if current is None:
                continue  # comment/ID/error lines before we've seen a type-1 (shouldn't happen)
            if len(line) >= 80 and line[79] in "3EI67":
                continue  # comment, error, ID, waveform-file, or phase-column-header line
            if len(line) >= 8 and line[6:8] in ("HZ", "H1", "H2") and line[10:11] in ("P", "S"):
                current["picks"].append(parse_pick_line(line))
    if current is not None:
        events.append(current)
    return events


def read_hypocenter_summary(path):
    """Read a *_headers_*.txt hypocenter summary: returns a list of dicts,
    one per event (no picks -- this format doesn't have them)."""
    events = []
    for line in open(path):
        parts = line.split()
        if len(parts) < 15:
            continue
        date, time = parts[0], parts[1]
        year, month, day = (int(x) for x in date.split("-"))
        hour, minute, sec_s = time.split(":")
        events.append(dict(
            year=year, month=month, day=day, hour=int(hour), minute=int(minute),
            second=float(sec_s), lat=float(parts[2]), lon=float(parts[3]),
            depth=float(parts[4]), err_h1=float(parts[5]), err_h2=float(parts[6]),
            err_z=float(parts[7]), mag1=float(parts[8]), n_stations=int(parts[9]),
            rms=float(parts[10]), agency=parts[11], mag1_agency=parts[12],
            gap=float(parts[13]),
            mag2=None if parts[14] in ("X.X",) else float(parts[14]),
        ))
    return events


def pick_travel_time(header, pick):
    """Seconds from the event's origin time to a pick, handling the pick's
    own hour:min:sec crossing midnight relative to the origin date (rare
    for a local network, but cheap to handle correctly): if the naive
    difference is very negative, assume the pick fell on the next day."""
    origin = header["hour"] * 3600 + header["minute"] * 60 + header["second"]
    pick_t = pick["hour"] * 3600 + pick["minute"] * 60 + pick["second"]
    tt = pick_t - origin
    if tt < -1800:  # more than 30 min "before" origin -> midnight rollover
        tt += 86400
    return tt


def event_datetime(header):
    return datetime.datetime(header["year"], header["month"], header["day"],
                              header["hour"], header["minute"],
                              int(header["second"]),
                              int((header["second"] % 1) * 1e6))
