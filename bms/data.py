"""
Loading the LG 18650HG2 test data (McMaster University, Kollmeyer 2020, CC BY 4.0).

Every test is returned with the same columns:
    t   time (s), starting at 0
    v   terminal voltage (V)
    i   current (A), POSITIVE = DISCHARGE (the raw files use negative = discharge)
    ah  charge removed since the start (Ah, positive = discharged)
    temp cell surface temperature (C)
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

TESTS = {
    "hppc": "HPPC_549.mat",
    "c20": "C20_549.mat",
    "cap1c": "Cap1C_551.mat",
    "udds": "UDDS_551.csv",
    "la92": "LA92_551.mat",
    "us06": "US06_551.mat",
    "mixed1": "Mixed1_551.mat",
}


def _from_mat(path: Path) -> dict:
    m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)["meas"]
    t = np.asarray(m.Time, float)
    i = -np.asarray(m.Current, float)
    return {"t": t - t[0], "v": np.asarray(m.Voltage, float), "i": i,
            "ah": -np.asarray(m.Ah, float), "temp": np.asarray(m.Battery_Temp_degC, float)}


def _hms(s: str) -> float:
    h, mnt, sec = s.split(":")
    return 3600 * float(h) + 60 * float(mnt) + float(sec)


def _from_csv(path: Path) -> dict:
    rows = []
    text = path.read_bytes().replace(b"\x00", b"").decode("latin-1")   # files contain stray NULs
    with io.StringIO(text, newline="") as f:
        reader = csv.reader(f)
        header = None
        for r in reader:
            if header is None:
                if r and r[0] == "Time Stamp":
                    header = r
                    next(reader)                      # units row
                continue
            if len(r) > 11 and r[8]:
                rows.append((_hms(r[3]), float(r[8]), float(r[9]), float(r[11]), float(r[10])))
    a = np.array(rows)
    t = a[:, 0] - a[0, 0]
    return {"t": t, "v": a[:, 1], "i": -a[:, 2], "ah": -a[:, 3] + a[0, 3], "temp": a[:, 4]}


def load(name: str) -> dict:
    """Load one test by short name (see TESTS)."""
    path = RAW / TESTS[name]
    if not path.exists():
        raise FileNotFoundError(f"{path} missing - run: python data/download.py")
    d = _from_mat(path) if path.suffix == ".mat" else _from_csv(path)
    d["name"] = name
    return d


def resample(d: dict, dt: float) -> dict:
    """Uniform time grid (BMS firmware samples at a fixed rate). Current is
    averaged over each interval so the charge is conserved."""
    t_new = np.arange(0.0, d["t"][-1], dt)
    q = np.concatenate([[0.0], np.cumsum(0.5 * (d["i"][1:] + d["i"][:-1]) * np.diff(d["t"]))])
    q_new = np.interp(np.append(t_new, t_new[-1] + dt), d["t"], q)
    out = {"t": t_new, "i": np.diff(q_new) / dt, "v": np.interp(t_new, d["t"], d["v"]),
           "temp": np.interp(t_new, d["t"], d["temp"]), "name": d.get("name", "")}
    out["ah"] = np.interp(t_new, d["t"], d["ah"])
    return out
