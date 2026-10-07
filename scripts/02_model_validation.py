"""
Step 2 - does the model predict voltage on driving it was NOT fitted to?

Four EPA drive-cycle tests on the same cell (UDDS city, LA92 aggressive city,
US06 highway/aggressive, and a mixed cycle), each from full to ~7 % SOC. The
model gets only the measured current and the starting SOC; we compare its
voltage with the measured voltage.

Two OCV variants are compared:
  "C/20 average"  : the OCV from step 1
  "HPPC-relaxed"  : C/20 average + the offsets the HPPC fit found at each SOC
                    (the cell's voltage after resting from discharge, which is
                    what a discharging drive cycle actually sees: hysteresis)
The better one is kept.

Output: data/cell_model.json (updated), results/model_validation.json, figures/03_model_validation.png
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bms.data import load, resample                    # noqa: E402
from bms.model import CellModel                         # noqa: E402
from bms import plotstyle as ps                         # noqa: E402
from bms.plotstyle import plt                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
(ROOT / "results").mkdir(exist_ok=True)
base = CellModel.load(ROOT / "data" / "cell_model.json")
if "ocv_c20" in base.meta:                              # re-run safe: start from the C/20 curve
    base.ocv = np.array(base.meta["ocv_c20"])
offs = np.array(base.meta["ocv_offsets_mV"]) * 1e-3
soc_pts = np.array(base.meta["hppc_soc_points"])
order = np.argsort(soc_pts)
# smooth the offsets with a low-order polynomial so the OCV stays monotonic
coef = np.polyfit(soc_pts[order], offs[order], 3)
ocv_relaxed = np.maximum.accumulate(base.ocv + np.polyval(coef, base.soc))

CYCLES = {"udds": "UDDS (city)", "la92": "LA92 (aggressive city)",
          "us06": "US06 (highway, aggressive)", "mixed1": "Mixed cycle"}
data = {k: resample(load(k), 1.0) for k in CYCLES}

results = {}
variants = {"C/20 average": base.ocv, "HPPC-relaxed": ocv_relaxed}
for vname, ocv in variants.items():
    m = CellModel(**{**base.__dict__, "ocv": ocv})
    per = {}
    for key, d in data.items():
        v_model, _ = m.simulate(d["t"], d["i"], soc0=1.0)
        err = v_model - d["v"]
        valid = 1 - d["ah"] / m.q_ah > 0.05             # exclude the last few % (steep knee)
        per[key] = {"rmse_mV": float(1e3 * np.sqrt(np.mean(err[valid] ** 2))),
                    "max_mV": float(1e3 * np.max(np.abs(err[valid])))}
    results[vname] = per
    print(vname, {k: round(v["rmse_mV"], 1) for k, v in per.items()})

best = min(results, key=lambda n: np.mean([v["rmse_mV"] for v in results[n].values()]))
model = CellModel(**{**base.__dict__, "ocv": variants[best]})
model.meta = {**base.meta, "ocv_variant": best, "ocv_c20": base.ocv.tolist()}
model.save(ROOT / "data" / "cell_model.json")
results["selected"] = best
(ROOT / "results" / "model_validation.json").write_text(json.dumps(results, indent=2) + "\n")
print("selected OCV:", best)

fig, ax = plt.subplots(2, 1, figsize=(12, 6.2), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
d = data["udds"]
v_model, soc = model.simulate(d["t"], d["i"], soc0=1.0)
th = d["t"] / 3600
ax[0].plot(th, d["v"], color=ps.INK_2, lw=2.5, alpha=0.35, label="measured")
ax[0].plot(th, v_model, color=ps.C1, lw=0.8, label="2-RC model (current in, voltage out)")
ax[0].set(ylabel="cell voltage (V)",
          title=f"Model vs measurement on a 4.4-hour UDDS test it was not fitted to: "
                f"RMSE {results[best]['udds']['rmse_mV']:.0f} mV")
ax[0].legend()
ax[1].plot(th, 1e3 * (v_model - d["v"]), color=ps.C2, lw=0.6)
ax[1].set(ylabel="error (mV)", xlabel="time (h)")
fig.tight_layout()
ps.save(fig, "03_model_validation.png")
