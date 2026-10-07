"""
Step 6 - software-in-the-loop (SIL): run the actual C firmware on a PC.

  A. Equivalence: the C EKF (single precision, firmware tables) and the
     Python reference EKF see the same real UDDS sensor data. They must agree.
     (Open loop: in the recorded test, regen pulses up to 5.3 A did flow, even
     though the firmware's 4.5 A charge limit would have blocked them.)
  B. Fault injection: a scripted test with an over-current, an over-voltage
     during charge and an over-temperature. The firmware's closed-loop plant
     (sil_main.c) stops the current when a FET opens. Every trip and recovery
     time is checked against the requirement in bms_config.h.

Output: results/sil.json, figures/07_sil_faults.png
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bms.data import load, resample                    # noqa: E402
from bms.model import CellModel                         # noqa: E402
from bms.ekf import SocEKF, EKFTuning                   # noqa: E402
from bms import plotstyle as ps                         # noqa: E402
from bms.plotstyle import plt                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FW = ROOT / "firmware"
EXE = FW / ("bms_sil.exe" if os.name == "nt" else "bms_sil")
make = "mingw32-make" if os.name == "nt" else "make"
subprocess.run([make, "-C", str(FW), "all"], check=True, capture_output=True)


def sil(rows, soc0, open_loop=False):
    with tempfile.TemporaryDirectory() as d:
        fin, fout = Path(d) / "in.csv", Path(d) / "out.csv"
        np.savetxt(fin, rows, delimiter=",", fmt="%.6f",
                   header="t_s,i_request_a,v1,v2,v3,v4,temp1_c,temp2_c", comments="")
        args = [str(EXE), str(fin), str(fout), str(soc0)] + (["--open-loop"] if open_loop else [])
        subprocess.run(args, check=True, capture_output=True)
        out = np.genfromtxt(fout, delimiter=",", names=True)
    return out


results = {}

# ------------------------------------------------------------ A. equivalence --
rng = np.random.default_rng(11)
d = resample(load("udds"), 1.0)
t = d["t"]
# what the firmware measures at sample k: the average current over the previous period
i_meas = np.concatenate([[0.0], d["i"][:-1]]) + 0.050 + rng.normal(0, 0.020, len(t))
volts = np.stack([d["v"] + rng.normal(0, 0.002, len(t)) for _ in range(4)], axis=1)
rows = np.column_stack([t, i_meas, volts, np.full(len(t), 25.0), np.full(len(t), 25.0)])
out = sil(rows, 0.5, open_loop=True)   # recorded data: the current did flow

fw_model = CellModel.load(ROOT / "data" / "cell_model_firmware.json")
ref = SocEKF(fw_model, EKFTuning(soc0=0.5, sigma_soc0=0.3, sigma_i=0.05, sigma_v1=1e-4, sigma_v=0.02))
soc_ref = np.empty(len(t))
for k in range(len(t)):
    dt = t[k] - t[k - 1] if k else 0.0
    soc_ref[k] = ref.step(i_meas[k], volts[k, 0], dt)
truth = 1 - d["ah"] / fw_model.q_ah
diff = np.abs(out["soc1"] - soc_ref)
s = t >= 600
results["equivalence"] = {
    "samples": int(len(t)),
    "max_abs_diff_C_vs_python_pct": float(100 * diff.max()),
    "c_firmware_rmse_vs_truth_pct": float(100 * np.sqrt(np.mean((out["soc1"][s] - truth[s]) ** 2))),
    "python_rmse_vs_truth_pct": float(100 * np.sqrt(np.mean((soc_ref[s] - truth[s]) ** 2))),
}
print("A. C firmware vs Python reference:", {k: round(v, 4) for k, v in results["equivalence"].items()})

# -------------------------------------------------------- B. fault injection --
dt = 0.1
tt = np.arange(0, 200, dt)
ireq = np.zeros_like(tt)
ireq[(tt >= 10) & (tt < 30)] = 10.0
ireq[(tt >= 30) & (tt < 33)] = 30.0          # over-current (limit 22 A for 0.5 s)
ireq[(tt >= 33) & (tt < 70)] = 8.0
ireq[(tt >= 80) & (tt < 125)] = -2.0         # charging
ireq[(tt >= 130) & (tt < 200)] = 5.0
v = np.full((len(tt), 4), 3.90)
ramp = (tt >= 95) & (tt < 110)
v[ramp, 2] = np.interp(tt[ramp], [95, 102, 110], [4.15, 4.27, 4.27])   # cell 3 over-charges
v[(tt >= 110) & (tt < 125), 2] = 4.12                                # ...and relaxes
temp = np.full(len(tt), 25.0)
temp = np.where(tt >= 130, np.interp(tt, [130, 165, 190, 200], [25, 66, 45, 40]), temp)
rows = np.column_stack([tt, ireq, v, temp, np.full(len(tt), 25.0)])
o = sil(rows, 0.6)


def first(cond, after=0.0):
    idx = np.where(cond & (o["t_s"] >= after))[0]
    return float(o["t_s"][idx[0]]) if len(idx) else None


t_ocd = first(o["dsg_fet"] == 0, 29)
t_ocd_retry = first(o["dsg_fet"] == 1, t_ocd + 1) if t_ocd else None
t_ov_cross = float(np.interp(4.25, v[ramp, 2][tt[ramp] < 102], tt[ramp][tt[ramp] < 102]))
t_ov = first(o["chg_fet"] == 0, 95)
t_ov_clear = first(o["chg_fet"] == 1, t_ov + 0.1) if t_ov else None
t_ot_cross = float(np.interp(60.0, [25, 66], [130, 165]))
t_ot = first(o["dsg_fet"] == 0, 130)
t_ot_clear_cross = 165 + (66 - 50) / (66 - 45) * (190 - 165)      # T falls through 50 C
t_ot_clear = first(o["dsg_fet"] == 1, t_ot + 0.1) if t_ot else None
checks = [
    ("discharge over-current trips 0.5 s after I > 22 A", 30.0 + 0.5, t_ocd),
    ("...and retries 30 s later", (t_ocd or 0) + 30.0, t_ocd_retry),
    ("cell over-voltage trips 1 s after V > 4.25 V", t_ov_cross + 1.0, t_ov),
    ("...clears 1 s after V < 4.15 V", 110.0 + 1.0, t_ov_clear),
    ("over-temperature trips 2 s after T > 60 C", t_ot_cross + 2.0, t_ot),
    ("...clears 2 s after T < 50 C", t_ot_clear_cross + 2.0, t_ot_clear),
]
ver = []
for name, expect, got in checks:
    ok = got is not None and abs(got - expect) <= 0.15
    ver.append({"requirement": name, "expected_s": round(expect, 2),
                "measured_s": None if got is None else round(got, 2), "pass": bool(ok)})
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: expected {expect:.2f} s, firmware {got}")
results["fault_injection"] = ver
(ROOT / "results" / "sil.json").write_text(json.dumps(results, indent=2) + "\n")

fig, ax = plt.subplots(4, 1, figsize=(12, 8.4), sharex=True,
                       gridspec_kw={"height_ratios": [1.4, 1, 1, 0.9]})
ax[0].plot(tt, ireq, color=ps.INK_2, lw=3, alpha=0.35, label="current the load / charger asks for")
ax[0].plot(o["t_s"], o["i_actual_a"], color=ps.C1, lw=1.4, label="current that actually flows")
ax[0].set(ylabel="A", title="Software-in-the-loop: the real C firmware decides when the FETs open")
ax[0].legend(loc="upper right")
ax[1].plot(tt, v.max(axis=1), color=ps.C2, lw=1.4)
ax[1].axhline(4.25, color=ps.INK_2, lw=1, ls=":")
ax[1].set(ylabel="highest cell (V)")
ax[2].plot(tt, temp, color=ps.C8, lw=1.4)
ax[2].axhline(60, color=ps.INK_2, lw=1, ls=":")
ax[2].set(ylabel="temperature (C)")
ax[3].step(o["t_s"], o["dsg_fet"] + 1.2, where="post", color=ps.C1, lw=1.6, label="discharge FET")
ax[3].step(o["t_s"], o["chg_fet"], where="post", color=ps.C3, lw=1.6, label="charge FET")
ax[3].set(yticks=[0, 1, 1.2, 2.2], yticklabels=["off", "on", "off", "on"], xlabel="time (s)")
ax[3].legend(loc="center right", ncol=2)
for x, lab in ((30, "over-current"), (100, "over-voltage"), (155, "over-temperature")):
    ax[0].annotate(lab, (x, ax[0].get_ylim()[1] * 0.85), ha="center", fontsize=8.5, color=ps.INK_2)
fig.tight_layout()
ps.save(fig, "07_sil_faults.png")
