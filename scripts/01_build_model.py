"""
Step 1 - build the cell model from lab data.

  OCV(SOC): average of a C/20 discharge and a C/20 charge. At 1/20 of the
            capacity per hour the cell is almost at rest, and averaging the two
            directions cancels the remaining IR drop and most of the hysteresis.
  R0, R1, tau1, R2, tau2 vs SOC: fitted to the HPPC test, which applies 10 s
            discharge and charge pulses (1C-6C) at SOC levels from 100 % down to 0 %.

Output: data/cell_model.json, figures/01_ocv.png, figures/02_hppc_fit.png
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bms.data import load                               # noqa: E402
from bms.model import CellModel                         # noqa: E402
from bms import plotstyle as ps                         # noqa: E402
from bms.plotstyle import plt                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

# ------------------------------------------------------------------ 1. OCV --
c20 = load("c20")
t, v, i = c20["t"], c20["v"], c20["i"]
dis = i > 0.05
chg = i < -0.05
q_dis = np.cumsum(np.where(dis, i, 0.0) * np.gradient(t)) / 3600
q_chg = np.cumsum(np.where(chg, -i, 0.0) * np.gradient(t)) / 3600
Q = float(q_dis[dis][-1])                     # capacity: C/20 discharge 4.2 V -> 2.8 V
soc_dis = 1 - q_dis[dis] / Q
soc_chg = (q_chg[chg] - q_chg[chg][0]) / (q_chg[chg][-1] - q_chg[chg][0])
grid = np.linspace(0, 1, 201)
v_dis = np.interp(grid, soc_dis[::-1], v[dis][::-1])
v_chg = np.interp(grid, soc_chg, v[chg])
ocv = 0.5 * (v_dis + v_chg)
ocv = np.maximum.accumulate(ocv)              # an OCV curve must rise with SOC
print(f"capacity Q = {Q:.3f} Ah; OCV {ocv[0]:.3f} V (0 %) .. {ocv[-1]:.3f} V (100 %); "
      f"mean charge/discharge gap {1e3 * np.mean(v_chg - v_dis):.1f} mV")

# ------------------------------------------------------------------ 2. HPPC --
hp = load("hppc")
th, vh, ih = hp["t"], hp["v"], hp["i"]
soc_h = 1 - np.concatenate([[0.0], np.cumsum(0.5 * (ih[1:] + ih[:-1]) * np.diff(th))]) / 3600 / Q
# A test block = the 4 discharge / 4 charge pulse pairs at one SOC. Blocks are
# separated by slow "SOC step" discharges (|I| < 3 A for more than 100 s).
moving = np.abs(ih) > 0.05
slow = moving & (np.abs(ih) < 2.0)
starts = []
k = 0
while k < len(ih):
    if slow[k]:
        j = k
        while j < len(ih) and slow[j]:
            j += 1
        if th[j - 1] - th[k] > 100:
            starts.append((k, j))
        k = j
    k += 1
edges = [0] + [s for s, _ in starts] + [len(ih)]
blocks = []
for a, b in zip(edges[:-1], edges[1:]):
    seg = slice(a, b)
    big = np.where(np.abs(ih[seg]) > 2.5)[0]
    if len(big) < 20:
        continue
    p0 = a + max(big[0] - 50, 0)              # 5 s before the first pulse
    p1 = b - 1                                # through the long rest that follows
    blocks.append((p0, p1))
print(f"HPPC: {len(blocks)} pulse blocks found")


def sim_block(p, sl, soc0):
    r0, r1, t1, r2, t2, dv = p
    tt, ii = th[sl], ih[sl]
    v1 = v2 = 0.0
    out = np.empty(len(tt))
    s = soc0
    for k in range(len(tt)):
        if k:
            dt = tt[k] - tt[k - 1]
            a1, a2 = np.exp(-dt / t1), np.exp(-dt / t2)
            v1 = a1 * v1 + r1 * (1 - a1) * ii[k - 1]
            v2 = a2 * v2 + r2 * (1 - a2) * ii[k - 1]
            s -= ii[k - 1] * dt / 3600 / Q
        out[k] = np.interp(s, grid, ocv) + dv - r0 * ii[k] - v1 - v2
    return out


rows = []
fits = []
for p0, p1 in blocks:
    # use the first ~15 minutes of the block: all pulses and most of the relaxation
    stop = min(p1, np.searchsorted(th, th[p0] + 900))
    sl = slice(p0, stop)
    soc0 = soc_h[p0]
    # decimate the long rests to keep the fit fast: keep every sample while current flows
    keep = np.unique(np.concatenate([np.where(np.abs(ih[sl]) > 0.05)[0],
                                     np.arange(0, stop - p0, 10)]))

    def res(p):
        return (sim_block(p, sl, soc0) - vh[sl])[keep]

    fit = least_squares(res, [0.02, 0.01, 3.0, 0.01, 60.0, 0.0],
                        bounds=([1e-3, 1e-4, 0.3, 1e-4, 15.0, -0.05], [0.1, 0.1, 15.0, 0.2, 900.0, 0.05]),
                        x_scale=[0.01, 0.01, 1.0, 0.01, 50.0, 0.01])
    rmse = float(np.sqrt(np.mean(fit.fun ** 2)))
    rows.append([soc0, *fit.x, rmse])
    fits.append((sl, soc0, fit.x))
    print(f"  SOC {100 * soc0:5.1f} %: R0 {1e3 * fit.x[0]:5.1f} mOhm  R1 {1e3 * fit.x[1]:5.1f} mOhm  "
          f"tau1 {fit.x[2]:5.1f} s  R2 {1e3 * fit.x[3]:5.1f} mOhm  tau2 {fit.x[4]:6.1f} s  "
          f"OCV offset {1e3 * fit.x[5]:+5.1f} mV  fit RMSE {1e3 * rmse:4.1f} mV", flush=True)

rows = np.array(sorted(rows))
# Near empty the diffusion response no longer fits a 2-RC shape and the fit
# runs into its parameter bounds (R2 = 200 mOhm). Those blocks are not trusted:
# the tables hold the nearest valid fit instead.
valid = (rows[:, 4] < 0.99 * 0.2) & (rows[:, 5] < 0.99 * 900.0)
print(f"valid HPPC fits: {valid.sum()} of {len(rows)} "
      f"(dropped SOC: {', '.join(f'{100 * s:.0f} %' for s in rows[~valid, 0])})")
all_rows, rows = rows, rows[valid]
soc_p = rows[:, 0]


def table(col):
    """Parameter vs SOC on the model grid (held constant beyond the trusted range)."""
    return np.interp(grid, soc_p, rows[:, col])


model = CellModel(q_ah=Q, soc=grid, ocv=ocv, r0=table(1), r1=table(2), tau1=table(3),
                  r2=table(4), tau2=table(5),
                  meta={"source": "LG HG2 cell, McMaster University (Kollmeyer 2020), CC BY 4.0",
                        "hppc_fit_rmse_mV": float(1e3 * rows[:, 7].mean()),
                        "ocv_offsets_mV": (1e3 * rows[:, 6]).round(1).tolist(),
                        "hppc_soc_points": soc_p.round(4).tolist(),
                        "hppc_dropped_soc_points": all_rows[~valid, 0].round(4).tolist()})
model.save(ROOT / "data" / "cell_model.json")
print(f"mean HPPC fit RMSE {1e3 * rows[:, 7].mean():.1f} mV; saved data/cell_model.json")

# ---------------------------------------------------------------- figures --
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
a = ax[0]
a.plot(100 * grid, v_dis, color=ps.C2, lw=1.4, label="C/20 discharge")
a.plot(100 * grid, v_chg, color=ps.C3, lw=1.4, label="C/20 charge")
a.plot(100 * grid, ocv, color=ps.C1, label="OCV (average)")
a.set(title=f"Open-circuit voltage, LG HG2 (Q = {Q:.2f} Ah)", xlabel="state of charge (%)", ylabel="V")
a.legend()
a = ax[1]
a.plot(100 * grid, 1e3 * model.docv_at(grid), color=ps.C1)
a.set(title="OCV slope: how much voltage tells you about SOC", xlabel="state of charge (%)",
      ylabel="dOCV/dSOC (mV per 1 % SOC)")
a.set_ylim(0, None)
ax[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y / 100:.0f}"))
fig.tight_layout()
ps.save(fig, "01_ocv.png")

fig, ax = plt.subplots(1, 3, figsize=(14, 4.0))
a = ax[0]
a.plot(100 * grid, 1e3 * model.r0, "-", color=ps.C1, lw=1.2)
a.plot(100 * grid, 1e3 * model.r1, "-", color=ps.C2, lw=1.2)
a.plot(100 * grid, 1e3 * model.r2, "-", color=ps.C3, lw=1.2)
a.plot(100 * soc_p, 1e3 * rows[:, 1], "o", color=ps.C1, ms=5, label="R0 (instant)")
a.plot(100 * soc_p, 1e3 * rows[:, 2], "s", color=ps.C2, ms=5, label="R1 (fast RC)")
a.plot(100 * soc_p, 1e3 * rows[:, 4], "^", color=ps.C3, ms=5, label="R2 (slow RC)")
a.set(title="Fitted resistances (markers) and tables (lines)", xlabel="SOC (%)", ylabel="mΩ")
a.legend()
a = ax[1]
a.semilogy(100 * grid, model.tau1, "-", color=ps.C2, lw=1.2)
a.semilogy(100 * grid, model.tau2, "-", color=ps.C3, lw=1.2)
a.semilogy(100 * soc_p, rows[:, 3], "s", color=ps.C2, ms=5, label="τ1")
a.semilogy(100 * soc_p, rows[:, 5], "^", color=ps.C3, ms=5, label="τ2")
a.set(title="Time constants", xlabel="SOC (%)", ylabel="seconds")
a.legend()
a = ax[2]
sl, soc0, p = fits[len(fits) // 2]
vm = sim_block(p, sl, soc0)
tt = th[sl] - th[sl][0]
a.plot(tt, vh[sl], color=ps.INK_2, lw=3.5, alpha=0.35, label="measured")
a.plot(tt, vm, color=ps.C1, lw=1.3, label="2-RC model")
a.set(title=f"HPPC pulses at {100 * soc0:.0f} % SOC", xlabel="time (s)", ylabel="V", xlim=(0, 300))
a.legend()
fig.tight_layout()
ps.save(fig, "02_hppc_fit.png")
