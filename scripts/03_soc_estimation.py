"""
Step 3 - state-of-charge estimation on real drive-cycle data.

Ground truth: lab-grade Coulomb counting (0.1 % current accuracy) from a known
full charge. The BMS, however, sees realistic sensors:
    current : +50 mA offset (a cheap shunt amplifier) + 20 mA noise
    voltage : 2 mV noise, quantized to the BQ76920 ADC step (0.38 mV)
and it starts with NO idea of the SOC (guess 50 %, truth 100 %), as after a
long storage or a controller reset.

Compared:
    EKF                        cell model + voltage feedback
    Coulomb counting (CC)      same wrong start and same biased current
    CC, perfect start          the best case for CC: only the sensor offset hurts

Output: results/soc_estimation.json, figures/04_soc_udds.png, figures/05_soc_summary.png
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bms.data import load, resample                    # noqa: E402
from bms.model import CellModel                         # noqa: E402
from bms.ekf import SocEKF, EKFTuning, coulomb_count    # noqa: E402
from bms import plotstyle as ps                         # noqa: E402
from bms.plotstyle import plt                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
model = CellModel.load(ROOT / "data" / "cell_model.json")
rng = np.random.default_rng(7)
I_BIAS, I_NOISE, V_NOISE, V_LSB = 0.050, 0.020, 0.002, 382e-6
TUNE = EKFTuning(soc0=0.5, sigma_soc0=0.3, sigma_i=0.05, sigma_v1=1e-4, sigma_v=0.02)
CYCLES = {"udds": "UDDS", "la92": "LA92", "us06": "US06", "mixed1": "Mixed"}
SKIP = 600  # s: score after the first 10 minutes (the EKF's convergence window)

results, traces = {}, {}
for key, label in CYCLES.items():
    d = resample(load(key), 1.0)
    t = d["t"]
    truth = 1 - d["ah"] / model.q_ah
    i_meas = d["i"] + I_BIAS + rng.normal(0, I_NOISE, len(t))
    v_meas = np.round((d["v"] + rng.normal(0, V_NOISE, len(t))) / V_LSB) * V_LSB
    ekf = SocEKF(model, TUNE).run(t, i_meas, v_meas)
    cc_bad = coulomb_count(t, i_meas, TUNE.soc0, model.q_ah)
    cc_good = coulomb_count(t, i_meas, 1.0, model.q_ah)
    s = t >= SKIP
    r = {}
    for name, est in (("EKF", ekf), ("CC, wrong start", cc_bad), ("CC, perfect start", cc_good)):
        e = 100 * (est - truth)
        within = np.where(np.abs(e) < 2.0)[0]
        r[name] = {"rmse_pct": float(np.sqrt(np.mean(e[s] ** 2))), "max_pct": float(np.max(np.abs(e[s]))),
                   "final_pct": float(e[-1]),
                   "time_to_2pct_s": float(t[within[0]]) if len(within) else None}
    results[label] = r
    traces[label] = (t, truth, ekf, cc_bad, cc_good)
    print(f"{label:6s} EKF RMSE {r['EKF']['rmse_pct']:.2f} % (max {r['EKF']['max_pct']:.2f} %, "
          f"within 2 % after {r['EKF']['time_to_2pct_s']:.0f} s) | CC wrong start "
          f"{r['CC, wrong start']['rmse_pct']:.1f} % | CC perfect start {r['CC, perfect start']['rmse_pct']:.2f} % "
          f"(final {r['CC, perfect start']['final_pct']:+.2f} %)", flush=True)

# Harder: switch the BMS on in the MIDDLE of a drive, in the flat part of the
# OCV curve, with a guess that is 40 % too low.
mid = {}
for key, label in CYCLES.items():
    d = resample(load(key), 1.0)
    k0 = int(np.searchsorted(d["t"], 3600.0 if d["t"][-1] > 5400 else 1200.0))
    t = d["t"][k0:] - d["t"][k0]
    truth = 1 - d["ah"][k0:] / model.q_ah
    i_meas = d["i"][k0:] + I_BIAS + rng.normal(0, I_NOISE, len(t))
    v_meas = np.round((d["v"][k0:] + rng.normal(0, V_NOISE, len(t))) / V_LSB) * V_LSB
    tune = EKFTuning(**{**TUNE.__dict__, "soc0": float(truth[0] - 0.4)})
    est = SocEKF(model, tune).run(t, i_meas, v_meas)
    err = 100 * np.abs(est - truth)
    settled = np.where(np.convolve(err < 2.0, np.ones(60), "valid") == 60)[0]   # inside 2 % for 1 min
    mid[label] = {"true_soc_at_start_pct": float(100 * truth[0]), "guess_pct": float(100 * tune.soc0),
                  "time_to_within_2pct_s": float(t[settled[0]]) if len(settled) else None,
                  "rmse_after_10min_pct": float(np.sqrt(np.mean(err[t >= 600] ** 2)))}
    print(f"{label:6s} mid-drive start: true {100 * truth[0]:.0f} %, guess {100 * tune.soc0:.0f} % -> "
          f"within 2 % after {mid[label]['time_to_within_2pct_s']} s, RMSE {mid[label]['rmse_after_10min_pct']:.2f} %")
results["mid_drive_start"] = mid

results["sensor_model"] = {"current_bias_A": I_BIAS, "current_noise_A": I_NOISE,
                           "voltage_noise_V": V_NOISE, "voltage_lsb_V": V_LSB}
results["tuning"] = TUNE.__dict__
(ROOT / "results" / "soc_estimation.json").write_text(json.dumps(results, indent=2) + "\n")

t, truth, ekf, cc_bad, cc_good = traces["UDDS"]
th = t / 3600
fig, ax = plt.subplots(2, 1, figsize=(12, 6.4), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
ax[0].plot(th, 100 * truth, color=ps.INK_2, lw=3.5, alpha=0.35, label="true SOC (lab)")
ax[0].plot(th, 100 * cc_bad, color=ps.C2, lw=1.4, label="Coulomb counting, wrong start")
ax[0].plot(th, 100 * cc_good, color=ps.C4, lw=1.4, label="Coulomb counting, perfect start")
ax[0].plot(th, 100 * ekf, color=ps.C1, lw=1.4, label="EKF (starts at 50 %)")
ax[0].set(ylabel="SOC (%)", title="UDDS city cycle, real cell data, realistic sensor errors")
ax[0].legend(ncol=2)
ax[1].plot(th, 100 * (cc_good - truth), color=ps.C4, lw=1.2, label="CC, perfect start")
ax[1].plot(th, 100 * (ekf - truth), color=ps.C1, lw=1.2, label="EKF")
ax[1].axhline(0, color=ps.INK_2, lw=0.8)
ax[1].set(ylabel="error (% SOC)", xlabel="time (h)", ylim=(-4, 4))
ax[1].legend(ncol=2)
fig.tight_layout()
ps.save(fig, "04_soc_udds.png")

fig, a = plt.subplots(figsize=(8, 3.8))
labels = list(CYCLES.values())
x = np.arange(len(labels))
w = 0.38
for k, (name, col) in enumerate((("CC, perfect start", ps.C4), ("EKF", ps.C1))):
    vals = [results[c][name]["rmse_pct"] for c in labels]
    bars = a.bar(x + (k - 0.5) * w, vals, w, color=col, edgecolor=ps.SURFACE, linewidth=2, label=name)
    a.bar_label(bars, fmt="%.2f", fontsize=8.5, color=ps.INK, padding=2)
a.set(xticks=x, xticklabels=labels, ylabel="SOC RMSE (%)",
      title="SOC error after 10 min: EKF starting blind beats a perfectly initialized counter")
a.grid(axis="x", visible=False)
a.legend()
fig.tight_layout()
ps.save(fig, "05_soc_summary.png")
