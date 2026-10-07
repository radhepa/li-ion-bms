"""
Step 4 - why a pack needs balancing: 4 mismatched cells in series, 8 charge /
drive cycles, with and without passive balancing.

The cells start 13 % apart in SOC (e.g. after storage with different
self-discharge) and differ by a few % in capacity and resistance. Each cycle:
CC-CV charge (+30 min rest), then the real UDDS current profile until the
first cell reaches 2.8 V.

Output: results/pack_balancing.json, figures/06_balancing.png
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bms.data import load, resample                    # noqa: E402
from bms.model import CellModel                         # noqa: E402
from bms.pack import Pack, PackConfig                   # noqa: E402
from bms import plotstyle as ps                         # noqa: E402
from bms.plotstyle import plt                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
model = CellModel.load(ROOT / "data" / "cell_model.json")
udds = resample(load("udds"), 1.0)
profile = udds["i"][udds["t"] < 1400]                   # one UDDS schedule (~23 min), repeated
CYCLES = 8
out = {}
for balancing in (False, True):
    pack = Pack(model, PackConfig())
    rows = []
    for c in range(CYCLES):
        pack.charge(balancing)
        spread = 100 * (pack.soc.max() - pack.soc.min())
        ah = pack.discharge(profile)
        rows.append({"cycle": c + 1, "soc_spread_after_charge_pct": spread, "usable_ah": ah,
                     "soc_after_discharge": pack.soc.round(3).tolist()})
        print(f"  balancing={balancing!s:5}  cycle {c + 1}: spread after charge {spread:5.2f} %  "
              f"usable {ah:.3f} Ah", flush=True)
    out["with_balancing" if balancing else "no_balancing"] = rows
out["cells"] = {"capacity_ah": Pack(model).q.round(3).tolist(), "r_scale": Pack(model).rs.round(3).tolist()}
out["ideal_usable_ah"] = float(min(Pack(model).q))
(ROOT / "results" / "pack_balancing.json").write_text(json.dumps(out, indent=2) + "\n")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.0))
cyc = np.arange(1, CYCLES + 1)
for key, col, label in (("no_balancing", ps.C2, "no balancing"), ("with_balancing", ps.C1, "passive balancing, 45 mA")):
    r = out[key]
    ax[0].plot(cyc, [x["soc_spread_after_charge_pct"] for x in r], "o-", color=col, label=label)
    ax[1].plot(cyc, [x["usable_ah"] for x in r], "o-", color=col, label=label)
ax[1].axhline(out["ideal_usable_ah"], color=ps.INK_2, lw=1, ls=":")
ax[1].text(CYCLES, out["ideal_usable_ah"], "limit: smallest cell", ha="right", va="bottom", fontsize=8,
           color=ps.INK_2)
ax[0].set(title="Cell SOC spread at the end of each charge", xlabel="cycle", ylabel="max - min SOC (%)")
ax[1].set(title="Charge the pack can deliver (UDDS current)", xlabel="cycle", ylabel="Ah")
for a in ax:
    a.legend()
    a.set_xticks(cyc)
fig.tight_layout()
ps.save(fig, "06_balancing.png")
