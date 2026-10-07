"""
Cell model: open-circuit voltage curve + 2-RC Thevenin equivalent circuit.

            R0        R1             R2
    +---/\/\/\---+--/\/\/\--+---+--/\/\/\--+---+
    |            |          |   |          |   |
  OCV(SOC)       +---||-----+   +---||-----+   v (terminal)
    |                C1             C2
    -
    v = OCV(SOC) - R0*i - V1 - V2          (i > 0 = discharge)
    dV1/dt = -V1/(R1 C1) + i/C1            tau1 = R1*C1  (seconds: charge transfer)
    dV2/dt = -V2/(R2 C2) + i/C2            tau2 = R2*C2  (minutes: diffusion)
    dSOC/dt = -i / (3600 * Q)

All parameters depend on SOC and are stored as tables on a common SOC grid,
which is exactly how they are stored in the firmware (cell_tables.h).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class CellModel:
    q_ah: float                       # capacity, Ah
    soc: np.ndarray                   # grid, 0..1
    ocv: np.ndarray                   # V
    r0: np.ndarray                    # ohm
    r1: np.ndarray
    tau1: np.ndarray                  # s
    r2: np.ndarray
    tau2: np.ndarray
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------- lookups --
    def _lut(self, arr, soc):
        return np.interp(np.clip(soc, 0.0, 1.0), self.soc, arr)

    def ocv_at(self, soc):
        """OCV, extended linearly past 0 % and 100 % with the end slopes, so an
        over-charged or over-discharged cell still moves its voltage."""
        soc = np.asarray(soc, float)
        v = self._lut(self.ocv, soc)
        top = (self.ocv[-1] - self.ocv[-3]) / (self.soc[-1] - self.soc[-3])
        bot = (self.ocv[2] - self.ocv[0]) / (self.soc[2] - self.soc[0])
        v = np.where(soc > 1.0, self.ocv[-1] + top * (soc - 1.0), v)
        return np.where(soc < 0.0, self.ocv[0] + bot * soc, v)

    def docv_at(self, soc, h=0.005):
        return (self.ocv_at(soc + h) - self.ocv_at(soc - h)) / (2 * h)

    def params_at(self, soc):
        return (self._lut(self.r0, soc), self._lut(self.r1, soc), self._lut(self.tau1, soc),
                self._lut(self.r2, soc), self._lut(self.tau2, soc))

    # ---------------------------------------------------------- simulation --
    def simulate(self, t, i, soc0=1.0, v1=0.0, v2=0.0):
        """Terminal voltage for a current profile (zero-order hold between samples)."""
        n = len(t)
        v = np.empty(n)
        soc = np.empty(n)
        s = soc0
        for k in range(n):
            if k > 0:
                dt = t[k] - t[k - 1]
                ik = i[k - 1]
                r0, r1, t1, r2, t2 = self.params_at(s)
                a1, a2 = np.exp(-dt / t1), np.exp(-dt / t2)
                v1 = a1 * v1 + r1 * (1 - a1) * ik
                v2 = a2 * v2 + r2 * (1 - a2) * ik
                s -= ik * dt / (3600.0 * self.q_ah)
            r0 = self._lut(self.r0, s)
            soc[k] = s
            v[k] = self.ocv_at(s) - r0 * i[k] - v1 - v2
        return v, soc

    # ------------------------------------------------------------------ io --
    def save(self, path):
        d = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in self.__dict__.items()}
        Path(path).write_text(json.dumps(d, indent=1) + "\n")

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text())
        for k in ("soc", "ocv", "r0", "r1", "tau1", "r2", "tau2"):
            d[k] = np.array(d[k])
        return cls(**d)
