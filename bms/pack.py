"""
4S pack simulation: four cells in series that are NOT identical, and passive
balancing.

In a series string every cell carries the same current, so the pack must stop
discharging when the WEAKEST/emptiest cell hits its minimum voltage and stop
charging when the FULLEST cell hits its maximum. Any SOC mismatch therefore
wastes capacity at both ends. Passive balancing bleeds a little charge from
the higher cells through a resistor so the string lines up again.

Balancing rule (also in firmware/src/balancing.c):
    a cell bleeds when it is more than 10 mV above the lowest cell and above
    3.9 V, only while charging or resting, and never two neighbours at once
    (a BQ76920 rule) - odd and even cells take turns every minute.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import CellModel


@dataclass
class PackConfig:
    n: int = 4
    q_spread: float = 0.03          # 1-sigma capacity spread
    r_spread: float = 0.08          # 1-sigma resistance spread
    soc0: tuple = (0.97, 0.90, 0.84, 0.93)
    i_bal: float = 0.045            # A, bleed current per cell: 4.1 V / (2 x 43 + 4) ohm (BQ76920)
    v_cut: float = 2.8              # V, discharge cutoff per cell
    v_max: float = 4.2              # V, charge limit per cell
    i_chg: float = 1.5              # A, constant-current charge
    i_end: float = 0.10             # A, end of the constant-voltage phase
    bal_window_mv: float = 10.0
    bal_min_v: float = 3.9
    rest_s: float = 1800.0          # rest after charging (balancing continues)
    dt: float = 2.0
    seed: int = 3


class Pack:
    def __init__(self, model: CellModel, cfg: PackConfig = PackConfig()):
        self.m, self.cfg = model, cfg
        rng = np.random.default_rng(cfg.seed)
        self.q = model.q_ah * (1 + cfg.q_spread * rng.standard_normal(cfg.n))
        self.rs = 1 + cfg.r_spread * rng.standard_normal(cfg.n)
        self.soc = np.array(cfg.soc0[: cfg.n], float)
        self.v1 = np.zeros(cfg.n)
        self.v2 = np.zeros(cfg.n)
        self.t = 0.0
        self.bleed_ah = np.zeros(cfg.n)

    def voltages(self, i):
        r0 = np.array([self.m.params_at(s)[0] for s in self.soc]) * self.rs
        return self.m.ocv_at(self.soc) - r0 * i - self.v1 - self.v2

    def step(self, i, balance):
        """Advance dt with string current i (+ = discharge); balance = cells to bleed."""
        dt = self.cfg.dt
        for k in range(self.cfg.n):
            _, r1, t1, r2, t2 = self.m.params_at(self.soc[k])
            ik = i + (self.cfg.i_bal if balance[k] else 0.0)        # bleed adds to the cell's discharge
            a1, a2 = np.exp(-dt / t1), np.exp(-dt / t2)
            self.v1[k] = a1 * self.v1[k] + r1 * self.rs[k] * (1 - a1) * ik
            self.v2[k] = a2 * self.v2[k] + r2 * self.rs[k] * (1 - a2) * ik
            self.soc[k] -= ik * dt / 3600 / self.q[k]
            if balance[k]:
                self.bleed_ah[k] += self.cfg.i_bal * dt / 3600
        self.t += dt

    def balance_mask(self, v, enabled):
        c = self.cfg
        if not enabled:
            return np.zeros(c.n, bool)
        want = (v > v.min() + c.bal_window_mv * 1e-3) & (v > c.bal_min_v)
        phase = int(self.t // 60) % 2                    # odd/even cells take turns
        allowed = (np.arange(c.n) % 2) == phase
        return want & allowed

    def discharge(self, current_profile, dt_profile=1.0):
        """Repeat a drive-cycle current profile until the first cell hits v_cut.
        Returns Ah delivered by the pack."""
        ah, k = 0.0, 0
        n = len(current_profile)
        while True:
            i = current_profile[int(k * self.cfg.dt / dt_profile) % n]
            v = self.voltages(i)
            if v.min() <= self.cfg.v_cut:
                return ah
            self.step(i, np.zeros(self.cfg.n, bool))
            ah += i * self.cfg.dt / 3600
            k += 1

    def charge(self, balancing: bool):
        """CC-CV on the string, limited by the highest cell; then rest. Balancing
        (if enabled) runs during both. Returns Ah put in."""
        c = self.cfg
        ah = 0.0
        while True:
            # constant current, unless that would push any cell above v_max:
            # then the largest current that holds the highest cell AT v_max (CV)
            r0 = np.array([self.m.params_at(s)[0] for s in self.soc]) * self.rs
            emf = self.m.ocv_at(self.soc) - self.v1 - self.v2
            i_cv = np.min((c.v_max - emf) / r0)
            i = -min(c.i_chg, max(i_cv, 0.0))
            if -i < c.i_end:
                break
            v = self.voltages(i)
            self.step(i, self.balance_mask(v, balancing))
            ah -= i * c.dt / 3600
        t_end = self.t + c.rest_s
        while self.t < t_end:
            v = self.voltages(0.0)
            self.step(0.0, self.balance_mask(v, balancing))
        return ah
