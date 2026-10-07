"""
Extended Kalman Filter for state of charge.

Two imperfect sources of information, combined optimally:

  * Coulomb counting (integrate current): smooth and precise over minutes,
    but it needs the right starting SOC and every current-sensor offset adds
    up forever.
  * Voltage: tied to SOC through the OCV curve, so it never drifts, but it is
    noisy and only weakly related to SOC in the flat middle of the curve, and
    under load it is shifted by R0*i and the RC voltages.

The EKF runs the cell model forward with the measured current (predict), then
compares the model's voltage with the measured voltage and nudges the state by
an amount set by how much it trusts each source (update). Its uncertainty P
grows during prediction and shrinks with every voltage measurement.

    state  x = [SOC, V1, V2]
    predict:  SOC -= i dt / (3600 Q);  V_k = a_k V_k + R_k (1 - a_k) i,  a_k = exp(-dt / tau_k)
    measure:  v = OCV(SOC) - R0 i - V1 - V2,   H = [dOCV/dSOC, -1, -1]
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import CellModel


@dataclass
class EKFTuning:
    soc0: float = 0.5          # initial guess when the true SOC is unknown
    sigma_soc0: float = 0.3    # how unsure that guess is
    sigma_i: float = 0.05      # current-sensor error (A) -> SOC process noise
    sigma_v1: float = 1e-4     # RC-voltage process noise per sqrt(s), V
    sigma_v: float = 0.010     # voltage measurement + model error, V


class SocEKF:
    def __init__(self, model: CellModel, tune: EKFTuning = EKFTuning()):
        self.m = model
        self.tune = tune
        self.x = np.array([tune.soc0, 0.0, 0.0])
        self.P = np.diag([tune.sigma_soc0 ** 2, 1e-4, 1e-4])

    def step(self, i: float, v: float, dt: float) -> float:
        """One sample: current i (A, + = discharge) applied for dt, then voltage v."""
        m, tune = self.m, self.tune
        soc, v1, v2 = self.x
        r0, r1, t1, r2, t2 = m.params_at(soc)
        a1, a2 = np.exp(-dt / t1), np.exp(-dt / t2)
        # predict
        self.x = np.array([soc - i * dt / (3600.0 * m.q_ah),
                           a1 * v1 + r1 * (1 - a1) * i,
                           a2 * v2 + r2 * (1 - a2) * i])
        F = np.diag([1.0, a1, a2])
        qs = (tune.sigma_i * dt / (3600.0 * m.q_ah)) ** 2
        Q = np.diag([qs, tune.sigma_v1 ** 2 * dt, tune.sigma_v1 ** 2 * dt])
        self.P = F @ self.P @ F.T + Q
        # update with the voltage measurement
        soc = self.x[0]
        y_hat = m.ocv_at(soc) - m.params_at(soc)[0] * i - self.x[1] - self.x[2]
        H = np.array([m.docv_at(soc), -1.0, -1.0])
        S = H @ self.P @ H + tune.sigma_v ** 2
        K = self.P @ H / S
        self.x = self.x + K * (v - y_hat)
        I_KH = np.eye(3) - np.outer(K, H)
        self.P = I_KH @ self.P @ I_KH.T + np.outer(K, K) * tune.sigma_v ** 2   # Joseph form
        self.x[0] = np.clip(self.x[0], -0.05, 1.05)
        return float(self.x[0])

    def run(self, t, i, v):
        out = np.empty(len(t))
        out[0] = self.x[0]
        for k in range(1, len(t)):
            out[k] = self.step(i[k - 1], v[k], t[k] - t[k - 1])
        return out


def coulomb_count(t, i, soc0, q_ah):
    q = np.concatenate([[0.0], np.cumsum(i[:-1] * np.diff(t))])
    return soc0 - q / 3600.0 / q_ah
