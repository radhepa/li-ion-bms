/*
 * soc_ekf.c - EKF for state of charge (see bms/ekf.py for the derivation;
 * this is the same algorithm in single-precision C, as it would run on a
 * Cortex-M4F microcontroller).
 */
#include <math.h>
#include "bms_config.h"
#include "cell_tables.h"
#include "soc_ekf.h"

/* linear interpolation in a table on the uniform SOC grid, clamped at the ends */
static float lut(const float *tab, float soc)
{
    if (soc <= 0.0f)
        return tab[0];
    if (soc >= 1.0f)
        return tab[CELL_TABLE_N - 1];
    float pos = soc * (float)(CELL_TABLE_N - 1);
    int k = (int)pos;
    float f = pos - (float)k;
    return tab[k] + f * (tab[k + 1] - tab[k]);
}

float cell_ocv(float soc)
{
    const float step = 1.0f / (float)(CELL_TABLE_N - 1);
    /* past the ends: continue with the end slope (an over-charged cell's voltage keeps rising) */
    if (soc > 1.0f) {
        float top = (CELL_OCV[CELL_TABLE_N - 1] - CELL_OCV[CELL_TABLE_N - 3]) / (2.0f * step);
        return CELL_OCV[CELL_TABLE_N - 1] + top * (soc - 1.0f);
    }
    if (soc < 0.0f) {
        float bot = (CELL_OCV[2] - CELL_OCV[0]) / (2.0f * step);
        return CELL_OCV[0] + bot * soc;
    }
    return lut(CELL_OCV, soc);
}

float cell_docv(float soc)
{
    const float h = 0.005f;
    return (cell_ocv(soc + h) - cell_ocv(soc - h)) / (2.0f * h);
}

void cell_params(float soc, float *r0, float *r1, float *tau1, float *r2, float *tau2)
{
    *r0 = lut(CELL_R0, soc);
    *r1 = lut(CELL_R1, soc);
    *tau1 = lut(CELL_TAU1, soc);
    *r2 = lut(CELL_R2, soc);
    *tau2 = lut(CELL_TAU2, soc);
}

void ekf_init(SocEkf *e, float soc0)
{
    e->x[0] = soc0;
    e->x[1] = 0.0f;
    e->x[2] = 0.0f;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            e->P[r][c] = 0.0f;
    e->P[0][0] = EKF_SIGMA_SOC0 * EKF_SIGMA_SOC0;
    e->P[1][1] = 1e-4f;
    e->P[2][2] = 1e-4f;
}

float ekf_step(SocEkf *e, float i, float v, float dt)
{
    float r0, r1, t1, r2, t2;
    cell_params(e->x[0], &r0, &r1, &t1, &r2, &t2);
    float a1 = expf(-dt / t1), a2 = expf(-dt / t2);

    /* ---- predict: run the model forward with the measured current ---- */
    e->x[0] -= i * dt / (3600.0f * CELL_CAPACITY_AH);
    e->x[1] = a1 * e->x[1] + r1 * (1.0f - a1) * i;
    e->x[2] = a2 * e->x[2] + r2 * (1.0f - a2) * i;
    /* P = F P F' + Q, with F = diag(1, a1, a2) */
    float f[3] = {1.0f, a1, a2};
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            e->P[r][c] *= f[r] * f[c];
    float qs = EKF_SIGMA_I * dt / (3600.0f * CELL_CAPACITY_AH);
    e->P[0][0] += qs * qs;
    e->P[1][1] += EKF_SIGMA_V1 * EKF_SIGMA_V1 * dt;
    e->P[2][2] += EKF_SIGMA_V1 * EKF_SIGMA_V1 * dt;

    /* ---- update: compare the model's voltage with the measurement ---- */
    float soc = e->x[0];
    cell_params(soc, &r0, &r1, &t1, &r2, &t2);
    float y_hat = cell_ocv(soc) - r0 * i - e->x[1] - e->x[2];
    float H[3] = {cell_docv(soc), -1.0f, -1.0f};
    float PH[3];
    for (int r = 0; r < 3; r++)
        PH[r] = e->P[r][0] * H[0] + e->P[r][1] * H[1] + e->P[r][2] * H[2];
    float S = H[0] * PH[0] + H[1] * PH[1] + H[2] * PH[2] + EKF_SIGMA_V * EKF_SIGMA_V;
    float K[3] = {PH[0] / S, PH[1] / S, PH[2] / S};
    float innov = v - y_hat;
    for (int r = 0; r < 3; r++)
        e->x[r] += K[r] * innov;

    /* Joseph form: P = (I - K H) P (I - K H)' + K R K'  (stays symmetric, positive) */
    float A[3][3], T[3][3];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            A[r][c] = (r == c ? 1.0f : 0.0f) - K[r] * H[c];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            T[r][c] = A[r][0] * e->P[0][c] + A[r][1] * e->P[1][c] + A[r][2] * e->P[2][c];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            e->P[r][c] = T[r][0] * A[c][0] + T[r][1] * A[c][1] + T[r][2] * A[c][2]
                         + K[r] * K[c] * EKF_SIGMA_V * EKF_SIGMA_V;

    if (e->x[0] < -0.05f) e->x[0] = -0.05f;
    if (e->x[0] > 1.05f) e->x[0] = 1.05f;
    return e->x[0];
}
