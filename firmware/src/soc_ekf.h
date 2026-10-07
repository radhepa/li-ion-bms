/*
 * soc_ekf.h - Extended Kalman Filter state-of-charge estimator for one cell.
 * State: x = [SOC, V1, V2] (2-RC equivalent circuit model).
 */
#ifndef SOC_EKF_H
#define SOC_EKF_H

typedef struct {
    float x[3];        /* SOC (0..1), V1 (V), V2 (V) */
    float P[3][3];     /* state covariance           */
} SocEkf;

void  ekf_init(SocEkf *e, float soc0);
/* current i (A, + = discharge) applied for dt seconds, then voltage v measured */
float ekf_step(SocEkf *e, float i, float v, float dt);

/* model lookups, shared with the unit tests */
float cell_ocv(float soc);
float cell_docv(float soc);
void  cell_params(float soc, float *r0, float *r1, float *tau1, float *r2, float *tau2);

#endif
