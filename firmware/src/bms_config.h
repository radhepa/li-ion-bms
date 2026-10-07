/*
 * bms_config.h - limits and tuning for a 4S pack of LG HG2 18650 cells.
 *
 * Cell ratings (LG HG2): 4.20 V max, 2.5 V min, 20 A continuous discharge,
 * 4 A max charge, charge 0..50 C. Every limit below keeps margin to those.
 */
#ifndef BMS_CONFIG_H
#define BMS_CONFIG_H

#define BMS_N_CELLS            4

/* ---- cell voltage ---------------------------------------------------- */
#define CELL_OV_V              4.25f   /* over-voltage trip                */
#define CELL_OV_RECOVER_V      4.15f   /* ...cleared below this            */
#define CELL_UV_V              2.70f   /* under-voltage trip               */
#define CELL_UV_RECOVER_V      3.00f   /* ...cleared above this            */
#define CELL_DEBOUNCE_MS       1000u

/* ---- current (A, + = discharge) ---------------------------------------- */
#define OCD_A                  22.0f   /* discharge over-current           */
#define OCD_DEBOUNCE_MS        500u
#define OCC_A                  4.5f    /* charge over-current (magnitude)  */
#define OCC_DEBOUNCE_MS        1000u
#define SCD_A                  60.0f   /* short circuit: the AFE trips in ~200 us in hardware;
                                          firmware sees it within one sample */
#define CURRENT_RETRY_MS       30000u  /* OCD / OCC / SCD: retry after 30 s */
#define CHARGING_A            (-0.05f) /* below this the pack is charging   */
#define RESTING_A              0.05f   /* |I| below this = resting          */

/* ---- temperature (C) ---------------------------------------------------- */
#define OT_DSG_C               60.0f
#define OT_DSG_RECOVER_C       50.0f
#define OT_CHG_C               45.0f
#define OT_CHG_RECOVER_C       40.0f
#define UT_CHG_C               0.0f    /* no charging below 0 C: lithium plating */
#define UT_CHG_RECOVER_C       3.0f
#define TEMP_DEBOUNCE_MS       2000u

/* ---- passive balancing --------------------------------------------------- */
#define BAL_WINDOW_V           0.010f  /* bleed cells > min + 10 mV        */
#define BAL_MIN_V              3.90f   /* only near the top of charge      */
#define BAL_SWAP_MS            60000u  /* odd / even cells take turns      */
#define BAL_RESISTOR_OHM       90.0f   /* 2 x 43 ohm input resistors + switch: ~45 mA */

/* ---- SOC estimator (EKF) -------------------------------------------------- */
#define EKF_SIGMA_SOC0         0.30f
#define EKF_SIGMA_I            0.05f   /* A   */
#define EKF_SIGMA_V1           1e-4f   /* V/sqrt(s) */
#define EKF_SIGMA_V            0.020f  /* V   */

#endif
