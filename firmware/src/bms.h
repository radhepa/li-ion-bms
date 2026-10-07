/*
 * bms.h - top level: one call per sample period does everything.
 *
 *   measurements in  ->  protection  ->  FET commands
 *                    ->  SOC (one EKF per cell)
 *                    ->  balancing mask
 *
 * Hardware-independent on purpose: on the real board an I2C driver fills
 * BmsInputs from the BQ76920 and applies BmsOutputs; on a PC the
 * software-in-the-loop runner (sil/sil_main.c) does the same from a file.
 */
#ifndef BMS_H
#define BMS_H

#include <stdint.h>
#include "bms_config.h"
#include "protection.h"
#include "soc_ekf.h"

typedef struct {
    float cell_v[BMS_N_CELLS];   /* V                         */
    float current_a;             /* A, + = discharge          */
    float temp_c[2];             /* two NTCs on the cells, C  */
} BmsInputs;

typedef struct {
    uint8_t  chg_fet;            /* 1 = charge FET on          */
    uint8_t  dsg_fet;            /* 1 = discharge FET on       */
    uint8_t  balance_mask;       /* bit c = bleed cell c+1     */
    uint16_t faults;             /* FAULT_* bits               */
    float    soc[BMS_N_CELLS];   /* 0..1                       */
    float    pack_soc;           /* the lowest cell limits the pack */
} BmsOutputs;

typedef struct {
    Protection prot;
    SocEkf ekf[BMS_N_CELLS];
    uint32_t t_ms;
    uint8_t last_mask;
} Bms;

void bms_init(Bms *b, float soc_guess);
void bms_step(Bms *b, const BmsInputs *in, uint32_t dt_ms, BmsOutputs *out);

#endif
