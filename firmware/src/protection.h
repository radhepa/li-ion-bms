/*
 * protection.h - fault detection with debounce, hysteresis and retry.
 */
#ifndef PROTECTION_H
#define PROTECTION_H

#include <stdint.h>
#include "bms_config.h"

enum {
    FAULT_OV     = 1u << 0,   /* a cell above CELL_OV_V              */
    FAULT_UV     = 1u << 1,   /* a cell below CELL_UV_V              */
    FAULT_OCD    = 1u << 2,   /* discharge over-current               */
    FAULT_OCC    = 1u << 3,   /* charge over-current                  */
    FAULT_SCD    = 1u << 4,   /* short circuit                        */
    FAULT_OT_DSG = 1u << 5,   /* too hot to discharge (stops all)     */
    FAULT_OT_CHG = 1u << 6,   /* too hot to charge                    */
    FAULT_UT_CHG = 1u << 7,   /* too cold to charge                   */
    FAULT_COUNT  = 8
};

typedef struct {
    float cell_v[BMS_N_CELLS];
    float current_a;          /* + = discharge */
    float temp_c;             /* hottest (or, for UT, coldest) sensor handled inside */
    float temp_min_c;
} ProtInputs;

typedef struct {
    uint16_t active;                      /* FAULT_* bits currently latched      */
    uint32_t set_ms[FAULT_COUNT];         /* how long the trip condition has held */
    uint32_t clear_ms[FAULT_COUNT];       /* how long the clear condition has held */
} Protection;

void protection_init(Protection *p);
void protection_update(Protection *p, const ProtInputs *in, uint32_t dt_ms);
int  protection_charge_allowed(const Protection *p);
int  protection_discharge_allowed(const Protection *p);

#endif
