/*
 * protection.c - every fault follows the same pattern:
 *
 *   trip  : the bad condition must hold continuously for a debounce time
 *           (so one noisy sample can't open the FETs);
 *   clear : a DIFFERENT, safer condition must hold for a while (hysteresis),
 *           or for current faults a fixed retry delay passes, so the pack
 *           doesn't chatter on and off at the threshold.
 */
#include <string.h>
#include "protection.h"

static int bit_index(uint16_t bit)
{
    int k = 0;
    while (bit > 1u) {
        bit >>= 1;
        k++;
    }
    return k;
}

/* advance one fault's state machine */
static void check(Protection *p, uint16_t bit, int trip_cond, uint32_t trip_ms,
                  int clear_cond, uint32_t clear_ms, uint32_t dt)
{
    int k = bit_index(bit);
    if (!(p->active & bit)) {
        p->set_ms[k] = trip_cond ? p->set_ms[k] + dt : 0u;
        if (trip_cond && p->set_ms[k] >= trip_ms) {
            p->active |= bit;
            p->clear_ms[k] = 0u;
        }
    } else {
        p->clear_ms[k] = clear_cond ? p->clear_ms[k] + dt : 0u;
        if (p->clear_ms[k] >= clear_ms) {
            p->active &= (uint16_t)~bit;
            p->set_ms[k] = 0u;
        }
    }
}

void protection_init(Protection *p)
{
    memset(p, 0, sizeof *p);
}

void protection_update(Protection *p, const ProtInputs *in, uint32_t dt)
{
    float vmax = in->cell_v[0], vmin = in->cell_v[0];
    for (int c = 1; c < BMS_N_CELLS; c++) {
        if (in->cell_v[c] > vmax) vmax = in->cell_v[c];
        if (in->cell_v[c] < vmin) vmin = in->cell_v[c];
    }
    float i = in->current_a;
    int charging = i < CHARGING_A;

    check(p, FAULT_OV, vmax > CELL_OV_V, CELL_DEBOUNCE_MS, vmax < CELL_OV_RECOVER_V, CELL_DEBOUNCE_MS, dt);
    /* under-voltage clears once the cells recover OR a charger pushes current in */
    check(p, FAULT_UV, vmin < CELL_UV_V, CELL_DEBOUNCE_MS,
          vmin > CELL_UV_RECOVER_V || charging, CELL_DEBOUNCE_MS, dt);
    /* current faults: the FET is open, so "clear" is simply a retry timer */
    check(p, FAULT_OCD, i > OCD_A, OCD_DEBOUNCE_MS, 1, CURRENT_RETRY_MS, dt);
    check(p, FAULT_OCC, -i > OCC_A, OCC_DEBOUNCE_MS, 1, CURRENT_RETRY_MS, dt);
    check(p, FAULT_SCD, i > SCD_A, 0u, 1, CURRENT_RETRY_MS, dt);
    check(p, FAULT_OT_DSG, in->temp_c > OT_DSG_C, TEMP_DEBOUNCE_MS,
          in->temp_c < OT_DSG_RECOVER_C, TEMP_DEBOUNCE_MS, dt);
    check(p, FAULT_OT_CHG, in->temp_c > OT_CHG_C, TEMP_DEBOUNCE_MS,
          in->temp_c < OT_CHG_RECOVER_C, TEMP_DEBOUNCE_MS, dt);
    check(p, FAULT_UT_CHG, in->temp_min_c < UT_CHG_C, TEMP_DEBOUNCE_MS,
          in->temp_min_c > UT_CHG_RECOVER_C, TEMP_DEBOUNCE_MS, dt);
}

int protection_charge_allowed(const Protection *p)
{
    return !(p->active & (FAULT_OV | FAULT_OCC | FAULT_SCD | FAULT_OT_DSG | FAULT_OT_CHG | FAULT_UT_CHG));
}

int protection_discharge_allowed(const Protection *p)
{
    return !(p->active & (FAULT_UV | FAULT_OCD | FAULT_SCD | FAULT_OT_DSG));
}
