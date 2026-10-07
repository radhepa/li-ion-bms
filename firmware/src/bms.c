#include "bms.h"
#include "balancing.h"

void bms_init(Bms *b, float soc_guess)
{
    protection_init(&b->prot);
    for (int c = 0; c < BMS_N_CELLS; c++)
        ekf_init(&b->ekf[c], soc_guess);
    b->t_ms = 0u;
    b->last_mask = 0u;
}

void bms_step(Bms *b, const BmsInputs *in, uint32_t dt_ms, BmsOutputs *out)
{
    b->t_ms += dt_ms;

    /* 1. protection */
    ProtInputs p;
    for (int c = 0; c < BMS_N_CELLS; c++)
        p.cell_v[c] = in->cell_v[c];
    p.current_a = in->current_a;
    p.temp_c = in->temp_c[0] > in->temp_c[1] ? in->temp_c[0] : in->temp_c[1];
    p.temp_min_c = in->temp_c[0] < in->temp_c[1] ? in->temp_c[0] : in->temp_c[1];
    protection_update(&b->prot, &p, dt_ms);
    out->faults = b->prot.active;
    out->chg_fet = (uint8_t)protection_charge_allowed(&b->prot);
    out->dsg_fet = (uint8_t)protection_discharge_allowed(&b->prot);

    /* 2. SOC: a bleeding cell also supplies its balancing current */
    float dt = (float)dt_ms * 1e-3f;
    float pack = 1.0f;
    for (int c = 0; c < BMS_N_CELLS; c++) {
        float i_cell = in->current_a;
        if (b->last_mask & (1u << c))
            i_cell += in->cell_v[c] / BAL_RESISTOR_OHM;
        out->soc[c] = ekf_step(&b->ekf[c], i_cell, in->cell_v[c], dt);
        if (out->soc[c] < pack)
            pack = out->soc[c];
    }
    out->pack_soc = pack;

    /* 3. balancing (not while any fault is active) */
    out->balance_mask = out->faults ? 0u : balancing_mask(in->cell_v, in->current_a, b->t_ms);
    b->last_mask = out->balance_mask;
}
