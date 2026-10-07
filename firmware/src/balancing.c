#include "balancing.h"

uint8_t balancing_mask(const float cell_v[BMS_N_CELLS], float current_a, uint32_t t_ms)
{
    if (current_a > RESTING_A)                     /* discharging: never bleed */
        return 0u;
    float vmin = cell_v[0];
    for (int c = 1; c < BMS_N_CELLS; c++)
        if (cell_v[c] < vmin)
            vmin = cell_v[c];
    int phase = (int)((t_ms / BAL_SWAP_MS) % 2u);  /* 0: cells 1,3   1: cells 2,4 */
    uint8_t mask = 0u;
    for (int c = 0; c < BMS_N_CELLS; c++) {
        if (c % 2 != phase)
            continue;
        if (cell_v[c] > BAL_MIN_V && cell_v[c] > vmin + BAL_WINDOW_V)
            mask |= (uint8_t)(1u << c);
    }
    return mask;
}
