/*
 * balancing.h - passive cell balancing decisions.
 */
#ifndef BALANCING_H
#define BALANCING_H

#include <stdint.h>
#include "bms_config.h"

/* Returns a bit mask of cells to bleed (bit c = cell c+1).
 * Rules: only while charging or resting, only cells above BAL_MIN_V and more
 * than BAL_WINDOW_V above the lowest cell, and never two neighbouring cells at
 * once (the BQ76920 forbids it): odd and even cells alternate every BAL_SWAP_MS. */
uint8_t balancing_mask(const float cell_v[BMS_N_CELLS], float current_a, uint32_t t_ms);

#endif
