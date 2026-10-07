/*
 * test_bms.c - unit tests for the BMS firmware. Run with:  make test
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "balancing.h"
#include "bms.h"
#include "cell_tables.h"
#include "protection.h"
#include "soc_ekf.h"

static int checks = 0, failures = 0;
#define CHECK(cond, msg)                                              \
    do {                                                              \
        checks++;                                                     \
        if (!(cond)) {                                                \
            failures++;                                               \
            printf("FAIL line %d: %s\n", __LINE__, msg);              \
        }                                                             \
    } while (0)

static ProtInputs nominal(void)
{
    ProtInputs p = {{3.7f, 3.7f, 3.7f, 3.7f}, 0.0f, 25.0f, 25.0f};
    return p;
}

/* run the protection logic for `ms` milliseconds in 10 ms steps */
static void run(Protection *p, const ProtInputs *in, int ms)
{
    for (int t = 0; t < ms; t += 10)
        protection_update(p, in, 10u);
}

static void test_overvoltage(void)
{
    Protection p;
    protection_init(&p);
    ProtInputs in = nominal();
    in.cell_v[2] = 4.26f;
    run(&p, &in, 900);
    CHECK(!(p.active & FAULT_OV), "OV must not trip before the 1 s debounce");
    run(&p, &in, 200);
    CHECK(p.active & FAULT_OV, "OV trips after 1 s");
    CHECK(!protection_charge_allowed(&p) && protection_discharge_allowed(&p), "OV blocks charge only");
    in.cell_v[2] = 4.20f;                       /* below the trip, above the recovery level */
    run(&p, &in, 3000);
    CHECK(p.active & FAULT_OV, "OV holds until the cell drops below 4.15 V (hysteresis)");
    in.cell_v[2] = 4.10f;
    run(&p, &in, 1100);
    CHECK(!(p.active & FAULT_OV), "OV clears after 1 s below 4.15 V");
}

static void test_noise_rejection(void)
{
    Protection p;
    protection_init(&p);
    ProtInputs in = nominal();
    for (int k = 0; k < 50; k++) {              /* 500 ms of alternating spikes */
        in.cell_v[0] = (k % 2) ? 4.30f : 4.10f;
        protection_update(&p, &in, 10u);
    }
    CHECK(!(p.active & FAULT_OV), "single-sample spikes must not trip OV");
}

static void test_overcurrent(void)
{
    Protection p;
    protection_init(&p);
    ProtInputs in = nominal();
    in.current_a = 25.0f;
    run(&p, &in, 400);
    CHECK(!(p.active & FAULT_OCD), "OCD waits for 500 ms");
    run(&p, &in, 200);
    CHECK(p.active & FAULT_OCD, "OCD trips");
    CHECK(!protection_discharge_allowed(&p), "OCD opens the discharge FET");
    in.current_a = 0.0f;
    run(&p, &in, 29000);
    CHECK(p.active & FAULT_OCD, "OCD retries only after 30 s");
    run(&p, &in, 1100);
    CHECK(!(p.active & FAULT_OCD), "OCD retry");

    protection_init(&p);
    in.current_a = 80.0f;
    protection_update(&p, &in, 10u);
    CHECK(p.active & FAULT_SCD, "short circuit trips on the first sample");
}

static void test_temperature(void)
{
    Protection p;
    protection_init(&p);
    ProtInputs in = nominal();
    in.temp_min_c = -5.0f;
    run(&p, &in, 2100);
    CHECK(!protection_charge_allowed(&p), "no charging below 0 C");
    CHECK(protection_discharge_allowed(&p), "discharging below 0 C is fine");
    in = nominal();
    protection_init(&p);
    in.temp_c = 65.0f;
    run(&p, &in, 2100);
    CHECK(!protection_charge_allowed(&p) && !protection_discharge_allowed(&p), "over-temperature stops both");
}

static void test_balancing(void)
{
    srand(1);
    int adjacent = 0, during_discharge = 0;
    for (int k = 0; k < 20000; k++) {
        float v[4];
        for (int c = 0; c < 4; c++)
            v[c] = 3.85f + 0.3f * (float)rand() / (float)RAND_MAX;
        uint8_t m = balancing_mask(v, -1.0f, (uint32_t)(k * 997));
        if (m & (m >> 1))
            adjacent++;
        if (balancing_mask(v, 5.0f, (uint32_t)k))
            during_discharge++;
    }
    CHECK(adjacent == 0, "never bleed two neighbouring cells at once");
    CHECK(during_discharge == 0, "never bleed while discharging");
    float v[4] = {4.10f, 4.05f, 4.10f, 4.12f};
    CHECK(balancing_mask(v, 0.0f, 0u) == 0x5u, "phase 0: cells 1 and 3 above min+10 mV");
    CHECK(balancing_mask(v, 0.0f, BAL_SWAP_MS) == 0x8u, "phase 1: cell 4 (cell 2 is the lowest)");
}

static void test_ocv_table(void)
{
    int monotonic = 1;
    for (int k = 1; k < CELL_TABLE_N; k++)
        if (CELL_OCV[k] < CELL_OCV[k - 1])
            monotonic = 0;
    CHECK(monotonic, "OCV table rises with SOC");
    CHECK(fabsf(cell_ocv(1.0f) - CELL_OCV[CELL_TABLE_N - 1]) < 1e-6f, "OCV at 100 %");
    CHECK(cell_ocv(1.02f) > cell_ocv(1.0f), "OCV keeps rising past 100 %");
    CHECK(cell_docv(0.5f) > 0.0f, "positive slope");
}

static void test_ekf_converges(void)
{
    /* simulate a cell with the firmware's own model: 1 A discharge from 80 % */
    float soc = 0.80f, v1 = 0.0f, v2 = 0.0f, dt = 1.0f;
    SocEkf e;
    ekf_init(&e, 0.40f);                        /* starts 40 % wrong */
    float est = 0.0f;
    for (int k = 0; k < 1800; k++) {
        float i = (k / 60) % 2 ? 2.0f : 0.5f;   /* changing load */
        float r0, r1, t1, r2, t2;
        cell_params(soc, &r0, &r1, &t1, &r2, &t2);
        float a1 = expf(-dt / t1), a2 = expf(-dt / t2);
        v1 = a1 * v1 + r1 * (1 - a1) * i;
        v2 = a2 * v2 + r2 * (1 - a2) * i;
        soc -= i * dt / (3600.0f * CELL_CAPACITY_AH);
        cell_params(soc, &r0, &r1, &t1, &r2, &t2);
        float v = cell_ocv(soc) - r0 * i - v1 - v2;
        est = ekf_step(&e, i, v, dt);
    }
    CHECK(fabsf(est - soc) < 0.01f, "EKF converges to within 1 % SOC in 30 min");
}

int main(void)
{
    test_overvoltage();
    test_noise_rejection();
    test_overcurrent();
    test_temperature();
    test_balancing();
    test_ocv_table();
    test_ekf_converges();
    printf("%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
