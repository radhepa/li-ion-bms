/*
 * sil_main.c - software-in-the-loop runner: the real firmware (bms_step), fed
 * from a CSV file instead of the hardware.
 *
 * input  columns: t_s, i_request_a, v1, v2, v3, v4, temp1_c, temp2_c
 * output columns: t_s, i_actual_a, chg_fet, dsg_fet, balance_mask, faults,
 *                 soc1, soc2, soc3, soc4, pack_soc
 *
 * A tiny "plant" closes the loop on current: a discharge request only flows
 * while the discharge FET is on, a charge request only while the charge FET
 * is on - so when protection trips, the current really stops.
 *
 * With --open-loop the requested current always flows (used to compare the
 * estimator against a reference on recorded data).
 *
 *   usage: bms_sil in.csv out.csv [soc_guess] [--open-loop]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bms.h"

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: bms_sil in.csv out.csv [soc_guess]\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "r"), *out = fopen(argv[2], "w");
    if (!in || !out) {
        fprintf(stderr, "bms_sil: cannot open files\n");
        return 1;
    }
    float soc0 = argc > 3 ? (float)atof(argv[3]) : 0.5f;
    int open_loop = argc > 4 && strcmp(argv[4], "--open-loop") == 0;
    Bms bms;
    bms_init(&bms, soc0);
    BmsOutputs o = {1, 1, 0, 0, {0}, 0};
    o.chg_fet = o.dsg_fet = 1;

    char line[512];
    if (!fgets(line, sizeof line, in)) {           /* header */
        fprintf(stderr, "bms_sil: empty input\n");
        return 1;
    }
    fprintf(out, "t_s,i_actual_a,chg_fet,dsg_fet,balance_mask,faults,soc1,soc2,soc3,soc4,pack_soc\n");
    double t_prev = -1.0;
    long rows = 0;
    while (fgets(line, sizeof line, in)) {
        double t;
        float ireq, v[4], t1, t2;
        if (sscanf(line, "%lf,%f,%f,%f,%f,%f,%f,%f", &t, &ireq, &v[0], &v[1], &v[2], &v[3], &t1, &t2) != 8)
            continue;
        /* plant: the FETs decide whether the requested current can flow */
        float i = ireq;
        if (!open_loop && ireq > 0.0f && !o.dsg_fet) i = 0.0f;
        if (!open_loop && ireq < 0.0f && !o.chg_fet) i = 0.0f;
        BmsInputs bi = {{v[0], v[1], v[2], v[3]}, i, {t1, t2}};
        uint32_t dt_ms = t_prev < 0 ? 0u : (uint32_t)((t - t_prev) * 1000.0 + 0.5);
        t_prev = t;
        bms_step(&bms, &bi, dt_ms, &o);
        fprintf(out, "%.3f,%.4f,%u,%u,%u,%u,%.6f,%.6f,%.6f,%.6f,%.6f\n", t, (double)i, o.chg_fet, o.dsg_fet,
                o.balance_mask, o.faults, (double)o.soc[0], (double)o.soc[1], (double)o.soc[2],
                (double)o.soc[3], (double)o.pack_soc);
        rows++;
    }
    fclose(in);
    fclose(out);
    fprintf(stderr, "bms_sil: %ld samples processed\n", rows);
    return 0;
}
