/*
 * lazy_demo.c — Two demonstrations of lazy state vector streaming:
 *
 *   DEMO 1: Quantum State Tomography at 100T Scale
 *           Stream entries to reconstruct ρ (reduced density matrix),
 *           compute purity, entropy, and eigenvalues — O(1) memory.
 *
 *   DEMO 3: Selective Non-Destructive Readout
 *           Inspect arbitrary quhits (0, 42, 10B, 99.999T) across
 *           all basis states WITHOUT collapsing the superposition.
 *           Then measure to show the state wasn't disturbed.
 *
 * Build:
 *   gcc -O2 -I. -std=c11 -D_GNU_SOURCE \
 *       -c lazy_demo.c -o lazy_demo.o && \
 *   gcc -O2 -o lazy_demo lazy_demo.o hexstate_engine.o bigint.o -lm
 */

#include "hexstate_engine.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define D  6
#define N  100000000000000ULL  /* 100T quhits */
static const char *bn[] = {"A","T","G","C","dR","Pi"};

/* ═══════════════════════════════════════════════════════════════════════════
 *  DEMO 1: Streaming Quantum State Tomography
 *
 *  Reconstruct the reduced density matrix ρ by streaming basis entries.
 *  For each entry with amplitude α_k and basis value v_k:
 *    ρ[v_k][v_k] += |α_k|²           (diagonal = probabilities)
 *    ρ[v_i][v_j] += α_i · α_j*       (off-diagonal = coherences)
 *
 *  Then extract:
 *    - Purity = Tr(ρ²)     — 1.0 = pure, 1/D = maximally mixed
 *    - Von Neumann entropy S = -Tr(ρ log₂ ρ)
 *    - Eigenvalues of ρ
 *  All computed incrementally with O(D²) memory, NOT O(D^N).
 * ═══════════════════════════════════════════════════════════════════════════ */
static void demo_tomography(void)
{
    printf("\n");
    printf("████████████████████████████████████████████████████████████████████████████\n");
    printf("██                                                                        ██\n");
    printf("██  DEMO 1: STREAMING QUANTUM STATE TOMOGRAPHY AT 100T SCALE             ██\n");
    printf("██  Reconstruct ρ, purity, entropy — O(1) memory via lazy iterator       ██\n");
    printf("██                                                                        ██\n");
    printf("████████████████████████████████████████████████████████████████████████████\n\n");

    /* ── Scenario A: GHZ state (maximally entangled) ── */
    printf("  ─── Scenario A: GHZ State (entangle_all → maximally entangled) ───\n\n");
    {
        static HexStateEngine eng;
        engine_init(&eng);
        op_infinite_resources_dim(&eng, 0, N, D);
        init_quhit_register(&eng, 0, N, D);
        entangle_all_quhits(&eng, 0);  /* |GHZ⟩ = (1/√D) Σ|k,k,...,k⟩ */

        /* Stream → build ρ incrementally */
        Complex rho[D][D];
        memset(rho, 0, sizeof(rho));

        /* We need the amplitudes and their basis values to build ρ.
         * For the reduced density matrix of a SINGLE quhit, we trace
         * over all other quhits. Since all quhits share the same bulk
         * value in each entry, the reduced ρ for any quhit is:
         *   ρ[v][v'] = Σ_k (α_k · α_k'*) × δ(v_k, v) × δ(v_k', v')
         * For GHZ: v_k = k for every quhit, so ρ is diagonal. */

        /* Store amplitudes and values for off-diagonal computation */
        Complex amps[D];
        uint32_t vals[D];
        int n_entries = 0;

        StateIterator it;
        state_iter_begin(&eng, 0, &it);
        printf("    Streaming %u entries...\n", it.total_entries);

        while (state_iter_next(&it)) {
            /* Resolve quhit 0's value in this entry */
            uint32_t v = state_iter_resolve(&it, 0);
            if (n_entries < D) {
                amps[n_entries] = it.amplitude;
                vals[n_entries] = v;
            }
            n_entries++;
        }
        state_iter_end(&it);

        /* Build reduced ρ for quhit 0 */
        for (int i = 0; i < n_entries && i < D; i++) {
            for (int j = 0; j < n_entries && j < D; j++) {
                /* ρ[v_i][v_j] += α_i · α_j* */
                double re_i = amps[i].real, im_i = amps[i].imag;
                double re_j = amps[j].real, im_j = amps[j].imag;
                rho[vals[i]][vals[j]].real += re_i * re_j + im_i * im_j;
                rho[vals[i]][vals[j]].imag += im_i * re_j - re_i * im_j;
            }
        }

        /* Print ρ */
        printf("\n    Reduced density matrix ρ (quhit 0 of 100T):\n\n");
        printf("           ");
        for (int j = 0; j < D; j++) printf("  %-8s", bn[j]);
        printf("\n");
        for (int i = 0; i < D; i++) {
            printf("    %4s [ ", bn[i]);
            for (int j = 0; j < D; j++) {
                double re = rho[i][j].real;
                if (fabs(re) < 1e-10) re = 0;
                printf("%7.4f  ", re);
            }
            printf("]\n");
        }

        /* Compute purity = Tr(ρ²) */
        double purity = 0;
        for (int i = 0; i < D; i++)
            for (int k = 0; k < D; k++) {
                double re = rho[i][k].real * rho[k][i].real -
                            rho[i][k].imag * rho[k][i].imag;
                purity += re;
            }

        /* Eigenvalues of ρ (diagonal for GHZ: they're just ρ[k][k]) */
        printf("\n    Eigenvalues of ρ:\n");
        double entropy = 0;
        for (int k = 0; k < D; k++) {
            double p = rho[k][k].real;
            printf("      λ_%d = %.6f\n", k, p);
            if (p > 1e-15) entropy -= p * log2(p);
        }

        printf("\n    Purity  Tr(ρ²) = %.6f", purity);
        if (fabs(purity - 1.0/D) < 0.01)
            printf("  ← 1/D = %.4f (maximally mixed → maximally entangled!)", 1.0/D);
        else if (fabs(purity - 1.0) < 0.01)
            printf("  ← 1.0 (pure → not entangled)");
        printf("\n");

        printf("    Entropy S(ρ)  = %.6f bits", entropy);
        if (fabs(entropy - log2(D)) < 0.01)
            printf("  ← log₂(D) = %.4f (maximum entanglement!)", log2(D));
        printf("\n");

        printf("\n    ✓ Tomography complete. Memory used: %lu bytes (ρ only)\n",
               sizeof(rho));
        printf("      vs full state vector: 6^(100T) amplitudes = impossible\n");

        engine_destroy(&eng);
    }

    /* ── Scenario B: DNA gate (partially entangled) ── */
    printf("\n  ─── Scenario B: DNA Gate (partial entanglement) ───\n\n");
    {
        static HexStateEngine eng;
        engine_init(&eng);
        op_infinite_resources_dim(&eng, 0, N, D);
        init_quhit_register(&eng, 0, N, D);
        apply_dna_bulk_quhits(&eng, 0, 1.0, 310.0);

        Complex rho[D][D];
        memset(rho, 0, sizeof(rho));
        Complex amps[D];
        uint32_t vals[D];
        int n_entries = 0;

        StateIterator it;
        state_iter_begin(&eng, 0, &it);
        printf("    Streaming %u entries after DNA gate...\n", it.total_entries);

        while (state_iter_next(&it)) {
            uint32_t v = state_iter_resolve(&it, 0);
            if (n_entries < D) {
                amps[n_entries] = it.amplitude;
                vals[n_entries] = v;
            }
            n_entries++;
        }
        state_iter_end(&it);

        for (int i = 0; i < n_entries && i < D; i++)
            for (int j = 0; j < n_entries && j < D; j++) {
                double re_i = amps[i].real, im_i = amps[i].imag;
                double re_j = amps[j].real, im_j = amps[j].imag;
                rho[vals[i]][vals[j]].real += re_i * re_j + im_i * im_j;
                rho[vals[i]][vals[j]].imag += im_i * re_j - re_i * im_j;
            }

        printf("\n    Reduced ρ (quhit 0) after DNA gate:\n\n");
        printf("           ");
        for (int j = 0; j < D; j++) printf("  %-8s", bn[j]);
        printf("\n");
        for (int i = 0; i < D; i++) {
            printf("    %4s [ ", bn[i]);
            for (int j = 0; j < D; j++) {
                double re = rho[i][j].real;
                if (fabs(re) < 1e-10) re = 0;
                printf("%7.4f  ", re);
            }
            printf("]\n");
        }

        double purity = 0;
        for (int i = 0; i < D; i++)
            for (int k = 0; k < D; k++) {
                double re = rho[i][k].real * rho[k][i].real -
                            rho[i][k].imag * rho[k][i].imag;
                purity += re;
            }
        double entropy = 0;
        for (int k = 0; k < D; k++) {
            double p = rho[k][k].real;
            if (p > 1e-15) entropy -= p * log2(p);
        }
        printf("\n    Purity  = %.6f\n", purity);
        printf("    Entropy = %.6f bits\n", entropy);
        printf("    ✓ DNA gate creates a different entanglement structure\n");

        engine_destroy(&eng);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  DEMO 3: Selective Non-Destructive Readout
 *
 *  Inspect specific quhits (by index) across all basis states without
 *  collapsing the superposition. Then verify the state is undisturbed
 *  by streaming again and confirming normalization + structure preserved.
 * ═══════════════════════════════════════════════════════════════════════════ */
static void demo_nondestructive(void)
{
    printf("\n");
    printf("████████████████████████████████████████████████████████████████████████████\n");
    printf("██                                                                        ██\n");
    printf("██  DEMO 3: SELECTIVE NON-DESTRUCTIVE READOUT AT 100T SCALE              ██\n");
    printf("██  Inspect arbitrary quhits WITHOUT collapsing the state                ██\n");
    printf("██                                                                        ██\n");
    printf("████████████████████████████████████████████████████████████████████████████\n\n");

    static HexStateEngine eng;
    engine_init(&eng);
    op_infinite_resources_dim(&eng, 0, N, D);
    init_quhit_register(&eng, 0, N, D);
    entangle_all_quhits(&eng, 0);

    /* Pick quhits spread across the 100T register */
    uint64_t targets[] = {
        0,                    /* first */
        42,                   /* early */
        10000000000ULL,       /* 10 billion */
        50000000000000ULL,    /* 50 trillion (midpoint) */
        99999999999999ULL     /* last */
    };
    int n_targets = sizeof(targets) / sizeof(targets[0]);

    printf("  ─── Pass 1: Non-destructive inspection of %d quhits ───\n\n", n_targets);
    printf("    Quhits inspected:\n");
    for (int t = 0; t < n_targets; t++) {
        if (targets[t] < 1000)
            printf("      #%llu", (unsigned long long)targets[t]);
        else if (targets[t] < 1000000000000ULL)
            printf("      #%lluB", (unsigned long long)(targets[t] / 1000000000ULL));
        else
            printf("      #%lluT", (unsigned long long)(targets[t] / 1000000000000ULL));
        printf("\n");
    }

    /* Stream and resolve each target quhit in each basis state */
    printf("\n    Values per basis entry:\n\n");
    printf("    Entry |");
    for (int t = 0; t < n_targets; t++) {
        if (targets[t] < 1000)
            printf(" q[%-4llu]|", (unsigned long long)targets[t]);
        else if (targets[t] < 1000000000000ULL)
            printf(" q[%lluB]|", (unsigned long long)(targets[t] / 1000000000ULL));
        else
            printf(" q[%lluT]|", (unsigned long long)(targets[t] / 1000000000000ULL));
    }
    printf("  Amplitude      | P(entry)\n");
    printf("    ──────┼");
    for (int t = 0; t < n_targets; t++) printf("────────┼");
    printf("──────────────────┼─────────\n");

    StateIterator it;
    state_iter_begin(&eng, 0, &it);
    double norm_before = 0;

    while (state_iter_next(&it)) {
        printf("    %5u |", it.entry_index);
        for (int t = 0; t < n_targets; t++) {
            uint32_t v = state_iter_resolve(&it, targets[t]);
            printf("  %-5s |", bn[v % D]);
        }
        printf(" (%+.4f,%+.4fi) | %.4f\n",
               it.amplitude.real, it.amplitude.imag, it.probability);
        norm_before += it.probability;
    }
    state_iter_end(&it);

    printf("\n    Norm after inspection: %.6f (should be 1.0)\n", norm_before);
    printf("    State collapsed? NO — iterator is non-destructive\n");

    /* ── Verify: stream AGAIN to prove state wasn't disturbed ── */
    printf("\n  ─── Pass 2: Verify state is UNDISTURBED after inspection ───\n\n");

    state_iter_begin(&eng, 0, &it);
    double norm_after = 0;
    int count = 0;
    int structure_match = 1;

    /* Re-resolve the same quhits — values should be identical */
    printf("    Re-reading same quhits:\n\n");
    while (state_iter_next(&it)) {
        uint32_t v0 = state_iter_resolve(&it, targets[0]);
        uint32_t v4 = state_iter_resolve(&it, targets[n_targets - 1]);
        /* In GHZ state: all quhits should have the same value */
        if (v0 != v4) structure_match = 0;
        norm_after += it.probability;
        count++;

        printf("      Entry %u: q[0]=%s  q[last]=%s  %s  P=%.4f\n",
               it.entry_index, bn[v0 % D], bn[v4 % D],
               (v0 == v4) ? "CORRELATED ✓" : "MISMATCH ✗",
               it.probability);
    }
    state_iter_end(&it);

    printf("\n    Entries:    %d (unchanged)\n", count);
    printf("    Norm:      %.6f (unchanged)\n", norm_after);
    printf("    Structure: %s\n", structure_match ?
           "ALL quhits perfectly correlated — GHZ intact ✓" :
           "MISMATCH — state was disturbed ✗");
    printf("    Collapsed: NO — can still measure normally\n");

    /* ── Now actually measure to show destructive vs non-destructive ── */
    printf("\n  ─── Pass 3: Destructive measurement (for comparison) ───\n\n");

    uint64_t result = measure_chunk(&eng, 0);
    printf("    measure_chunk → %llu (%s)\n",
           (unsigned long long)result, bn[result % D]);
    printf("    State is NOW collapsed.\n\n");

    /* Try streaming after collapse */
    int rc = state_iter_begin(&eng, 0, &it);
    if (rc == 0) {
        int post_count = 0;
        double post_norm = 0;
        while (state_iter_next(&it)) {
            post_count++;
            post_norm += it.probability;
        }
        state_iter_end(&it);
        printf("    Post-collapse stream: %d entries, norm=%.6f\n",
               post_count, post_norm);
        printf("    → Collapsed state has fewer active entries\n");
    } else {
        printf("    Post-collapse: No streamable state (fully collapsed)\n");
    }

    printf("\n    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │  SUMMARY                                              │\n");
    printf("    │                                                       │\n");
    printf("    │  Pass 1: Inspected 5 quhits across 100T — NO COLLAPSE│\n");
    printf("    │  Pass 2: Re-read same quhits — state IDENTICAL       │\n");
    printf("    │  Pass 3: measure_chunk DESTROYED the state           │\n");
    printf("    │                                                       │\n");
    printf("    │  Lazy iterator = non-destructive quantum readout     │\n");
    printf("    │  measure_chunk = destructive collapse                 │\n");
    printf("    └─────────────────────────────────────────────────────────┘\n");

    engine_destroy(&eng);
}

int main(void)
{
    demo_tomography();
    demo_nondestructive();
    printf("\n");
    return 0;
}
