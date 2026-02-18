/*
 * hilbert_memory_test.c — Does the Hilbert Space Remember?
 *
 * PROTOCOL:
 *   1. Create a 100T quhit register
 *   2. Put all quhits in superposition (GHZ via DFT-bulk)
 *   3. Promote a specific quhit (e.g. q42) by applying DFT to it
 *      → This promotes q42 into addr[] entries, giving it its own basis
 *   4. Inspect the state: how many entries, what does q42 look like?
 *   5. "Release" q42 by manually removing it from addr[] in all entries
 *      → q42 now falls back to lazy_resolve: (bulk + 42) % 6
 *   6. Re-grab q42: apply DFT again to promote it back
 *   7. Inspect: does the Hilbert space remember the first gate?
 *      → If entries differ from step 4, the release destroyed information
 *      → If entries match, the Hilbert space retained the state
 *
 *   We also test WITHOUT release as a control:
 *   8. Apply DFT twice to q42 (no release in between)
 *   9. Compare results
 *
 * This probes whether the quhit register's sparse entry system
 * genuinely stores quantum state that persists independently of
 * which quhits are currently promoted.
 *
 * Build:
 *   gcc -O2 -I. -o hilbert_memory_test hilbert_memory_test.c \
 *       hexstate_engine.o bigint.o -lm
 */
#include "hexstate_engine.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define D       6
#define N       100000000000000ULL  /* 100T quhits */

static HexStateEngine eng;

/* Print the state of the register: entries and a specific quhit's resolved values */
static void print_state(const char *label, uint64_t quhit_idx)
{
    int r = find_quhit_reg(&eng, 0);
    if (r < 0) { printf("  [ERROR] No register found\n"); return; }

    uint32_t nz = eng.quhit_regs[r].num_nonzero;
    printf("\n  ── %s ──\n", label);
    printf("  Entries: %u\n", nz);

    for (uint32_t e = 0; e < nz && e < 20; e++) {
        QuhitBasisEntry *ent = &eng.quhit_regs[r].entries[e];
        double p = ent->amplitude.real * ent->amplitude.real +
                   ent->amplitude.imag * ent->amplitude.imag;

        printf("    [%2u] bulk=%u, amp=(%.4f, %.4f), |amp|²=%.6f",
               e, ent->bulk_value, ent->amplitude.real, ent->amplitude.imag, p);

        /* Show addr[] contents */
        if (ent->num_addr > 0) {
            printf(", addr={");
            for (uint8_t i = 0; i < ent->num_addr; i++) {
                printf("q%lu=%u", (unsigned long)ent->addr[i].quhit_idx, ent->addr[i].value);
                if (i + 1 < ent->num_addr) printf(", ");
            }
            printf("}");
        }

        /* Show lazy_resolve for the target quhit */
        uint32_t v = 0;
        /* inline lazy_resolve */
        int found = 0;
        for (uint8_t i = 0; i < ent->num_addr; i++) {
            if (ent->addr[i].quhit_idx == quhit_idx) {
                v = ent->addr[i].value;
                found = 1;
                break;
            }
        }
        if (!found) {
            if (eng.quhit_regs[r].bulk_rule == 1)
                v = (uint32_t)((ent->bulk_value + quhit_idx) % D);
            else
                v = ent->bulk_value;
        }
        printf(" → q%lu resolves to %u%s",
               (unsigned long)quhit_idx, v, found ? " (PROMOTED)" : " (bulk-derived)");
        printf("\n");
    }
    if (nz > 20) printf("    ... (%u more entries)\n", nz - 20);
}

/* Remove a quhit from addr[] in all entries (simulate "release from memory") */
static void release_quhit(uint64_t quhit_idx)
{
    int r = find_quhit_reg(&eng, 0);
    if (r < 0) return;

    uint32_t nz = eng.quhit_regs[r].num_nonzero;
    for (uint32_t e = 0; e < nz; e++) {
        QuhitBasisEntry *ent = &eng.quhit_regs[r].entries[e];
        for (uint8_t i = 0; i < ent->num_addr; i++) {
            if (ent->addr[i].quhit_idx == quhit_idx) {
                /* Remove by shifting remaining entries */
                for (uint8_t j = i; j + 1 < ent->num_addr; j++)
                    ent->addr[j] = ent->addr[j + 1];
                ent->num_addr--;
                break;
            }
        }
    }
}

/* Compute the marginal probability distribution for a quhit (without measuring) */
static void marginal_probs(uint64_t quhit_idx, double *probs_out)
{
    int r = find_quhit_reg(&eng, 0);
    if (r < 0) return;

    memset(probs_out, 0, D * sizeof(double));
    uint32_t nz = eng.quhit_regs[r].num_nonzero;
    for (uint32_t e = 0; e < nz; e++) {
        QuhitBasisEntry *ent = &eng.quhit_regs[r].entries[e];

        /* lazy_resolve inline */
        uint32_t v = eng.quhit_regs[r].bulk_rule == 1 ?
            (uint32_t)((ent->bulk_value + quhit_idx) % D) : ent->bulk_value;
        for (uint8_t i = 0; i < ent->num_addr; i++) {
            if (ent->addr[i].quhit_idx == quhit_idx) {
                v = ent->addr[i].value;
                break;
            }
        }

        double p = ent->amplitude.real * ent->amplitude.real +
                   ent->amplitude.imag * ent->amplitude.imag;
        if (v < D) probs_out[v] += p;
    }
}

static void print_marginal(const char *label, uint64_t quhit_idx)
{
    double probs[D];
    marginal_probs(quhit_idx, probs);
    printf("  Marginal P(q%lu): [", (unsigned long)quhit_idx);
    for (int i = 0; i < D; i++) {
        printf("%.4f", probs[i]);
        if (i + 1 < D) printf(", ");
    }
    printf("]\n");
}

int main(void)
{
    setbuf(stdout, NULL);
    engine_init(&eng);

    printf("\n");
    printf("  ╔═══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║  HILBERT SPACE MEMORY TEST                                       ║\n");
    printf("  ║  Does the Hilbert space remember what we did to a quhit?         ║\n");
    printf("  ╚═══════════════════════════════════════════════════════════════════╝\n\n");

    uint64_t TARGET = 42;  /* the quhit we'll promote, gate, release, re-grab */

    /* ═══════════════════════════════════════════════════════════════════
     *  TEST 1: BASELINE — What does an untouched quhit look like?
     * ═══════════════════════════════════════════════════════════════════ */
    printf("  ═══ TEST 1: BASELINE (untouched quhit) ═══\n");
    eng.num_quhit_regs = 0;
    init_quhit_register(&eng, 0, N, D);
    eng.quhit_regs[0].bulk_rule = 1;  /* cyclic: V(k) = (bulk + k) % D */

    print_state("Initial state (|0⟩ bulk)", TARGET);
    print_marginal("  Before superposition", TARGET);

    /* Create GHZ superposition */
    entangle_all_quhits(&eng, 0);
    print_state("After DFT-bulk (GHZ superposition)", TARGET);
    print_marginal("  After DFT-bulk", TARGET);

    /* ═══════════════════════════════════════════════════════════════════
     *  TEST 2: PROMOTE + GATE — Apply DFT to q42 specifically
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n\n  ═══ TEST 2: PROMOTE + GATE quhit %lu ═══\n", (unsigned long)TARGET);

    /* Save the marginals of another quhit (q7) as a control */
    double control_before[D];
    marginal_probs(7, control_before);

    /* Apply DFT to q42 — this promotes it into addr[] */
    apply_dft_quhit(&eng, 0, TARGET, D);
    print_state("After DFT on q42 (promoted)", TARGET);
    print_marginal("  q42 after gate", TARGET);
    print_marginal("  q7 (control, should be unchanged)", 7);

    double control_after[D];
    marginal_probs(7, control_after);
    int control_ok = 1;
    for (int i = 0; i < D; i++) {
        if (fabs(control_before[i] - control_after[i]) > 1e-6)
            control_ok = 0;
    }
    printf("  Control quhit (q7) unchanged? %s\n",
           control_ok ? "✓ YES" : "✗ NO — gate on q42 affected q7!");

    /* Save q42's marginal after gating */
    double q42_after_gate[D];
    marginal_probs(TARGET, q42_after_gate);

    /* ═══════════════════════════════════════════════════════════════════
     *  TEST 3: RELEASE — Remove q42 from addr[] (forget the promotion)
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n\n  ═══ TEST 3: RELEASE q42 from addr[] ═══\n");
    printf("  Removing q42 from all entry addr[] arrays...\n");
    release_quhit(TARGET);
    print_state("After release (q42 removed from addr[])", TARGET);
    print_marginal("  q42 after release (now bulk-derived)", TARGET);

    double q42_after_release[D];
    marginal_probs(TARGET, q42_after_release);

    int release_changed = 0;
    for (int i = 0; i < D; i++) {
        if (fabs(q42_after_gate[i] - q42_after_release[i]) > 1e-6)
            release_changed = 1;
    }
    printf("\n  Did release change q42's marginal? %s\n",
           release_changed ? "✓ YES — release destroyed the gate information"
                          : "✗ NO — marginal survived (but may be coincidental)");

    /* ═══════════════════════════════════════════════════════════════════
     *  TEST 4: RE-GRAB — Apply DFT to q42 again after release
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n\n  ═══ TEST 4: RE-GRAB q42 (apply DFT again) ═══\n");
    apply_dft_quhit(&eng, 0, TARGET, D);
    print_state("After re-DFT on released q42", TARGET);
    print_marginal("  q42 after re-grab + gate", TARGET);

    double q42_after_regrab[D];
    marginal_probs(TARGET, q42_after_regrab);

    /* Compare: is re-grab(release(gate(q42))) the same as gate(gate(q42))? */

    /* ═══════════════════════════════════════════════════════════════════
     *  TEST 5: CONTROL — Apply DFT twice WITHOUT release
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n\n  ═══ TEST 5: CONTROL — DFT twice WITHOUT release ═══\n");
    eng.num_quhit_regs = 0;
    init_quhit_register(&eng, 0, N, D);
    eng.quhit_regs[0].bulk_rule = 1;
    entangle_all_quhits(&eng, 0);
    apply_dft_quhit(&eng, 0, TARGET, D);
    apply_dft_quhit(&eng, 0, TARGET, D);
    print_state("DFT² on q42 (no release)", TARGET);
    print_marginal("  q42 after DFT²", TARGET);

    double q42_dft_squared[D];
    marginal_probs(TARGET, q42_dft_squared);

    /* Compare DFT²(q42) vs DFT(release(DFT(q42))) */
    int same_as_dft2 = 1;
    for (int i = 0; i < D; i++) {
        if (fabs(q42_after_regrab[i] - q42_dft_squared[i]) > 1e-6)
            same_as_dft2 = 0;
    }

    /* ═══════════════════════════════════════════════════════════════════
     *  TEST 6: RELEASE + DIFFERENT QUHIT — Does releasing q42 affect q99?
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n\n  ═══ TEST 6: Does releasing q42 affect an unrelated quhit (q99)? ═══\n");
    eng.num_quhit_regs = 0;
    init_quhit_register(&eng, 0, N, D);
    eng.quhit_regs[0].bulk_rule = 1;
    entangle_all_quhits(&eng, 0);

    /* Gate both q42 and q99 */
    apply_dft_quhit(&eng, 0, TARGET, D);
    apply_dft_quhit(&eng, 0, 99, D);

    double q99_before_release[D];
    marginal_probs(99, q99_before_release);
    print_marginal("  q99 before releasing q42", 99);

    /* Release q42 only */
    release_quhit(TARGET);

    double q99_after_release[D];
    marginal_probs(99, q99_after_release);
    print_marginal("  q99 after releasing q42", 99);

    int q99_changed = 0;
    for (int i = 0; i < D; i++) {
        if (fabs(q99_before_release[i] - q99_after_release[i]) > 1e-6)
            q99_changed = 1;
    }
    printf("  Releasing q42 changed q99? %s\n",
           q99_changed ? "✓ YES — entanglement detected: releasing one affects the other!"
                      : "✗ NO — q99 is independent of q42's addr[] status");

    /* ═══════════════════════════════════════════════════════════════════
     *  SUMMARY
     * ═══════════════════════════════════════════════════════════════════ */
    printf("\n\n");
    printf("  ╔═══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║  RESULTS                                                         ║\n");
    printf("  ║                                                                   ║\n");
    printf("  ║  Test 2: Gate q42 → promoted to addr[]                           ║\n");
    printf("  ║  Test 3: Release q42 → marginal %s                       ║\n",
           release_changed ? "CHANGED (info lost)" : "SAME   (survived) ");
    printf("  ║  Test 4: Re-grab → DFT(release(DFT)) %s DFT²              ║\n",
           same_as_dft2 ? "==" : "!=");
    printf("  ║  Test 6: Release q42 → q99 %s                            ║\n",
           q99_changed ? "CHANGED" : "UNCHANGED");
    printf("  ║                                                                   ║\n");

    if (release_changed) {
        printf("  ║  CONCLUSION: The Hilbert space stores gate info in addr[].       ║\n");
        printf("  ║  Releasing a quhit from addr[] erases its individual state.      ║\n");
        printf("  ║  It reverts to bulk-derived: V(k) = (bulk + k) %% D.             ║\n");
        printf("  ║  The Hilbert space DOES NOT remember the gate after release.     ║\n");
    } else {
        printf("  ║  CONCLUSION: The Hilbert space retained the gate info despite    ║\n");
        printf("  ║  the quhit being removed from addr[]. The entries themselves     ║\n");
        printf("  ║  encode the quantum state through amplitude structure.           ║\n");
    }

    printf("  ╚═══════════════════════════════════════════════════════════════════╝\n\n");

    engine_destroy(&eng);
    return 0;
}
