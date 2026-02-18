/*
 * mirror_supremacy.c — Entanglement Persistence Verification
 *
 * EXPERIMENT:
 *   Run two versions of the same 5-qudit sliding window circuit:
 *
 *   (A) QUANTUM carry: inject the full 6-amplitude superposition
 *       |ψ_carry⟩ = Σ_v α_v |v⟩  (phases and magnitudes preserved)
 *
 *   (B) CLASSICAL carry: measure the carry qudit, inject a single
 *       definite value |v_measured⟩ (no superposition, no phases)
 *
 *   If entanglement persists:
 *     • The marginal distributions will DIFFER
 *     • The pair correlations will DIFFER
 *     • The mutual information will be HIGHER in (A)
 *
 *   If entanglement does NOT persist:
 *     • Both versions produce identical statistics
 *     • The carry was effectively classical all along
 *
 * BUILD:
 *   gcc -O2 -I. -std=c11 -D_GNU_SOURCE \
 *       -o mirror_supremacy mirror_supremacy.c hexstate_engine.o bigint.o -lm
 */

#include "hexstate_engine.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>

#define D         6
#define W         5
#define N_QUHITS  100000000000000ULL

/* ── Silence engine prints ────────────────────────────────────── */
static int saved_fd = -1;
static void hush(void) {
    fflush(stdout);
    saved_fd = dup(STDOUT_FILENO);
    int devnull = open("/dev/null", O_WRONLY);
    dup2(devnull, STDOUT_FILENO);
    close(devnull);
}
static void unhush(void) {
    if (saved_fd >= 0) {
        fflush(stdout);
        dup2(saved_fd, STDOUT_FILENO);
        close(saved_fd);
        saved_fd = -1;
    }
}

/* ═══════════════════════════════════════════════════════════════════
 *  CARRY STATE
 * ═══════════════════════════════════════════════════════════════════ */

typedef struct {
    Complex amp[D];
    int live;
} CarryState;

static CarryState extract_carry(HexStateEngine *eng, uint64_t q_carry)
{
    CarryState cs;
    memset(&cs, 0, sizeof(cs));

    int r = find_quhit_reg(eng, 0);
    if (r < 0) return cs;

    uint32_t nz = eng->quhit_regs[r].num_nonzero;
    uint32_t dim = eng->quhit_regs[r].dim;

    for (uint32_t e = 0; e < nz; e++) {
        QuhitBasisEntry *entry = &eng->quhit_regs[r].entries[e];
        uint32_t v = UINT32_MAX;
        for (uint8_t i = 0; i < entry->num_addr; i++) {
            if (entry->addr[i].quhit_idx == q_carry) {
                v = entry->addr[i].value;
                break;
            }
        }
        if (v == UINT32_MAX)
            v = (uint32_t)((entry->bulk_value + q_carry) % dim);
        if (v < D) {
            cs.amp[v].real += entry->amplitude.real;
            cs.amp[v].imag += entry->amplitude.imag;
        }
    }

    int nonzero = 0;
    for (int v = 0; v < D; v++) {
        double p = cs.amp[v].real * cs.amp[v].real +
                   cs.amp[v].imag * cs.amp[v].imag;
        if (p > 1e-15) nonzero++;
    }
    cs.live = (nonzero > 1) ? 1 : 0;

    return cs;
}

/* Classical carry: collapse to a definite value via Born rule */
static CarryState collapse_carry(const CarryState *quantum)
{
    CarryState cs;
    memset(&cs, 0, sizeof(cs));

    /* Compute probabilities */
    double probs[D];
    double total = 0.0;
    for (int v = 0; v < D; v++) {
        probs[v] = quantum->amp[v].real * quantum->amp[v].real +
                   quantum->amp[v].imag * quantum->amp[v].imag;
        total += probs[v];
    }
    if (total < 1e-15) { cs.amp[0].real = 1.0; cs.live = 0; return cs; }

    /* Born rule sampling */
    double r = (double)(rand() % 1000000) / 1000000.0;
    double cumul = 0.0;
    int result = D - 1;
    for (int v = 0; v < D; v++) {
        cumul += probs[v] / total;
        if (cumul >= r) { result = v; break; }
    }

    cs.amp[result].real = 1.0;
    cs.live = 0;
    return cs;
}

static int inject_carry(HexStateEngine *eng, uint64_t q_carry,
                         const CarryState *cs)
{
    int r = find_quhit_reg(eng, 0);
    if (r < 0) return 0;

    eng->quhit_regs[r].num_nonzero = 0;
    eng->quhit_regs[r].collapsed = 0;

    int injected = 0;
    for (uint32_t v = 0; v < D; v++) {
        double p = cs->amp[v].real * cs->amp[v].real +
                   cs->amp[v].imag * cs->amp[v].imag;
        if (p < 1e-15) continue;

        QuhitBasisEntry e;
        memset(&e, 0, sizeof(e));
        e.bulk_value = v;
        e.num_addr = 1;
        e.addr[0].quhit_idx = q_carry;
        e.addr[0].value = v;
        e.amplitude = cs->amp[v];

        uint32_t nz = eng->quhit_regs[r].num_nonzero;
        if (nz < MAX_QUHIT_HILBERT_ENTRIES) {
            eng->quhit_regs[r].entries[nz] = e;
            eng->quhit_regs[r].num_nonzero = nz + 1;
            injected++;
        }
    }
    return injected;
}

/* ═══════════════════════════════════════════════════════════════════
 *  WINDOW EXECUTION
 * ═══════════════════════════════════════════════════════════════════ */

static void run_window_circuit(HexStateEngine *eng, int window_idx,
                                const CarryState *carry_in,
                                int out_digits[4],
                                CarryState *carry_out,
                                int *peak_entries, int *carry_entries)
{
    int r = find_quhit_reg(eng, 0);
    if (r < 0) return;

    uint64_t q[W] = {0, 1, 2, 3, 4};
    *peak_entries = (int)eng->quhit_regs[r].num_nonzero;

    /* Step 1: Inject carry or DFT */
    if (window_idx > 0 && carry_in->live) {
        inject_carry(eng, q[0], carry_in);
    } else if (window_idx > 0) {
        inject_carry(eng, q[0], carry_in);  /* single definite value */
    } else {
        hush(); apply_dft_quhit(eng, 0, q[0], D); unhush();
    }

    /* Step 2: DFT on q1-q4 */
    for (int i = 1; i < W; i++) {
        hush(); apply_dft_quhit(eng, 0, q[i], D); unhush();
        int ent = (int)eng->quhit_regs[r].num_nonzero;
        if (ent > *peak_entries) *peak_entries = ent;
    }

    /* Step 3: CZ chain */
    for (int i = 0; i < W - 1; i++) {
        hush(); apply_cz_quhits(eng, 0, q[i], 0, q[i + 1]); unhush();
    }

    /* Step 4: DNA on all */
    for (int i = 0; i < W; i++) {
        hush(); apply_dna_quhit(eng, 0, q[i], 1.0, 310.0); unhush();
        int ent = (int)eng->quhit_regs[r].num_nonzero;
        if (ent > *peak_entries) *peak_entries = ent;
    }

    /* Step 5: Second DFT */
    for (int i = 0; i < W; i++) {
        hush(); apply_dft_quhit(eng, 0, q[i], D); unhush();
        int ent = (int)eng->quhit_regs[r].num_nonzero;
        if (ent > *peak_entries) *peak_entries = ent;
    }

    /* Step 6: Measure q0-q3 */
    for (int i = 0; i < 4; i++) {
        hush(); uint64_t m = measure_quhit(eng, 0, q[i]); unhush();
        out_digits[i] = (int)(m % D);
    }
    *carry_entries = (int)eng->quhit_regs[r].num_nonzero;

    /* Step 7: Extract carry */
    *carry_out = extract_carry(eng, q[4]);
}

/* ═══════════════════════════════════════════════════════════════════
 *  Run a chain of windows (quantum or classical carry)
 * ═══════════════════════════════════════════════════════════════════ */

typedef struct {
    int n_windows;
    int marginal[D];
    int pair_hist[D][D];
    double mutual_info;
    int peak_entries;
    int carry_live_count;
    double carry_entries_avg;
    double total_ms;
    /* Carry amplitude profile: average |α_v|² across all windows */
    double carry_prob_avg[D];
    /* Shannon entropy of carry state */
    double carry_entropy_avg;
} RunResult;

static RunResult run_chain(int n_windows, int classical_carry,
                            const char *label)
{
    RunResult rr;
    memset(&rr, 0, sizeof(rr));
    rr.n_windows = n_windows;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    CarryState carry;
    memset(&carry, 0, sizeof(carry));
    carry.live = 0;

    int prev_digit = -1;

    for (int w = 0; w < n_windows; w++) {
        static HexStateEngine eng;
        hush();
        engine_init(&eng);
        init_quhit_register(&eng, 0, N_QUHITS, D);
        unhush();

        eng.quhit_regs[0].bulk_rule = 1;
        hush(); entangle_all_quhits(&eng, 0); unhush();

        CarryState effective_carry = carry;
        if (classical_carry && w > 0 && carry.live) {
            effective_carry = collapse_carry(&carry);
        }

        int digits[4];
        CarryState next_carry;
        int peak = 0, carry_ent = 0;

        run_window_circuit(&eng, w, &effective_carry, digits,
                           &next_carry, &peak, &carry_ent);

        hush(); engine_destroy(&eng); unhush();

        carry = next_carry;

        if (peak > rr.peak_entries) rr.peak_entries = peak;
        if (next_carry.live) rr.carry_live_count++;
        rr.carry_entries_avg += carry_ent;

        /* Carry amplitude profile */
        double carry_total = 0.0;
        for (int v = 0; v < D; v++) {
            double p = next_carry.amp[v].real * next_carry.amp[v].real +
                       next_carry.amp[v].imag * next_carry.amp[v].imag;
            rr.carry_prob_avg[v] += p;
            carry_total += p;
        }
        /* Carry entropy */
        if (carry_total > 1e-15) {
            double entropy = 0.0;
            for (int v = 0; v < D; v++) {
                double p = (next_carry.amp[v].real * next_carry.amp[v].real +
                            next_carry.amp[v].imag * next_carry.amp[v].imag) / carry_total;
                if (p > 1e-15) entropy -= p * log2(p);
            }
            rr.carry_entropy_avg += entropy;
        }

        for (int i = 0; i < 4; i++) {
            rr.marginal[digits[i] % D]++;
            if (prev_digit >= 0)
                rr.pair_hist[prev_digit][digits[i] % D]++;
            prev_digit = digits[i] % D;
        }

        if (w < 5 || w == n_windows - 1) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double ms = (t1.tv_sec - t0.tv_sec) * 1000.0 +
                        (t1.tv_nsec - t0.tv_nsec) / 1e6;
            printf("    %s W%2d: [%d,%d,%d,%d] carry=%d ent %s  peak=%d  [%.0fms]\n",
                   label, w,
                   digits[0], digits[1], digits[2], digits[3],
                   carry_ent, next_carry.live ? "LIVE" : "dead",
                   peak, ms);
        } else if (w == 5) {
            printf("    %s ...\n", label);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    rr.total_ms = (t1.tv_sec - t0.tv_sec) * 1000.0 +
                  (t1.tv_nsec - t0.tv_nsec) / 1e6;
    rr.carry_entries_avg /= n_windows;
    rr.carry_entropy_avg /= n_windows;
    for (int v = 0; v < D; v++)
        rr.carry_prob_avg[v] /= n_windows;

    /* Mutual information */
    int n_pairs = n_windows * 4 - 1;
    double m0[D], m1[D];
    memset(m0, 0, sizeof(m0));
    memset(m1, 0, sizeof(m1));
    for (int a = 0; a < D; a++)
        for (int b = 0; b < D; b++) {
            double p = (double)rr.pair_hist[a][b] / n_pairs;
            m0[a] += p;
            m1[b] += p;
        }

    double mi = 0.0;
    for (int a = 0; a < D; a++)
        for (int b = 0; b < D; b++) {
            double p_ab = (double)rr.pair_hist[a][b] / n_pairs;
            if (p_ab > 1e-10 && m0[a] > 1e-10 && m1[b] > 1e-10)
                mi += p_ab * log2(p_ab / (m0[a] * m1[b]));
        }
    rr.mutual_info = mi;

    return rr;
}

/* ═══════════════════════════════════════════════════════════════════
 *  COMPARISON PRINTER
 * ═══════════════════════════════════════════════════════════════════ */

static void print_comparison(const RunResult *quantum, const RunResult *classical,
                              int n_windows)
{
    int n_measured = n_windows * 4;

    printf("  ╔═══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ENTANGLEMENT PERSISTENCE TEST                                  ║\n");
    printf("  ║  %d windows × 4 measured = %d qudits                           ║\n",
           n_windows, n_measured);
    printf("  ╠═══════════════════════════════════════════════════════════════════╣\n");
    printf("  ║                                                                 ║\n");

    /* Marginal comparison */
    printf("  ║  MARGINAL DISTRIBUTION:                                         ║\n");
    printf("  ║  ─────────────────────────────────────────────────────────────   ║\n");
    printf("  ║  Digit │ Quantum  │ Classical │  Δ       │ Verdict              ║\n");
    printf("  ║  ──────┼──────────┼───────────┼──────────┼──────────            ║\n");

    double max_delta = 0.0;
    for (int k = 0; k < D; k++) {
        double qp = (double)quantum->marginal[k] / n_measured;
        double cp = (double)classical->marginal[k] / n_measured;
        double delta = fabs(qp - cp);
        if (delta > max_delta) max_delta = delta;
        const char *verdict = delta > 0.02 ? "DIFFERS" : "~same";
        printf("  ║  |%d⟩   │  %.4f  │  %.4f   │ %+.4f  │ %s            ║\n",
               k, qp, cp, qp - cp, verdict);
    }
    printf("  ║                                                                 ║\n");

    /* Mutual information */
    printf("  ║  MUTUAL INFORMATION:                                            ║\n");
    printf("  ║  ─────────────────────────────────────────────────────────────   ║\n");
    printf("  ║  Quantum:   I = %.6f bits                                    ║\n",
           quantum->mutual_info);
    printf("  ║  Classical: I = %.6f bits                                    ║\n",
           classical->mutual_info);
    printf("  ║  Δ(I) = %+.6f bits                                          ║\n",
           quantum->mutual_info - classical->mutual_info);
    printf("  ║                                                                 ║\n");

    /* Carry state analysis */
    printf("  ║  CARRY QUDIT STATE:                                             ║\n");
    printf("  ║  ─────────────────────────────────────────────────────────────   ║\n");
    printf("  ║  Quantum:   avg entropy = %.4f bits (max = %.3f)             ║\n",
           quantum->carry_entropy_avg, log2(D));
    printf("  ║  Classical: avg entropy = %.4f bits                           ║\n",
           classical->carry_entropy_avg);
    printf("  ║                                                                 ║\n");
    printf("  ║  Carry amplitude profile |α_v|² (quantum only):                ║\n");
    printf("  ║  ");
    for (int v = 0; v < D; v++)
        printf("v=%d:%.3f ", v, quantum->carry_prob_avg[v]);
    printf("  ║\n");
    printf("  ║                                                                 ║\n");

    /* Peak entries */
    printf("  ║  PEAK ENTRIES:                                                  ║\n");
    printf("  ║  Quantum:   %d / 7776                                         ║\n",
           quantum->peak_entries);
    printf("  ║  Classical: %d / 7776                                         ║\n",
           classical->peak_entries);
    printf("  ║  LIVE carry windows: Q=%d/%d  C=%d/%d                        ║\n",
           quantum->carry_live_count, n_windows,
           classical->carry_live_count, n_windows);
    printf("  ║                                                                 ║\n");

    /* Verdict */
    int entangled = (max_delta > 0.02) ||
                    (fabs(quantum->mutual_info - classical->mutual_info) > 0.01);
    printf("  ╠═══════════════════════════════════════════════════════════════════╣\n");
    if (entangled) {
        printf("  ║                                                                 ║\n");
        printf("  ║  ✓ ENTANGLEMENT PERSISTS ACROSS WINDOWS                        ║\n");
        printf("  ║                                                                 ║\n");
        printf("  ║  The quantum carry produces DIFFERENT statistics from the       ║\n");
        printf("  ║  classical carry. The carry qudit's phase information —         ║\n");
        printf("  ║  which only exists in superposition — affects the output.       ║\n");
        printf("  ║  This is impossible without genuine quantum state transfer.     ║\n");
    } else {
        printf("  ║                                                                 ║\n");
        printf("  ║  ✗ NO DIFFERENCE DETECTED                                      ║\n");
        printf("  ║                                                                 ║\n");
        printf("  ║  Quantum and classical carry produce identical statistics.      ║\n");
        printf("  ║  The carry superposition does not affect the output —           ║\n");
        printf("  ║  either decoherence or the circuit structure makes the          ║\n");
        printf("  ║  phases irrelevant.                                             ║\n");
    }
    printf("  ║                                                                 ║\n");
    printf("  ╚═══════════════════════════════════════════════════════════════════╝\n\n");
}

/* ═══════════════════════════════════════════════════════════════════ */

int main(void)
{
    srand(42);  /* reproducible Born sampling for classical carry */

    printf("\n");
    printf("  ╔═══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║                                                                 ║\n");
    printf("  ║  E N T A N G L E M E N T   P E R S I S T E N C E   T E S T   ║\n");
    printf("  ║                                                                 ║\n");
    printf("  ║  QUANTUM  carry: inject full amplitude vector (6 complex)      ║\n");
    printf("  ║  CLASSICAL carry: measure carry, inject definite |v⟩           ║\n");
    printf("  ║                                                                 ║\n");
    printf("  ║  If distributions DIFFER → entanglement persists               ║\n");
    printf("  ║  If distributions MATCH → carry was effectively classical      ║\n");
    printf("  ║                                                                 ║\n");
    printf("  ║  Circuit: H⊗5 → CZ-chain → DNA⊗5 → H⊗5 → Measure 4          ║\n");
    printf("  ║  Window stride: 4 measured + 1 carry                           ║\n");
    printf("  ╚═══════════════════════════════════════════════════════════════════╝\n\n");

    /* ══════ TEST 1: 20 windows ══════ */
    int N = 20;
    printf("  ════════════════════════════════════════════════════════════════════\n");
    printf("  ══  Running %d windows with QUANTUM carry...                    ══\n",N);
    printf("  ════════════════════════════════════════════════════════════════════\n\n");
    RunResult quantum = run_chain(N, 0, "[Q]");

    printf("\n");
    printf("  ════════════════════════════════════════════════════════════════════\n");
    printf("  ══  Running %d windows with CLASSICAL carry...                  ══\n",N);
    printf("  ════════════════════════════════════════════════════════════════════\n\n");
    RunResult classical = run_chain(N, 1, "[C]");

    printf("\n");
    print_comparison(&quantum, &classical, N);

    /* ══════ TEST 2: 50 windows (more statistics) ══════ */
    N = 50;
    printf("  ════════════════════════════════════════════════════════════════════\n");
    printf("  ══  Running %d windows with QUANTUM carry...                    ══\n",N);
    printf("  ════════════════════════════════════════════════════════════════════\n\n");
    RunResult q2 = run_chain(N, 0, "[Q]");

    printf("\n");
    printf("  ════════════════════════════════════════════════════════════════════\n");
    printf("  ══  Running %d windows with CLASSICAL carry...                  ══\n",N);
    printf("  ════════════════════════════════════════════════════════════════════\n\n");
    RunResult c2 = run_chain(N, 1, "[C]");

    printf("\n");
    print_comparison(&q2, &c2, N);

    return 0;
}
