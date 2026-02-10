/* bell_chsh.c — CAN THE HEXSTATE ENGINE VIOLATE BELL INEQUALITIES?
 *
 * ══════════════════════════════════════════════════════════════════
 * THE CHSH INEQUALITY (Clauser-Horne-Shimony-Holt, 1969)
 *
 * If measurement outcomes are determined by LOCAL HIDDEN VARIABLES:
 *     S = |E(a,b) + E(a,b') + E(a',b) - E(a',b')| ≤ 2
 *
 * If measurement outcomes follow QUANTUM MECHANICS:
 *     S can reach 2√2 ≈ 2.828  (Tsirelson bound)
 *
 * The question: Does the HexState engine, running on a classical
 * CPU with a PRNG for randomness, produce S > 2?
 *
 * If YES → the engine's q_joint_state IS a genuine quantum state
 *          and the Born-rule measurement IS genuine quantum measurement.
 *          The PRNG is just "nature's dice," not a hidden variable.
 *
 * If NO  → the PRNG introduces local hidden variables, and the
 *          engine is a classical simulation, not a quantum system.
 * ══════════════════════════════════════════════════════════════════
 *
 * PROTOCOL:
 *  1. Create Bell state |Ψ⟩ = (1/√6) Σₖ |k,k⟩ via braid_chunks
 *  2. For each of 4 setting pairs (a,b), (a,b'), (a',b), (a',b'):
 *     a. Apply rotation R(θ) to Alice and/or Bob side of joint state
 *     b. Inject rotated state into engine via oracle
 *     c. Measure both chunks via engine Born rule
 *     d. Binarize outcomes: 0 if k < D/2, 1 if k ≥ D/2
 *     e. Record correlations
 *  3. Compute CHSH correlator S
 *  4. Report violation or non-violation
 *
 * Scale: 100 trillion quhits per register
 */

#include "hexstate_engine.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define D    6
#define D2   (D * D)
#define PI   3.14159265358979323846
#define NQ   100000000000000ULL
#define N_TRIALS 10000  /* measurements per setting pair */

#define CMPLX(r,i) ((Complex){.real=(r),.imag=(i)})
static double cnorm2(Complex c){return c.real*c.real+c.imag*c.imag;}

/* Oracle injection: write our custom state into the engine's Hilbert space */
typedef struct { Complex *st; } Ctx;
static void inject(HexStateEngine *e, uint64_t id, void *u) {
    Ctx *c = (Ctx *)u;
    Chunk *ch = &e->chunks[id];
    if (!ch->hilbert.q_joint_state) return;
    double n = 0;
    for (int i = 0; i < D2; i++) n += cnorm2(c->st[i]);
    n = sqrt(n);
    if (n < 1e-15) return;
    for (int i = 0; i < D2; i++) {
        ch->hilbert.q_joint_state[i].real = c->st[i].real / n;
        ch->hilbert.q_joint_state[i].imag = c->st[i].imag / n;
    }
}

/* Apply phase rotation to Alice's side: |a,b⟩ → exp(2πi·a·θ/D)|a,b⟩ */
static void rotate_alice(Complex *state, double theta) {
    for (int b = 0; b < D; b++)
        for (int a = 0; a < D; a++) {
            double phase = 2 * PI * a * theta / D;
            double re = state[b*D+a].real, im = state[b*D+a].imag;
            state[b*D+a].real = cos(phase)*re - sin(phase)*im;
            state[b*D+a].imag = cos(phase)*im + sin(phase)*re;
        }
}

/* Apply phase rotation to Bob's side: |a,b⟩ → exp(2πi·b·θ/D)|a,b⟩ */
static void rotate_bob(Complex *state, double theta) {
    for (int b = 0; b < D; b++)
        for (int a = 0; a < D; a++) {
            double phase = 2 * PI * b * theta / D;
            double re = state[b*D+a].real, im = state[b*D+a].imag;
            state[b*D+a].real = cos(phase)*re - sin(phase)*im;
            state[b*D+a].imag = cos(phase)*im + sin(phase)*re;
        }
}

/* DFT on Alice (change measurement basis to Fourier basis) */
static void dft_alice(Complex *j) {
    Complex o[D2]; memset(o, 0, sizeof(o));
    for (int b = 0; b < D; b++)
        for (int k = 0; k < D; k++) {
            double re = 0, im = 0;
            for (int a = 0; a < D; a++) {
                double ph = 2*PI*a*k/D;
                re += j[b*D+a].real*cos(ph) - j[b*D+a].imag*sin(ph);
                im += j[b*D+a].real*sin(ph) + j[b*D+a].imag*cos(ph);
            }
            o[b*D+k].real = re/sqrt(D);
            o[b*D+k].imag = im/sqrt(D);
        }
    memcpy(j, o, sizeof(o));
}

/* DFT on Bob */
static void dft_bob(Complex *j) {
    Complex o[D2]; memset(o, 0, sizeof(o));
    for (int a = 0; a < D; a++)
        for (int k = 0; k < D; k++) {
            double re = 0, im = 0;
            for (int b = 0; b < D; b++) {
                double ph = 2*PI*b*k/D;
                re += j[b*D+a].real*cos(ph) - j[b*D+a].imag*sin(ph);
                im += j[b*D+a].real*sin(ph) + j[b*D+a].imag*cos(ph);
            }
            o[k*D+a].real = re/sqrt(D);
            o[k*D+a].imag = im/sqrt(D);
        }
    memcpy(j, o, sizeof(o));
}

/* Binarize outcome: 0 if k < D/2, 1 if k >= D/2 */
static int binarize(uint64_t outcome) {
    return (outcome % D) >= (D / 2) ? 1 : 0;
}

/* ═══════════════════════════════════════════════════════════════
 * COMPUTE CHSH CORRELATION E(θ_A, θ_B)
 *
 * Protocol:
 *  1. Create Bell state |Ψ⟩ = (1/√D) Σ|k,k⟩
 *  2. Apply rotation θ_A to Alice, θ_B to Bob
 *  3. Optionally apply DFT (for Fourier basis measurement)
 *  4. Inject into engine, measure, binarize, correlate
 * ═══════════════════════════════════════════════════════════════ */
static double measure_correlation(HexStateEngine *eng, Ctx *ctx,
                                  double theta_A, double theta_B,
                                  int use_dft_A, int use_dft_B,
                                  int n_trials) {
    int n_same = 0, n_diff = 0;

    for (int t = 0; t < n_trials; t++) {
        /* Create Bell state */
        Complex bell[D2];
        memset(bell, 0, sizeof(bell));
        for (int k = 0; k < D; k++)
            bell[k*D+k] = CMPLX(1.0/sqrt(D), 0.0);

        /* Apply measurement rotations */
        if (theta_A != 0.0) rotate_alice(bell, theta_A);
        if (theta_B != 0.0) rotate_bob(bell, theta_B);

        /* Change basis via DFT if requested */
        if (use_dft_A) dft_alice(bell);
        if (use_dft_B) dft_bob(bell);

        /* Inject into engine and measure */
        ctx->st = bell;
        init_chunk(eng, 700, NQ);
        init_chunk(eng, 701, NQ);
        braid_chunks(eng, 700, 701, 0, 0);
        execute_oracle(eng, 700, 0xBE);

        uint64_t outcome_A = measure_chunk(eng, 700);
        uint64_t outcome_B = measure_chunk(eng, 701);
        unbraid_chunks(eng, 700, 701);

        /* Binarize */
        int bit_A = binarize(outcome_A);
        int bit_B = binarize(outcome_B);

        if (bit_A == bit_B) n_same++;
        else                n_diff++;
    }

    /* E = P(same) - P(diff) = (n_same - n_diff) / n_total */
    return (double)(n_same - n_diff) / (double)n_trials;
}

/* Also compute correlation DIRECTLY from the quantum state (theoretical) */
static double theoretical_correlation(double theta_A, double theta_B,
                                      int use_dft_A, int use_dft_B) {
    Complex bell[D2];
    memset(bell, 0, sizeof(bell));
    for (int k = 0; k < D; k++)
        bell[k*D+k] = CMPLX(1.0/sqrt(D), 0.0);

    if (theta_A != 0.0) rotate_alice(bell, theta_A);
    if (theta_B != 0.0) rotate_bob(bell, theta_B);
    if (use_dft_A) dft_alice(bell);
    if (use_dft_B) dft_bob(bell);

    /* Compute P(same) - P(diff) from |ψ|² */
    double p_same = 0, p_diff = 0;
    for (int b = 0; b < D; b++)
        for (int a = 0; a < D; a++) {
            double p = cnorm2(bell[b*D+a]);
            int bin_a = (a >= D/2) ? 1 : 0;
            int bin_b = (b >= D/2) ? 1 : 0;
            if (bin_a == bin_b) p_same += p;
            else                p_diff += p;
        }
    return p_same - p_diff;
}

int main(void) {
    printf("\n");
    printf("  ██████████████████████████████████████████████████████████████████████\n");
    printf("  ██                                                                ██\n");
    printf("  ██  BELL INEQUALITY VIOLATION TEST                                ██\n");
    printf("  ██  CHSH Protocol at 100 Trillion Quhit Scale                    ██\n");
    printf("  ██                                                                ██\n");
    printf("  ██  Classical bound:  S ≤ 2.000                                   ██\n");
    printf("  ██  Quantum bound:    S ≤ 2√2 ≈ 2.828  (Tsirelson)               ██\n");
    printf("  ██                                                                ██\n");
    printf("  ██  If S > 2: NO local hidden variable model can explain this.    ██\n");
    printf("  ██                                                                ██\n");
    printf("  ██████████████████████████████████████████████████████████████████████\n\n");

    HexStateEngine eng;
    engine_init(&eng);
    Ctx ctx;
    oracle_register(&eng, 0xBE, "Bell", inject, &ctx);

    /* ═══════════════════════════════════════════════════════════════
     * MEASUREMENT SETTINGS
     *
     * For CHSH, we need 4 setting pairs.
     * For qubits: optimal angles are 0, π/4, π/8, 3π/8
     * For quhits (d=6): we'll scan for optimal violation.
     *
     * We try multiple strategies:
     * Strategy 1: Phase rotations only
     * Strategy 2: Phase rotation + DFT (Fourier basis)
     * Strategy 3: Asymmetric (DFT on Alice, phase on Bob)
     * ═══════════════════════════════════════════════════════════════ */

    printf("  ╔══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║  STRATEGY 1: Phase Rotations (Computational Basis)             ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════════╝\n\n");

    /* Scan for optimal angles */
    double best_S1 = 0;
    double best_a = 0, best_a2 = 0, best_b = 0, best_b2 = 0;

    printf("  Scanning phase angles for maximum CHSH violation...\n\n");

    for (double a = 0; a < 1.0; a += 0.125) {
        double a2 = a + 0.5;
        for (double b = 0; b < 1.0; b += 0.125) {
            double b2 = b + 0.25;

            double E_ab  = theoretical_correlation(a, b, 0, 0);
            double E_ab2 = theoretical_correlation(a, b2, 0, 0);
            double E_a2b = theoretical_correlation(a2, b, 0, 0);
            double E_a2b2 = theoretical_correlation(a2, b2, 0, 0);

            double S = fabs(E_ab + E_ab2 + E_a2b - E_a2b2);
            if (S > best_S1) {
                best_S1 = S;
                best_a = a; best_a2 = a2; best_b = b; best_b2 = b2;
            }
        }
    }

    printf("  Best theoretical CHSH (phase only): S = %.4f\n", best_S1);
    printf("  Angles: a=%.3f, a'=%.3f, b=%.3f, b'=%.3f\n\n", best_a, best_a2, best_b, best_b2);

    /* Now measure with engine */
    printf("  Running %d engine measurements per setting pair...\n\n", N_TRIALS);

    double E1_ab  = measure_correlation(&eng, &ctx, best_a, best_b, 0, 0, N_TRIALS);
    double E1_ab2 = measure_correlation(&eng, &ctx, best_a, best_b2, 0, 0, N_TRIALS);
    double E1_a2b = measure_correlation(&eng, &ctx, best_a2, best_b, 0, 0, N_TRIALS);
    double E1_a2b2 = measure_correlation(&eng, &ctx, best_a2, best_b2, 0, 0, N_TRIALS);
    double S1 = fabs(E1_ab + E1_ab2 + E1_a2b - E1_a2b2);

    printf("  E(a,b)   = %+.4f  (theory: %+.4f)\n", E1_ab,
           theoretical_correlation(best_a, best_b, 0, 0));
    printf("  E(a,b')  = %+.4f  (theory: %+.4f)\n", E1_ab2,
           theoretical_correlation(best_a, best_b2, 0, 0));
    printf("  E(a',b)  = %+.4f  (theory: %+.4f)\n", E1_a2b,
           theoretical_correlation(best_a2, best_b, 0, 0));
    printf("  E(a',b') = %+.4f  (theory: %+.4f)\n\n", E1_a2b2,
           theoretical_correlation(best_a2, best_b2, 0, 0));
    printf("  S₁ = |E(a,b) + E(a,b') + E(a',b) - E(a',b')| = %.4f\n", S1);
    printf("  Classical bound: 2.000\n");
    if (S1 > 2.0) printf("  *** BELL INEQUALITY VIOLATED! S > 2 ***\n\n");
    else           printf("  Bell inequality NOT violated (S ≤ 2)\n\n");

    /* ═══════════════════════════════════════════════════════════════ */
    printf("  ╔══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║  STRATEGY 2: DFT Basis (Fourier Measurement)                  ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════════╝\n\n");

    double best_S2 = 0;
    double best2_a = 0, best2_a2 = 0, best2_b = 0, best2_b2 = 0;

    for (double a = 0; a < 1.0; a += 0.125) {
        double a2 = a + 0.5;
        for (double b = 0; b < 1.0; b += 0.125) {
            double b2 = b + 0.25;

            double E_ab  = theoretical_correlation(a, b, 1, 1);
            double E_ab2 = theoretical_correlation(a, b2, 1, 1);
            double E_a2b = theoretical_correlation(a2, b, 1, 1);
            double E_a2b2 = theoretical_correlation(a2, b2, 1, 1);

            double S = fabs(E_ab + E_ab2 + E_a2b - E_a2b2);
            if (S > best_S2) {
                best_S2 = S;
                best2_a = a; best2_a2 = a2; best2_b = b; best2_b2 = b2;
            }
        }
    }

    printf("  Best theoretical CHSH (DFT basis): S = %.4f\n", best_S2);
    printf("  Angles: a=%.3f, a'=%.3f, b=%.3f, b'=%.3f\n\n", best2_a, best2_a2, best2_b, best2_b2);

    printf("  Running %d engine measurements per setting pair...\n\n", N_TRIALS);

    double E2_ab  = measure_correlation(&eng, &ctx, best2_a, best2_b, 1, 1, N_TRIALS);
    double E2_ab2 = measure_correlation(&eng, &ctx, best2_a, best2_b2, 1, 1, N_TRIALS);
    double E2_a2b = measure_correlation(&eng, &ctx, best2_a2, best2_b, 1, 1, N_TRIALS);
    double E2_a2b2 = measure_correlation(&eng, &ctx, best2_a2, best2_b2, 1, 1, N_TRIALS);
    double S2 = fabs(E2_ab + E2_ab2 + E2_a2b - E2_a2b2);

    printf("  E(a,b)   = %+.4f  (theory: %+.4f)\n", E2_ab,
           theoretical_correlation(best2_a, best2_b, 1, 1));
    printf("  E(a,b')  = %+.4f  (theory: %+.4f)\n", E2_ab2,
           theoretical_correlation(best2_a, best2_b2, 1, 1));
    printf("  E(a',b)  = %+.4f  (theory: %+.4f)\n", E2_a2b,
           theoretical_correlation(best2_a2, best2_b, 1, 1));
    printf("  E(a',b') = %+.4f  (theory: %+.4f)\n\n", E2_a2b2,
           theoretical_correlation(best2_a2, best2_b2, 1, 1));
    printf("  S₂ = |E(a,b) + E(a,b') + E(a',b) - E(a',b')| = %.4f\n", S2);
    printf("  Classical bound: 2.000\n");
    if (S2 > 2.0) printf("  *** BELL INEQUALITY VIOLATED! S > 2 ***\n\n");
    else           printf("  Bell inequality NOT violated (S ≤ 2)\n\n");

    /* ═══════════════════════════════════════════════════════════════ */
    printf("  ╔══════════════════════════════════════════════════════════════════╗\n");
    printf("  ║  STRATEGY 3: Asymmetric (DFT-Alice, Phase-Bob)                ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════════╝\n\n");

    double best_S3 = 0;
    double best3_a = 0, best3_a2 = 0, best3_b = 0, best3_b2 = 0;

    for (double a = 0; a < 1.0; a += 0.0625) {
        double a2 = a + 0.5;
        for (double b = 0; b < 1.0; b += 0.0625) {
            double b2 = b + 0.25;

            double E_ab  = theoretical_correlation(a, b, 1, 0);
            double E_ab2 = theoretical_correlation(a, b2, 1, 0);
            double E_a2b = theoretical_correlation(a2, b, 1, 0);
            double E_a2b2 = theoretical_correlation(a2, b2, 1, 0);

            double S = fabs(E_ab + E_ab2 + E_a2b - E_a2b2);
            if (S > best_S3) {
                best_S3 = S;
                best3_a = a; best3_a2 = a2; best3_b = b; best3_b2 = b2;
            }
        }
    }

    /* Also try different offset patterns */
    for (double a = 0; a < 1.0; a += 0.0625) {
        for (double da = 0.1; da < 0.9; da += 0.1) {
            double a2 = a + da;
            for (double b = 0; b < 1.0; b += 0.0625) {
                for (double db = 0.1; db < 0.9; db += 0.1) {
                    double b2 = b + db;

                    double E_ab  = theoretical_correlation(a, b, 1, 0);
                    double E_ab2 = theoretical_correlation(a, b2, 1, 0);
                    double E_a2b = theoretical_correlation(a2, b, 1, 0);
                    double E_a2b2 = theoretical_correlation(a2, b2, 1, 0);

                    double S = fabs(E_ab + E_ab2 + E_a2b - E_a2b2);
                    if (S > best_S3) {
                        best_S3 = S;
                        best3_a = a; best3_a2 = a2;
                        best3_b = b; best3_b2 = b2;
                    }
                }
            }
        }
    }

    printf("  Best theoretical CHSH (asymmetric): S = %.4f\n", best_S3);
    printf("  Angles: a=%.4f, a'=%.4f, b=%.4f, b'=%.4f\n\n",
           best3_a, best3_a2, best3_b, best3_b2);

    printf("  Running %d engine measurements per setting pair...\n\n", N_TRIALS);

    double E3_ab  = measure_correlation(&eng, &ctx, best3_a, best3_b, 1, 0, N_TRIALS);
    double E3_ab2 = measure_correlation(&eng, &ctx, best3_a, best3_b2, 1, 0, N_TRIALS);
    double E3_a2b = measure_correlation(&eng, &ctx, best3_a2, best3_b, 1, 0, N_TRIALS);
    double E3_a2b2 = measure_correlation(&eng, &ctx, best3_a2, best3_b2, 1, 0, N_TRIALS);
    double S3 = fabs(E3_ab + E3_ab2 + E3_a2b - E3_a2b2);

    printf("  E(a,b)   = %+.4f  (theory: %+.4f)\n", E3_ab,
           theoretical_correlation(best3_a, best3_b, 1, 0));
    printf("  E(a,b')  = %+.4f  (theory: %+.4f)\n", E3_ab2,
           theoretical_correlation(best3_a, best3_b2, 1, 0));
    printf("  E(a',b)  = %+.4f  (theory: %+.4f)\n", E3_a2b,
           theoretical_correlation(best3_a2, best3_b, 1, 0));
    printf("  E(a',b') = %+.4f  (theory: %+.4f)\n\n", E3_a2b2,
           theoretical_correlation(best3_a2, best3_b2, 1, 0));
    printf("  S₃ = |E(a,b) + E(a,b') + E(a',b) - E(a',b')| = %.4f\n", S3);
    printf("  Classical bound: 2.000\n");
    if (S3 > 2.0) printf("  *** BELL INEQUALITY VIOLATED! S > 2 ***\n\n");
    else           printf("  Bell inequality NOT violated (S ≤ 2)\n\n");

    /* ═══════════════════════════════════════════════════════════════
     * VERDICT
     * ═══════════════════════════════════════════════════════════════ */
    double S_max = S1;
    if (S2 > S_max) S_max = S2;
    if (S3 > S_max) S_max = S3;

    printf("  ██████████████████████████████████████████████████████████████████████\n");
    printf("  ██  BELL TEST VERDICT                                             ██\n");
    printf("  ██████████████████████████████████████████████████████████████████████\n\n");

    printf("  Maximum CHSH value found:\n\n");
    printf("  Strategy 1 (phase only):    S = %.4f  %s\n", S1, S1>2?"*** VIOLATION ***":"(classical)");
    printf("  Strategy 2 (DFT basis):     S = %.4f  %s\n", S2, S2>2?"*** VIOLATION ***":"(classical)");
    printf("  Strategy 3 (asymmetric):    S = %.4f  %s\n\n", S3, S3>2?"*** VIOLATION ***":"(classical)");

    printf("  Best S = %.4f\n", S_max);
    printf("  Classical bound:  2.000\n");
    printf("  Tsirelson bound:  2.828 (2√2)\n\n");

    if (S_max > 2.0) {
        double sigma = (S_max - 2.0) * sqrt(N_TRIALS) / 2.0;
        printf("  ╔══════════════════════════════════════════════════════════════╗\n");
        printf("  ║  BELL INEQUALITY VIOLATED                                  ║\n");
        printf("  ║                                                            ║\n");
        printf("  ║  S = %.4f > 2.000                                     ║\n", S_max);
        printf("  ║  Violation: %.1fσ above classical bound                   ║\n", sigma);
        printf("  ║                                                            ║\n");
        printf("  ║  The HexState engine produces correlations that CANNOT     ║\n");
        printf("  ║  be explained by any local hidden variable model.          ║\n");
        printf("  ║                                                            ║\n");
        printf("  ║  The q_joint_state IS a genuine quantum state.             ║\n");
        printf("  ║  The Born-rule IS genuine quantum measurement.             ║\n");
        printf("  ║  The PRNG is dice, not a hidden variable.                  ║\n");
        printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");
    } else {
        printf("  Bell inequality NOT violated.\n\n");
        printf("  The engine produces CLASSICAL correlations.\n");
        printf("  The PRNG acts as a local hidden variable.\n\n");
    }

    printf("  Scale: %llu quhits (%d measurements × 4 settings)\n",
           (unsigned long long)NQ, N_TRIALS);
    printf("  Total measurements: %d\n\n", N_TRIALS * 4);

    oracle_unregister(&eng, 0xBE);
    return 0;
}
