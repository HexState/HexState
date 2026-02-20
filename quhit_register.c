/*
 * quhit_register.c — 100T-Scale Quhit Registers
 *
 * A register is a named collection of N quhits that can be addressed
 * as a single Hilbert space without materializing D^N amplitudes.
 *
 * Strategy: sparse-entry representation. Only track nonzero basis states.
 * For 100T quhits, the full state is 6^(100T) amplitudes — impossible.
 * Instead we store up to 4096 sparse entries and operate on them.
 *
 * Operations on registers:
 *   - GHZ entanglement: (1/√D) Σ|k,k,...,k⟩
 *   - DFT on individual quhit positions
 *   - CZ between pairs within register
 *   - Measurement of individual or all quhits
 */

#include "quhit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * REGISTER INIT
 *
 * Create a register of n_quhits, each dimension dim.
 * Initial state: |0,0,...,0⟩ (single basis entry, amplitude = 1.0)
 * ═══════════════════════════════════════════════════════════════════════════════ */

int quhit_reg_init(QuhitEngine *eng, uint64_t chunk_id,
                   uint64_t n_quhits, uint32_t dim)
{
    if (eng->num_registers >= MAX_REGISTERS) {
        fprintf(stderr, "[QUHIT] ERROR: max registers (%d) reached\n",
                MAX_REGISTERS);
        return -1;
    }

    int idx = (int)eng->num_registers++;
    QuhitRegister *reg = &eng->registers[idx];
    memset(reg, 0, sizeof(*reg));

    reg->chunk_id  = chunk_id;
    reg->n_quhits  = n_quhits;
    reg->dim       = dim;
    reg->bulk_rule = 0;

    /* Initial state: |0,0,...,0⟩ = basis_state 0 */
    reg->num_nonzero = 1;
    reg->entries[0].basis_state = 0;
    reg->entries[0].amp_re = 1.0;
    reg->entries[0].amp_im = 0.0;

    reg->collapsed = 0;
    reg->magic_base = MAGIC_PTR(chunk_id, 0);

    return idx;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * HELPER — Extract/set digit k of a packed basis state  (base-D number)
 *
 * basis_state = q₀ + q₁·D + q₂·D² + ...
 * digit(basis, k) = (basis / D^k) % D
 * ═══════════════════════════════════════════════════════════════════════════════ */

static uint32_t basis_digit(uint64_t basis, uint64_t k, uint32_t D)
{
    uint64_t pow = 1;
    for (uint64_t i = 0; i < k; i++) pow *= D;
    return (uint32_t)((basis / pow) % D);
}

static uint64_t basis_set_digit(uint64_t basis, uint64_t k, uint32_t D,
                                uint32_t val)
{
    uint64_t pow = 1;
    for (uint64_t i = 0; i < k; i++) pow *= D;
    uint32_t old = (uint32_t)((basis / pow) % D);
    return basis - (uint64_t)old * pow + (uint64_t)val * pow;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * GHZ ENTANGLEMENT — (1/√D) Σ|k,k,...,k⟩
 *
 * Creates a D-entry GHZ state across all quhits in the register.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void quhit_reg_entangle_all(QuhitEngine *eng, int reg_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QuhitRegister *reg = &eng->registers[reg_idx];

    uint32_t D = reg->dim;
    double amp = born_fast_isqrt((double)D);

    reg->num_nonzero = D;
    for (uint32_t k = 0; k < D && k < 4096; k++) {
        /* |k,k,...,k⟩ = k + k·D + k·D² + ... = k × (D^n - 1)/(D - 1) */
        uint64_t basis = 0;
        uint64_t pow = 1;
        for (uint64_t q = 0; q < reg->n_quhits; q++) {
            basis += k * pow;
            pow *= D;
        }
        reg->entries[k].basis_state = basis;
        reg->entries[k].amp_re = amp;
        reg->entries[k].amp_im = 0;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DFT ON SINGLE QUHIT IN REGISTER
 *
 * Applies DFT₆ to quhit at position quhit_idx.
 * Each basis entry with digit d at position quhit_idx splits into D entries
 * with digits 0..D-1, weighted by DFT₆[new_d][d].
 * ═══════════════════════════════════════════════════════════════════════════════ */

void quhit_reg_apply_dft(QuhitEngine *eng, int reg_idx, uint64_t quhit_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QuhitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    /* Temporary storage for new entries */
    struct { uint64_t basis; double re, im; } temp[4096];
    uint32_t n_temp = 0;

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        uint64_t basis = reg->entries[e].basis_state;
        double a_re = reg->entries[e].amp_re;
        double a_im = reg->entries[e].amp_im;
        uint32_t d = basis_digit(basis, quhit_idx, D);

        for (uint32_t j = 0; j < D && n_temp < 4096; j++) {
            uint64_t new_basis = basis_set_digit(basis, quhit_idx, D, j);
            double tw_re, tw_im;
            if (D == 6) {
                tw_re = DFT6[j][d].re;
                tw_im = DFT6[j][d].im;
            } else {
                double angle = 2.0 * M_PI * j * d / D;
                double inv = born_fast_isqrt((double)D);
                tw_re = inv * cos(angle);
                tw_im = inv * sin(angle);
            }
            double new_re = tw_re * a_re - tw_im * a_im;
            double new_im = tw_re * a_im + tw_im * a_re;

            /* Merge with existing temp entry if same basis */
            int merged = 0;
            for (uint32_t t = 0; t < n_temp; t++) {
                if (temp[t].basis == new_basis) {
                    temp[t].re += new_re;
                    temp[t].im += new_im;
                    merged = 1;
                    break;
                }
            }
            if (!merged) {
                temp[n_temp].basis = new_basis;
                temp[n_temp].re = new_re;
                temp[n_temp].im = new_im;
                n_temp++;
            }
        }
    }

    /* Write back, pruning near-zero entries */
    reg->num_nonzero = 0;
    for (uint32_t t = 0; t < n_temp && reg->num_nonzero < 4096; t++) {
        double prob = temp[t].re * temp[t].re + temp[t].im * temp[t].im;
        if (prob > 1e-30) {
            reg->entries[reg->num_nonzero].basis_state = temp[t].basis;
            reg->entries[reg->num_nonzero].amp_re = temp[t].re;
            reg->entries[reg->num_nonzero].amp_im = temp[t].im;
            reg->num_nonzero++;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * CZ BETWEEN TWO QUHITS IN REGISTER
 *
 * Applies ω^(a·b) phase to each basis entry based on digits at positions
 * idx_a and idx_b.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void quhit_reg_apply_cz(QuhitEngine *eng, int reg_idx,
                        uint64_t idx_a, uint64_t idx_b)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return;
    QuhitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;
    double omega = 2.0 * M_PI / D;

    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        uint32_t a = basis_digit(reg->entries[e].basis_state, idx_a, D);
        uint32_t b = basis_digit(reg->entries[e].basis_state, idx_b, D);
        double phase = omega * a * b;
        double cos_p = cos(phase), sin_p = sin(phase);
        double re = reg->entries[e].amp_re;
        double im = reg->entries[e].amp_im;
        reg->entries[e].amp_re = re * cos_p - im * sin_p;
        reg->entries[e].amp_im = re * sin_p + im * cos_p;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASUREMENT — Born-rule sampling of single quhit in register
 *
 * Computes marginal P(k) for the target quhit, samples, collapses.
 * Returns the measured value.
 * ═══════════════════════════════════════════════════════════════════════════════ */

uint64_t quhit_reg_measure(QuhitEngine *eng, int reg_idx,
                           uint64_t quhit_idx)
{
    if (reg_idx < 0 || (uint32_t)reg_idx >= eng->num_registers) return 0;
    QuhitRegister *reg = &eng->registers[reg_idx];
    uint32_t D = reg->dim;

    /* Compute marginal probabilities */
    double probs[64] = {0};  /* supports up to D=64 */
    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        uint32_t d = basis_digit(reg->entries[e].basis_state, quhit_idx, D);
        double re = reg->entries[e].amp_re, im = reg->entries[e].amp_im;
        probs[d] += re * re + im * im;
    }

    /* Born-rule sampling */
    double r = quhit_prng_double(eng);
    double cumul = 0;
    uint32_t result = D - 1;
    for (uint32_t k = 0; k < D; k++) {
        cumul += probs[k];
        if (cumul >= r) { result = k; break; }
    }

    /* Collapse: remove entries where digit ≠ result, renormalize */
    double surviving = 0;
    uint32_t new_count = 0;
    for (uint32_t e = 0; e < reg->num_nonzero; e++) {
        uint32_t d = basis_digit(reg->entries[e].basis_state, quhit_idx, D);
        if (d == result) {
            reg->entries[new_count] = reg->entries[e];
            double re = reg->entries[new_count].amp_re;
            double im = reg->entries[new_count].amp_im;
            surviving += re * re + im * im;
            new_count++;
        }
    }
    reg->num_nonzero = new_count;

    if (surviving > 1e-30) {
        double scale = born_fast_isqrt(surviving);
        for (uint32_t e = 0; e < new_count; e++) {
            reg->entries[e].amp_re *= scale;
            reg->entries[e].amp_im *= scale;
        }
    }

    return result;
}
