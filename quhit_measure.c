/*
 * quhit_measure.c — Born-Rule Measurement, Collapse, Inspection
 *
 * Measurement samples from |ψ|² using Born rule.
 * Collapse zeros incompatible amplitudes and renormalizes.
 * Inspect reads state non-destructively (probabilities, entropy, purity).
 */

#include "quhit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * MEASUREMENT — Born rule sampling + collapse
 *
 * Local quhit: sample from 6 probabilities, collapse to |k⟩.
 * Entangled pair: compute marginal, sample, partial-collapse joint state,
 *                 then auto-determine partner's state.
 * ═══════════════════════════════════════════════════════════════════════════════ */

uint32_t quhit_measure(QuhitEngine *eng, uint32_t id)
{
    if (id >= eng->num_quhits) return 0;
    Quhit *q = &eng->quhits[id];

    /* Already collapsed? Return cached value */
    if (q->collapsed) return q->collapse_value;

    /* ── Entangled measurement ── */
    if (q->pair_id >= 0) {
        QuhitPair *p = &eng->pairs[q->pair_id];
        uint8_t side = q->pair_side;

        /* Compute marginal probabilities for our side */
        double probs[QUHIT_D] = {0};
        for (int a = 0; a < QUHIT_D; a++) {
            for (int b = 0; b < QUHIT_D; b++) {
                int idx = a * QUHIT_D + b;
                double prob = p->joint.re[idx] * p->joint.re[idx]
                            + p->joint.im[idx] * p->joint.im[idx];
                int my_val = (side == 0) ? a : b;
                probs[my_val] += prob;
            }
        }

        /* Born-rule sampling */
        double r = quhit_prng_double(eng);
        double cumul = 0;
        uint32_t result = QUHIT_D - 1;
        for (int k = 0; k < QUHIT_D; k++) {
            cumul += probs[k];
            if (cumul >= r) { result = k; break; }
        }

        /* Partial collapse: zero all entries where our index ≠ result */
        double surviving = 0;
        for (int a = 0; a < QUHIT_D; a++) {
            for (int b = 0; b < QUHIT_D; b++) {
                int idx = a * QUHIT_D + b;
                int my_val = (side == 0) ? a : b;
                if ((uint32_t)my_val != result) {
                    p->joint.re[idx] = 0;
                    p->joint.im[idx] = 0;
                } else {
                    surviving += p->joint.re[idx] * p->joint.re[idx]
                               + p->joint.im[idx] * p->joint.im[idx];
                }
            }
        }

        /* Renormalize using Quake fast inverse sqrt */
        if (surviving > 1e-30) {
            double scale = born_fast_isqrt(surviving);
            for (int i = 0; i < QUHIT_D2; i++) {
                p->joint.re[i] *= scale;
                p->joint.im[i] *= scale;
            }
        }

        /* Record outcome */
        q->collapsed = 1;
        q->collapse_value = result;
        eng->measured_values[id] = result;

        /* Update local state to |result⟩ */
        memset(&q->state, 0, sizeof(q->state));
        q->state.re[result] = 1.0;

        /* Partner's state is now determined — extract from collapsed row/col */
        uint32_t partner_id = (side == 0) ? p->id_b : p->id_a;
        Quhit *qp = &eng->quhits[partner_id];
        memset(&qp->state, 0, sizeof(qp->state));
        double pnorm = 0;
        for (int k = 0; k < QUHIT_D; k++) {
            int idx = (side == 0) ? result * QUHIT_D + k : k * QUHIT_D + result;
            qp->state.re[k] = p->joint.re[idx];
            qp->state.im[k] = p->joint.im[idx];
            pnorm += qp->state.re[k] * qp->state.re[k]
                   + qp->state.im[k] * qp->state.im[k];
        }
        if (pnorm > 1e-30) {
            double sc = born_fast_isqrt(pnorm);
            for (int k = 0; k < QUHIT_D; k++) {
                qp->state.re[k] *= sc;
                qp->state.im[k] *= sc;
            }
        }

        return result;
    }

    /* ── Local measurement ── */
    double r = quhit_prng_double(eng);
    int outcome = born_sample(q->state.re, q->state.im, QUHIT_D, r);

    /* Collapse to |outcome⟩ */
    born_collapse(q->state.re, q->state.im, QUHIT_D, outcome);

    q->collapsed = 1;
    q->collapse_value = (uint32_t)outcome;
    eng->measured_values[id] = (uint32_t)outcome;

    return (uint32_t)outcome;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PROBABILITY — P(outcome) for a specific outcome
 * ═══════════════════════════════════════════════════════════════════════════════ */

double quhit_prob(QuhitEngine *eng, uint32_t id, uint32_t outcome)
{
    if (id >= eng->num_quhits || outcome >= QUHIT_D) return 0;
    Quhit *q = &eng->quhits[id];

    if (q->pair_id >= 0) {
        /* Marginal probability from joint state */
        QuhitPair *p = &eng->pairs[q->pair_id];
        uint8_t side = q->pair_side;
        double prob = 0;
        for (int k = 0; k < QUHIT_D; k++) {
            int idx = (side == 0) ? outcome * QUHIT_D + k : k * QUHIT_D + outcome;
            prob += p->joint.re[idx] * p->joint.re[idx]
                  + p->joint.im[idx] * p->joint.im[idx];
        }
        return prob;
    }

    return q->state.re[outcome] * q->state.re[outcome]
         + q->state.im[outcome] * q->state.im[outcome];
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * INSPECT — Non-destructive state readout
 *
 * Fills a QuhitSnapshot with probabilities, entropy, purity, entanglement.
 * Does NOT collapse the state. This is the "God's eye view" of the substrate.
 * ═══════════════════════════════════════════════════════════════════════════════ */

void quhit_inspect(QuhitEngine *eng, uint32_t id, QuhitSnapshot *snap)
{
    memset(snap, 0, sizeof(*snap));
    if (id >= eng->num_quhits) return;
    Quhit *q = &eng->quhits[id];
    snap->quhit_id = id;
    snap->dim = QUHIT_D;

    /* ── Entangled: inspect joint state ── */
    if (q->pair_id >= 0) {
        QuhitPair *p = &eng->pairs[q->pair_id];
        uint8_t side = q->pair_side;
        snap->is_entangled = 1;

        /* Marginal probabilities for our side */
        uint32_t ne = 0;
        for (int a = 0; a < QUHIT_D; a++) {
            for (int b = 0; b < QUHIT_D; b++) {
                int idx = a * QUHIT_D + b;
                double prob = p->joint.re[idx] * p->joint.re[idx]
                            + p->joint.im[idx] * p->joint.im[idx];
                int my_val = (side == 0) ? a : b;
                snap->probs[my_val] += prob;
                snap->total_prob += prob;

                if (prob > 1e-30 && ne < QUHIT_D2) {
                    snap->entries[ne].idx_a = a;
                    snap->entries[ne].idx_b = b;
                    snap->entries[ne].amp_re = p->joint.re[idx];
                    snap->entries[ne].amp_im = p->joint.im[idx];
                    snap->entries[ne].probability = prob;
                    snap->entries[ne].phase_rad = atan2(p->joint.im[idx],
                                                        p->joint.re[idx]);
                    ne++;
                }
            }
        }
        snap->num_entries = ne;

        /* Reduced density matrix ρ = Tr_B(|ψ⟩⟨ψ|) */
        double rho_re[QUHIT_D * QUHIT_D] = {0};
        double rho_im[QUHIT_D * QUHIT_D] = {0};

        for (int i = 0; i < QUHIT_D; i++) {
            for (int j = 0; j < QUHIT_D; j++) {
                for (int k = 0; k < QUHIT_D; k++) {
                    int idx_i, idx_j;
                    if (side == 0) {
                        idx_i = i * QUHIT_D + k;
                        idx_j = j * QUHIT_D + k;
                    } else {
                        idx_i = k * QUHIT_D + i;
                        idx_j = k * QUHIT_D + j;
                    }
                    /* ρ[i,j] += ψ(i,k) × conj(ψ(j,k)) */
                    double ai_re = p->joint.re[idx_i], ai_im = p->joint.im[idx_i];
                    double aj_re = p->joint.re[idx_j], aj_im = p->joint.im[idx_j];
                    rho_re[i * QUHIT_D + j] += ai_re * aj_re + ai_im * aj_im;
                    rho_im[i * QUHIT_D + j] += ai_re * aj_im - ai_im * aj_re;
                }
            }
        }

        /* Purity = Tr(ρ²) */
        double purity = 0;
        for (int i = 0; i < QUHIT_D; i++)
            for (int j = 0; j < QUHIT_D; j++)
                purity += rho_re[i*QUHIT_D+j] * rho_re[i*QUHIT_D+j]
                        + rho_im[i*QUHIT_D+j] * rho_im[i*QUHIT_D+j];
        snap->purity = purity;

        /* Entropy = -Σ λ log₂(λ) from diagonal of ρ */
        double entropy = 0;
        for (int k = 0; k < QUHIT_D; k++) {
            double lam = rho_re[k * QUHIT_D + k];
            if (lam > 1e-14) entropy -= lam * log2(lam);
        }
        snap->entropy = entropy;

        /* Schmidt rank = number of nonzero eigenvalues */
        int rank = 0;
        for (int k = 0; k < QUHIT_D; k++)
            if (rho_re[k * QUHIT_D + k] > 1e-14) rank++;
        snap->schmidt_rank = rank;

        return;
    }

    /* ── Local: inspect local state ── */
    snap->is_entangled = 0;
    snap->schmidt_rank = 1;
    snap->total_prob = 0;
    double purity = 0;
    double entropy = 0;

    for (int k = 0; k < QUHIT_D; k++) {
        double p = q->state.re[k] * q->state.re[k]
                 + q->state.im[k] * q->state.im[k];
        snap->probs[k] = p;
        snap->total_prob += p;
        purity += p * p;
        if (p > 1e-14) entropy -= p * log2(p);
    }

    snap->purity  = purity;
    snap->entropy = entropy;
}
