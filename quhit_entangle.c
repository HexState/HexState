/*
 * quhit_entangle.c — Entanglement: Braid, Product, Unbraid
 *
 * Two modes of entanglement creation:
 *   Bell:    (1/√6) Σ|k,k⟩ — maximally entangled (576 bytes)
 *   Product: |ψ_a⟩ ⊗ |ψ_b⟩ — separable (576 bytes, becomes entangled via CZ)
 *
 * Disentangle extracts marginals back into local quhit states.
 */

#include "quhit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * HELPER — Allocate a pair slot
 * ═══════════════════════════════════════════════════════════════════════════════ */

static int alloc_pair(QuhitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (eng->num_pairs >= MAX_PAIRS) {
        fprintf(stderr, "[QUHIT] ERROR: max pairs (%d) reached\n", MAX_PAIRS);
        return -1;
    }

    int slot = (int)eng->num_pairs++;
    QuhitPair *p = &eng->pairs[slot];
    memset(p, 0, sizeof(*p));
    p->id_a   = id_a;
    p->id_b   = id_b;
    p->active = 1;

    Quhit *qa = &eng->quhits[id_a];
    Quhit *qb = &eng->quhits[id_b];
    qa->pair_id   = slot;
    qa->pair_side = 0;
    qb->pair_id   = slot;
    qb->pair_side = 1;

    return slot;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * BELL PAIR — (1/√6) Σ|k,k⟩
 *
 * Maximally entangled state. 576 bytes.
 * Uses quhit_management.h primitive.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int quhit_entangle_bell(QuhitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_quhits || id_b >= eng->num_quhits) return -1;

    int slot = alloc_pair(eng, id_a, id_b);
    if (slot < 0) return -1;

    /* Uses quhit_management.h: fills diagonal with 1/√6 */
    qm_entangle_bell(&eng->pairs[slot].joint);

    return slot;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRODUCT PAIR — |ψ_a⟩ ⊗ |ψ_b⟩
 *
 * Creates joint state from tensor product of current local states.
 * Initially separable. CZ gate then creates genuine entanglement.
 * Uses quhit_management.h primitive.
 * ═══════════════════════════════════════════════════════════════════════════════ */

int quhit_entangle_product(QuhitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_quhits || id_b >= eng->num_quhits) return -1;

    int slot = alloc_pair(eng, id_a, id_b);
    if (slot < 0) return -1;

    /* Uses quhit_management.h: tensor product */
    qm_entangle_product(&eng->pairs[slot].joint,
                        &eng->quhits[id_a].state,
                        &eng->quhits[id_b].state);

    return slot;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * DISENTANGLE — Break pair, extract marginals to local states
 *
 * After we're done with a pairwise interaction, we extract each quhit's
 * marginal state from the joint and return to independent local states.
 *
 * For side A: ψ_A(k) = √(Σ_b |ψ(k,b)|²) × phase_of_dominant_component
 * This preserves probabilities but not full phase coherence
 * (which is physically correct — tracing out the partner destroys coherence).
 * ═══════════════════════════════════════════════════════════════════════════════ */

void quhit_disentangle(QuhitEngine *eng, uint32_t id_a, uint32_t id_b)
{
    if (id_a >= eng->num_quhits || id_b >= eng->num_quhits) return;

    Quhit *qa = &eng->quhits[id_a];
    Quhit *qb = &eng->quhits[id_b];

    /* Verify they share a pair */
    if (qa->pair_id < 0 || qa->pair_id != qb->pair_id) return;

    QuhitPair *p = &eng->pairs[qa->pair_id];

    /* Extract marginal for A (rows) */
    for (int k = 0; k < QUHIT_D; k++) {
        double prob = 0;
        double dom_re = 0, dom_im = 0;
        double dom_p = 0;
        for (int b = 0; b < QUHIT_D; b++) {
            int idx = k * QUHIT_D + b;
            double r = p->joint.re[idx], i = p->joint.im[idx];
            double p2 = r * r + i * i;
            prob += p2;
            if (p2 > dom_p) { dom_p = p2; dom_re = r; dom_im = i; }
        }
        if (prob > 1e-30) {
            double amp = sqrt(prob);
            double phase = atan2(dom_im, dom_re);
            qa->state.re[k] = amp * cos(phase);
            qa->state.im[k] = amp * sin(phase);
        } else {
            qa->state.re[k] = 0;
            qa->state.im[k] = 0;
        }
    }

    /* Extract marginal for B (columns) */
    for (int k = 0; k < QUHIT_D; k++) {
        double prob = 0;
        double dom_re = 0, dom_im = 0;
        double dom_p = 0;
        for (int a = 0; a < QUHIT_D; a++) {
            int idx = a * QUHIT_D + k;
            double r = p->joint.re[idx], i = p->joint.im[idx];
            double p2 = r * r + i * i;
            prob += p2;
            if (p2 > dom_p) { dom_p = p2; dom_re = r; dom_im = i; }
        }
        if (prob > 1e-30) {
            double amp = sqrt(prob);
            double phase = atan2(dom_im, dom_re);
            qb->state.re[k] = amp * cos(phase);
            qb->state.im[k] = amp * sin(phase);
        } else {
            qb->state.re[k] = 0;
            qb->state.im[k] = 0;
        }
    }

    /* Renormalize both */
    sup_renormalize(qa->state.re, qa->state.im, QUHIT_D);
    sup_renormalize(qb->state.re, qb->state.im, QUHIT_D);

    /* Deactivate pair */
    p->active = 0;
    qa->pair_id = -1;
    qb->pair_id = -1;
}
