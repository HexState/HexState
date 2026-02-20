/*
 * quhit_core.c — Engine Lifecycle, PRNG, Quhit Init
 *
 * The foundation. Every quhit starts here.
 */

#include "quhit_engine.h"

/* ═══════════════════════════════════════════════════════════════════════════════
 * ENGINE LIFECYCLE
 * ═══════════════════════════════════════════════════════════════════════════════ */

void quhit_engine_init(QuhitEngine *eng)
{
    memset(eng, 0, sizeof(*eng));
    eng->prng_state = 0x5DEECE66DULL ^ 0xCAFEBABEULL;
}

void quhit_engine_destroy(QuhitEngine *eng)
{
    /* Deactivate all pairs (no heap allocs — QuhitJoint is inline) */
    for (uint32_t i = 0; i < eng->num_pairs; i++)
        eng->pairs[i].active = 0;

    eng->num_quhits    = 0;
    eng->num_pairs     = 0;
    eng->num_registers = 0;
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * PRNG — LCG (same constants as java.util.Random)
 *
 * Fast, deterministic, reproducible. Not cryptographic.
 * The quantum engine needs randomness for Born-rule sampling.
 * ═══════════════════════════════════════════════════════════════════════════════ */

uint64_t quhit_prng(QuhitEngine *eng)
{
    eng->prng_state = eng->prng_state * 6364136223846793005ULL
                    + 1442695040888963407ULL;
    return eng->prng_state;
}

double quhit_prng_double(QuhitEngine *eng)
{
    return (double)(quhit_prng(eng) >> 11) / (double)(1ULL << 53);
}

/* ═══════════════════════════════════════════════════════════════════════════════
 * QUHIT INITIALIZATION
 *
 * Each quhit = 6 complex amplitudes = 96 bytes.
 * No heap allocation — everything is inline in the engine struct.
 * ═══════════════════════════════════════════════════════════════════════════════ */

uint32_t quhit_init(QuhitEngine *eng)
{
    if (eng->num_quhits >= MAX_QUHITS) {
        fprintf(stderr, "[QUHIT] ERROR: max quhits (%d) reached\n", MAX_QUHITS);
        return UINT32_MAX;
    }

    uint32_t id = eng->num_quhits++;
    Quhit *q = &eng->quhits[id];

    q->id             = id;
    q->collapsed      = 0;
    q->collapse_value = 0;
    q->pair_id        = -1;
    q->pair_side      = 0;

    /* |0⟩ — uses quhit_management.h primitive */
    qm_init_zero(&q->state);

    return id;
}

uint32_t quhit_init_plus(QuhitEngine *eng)
{
    uint32_t id = quhit_init(eng);
    if (id == UINT32_MAX) return id;

    /* |+⟩ = (1/√6) Σ|k⟩ — uses quhit_management.h primitive */
    qm_init_plus(&eng->quhits[id].state);

    return id;
}

uint32_t quhit_init_basis(QuhitEngine *eng, uint32_t k)
{
    uint32_t id = quhit_init(eng);
    if (id == UINT32_MAX) return id;

    if (k < QUHIT_D) {
        QuhitState *s = &eng->quhits[id].state;
        memset(s, 0, sizeof(*s));
        s->re[k] = 1.0;
    }

    return id;
}

void quhit_reset(QuhitEngine *eng, uint32_t id)
{
    if (id >= eng->num_quhits) return;
    Quhit *q = &eng->quhits[id];

    /* If entangled, disentangle first */
    if (q->pair_id >= 0) {
        QuhitPair *p = &eng->pairs[q->pair_id];
        uint32_t partner = (q->pair_side == 0) ? p->id_b : p->id_a;
        quhit_disentangle(eng, id, partner);
    }

    q->collapsed      = 0;
    q->collapse_value = 0;
    q->pair_id        = -1;
    q->pair_side      = 0;
    qm_init_zero(&q->state);
}
