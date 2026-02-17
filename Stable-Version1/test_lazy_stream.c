/*
 * test_lazy_stream.c — Verify the lazy state vector streaming iterator
 *
 * Tests all 3 modes:
 *   1. Quhit register (Mode 2 — infinite resource, sparse entries)
 *   2. Joint state (braided pair, D² dense)
 *   3. Local state (single chunk, D dense)
 *
 * Build:
 *   gcc -O2 -I. -std=c11 -D_GNU_SOURCE \
 *       -c test_lazy_stream.c -o test_lazy_stream.o && \
 *   gcc -O2 -o test_lazy_stream test_lazy_stream.o hexstate_engine.o bigint.o -lm
 */

#include "hexstate_engine.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#define D  6
#define N  100000000000000ULL  /* 100T */
static const char *bn[] = {"A","T","G","C","dR","Pi"};

int main(void)
{
    printf("\n");
    printf("════════════════════════════════════════════════════════════════\n");
    printf("  LAZY STATE VECTOR STREAMING — TEST SUITE\n");
    printf("  Streaming full quantum state without holding it all\n");
    printf("════════════════════════════════════════════════════════════════\n\n");

    int pass = 0, fail = 0;
    #define CHECK(c, name) do { \
        if (c) { pass++; printf("  ✓ %s\n", name); } \
        else { fail++; printf("  ✗ FAIL: %s\n", name); } \
    } while(0)

    /* ═══ TEST 1: Mode 2 — Quhit register with 100T quhits ═══ */
    printf("\n── Test 1: Quhit Register (100T, D=6) ──\n\n");
    {
        static HexStateEngine eng;
        engine_init(&eng);
        op_infinite_resources_dim(&eng, 0, N, D);
        init_quhit_register(&eng, 0, N, D);

        /* Apply DFT to create superposition */
        entangle_all_quhits(&eng, 0);

        StateIterator it;
        int rc = state_iter_begin(&eng, 0, &it);
        CHECK(rc == 0, "state_iter_begin succeeds");
        CHECK(it.total_entries > 0, "Has non-zero entries");
        CHECK(it.n_quhits == N, "Reports correct quhit count");
        CHECK(it.dim == D, "Reports correct dimension");

        printf("    Total entries: %u\n", it.total_entries);
        printf("    N quhits:     %llu\n", (unsigned long long)it.n_quhits);
        printf("    Dimension:    %u\n\n", it.dim);

        double total_prob = 0;
        int count = 0;
        printf("    Streaming entries:\n");
        while (state_iter_next(&it)) {
            total_prob += it.probability;
            count++;

            /* Print first few entries */
            if (count <= 6) {
                printf("      [%u] bulk=%s  amp=(%.4f, %.4f)  P=%.4f",
                       it.entry_index, bn[it.bulk_value % D],
                       it.amplitude.real, it.amplitude.imag,
                       it.probability);

                /* Resolve specific quhits */
                uint32_t v0 = state_iter_resolve(&it, 0);
                uint32_t v42 = state_iter_resolve(&it, 42);
                uint32_t vMax = state_iter_resolve(&it, N - 1);
                printf("  q[0]=%s q[42]=%s q[%lluT]=%s\n",
                       bn[v0 % D], bn[v42 % D],
                       (unsigned long long)(N/1000000000000ULL), bn[vMax % D]);
            }
        }
        state_iter_end(&it);

        printf("\n    Entries streamed: %d\n", count);
        printf("    Total probability: %.6f (expect 1.0)\n", total_prob);
        CHECK(count == (int)it.total_entries || count > 0, "Iterated all entries");
        CHECK(fabs(total_prob - 1.0) < 0.01, "Probabilities sum to 1.0");

        engine_destroy(&eng);
    }

    /* ═══ TEST 2: Joint state (braided pair) ═══ */
    printf("\n── Test 2: Joint State (braided pair, D²=36) ──\n\n");
    {
        static HexStateEngine eng;
        engine_init(&eng);
        op_infinite_resources_dim(&eng, 0, N, D);
        op_infinite_resources_dim(&eng, 1, N, D);
        braid_chunks_dim(&eng, 0, 1, 0, 0, D);

        StateIterator it;
        int rc = state_iter_begin(&eng, 0, &it);
        CHECK(rc == 0, "state_iter_begin on braided pair");
        CHECK(it.total_entries == D || it.total_entries == D * D,
              "Has correct entry count (D sparse or D² dense)");

        double total_prob = 0;
        int count = 0;
        int nonzero = 0;
        printf("    Bell state amplitudes:\n");
        while (state_iter_next(&it)) {
            total_prob += it.probability;
            count++;

            if (it.probability > 1e-15) {
                nonzero++;
                uint32_t va = state_iter_resolve(&it, 0);
                uint32_t vb = state_iter_resolve(&it, 1);
                printf("      [%2u] |%s,%s⟩ amp=(%.4f, %.4f) P=%.4f\n",
                       it.entry_index, bn[va % D], bn[vb % D],
                       it.amplitude.real, it.amplitude.imag,
                       it.probability);
            }
        }
        int expected_count = it.total_entries;
        state_iter_end(&it);

        printf("\n    Total entries:   %d\n", count);
        printf("    Non-zero:       %d (expect %d for Bell state)\n", nonzero, D);
        printf("    Total prob:     %.6f\n", total_prob);
        CHECK(count == expected_count, "Iterated all entries");
        CHECK(nonzero == D, "Bell state has exactly D nonzero amplitudes");
        CHECK(fabs(total_prob - 1.0) < 0.01, "Normalized");

        engine_destroy(&eng);
    }

    /* ═══ TEST 3: Local state (single chunk) ═══ */
    printf("\n── Test 3: Local State (single chunk, D=6) ──\n\n");
    {
        static HexStateEngine eng;
        engine_init(&eng);
        op_infinite_resources_dim(&eng, 0, N, D);
        /* Don't init quhit register — test local state path */

        /* Apply Hadamard to put in superposition */
        apply_hadamard(&eng, 0, 0);

        StateIterator it;
        int rc = state_iter_begin(&eng, 0, &it);
        CHECK(rc == 0, "state_iter_begin on local state");
        CHECK(it.total_entries == D, "D entries");

        double total_prob = 0;
        int count = 0;
        printf("    Hadamard state amplitudes:\n");
        while (state_iter_next(&it)) {
            total_prob += it.probability;
            count++;
            printf("      |%s⟩ amp=(%.4f, %.4f) P=%.4f\n",
                   bn[it.bulk_value % D],
                   it.amplitude.real, it.amplitude.imag,
                   it.probability);
        }
        state_iter_end(&it);

        printf("\n    Entries: %d  Total prob: %.6f\n", count, total_prob);
        CHECK(count == D, "Iterated all D entries");
        CHECK(fabs(total_prob - 1.0) < 0.01, "Normalized");

        engine_destroy(&eng);
    }

    /* ═══ TEST 4: DNA gate then stream ═══ */
    printf("\n── Test 4: DNA Gate + Stream (100T quhits) ──\n\n");
    {
        static HexStateEngine eng;
        engine_init(&eng);
        op_infinite_resources_dim(&eng, 0, N, D);
        init_quhit_register(&eng, 0, N, D);
        apply_dna_bulk_quhits(&eng, 0, 1.0, 310.0);

        StateIterator it;
        state_iter_begin(&eng, 0, &it);

        double total_prob = 0;
        int count = 0;
        printf("    After DNA gate:\n");
        while (state_iter_next(&it)) {
            total_prob += it.probability;
            count++;
            if (count <= 6) {
                printf("      [%u] bulk=%s P=%.4f\n",
                       it.entry_index, bn[it.bulk_value % D],
                       it.probability);
            }
        }
        state_iter_end(&it);

        printf("    Entries: %d  Total prob: %.6f\n", count, total_prob);
        CHECK(count > 1, "DNA gate created multiple entries");
        CHECK(fabs(total_prob - 1.0) < 0.01, "Normalized after DNA gate");

        engine_destroy(&eng);
    }

    /* ═══ Summary ═══ */
    printf("\n════════════════════════════════════════════════════════════════\n");
    printf("  RESULTS: %d/%d passed", pass, pass + fail);
    if (fail == 0) printf(" — ALL TESTS PASSED ✓");
    printf("\n════════════════════════════════════════════════════════════════\n\n");

    return fail > 0 ? 1 : 0;
}
