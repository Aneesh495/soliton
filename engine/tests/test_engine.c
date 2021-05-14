#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "soliton/q.h"
#include "soliton/eval.h"
#include "soliton/ipc.h"

static void test_vector_arithmetic(void) {
    printf("[TEST] Vector arithmetic...\n");
    K a = ktn(KF_VEC, 4);
    K b = ktn(KF_VEC, 4);
    for (int i = 0; i < 4; i++) {
        kF(a)[i] = (double)(i + 1);
        kF(b)[i] = 10.0;
    }
    K c = v_add(a, b);
    assert(c->t == KF_VEC);
    assert(c->n == 4);
    assert(kF(c)[0] == 11.0);
    assert(kF(c)[3] == 14.0);
    r0(a);
    r0(b);
    r0(c);
    printf("  PASS: Vector arithmetic verified.\n");
}

static void test_reductions_and_stats(void) {
    printf("[TEST] Statistical reductions...\n");
    K x = ktn(KF_VEC, 5);
    for (int i = 0; i < 5; i++) {
        kF(x)[i] = (double)(i + 1); /* 1, 2, 3, 4, 5 */
    }
    K s = v_sum(x);
    assert(s->t == KF_ATOM);
    assert(s->f == 15.0);
    r0(s);
    
    K m = v_avg(x);
    assert(m->t == KF_ATOM);
    assert(m->f == 3.0);
    r0(m);
    
    K med = v_med(x);
    assert(med->t == KF_ATOM);
    assert(med->f == 3.0);
    r0(med);
    
    r0(x);
    printf("  PASS: Statistical reductions verified.\n");
}

static void test_ipc_roundtrip(void) {
    printf("[TEST] IPC wire protocol roundtrip...\n");
    K orig = ktn(KF_VEC, 3);
    kF(orig)[0] = 3.14159;
    kF(orig)[1] = 2.71828;
    kF(orig)[2] = 1.41421;
    
    size_t len = 0;
    uint8_t* wire = ipc_serialize(orig, IPC_SYNC, &len);
    assert(wire != NULL);
    assert(len > 8);
    
    uint8_t msg_type = 0;
    K decoded = ipc_deserialize(wire, len, &msg_type);
    assert(msg_type == IPC_SYNC);
    assert(decoded->t == KF_VEC);
    assert(decoded->n == 3);
    assert(kF(decoded)[0] == kF(orig)[0]);
    assert(kF(decoded)[1] == kF(orig)[1]);
    
    free(wire);
    r0(orig);
    r0(decoded);
    printf("  PASS: IPC roundtrip verified.\n");
}

static void test_eval_parser(void) {
    printf("[TEST] Expression parser and environment...\n");
    Env* env = env_new(NULL);
    
    K res1 = q_eval_string(env, "a: 10; b: 25; a + b");
    assert(res1 != NULL);
    assert(res1->j == 35);
    r0(res1);
    
    K res2 = q_eval_string(env, "2 * 3 + 4");
    /* Right to left: 3 + 4 = 7; 2 * 7 = 14 */
    assert(res2->j == 14);
    r0(res2);
    
    env_free(env);
    printf("  PASS: Parser and evaluation verified.\n");
}

int main(void) {
    k_mem_init();
    printf("=========================================\n");
    printf("Running Soliton Engine Verification Suite\n");
    printf("=========================================\n");
    
    test_vector_arithmetic();
    test_reductions_and_stats();
    test_ipc_roundtrip();
    test_eval_parser();
    
    printf("=========================================\n");
    printf("ALL ENGINE TESTS PASSED SUCCESSFULLY!\n");
    printf("Allocated bytes active: %zu\n", k_allocated_bytes());
    printf("=========================================\n");
    k_mem_cleanup();
    return 0;
}
