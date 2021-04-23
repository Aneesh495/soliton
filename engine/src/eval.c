#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "soliton/q.h"
#include "soliton/eval.h"

Env* env_new(Env* parent) {
    Env* env = (Env*)malloc(sizeof(Env));
    env->parent = parent;
    env->dict = xD(ktn(KS_VEC, 0), k_mixed_list(0));
    return env;
}

void env_free(Env* env) {
    if (!env) return;
    r0(env->dict);
    free(env);
}

K env_get(Env* env, S sym) {
    for (Env* e = env; e != NULL; e = e->parent) {
        K keys = e->dict->dict.k;
        K vals = e->dict->dict.v;
        for (int64_t i = 0; i < keys->n; i++) {
            if (kS(keys)[i] == sym) {
                return kK(vals)[i];
            }
        }
    }
    return NULL;
}

void env_set(Env* env, S sym, K val) {
    if (!env) return;
    K keys = env->dict->dict.k;
    K vals = env->dict->dict.v;
    for (int64_t i = 0; i < keys->n; i++) {
        if (kS(keys)[i] == sym) {
            r0(kK(vals)[i]);
            kK(vals)[i] = val;
            return;
        }
    }
    /* Append new entry */
    K new_keys = ktn(KS_VEC, keys->n + 1);
    memcpy(kS(new_keys), kS(keys), (size_t)keys->n * sizeof(S));
    kS(new_keys)[keys->n] = sym;
    
    K new_vals = k_mixed_list(vals->n + 1);
    for (int64_t i = 0; i < vals->n; i++) {
        kK(new_vals)[i] = r1(kK(vals)[i]);
    }
    kK(new_vals)[vals->n] = val;
    
    r0(env->dict->dict.k);
    r0(env->dict->dict.v);
    env->dict->dict.k = new_keys;
    env->dict->dict.v = new_vals;
}

K q_load_file(Env* env, const char* filepath) {
    FILE* f = fopen(filepath, "rb");
    if (!f) {
        fprintf(stderr, "Could not open file: %s\n", filepath);
        return kerror("file");
    }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    
    char* buf = (char*)malloc((size_t)sz + 1);
    size_t read_bytes = fread(buf, 1, (size_t)sz, f);
    buf[read_bytes] = '\0';
    fclose(f);
    
    K res = q_eval_string(env, buf);
    free(buf);
    return res;
}
