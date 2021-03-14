#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "soliton/q.h"

static size_t g_allocated_bytes = 0;

void k_mem_init(void) {
    g_allocated_bytes = 0;
    sym_init();
}

void k_mem_cleanup(void) {
    sym_cleanup();
}

size_t k_allocated_bytes(void) {
    return g_allocated_bytes;
}

static size_t type_width(int type) {
    switch (type) {
        case 1:  return 1;              /* bool */
        case 4:  return 1;              /* byte */
        case 5:  return sizeof(int16_t);/* short */
        case 6:  return sizeof(int32_t);/* int */
        case 7:  return sizeof(int64_t);/* long */
        case 8:  return sizeof(float);  /* real */
        case 9:  return sizeof(double); /* float */
        case 10: return 1;              /* char */
        case 11: return sizeof(S);      /* symbol */
        case 12: return sizeof(int64_t);/* timestamp */
        case 14: return sizeof(int32_t);/* date */
        case 15: return sizeof(int32_t);/* time */
        case 16: return sizeof(int64_t);/* timespan */
        case 0:  return sizeof(K);      /* mixed list */
        default: return sizeof(int64_t);
    }
}

K kb(bool b) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KB_ATOM;
    x->r = 1;
    x->g = b ? 1 : 0;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K kg(uint8_t b) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KG_ATOM;
    x->r = 1;
    x->g = b;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K kh(int16_t h) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KH_ATOM;
    x->r = 1;
    x->h = h;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K ki(int32_t i) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KI_ATOM;
    x->r = 1;
    x->i = i;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K kj(int64_t j) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KJ_ATOM;
    x->r = 1;
    x->j = j;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K ke(float e) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KE_ATOM;
    x->r = 1;
    x->e = e;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K kf(double f) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KF_ATOM;
    x->r = 1;
    x->f = f;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K kc(char c) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KC_ATOM;
    x->r = 1;
    x->c = c;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K ks(const char* s) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = KS_ATOM;
    x->r = 1;
    x->s = ss(s);
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K kerror(const char* msg) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = XE_ERR;
    x->r = 1;
    x->s = ss(msg ? msg : "err");
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K ktn(int type, int64_t length) {
    if (length < 0) length = 0;
    size_t elem_sz = type_width(type);
    size_t payload_sz = (length > 0) ? (size_t)length * elem_sz : 1;
    size_t total_sz = sizeof(struct k0) + payload_sz;
    
    K x = (K)calloc(1, total_sz);
    x->t = (int8_t)type;
    x->r = 1;
    x->n = length;
    g_allocated_bytes += total_sz;
    return x;
}

K kp(const char* s) {
    if (!s) return ktn(KC_VEC, 0);
    size_t len = strlen(s);
    K x = ktn(KC_VEC, (int64_t)len);
    memcpy(kC(x), s, len);
    return x;
}

K kpn(const char* s, int64_t n) {
    if (!s || n <= 0) return ktn(KC_VEC, 0);
    K x = ktn(KC_VEC, n);
    memcpy(kC(x), s, (size_t)n);
    return x;
}

K k_mixed_list(int64_t n) {
    return ktn(KM_LIST, n);
}

K xD(K keys, K vals) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = XD_DICT;
    x->r = 1;
    x->dict.k = keys;
    x->dict.v = vals;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

K xT(K dict) {
    K x = (K)calloc(1, sizeof(struct k0));
    x->t = XT_TABLE;
    x->r = 1;
    x->dict.k = dict; /* stores underlying dictionary */
    x->dict.v = NULL;
    g_allocated_bytes += sizeof(struct k0);
    return x;
}

void r0(K x) {
    if (!x) return;
    x->r--;
    if (x->r > 0) return;
    
    if (x->t < 0) {
        /* Atom */
        free(x);
        return;
    }
    
    if (x->t == KM_LIST) {
        /* Mixed list */
        for (int64_t i = 0; i < x->n; i++) {
            r0(kK(x)[i]);
        }
        free(x);
        return;
    }
    
    if (x->t == XD_DICT) {
        r0(x->dict.k);
        r0(x->dict.v);
        free(x);
        return;
    }
    
    if (x->t == XT_TABLE) {
        r0(x->dict.k);
        free(x);
        return;
    }
    
    if (x->t == XF_LAMBDA) {
        if (x->lambda.body) free(x->lambda.body);
        r0(x->lambda.params);
        r0(x->lambda.env);
        free(x);
        return;
    }
    
    if (x->t == XP_PROJ) {
        r0(x->proj.fn);
        r0(x->proj.args);
        free(x);
        return;
    }
    
    /* Homogeneous vector */
    free(x);
}
