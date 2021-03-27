#include <stdlib.h>
#include <string.h>
#include "soliton/q.h"
#include "soliton/eval.h"

static double get_val_f(K x, int64_t idx) {
    if (x->t < 0) {
        if (x->t == KF_ATOM) return x->f;
        if (x->t == KJ_ATOM) return (double)x->j;
        if (x->t == KI_ATOM) return (double)x->i;
        if (x->t == KH_ATOM) return (double)x->h;
        if (x->t == KB_ATOM || x->t == KG_ATOM) return (double)x->g;
        return 0.0;
    }
    if (x->t == KF_VEC) return kF(x)[idx];
    if (x->t == KJ_VEC) return (double)kJ(x)[idx];
    if (x->t == KI_VEC) return (double)kI(x)[idx];
    if (x->t == KH_VEC) return (double)kH(x)[idx];
    if (x->t == KB_VEC || x->t == KG_VEC) return (double)kG(x)[idx];
    return 0.0;
}

K v_eq(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t == KS_ATOM && y->t == KS_ATOM) {
        return kb(x->s == y->s);
    }
    if (x->t < 0 && y->t < 0) {
        return kb(get_val_f(x, 0) == get_val_f(y, 0));
    }
    if (x->t > 0 && y->t > 0) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KB_VEC, x->n);
        if (x->t == KS_VEC && y->t == KS_VEC) {
            for (int64_t i = 0; i < x->n; i++) kG(res)[i] = (kS(x)[i] == kS(y)[i]) ? 1 : 0;
            return res;
        }
        for (int64_t i = 0; i < x->n; i++) {
            kG(res)[i] = (get_val_f(x, i) == get_val_f(y, i)) ? 1 : 0;
        }
        return res;
    }
    if (x->t < 0 && y->t > 0) {
        K res = ktn(KB_VEC, y->n);
        double sx = get_val_f(x, 0);
        for (int64_t i = 0; i < y->n; i++) {
            kG(res)[i] = (sx == get_val_f(y, i)) ? 1 : 0;
        }
        return res;
    }
    if (x->t > 0 && y->t < 0) {
        K res = ktn(KB_VEC, x->n);
        double sy = get_val_f(y, 0);
        for (int64_t i = 0; i < x->n; i++) {
            kG(res)[i] = (get_val_f(x, i) == sy) ? 1 : 0;
        }
        return res;
    }
    return kerror("type");
}

K v_neq(K x, K y) {
    K eq = v_eq(x, y);
    if (eq->t == XE_ERR) return eq;
    if (eq->t == KB_ATOM) {
        eq->g = eq->g ? 0 : 1;
        return eq;
    }
    for (int64_t i = 0; i < eq->n; i++) {
        kG(eq)[i] = kG(eq)[i] ? 0 : 1;
    }
    return eq;
}

K v_lt(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t < 0 && y->t < 0) return kb(get_val_f(x, 0) < get_val_f(y, 0));
    if (x->t > 0 && y->t > 0) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KB_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) kG(res)[i] = (get_val_f(x, i) < get_val_f(y, i)) ? 1 : 0;
        return res;
    }
    if (x->t < 0 && y->t > 0) {
        K res = ktn(KB_VEC, y->n);
        double sx = get_val_f(x, 0);
        for (int64_t i = 0; i < y->n; i++) kG(res)[i] = (sx < get_val_f(y, i)) ? 1 : 0;
        return res;
    }
    if (x->t > 0 && y->t < 0) {
        K res = ktn(KB_VEC, x->n);
        double sy = get_val_f(y, 0);
        for (int64_t i = 0; i < x->n; i++) kG(res)[i] = (get_val_f(x, i) < sy) ? 1 : 0;
        return res;
    }
    return kerror("type");
}

K v_gt(K x, K y) {
    return v_lt(y, x);
}

K v_lte(K x, K y) {
    K g = v_gt(x, y);
    if (g->t == XE_ERR) return g;
    if (g->t == KB_ATOM) { g->g = g->g ? 0 : 1; return g; }
    for (int64_t i = 0; i < g->n; i++) kG(g)[i] = kG(g)[i] ? 0 : 1;
    return g;
}

K v_gte(K x, K y) {
    K l = v_lt(x, y);
    if (l->t == XE_ERR) return l;
    if (l->t == KB_ATOM) { l->g = l->g ? 0 : 1; return l; }
    for (int64_t i = 0; i < l->n; i++) kG(l)[i] = kG(l)[i] ? 0 : 1;
    return l;
}

K v_match(K x, K y) {
    if (x == y) return kb(true);
    if (!x || !y) return kb(false);
    if (x->t != y->t) return kb(false);
    if (x->t < 0) {
        if (x->t == KS_ATOM) return kb(x->s == y->s);
        if (x->t == KF_ATOM) return kb(x->f == y->f);
        if (x->t == KJ_ATOM) return kb(x->j == y->j);
        if (x->t == KC_ATOM) return kb(x->c == y->c);
        return kb(get_val_f(x, 0) == get_val_f(y, 0));
    }
    if (x->n != y->n) return kb(false);
    if (x->t == KM_LIST) {
        for (int64_t i = 0; i < x->n; i++) {
            K m = v_match(kK(x)[i], kK(y)[i]);
            bool b = m->g;
            r0(m);
            if (!b) return kb(false);
        }
        return kb(true);
    }
    if (x->t == KS_VEC) {
        for (int64_t i = 0; i < x->n; i++) {
            if (kS(x)[i] != kS(y)[i]) return kb(false);
        }
        return kb(true);
    }
    for (int64_t i = 0; i < x->n; i++) {
        if (get_val_f(x, i) != get_val_f(y, i)) return kb(false);
    }
    return kb(true);
}

K v_not(K x) {
    if (!x) return kerror("type");
    if (x->t == KB_ATOM) return kb(!x->g);
    if (x->t == KB_VEC) {
        K res = ktn(KB_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) kG(res)[i] = kG(x)[i] ? 0 : 1;
        return res;
    }
    return kerror("type");
}

K v_and(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t == KB_ATOM && y->t == KB_ATOM) return kb(x->g && y->g);
    if (x->t == KB_VEC && y->t == KB_VEC) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KB_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) kG(res)[i] = (kG(x)[i] && kG(y)[i]) ? 1 : 0;
        return res;
    }
    return v_min(x, y);
}

K v_or(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t == KB_ATOM && y->t == KB_ATOM) return kb(x->g || y->g);
    if (x->t == KB_VEC && y->t == KB_VEC) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KB_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) kG(res)[i] = (kG(x)[i] || kG(y)[i]) ? 1 : 0;
        return res;
    }
    return v_max(x, y);
}

K v_count(K x) {
    if (!x) return kj(0);
    if (x->t < 0) return kj(1);
    if (x->t == XT_TABLE) {
        K dict = x->dict.k;
        K vals = dict->dict.v;
        if (vals && vals->n > 0) return kj(kK(vals)[0]->n);
        return kj(0);
    }
    if (x->t == XD_DICT) {
        return kj(x->dict.k->n);
    }
    return kj(x->n);
}

K v_first(K x) {
    if (!x) return kerror("type");
    if (x->t < 0) return r1(x);
    if (x->n == 0) return kerror("length");
    if (x->t == KM_LIST) return r1(kK(x)[0]);
    if (x->t == KF_VEC) return kf(kF(x)[0]);
    if (x->t == KJ_VEC) return kj(kJ(x)[0]);
    if (x->t == KI_VEC) return ki(kI(x)[0]);
    if (x->t == KS_VEC) return ks(kS(x)[0]);
    if (x->t == KC_VEC) return kc(kC(x)[0]);
    if (x->t == KB_VEC) return kb(kG(x)[0] != 0);
    return kerror("type");
}

K v_last(K x) {
    if (!x) return kerror("type");
    if (x->t < 0) return r1(x);
    if (x->n == 0) return kerror("length");
    int64_t idx = x->n - 1;
    if (x->t == KM_LIST) return r1(kK(x)[idx]);
    if (x->t == KF_VEC) return kf(kF(x)[idx]);
    if (x->t == KJ_VEC) return kj(kJ(x)[idx]);
    if (x->t == KI_VEC) return ki(kI(x)[idx]);
    if (x->t == KS_VEC) return ks(kS(x)[idx]);
    if (x->t == KC_VEC) return kc(kC(x)[idx]);
    if (x->t == KB_VEC) return kb(kG(x)[idx] != 0);
    return kerror("type");
}

K v_reverse(K x) {
    if (!x) return kerror("type");
    if (x->t < 0) return r1(x);
    K res = ktn(x->t, x->n);
    for (int64_t i = 0; i < x->n; i++) {
        int64_t src = x->n - 1 - i;
        if (x->t == KM_LIST) kK(res)[i] = r1(kK(x)[src]);
        else if (x->t == KF_VEC) kF(res)[i] = kF(x)[src];
        else if (x->t == KJ_VEC) kJ(res)[i] = kJ(x)[src];
        else if (x->t == KI_VEC) kI(res)[i] = kI(x)[src];
        else if (x->t == KS_VEC) kS(res)[i] = kS(x)[src];
        else if (x->t == KC_VEC) kC(res)[i] = kC(x)[src];
        else if (x->t == KB_VEC) kG(res)[i] = kG(x)[src];
    }
    return res;
}

K v_take(K x, K y) {
    if (!x || !y) return kerror("type");
    int64_t n = (x->t == KJ_ATOM) ? x->j : (int64_t)get_val_f(x, 0);
    if (y->t < 0) {
        /* Repeat scalar atom */
        int64_t count = (n < 0) ? -n : n;
        int vec_type = (y->t == KF_ATOM) ? KF_VEC :
                       (y->t == KJ_ATOM) ? KJ_VEC :
                       (y->t == KI_ATOM) ? KI_VEC :
                       (y->t == KS_ATOM) ? KS_VEC :
                       (y->t == KC_ATOM) ? KC_VEC : KB_VEC;
        K res = ktn(vec_type, count);
        for (int64_t i = 0; i < count; i++) {
            if (vec_type == KF_VEC) kF(res)[i] = y->f;
            else if (vec_type == KJ_VEC) kJ(res)[i] = y->j;
            else if (vec_type == KI_VEC) kI(res)[i] = y->i;
            else if (vec_type == KS_VEC) kS(res)[i] = y->s;
            else if (vec_type == KC_VEC) kC(res)[i] = y->c;
            else if (vec_type == KB_VEC) kG(res)[i] = y->g;
        }
        return res;
    }
    
    if (n >= 0) {
        int64_t take_cnt = (n <= y->n) ? n : y->n;
        K res = ktn(y->t, take_cnt);
        for (int64_t i = 0; i < take_cnt; i++) {
            if (y->t == KM_LIST) kK(res)[i] = r1(kK(y)[i]);
            else if (y->t == KF_VEC) kF(res)[i] = kF(y)[i];
            else if (y->t == KJ_VEC) kJ(res)[i] = kJ(y)[i];
            else if (y->t == KI_VEC) kI(res)[i] = kI(y)[i];
            else if (y->t == KS_VEC) kS(res)[i] = kS(y)[i];
            else if (y->t == KC_VEC) kC(res)[i] = kC(y)[i];
            else if (y->t == KB_VEC) kG(res)[i] = kG(y)[i];
        }
        return res;
    } else {
        int64_t count = -n;
        if (count > y->n) count = y->n;
        int64_t start = y->n - count;
        K res = ktn(y->t, count);
        for (int64_t i = 0; i < count; i++) {
            int64_t src = start + i;
            if (y->t == KM_LIST) kK(res)[i] = r1(kK(y)[src]);
            else if (y->t == KF_VEC) kF(res)[i] = kF(y)[src];
            else if (y->t == KJ_VEC) kJ(res)[i] = kJ(y)[src];
            else if (y->t == KI_VEC) kI(res)[i] = kI(y)[src];
            else if (y->t == KS_VEC) kS(res)[i] = kS(y)[src];
            else if (y->t == KC_VEC) kC(res)[i] = kC(y)[src];
            else if (y->t == KB_VEC) kG(res)[i] = kG(y)[src];
        }
        return res;
    }
}

K v_drop(K x, K y) {
    if (!x || !y) return kerror("type");
    int64_t n = (x->t == KJ_ATOM) ? x->j : (int64_t)get_val_f(x, 0);
    if (y->t < 0) return kerror("type");
    if (n >= 0) {
        if (n >= y->n) return ktn(y->t, 0);
        int64_t rem = y->n - n;
        K res = ktn(y->t, rem);
        for (int64_t i = 0; i < rem; i++) {
            int64_t src = n + i;
            if (y->t == KM_LIST) kK(res)[i] = r1(kK(y)[src]);
            else if (y->t == KF_VEC) kF(res)[i] = kF(y)[src];
            else if (y->t == KJ_VEC) kJ(res)[i] = kJ(y)[src];
            else if (y->t == KI_VEC) kI(res)[i] = kI(y)[src];
            else if (y->t == KS_VEC) kS(res)[i] = kS(y)[src];
            else if (y->t == KC_VEC) kC(res)[i] = kC(y)[src];
            else if (y->t == KB_VEC) kG(res)[i] = kG(y)[src];
        }
        return res;
    } else {
        int64_t drop_cnt = -n;
        if (drop_cnt >= y->n) return ktn(y->t, 0);
        int64_t rem = y->n - drop_cnt;
        K res = ktn(y->t, rem);
        for (int64_t i = 0; i < rem; i++) {
            if (y->t == KM_LIST) kK(res)[i] = r1(kK(y)[i]);
            else if (y->t == KF_VEC) kF(res)[i] = kF(y)[i];
            else if (y->t == KJ_VEC) kJ(res)[i] = kJ(y)[i];
            else if (y->t == KI_VEC) kI(res)[i] = kI(y)[i];
            else if (y->t == KS_VEC) kS(res)[i] = kS(y)[i];
            else if (y->t == KC_VEC) kC(res)[i] = kC(y)[i];
            else if (y->t == KB_VEC) kG(res)[i] = kG(y)[i];
        }
        return res;
    }
}

K v_join(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t > 0 && y->t > 0 && x->t == y->t) {
        K res = ktn(x->t, x->n + y->n);
        if (x->t == KF_VEC) {
            memcpy(kF(res), kF(x), (size_t)x->n * sizeof(double));
            memcpy(kF(res) + x->n, kF(y), (size_t)y->n * sizeof(double));
        } else if (x->t == KJ_VEC) {
            memcpy(kJ(res), kJ(x), (size_t)x->n * sizeof(int64_t));
            memcpy(kJ(res) + x->n, kJ(y), (size_t)y->n * sizeof(int64_t));
        } else if (x->t == KS_VEC) {
            memcpy(kS(res), kS(x), (size_t)x->n * sizeof(S));
            memcpy(kS(res) + x->n, kS(y), (size_t)y->n * sizeof(S));
        } else if (x->t == KC_VEC) {
            memcpy(kC(res), kC(x), (size_t)x->n);
            memcpy(kC(res) + x->n, kC(y), (size_t)y->n);
        } else if (x->t == KM_LIST) {
            for (int64_t i = 0; i < x->n; i++) kK(res)[i] = r1(kK(x)[i]);
            for (int64_t i = 0; i < y->n; i++) kK(res)[x->n + i] = r1(kK(y)[i]);
        }
        return res;
    }
    /* Mixed join fallback */
    int64_t nx = (x->t < 0) ? 1 : x->n;
    int64_t ny = (y->t < 0) ? 1 : y->n;
    K res = k_mixed_list(nx + ny);
    if (x->t < 0) {
        kK(res)[0] = r1(x);
    } else {
        for (int64_t i = 0; i < x->n; i++) kK(res)[i] = (x->t == KM_LIST) ? r1(kK(x)[i]) : kf(get_val_f(x, i));
    }
    if (y->t < 0) {
        kK(res)[nx] = r1(y);
    } else {
        for (int64_t i = 0; i < y->n; i++) kK(res)[nx + i] = (y->t == KM_LIST) ? r1(kK(y)[i]) : kf(get_val_f(y, i));
    }
    return res;
}

K v_find(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t < 0) return kerror("type");
    for (int64_t i = 0; i < x->n; i++) {
        if (x->t == KS_VEC && y->t == KS_ATOM) {
            if (kS(x)[i] == y->s) return kj(i);
        } else if (x->t == KF_VEC && y->t == KF_ATOM) {
            if (kF(x)[i] == y->f) return kj(i);
        } else if (x->t == KJ_VEC && y->t == KJ_ATOM) {
            if (kJ(x)[i] == y->j) return kj(i);
        }
    }
    return kj(x->n);
}

K v_at(K x, K y) {
    if (!x || !y) return kerror("type");
    if (x->t == XD_DICT) {
        K keys = x->dict.k;
        K vals = x->dict.v;
        int64_t idx = -1;
        for (int64_t i = 0; i < keys->n; i++) {
            if (keys->t == KS_VEC && y->t == KS_ATOM && kS(keys)[i] == y->s) {
                idx = i;
                break;
            }
        }
        if (idx >= 0) return r1(kK(vals)[idx]);
        return kerror("find");
    }
    if (x->t == XT_TABLE) {
        /* Table indexed by column symbol or row int */
        if (y->t == KS_ATOM) {
            return v_at(x->dict.k, y);
        }
    }
    if (y->t < 0) {
        int64_t idx = (y->t == KJ_ATOM) ? y->j : (int64_t)get_val_f(y, 0);
        if (idx < 0 || idx >= x->n) return kerror("index");
        if (x->t == KM_LIST) return r1(kK(x)[idx]);
        if (x->t == KF_VEC) return kf(kF(x)[idx]);
        if (x->t == KJ_VEC) return kj(kJ(x)[idx]);
        if (x->t == KS_VEC) return ks(kS(x)[idx]);
        if (x->t == KC_VEC) return kc(kC(x)[idx]);
    }
    return kerror("type");
}

K v_dot(K x, K y) {
    if (!x || !y) return kerror("type");
    if (y->t == KM_LIST && y->n > 0) {
        K curr = r1(x);
        for (int64_t i = 0; i < y->n; i++) {
            K next = v_at(curr, kK(y)[i]);
            r0(curr);
            curr = next;
            if (curr->t == XE_ERR) return curr;
        }
        return curr;
    }
    return v_at(x, y);
}

K v_enlist(K x) {
    if (!x) return kerror("type");
    K res = k_mixed_list(1);
    kK(res)[0] = r1(x);
    return res;
}

K v_flip(K x) {
    if (!x) return kerror("type");
    if (x->t == XD_DICT) return xT(r1(x));
    if (x->t == XT_TABLE) return r1(x->dict.k);
    return kerror("type");
}

K v_keys(K x) {
    if (!x) return kerror("type");
    if (x->t == XD_DICT) return r1(x->dict.k);
    if (x->t == XT_TABLE) return v_keys(x->dict.k);
    return kerror("type");
}

K v_values(K x) {
    if (!x) return kerror("type");
    if (x->t == XD_DICT) return r1(x->dict.v);
    if (x->t == XT_TABLE) return v_values(x->dict.k);
    return kerror("type");
}
