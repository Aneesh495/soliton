#include <stdlib.h>
#include <string.h>
#include "soliton/q.h"
#include "soliton/eval.h"

static K apply_fn(Env* env, K fn, K arg1, K arg2) {
    if (!fn) return kerror("type");
    if (fn->t == XV_PRIM) {
        /* Builtin primitive verb */
        switch (fn->prim.id) {
            case '+': return v_add(arg1, arg2);
            case '-': return (arg2 ? v_sub(arg1, arg2) : v_neg(arg1));
            case '*': return v_mul(arg1, arg2);
            case '%': return v_div(arg1, arg2);
            case '&': return v_min(arg1, arg2);
            case '|': return v_max(arg1, arg2);
            case '=': return v_eq(arg1, arg2);
            case '<': return v_lt(arg1, arg2);
            case '>': return v_gt(arg1, arg2);
            case ',': return v_join(arg1, arg2);
            case '#': return v_take(arg1, arg2);
            case '_': return v_drop(arg1, arg2);
            default: break;
        }
    }
    if (fn->t == XF_LAMBDA) {
        /* Create closure environment */
        Env* local = env_new(env);
        if (fn->lambda.params && fn->lambda.params->n > 0) {
            S p0 = kS(fn->lambda.params)[0];
            env_set(local, p0, r1(arg1));
            if (fn->lambda.params->n > 1 && arg2) {
                S p1 = kS(fn->lambda.params)[1];
                env_set(local, p1, r1(arg2));
            }
        } else {
            /* Default implicit variables x, y */
            env_set(local, ss("x"), r1(arg1));
            if (arg2) env_set(local, ss("y"), r1(arg2));
        }
        K res = q_eval_string(local, fn->lambda.body);
        env_free(local);
        return res;
    }
    return kerror("eval");
}

K adv_each(Env* env, K fn, K list) {
    if (!fn || !list) return kerror("type");
    if (list->t < 0) {
        return apply_fn(env, fn, list, NULL);
    }
    K res = k_mixed_list(list->n);
    bool all_float = true;
    for (int64_t i = 0; i < list->n; i++) {
        K item = (list->t == KM_LIST) ? r1(kK(list)[i]) : v_at(list, kj(i));
        K out = apply_fn(env, fn, item, NULL);
        r0(item);
        if (out->t != KF_ATOM) all_float = false;
        kK(res)[i] = out;
    }
    if (all_float && list->n > 0) {
        K fvec = ktn(KF_VEC, list->n);
        for (int64_t i = 0; i < list->n; i++) {
            kF(fvec)[i] = kK(res)[i]->f;
        }
        r0(res);
        return fvec;
    }
    return res;
}

K adv_over(Env* env, K fn, K list) {
    if (!fn || !list) return kerror("type");
    if (list->t < 0) return r1(list);
    if (list->n == 0) return kerror("length");
    
    K acc = (list->t == KM_LIST) ? r1(kK(list)[0]) : v_at(list, kj(0));
    for (int64_t i = 1; i < list->n; i++) {
        K next = (list->t == KM_LIST) ? r1(kK(list)[i]) : v_at(list, kj(i));
        K step = apply_fn(env, fn, acc, next);
        r0(acc);
        r0(next);
        acc = step;
        if (acc->t == XE_ERR) break;
    }
    return acc;
}

K adv_scan(Env* env, K fn, K list) {
    if (!fn || !list) return kerror("type");
    if (list->t < 0) return r1(list);
    if (list->n == 0) return ktn(KF_VEC, 0);
    
    K res = k_mixed_list(list->n);
    K acc = (list->t == KM_LIST) ? r1(kK(list)[0]) : v_at(list, kj(0));
    kK(res)[0] = r1(acc);
    
    bool all_float = (acc->t == KF_ATOM);
    for (int64_t i = 1; i < list->n; i++) {
        K next = (list->t == KM_LIST) ? r1(kK(list)[i]) : v_at(list, kj(i));
        K step = apply_fn(env, fn, acc, next);
        r0(acc);
        r0(next);
        acc = step;
        if (acc->t != KF_ATOM) all_float = false;
        kK(res)[i] = r1(acc);
    }
    r0(acc);
    
    if (all_float) {
        K fvec = ktn(KF_VEC, list->n);
        for (int64_t i = 0; i < list->n; i++) {
            kF(fvec)[i] = kK(res)[i]->f;
        }
        r0(res);
        return fvec;
    }
    return res;
}

K adv_prior(Env* env, K fn, K list) {
    if (!fn || !list) return kerror("type");
    if (list->t < 0 || list->n < 2) return kerror("length");
    
    K res = k_mixed_list(list->n - 1);
    bool all_float = true;
    for (int64_t i = 1; i < list->n; i++) {
        K prev = (list->t == KM_LIST) ? r1(kK(list)[i - 1]) : v_at(list, kj(i - 1));
        K curr = (list->t == KM_LIST) ? r1(kK(list)[i]) : v_at(list, kj(i));
        K out = apply_fn(env, fn, curr, prev);
        r0(prev);
        r0(curr);
        if (out->t != KF_ATOM) all_float = false;
        kK(res)[i - 1] = out;
    }
    if (all_float && list->n > 1) {
        K fvec = ktn(KF_VEC, list->n - 1);
        for (int64_t i = 0; i < list->n - 1; i++) {
            kF(fvec)[i] = kK(res)[i]->f;
        }
        r0(res);
        return fvec;
    }
    return res;
}

K adv_each_right(Env* env, K fn, K x, K y) {
    if (!fn || !x || !y) return kerror("type");
    if (y->t < 0) return apply_fn(env, fn, x, y);
    
    K res = k_mixed_list(y->n);
    bool all_float = true;
    for (int64_t i = 0; i < y->n; i++) {
        K item = (y->t == KM_LIST) ? r1(kK(y)[i]) : v_at(y, kj(i));
        K out = apply_fn(env, fn, x, item);
        r0(item);
        if (out->t != KF_ATOM) all_float = false;
        kK(res)[i] = out;
    }
    if (all_float && y->n > 0) {
        K fvec = ktn(KF_VEC, y->n);
        for (int64_t i = 0; i < y->n; i++) kF(fvec)[i] = kK(res)[i]->f;
        r0(res);
        return fvec;
    }
    return res;
}

K adv_each_left(Env* env, K fn, K x, K y) {
    if (!fn || !x || !y) return kerror("type");
    if (x->t < 0) return apply_fn(env, fn, x, y);
    
    K res = k_mixed_list(x->n);
    bool all_float = true;
    for (int64_t i = 0; i < x->n; i++) {
        K item = (x->t == KM_LIST) ? r1(kK(x)[i]) : v_at(x, kj(i));
        K out = apply_fn(env, fn, item, y);
        r0(item);
        if (out->t != KF_ATOM) all_float = false;
        kK(res)[i] = out;
    }
    if (all_float && x->n > 0) {
        K fvec = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) kF(fvec)[i] = kK(res)[i]->f;
        r0(res);
        return fvec;
    }
    return res;
}
