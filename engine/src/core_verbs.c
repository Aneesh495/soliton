#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "soliton/q.h"
#include "soliton/eval.h"

static double get_scalar_float(K x) {
    if (x->t == KF_ATOM) return x->f;
    if (x->t == KJ_ATOM) return (double)x->j;
    if (x->t == KI_ATOM) return (double)x->i;
    if (x->t == KH_ATOM) return (double)x->h;
    if (x->t == KB_ATOM || x->t == KG_ATOM) return (double)x->g;
    if (x->t == KE_ATOM) return (double)x->e;
    return 0.0;
}

static int64_t get_scalar_int(K x) {
    if (x->t == KJ_ATOM) return x->j;
    if (x->t == KI_ATOM) return (int64_t)x->i;
    if (x->t == KH_ATOM) return (int64_t)x->h;
    if (x->t == KB_ATOM || x->t == KG_ATOM) return (int64_t)x->g;
    if (x->t == KF_ATOM) return (int64_t)x->f;
    return 0;
}

static bool is_numeric_atom(K x) {
    return (x->t == KB_ATOM || x->t == KG_ATOM || x->t == KH_ATOM || 
            x->t == KI_ATOM || x->t == KJ_ATOM || x->t == KE_ATOM || x->t == KF_ATOM);
}

static bool is_numeric_vec(K x) {
    return (x->t == KB_VEC || x->t == KG_VEC || x->t == KH_VEC || 
            x->t == KI_VEC || x->t == KJ_VEC || x->t == KE_VEC || x->t == KF_VEC);
}

static double get_vec_float(K x, int64_t i) {
    if (x->t == KF_VEC) return kF(x)[i];
    if (x->t == KJ_VEC) return (double)kJ(x)[i];
    if (x->t == KI_VEC) return (double)kI(x)[i];
    if (x->t == KH_VEC) return (double)kH(x)[i];
    if (x->t == KB_VEC || x->t == KG_VEC) return (double)kG(x)[i];
    if (x->t == KE_VEC) return (double)kE(x)[i];
    return 0.0;
}

K v_add(K x, K y) {
    if (!x || !y) return kerror("type");
    
    /* Atom + Atom */
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        if (x->t == KF_ATOM || y->t == KF_ATOM || x->t == KE_ATOM || y->t == KE_ATOM) {
            return kf(get_scalar_float(x) + get_scalar_float(y));
        }
        return kj(get_scalar_int(x) + get_scalar_int(y));
    }
    
    /* Atom + Vec */
    if (is_numeric_atom(x) && is_numeric_vec(y)) {
        K res = ktn(KF_VEC, y->n);
        double sx = get_scalar_float(x);
        for (int64_t i = 0; i < y->n; i++) {
            kF(res)[i] = sx + get_vec_float(y, i);
        }
        return res;
    }
    
    /* Vec + Atom */
    if (is_numeric_vec(x) && is_numeric_atom(y)) {
        K res = ktn(KF_VEC, x->n);
        double sy = get_scalar_float(y);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = get_vec_float(x, i) + sy;
        }
        return res;
    }
    
    /* Vec + Vec */
    if (is_numeric_vec(x) && is_numeric_vec(y)) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = get_vec_float(x, i) + get_vec_float(y, i);
        }
        return res;
    }
    
    return kerror("type");
}

K v_sub(K x, K y) {
    if (!x || !y) return kerror("type");
    
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        if (x->t == KF_ATOM || y->t == KF_ATOM || x->t == KE_ATOM || y->t == KE_ATOM) {
            return kf(get_scalar_float(x) - get_scalar_float(y));
        }
        return kj(get_scalar_int(x) - get_scalar_int(y));
    }
    if (is_numeric_atom(x) && is_numeric_vec(y)) {
        K res = ktn(KF_VEC, y->n);
        double sx = get_scalar_float(x);
        for (int64_t i = 0; i < y->n; i++) {
            kF(res)[i] = sx - get_vec_float(y, i);
        }
        return res;
    }
    if (is_numeric_vec(x) && is_numeric_atom(y)) {
        K res = ktn(KF_VEC, x->n);
        double sy = get_scalar_float(y);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = get_vec_float(x, i) - sy;
        }
        return res;
    }
    if (is_numeric_vec(x) && is_numeric_vec(y)) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = get_vec_float(x, i) - get_vec_float(y, i);
        }
        return res;
    }
    return kerror("type");
}

K v_mul(K x, K y) {
    if (!x || !y) return kerror("type");
    
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        if (x->t == KF_ATOM || y->t == KF_ATOM || x->t == KE_ATOM || y->t == KE_ATOM) {
            return kf(get_scalar_float(x) * get_scalar_float(y));
        }
        return kj(get_scalar_int(x) * get_scalar_int(y));
    }
    if (is_numeric_atom(x) && is_numeric_vec(y)) {
        K res = ktn(KF_VEC, y->n);
        double sx = get_scalar_float(x);
        for (int64_t i = 0; i < y->n; i++) {
            kF(res)[i] = sx * get_vec_float(y, i);
        }
        return res;
    }
    if (is_numeric_vec(x) && is_numeric_atom(y)) {
        K res = ktn(KF_VEC, x->n);
        double sy = get_scalar_float(y);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = get_vec_float(x, i) * sy;
        }
        return res;
    }
    if (is_numeric_vec(x) && is_numeric_vec(y)) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = get_vec_float(x, i) * get_vec_float(y, i);
        }
        return res;
    }
    return kerror("type");
}

K v_div(K x, K y) {
    if (!x || !y) return kerror("type");
    
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        double denom = get_scalar_float(y);
        if (denom == 0.0) return kf(INFINITY);
        return kf(get_scalar_float(x) / denom);
    }
    if (is_numeric_atom(x) && is_numeric_vec(y)) {
        K res = ktn(KF_VEC, y->n);
        double sx = get_scalar_float(x);
        for (int64_t i = 0; i < y->n; i++) {
            double d = get_vec_float(y, i);
            kF(res)[i] = (d == 0.0) ? INFINITY : (sx / d);
        }
        return res;
    }
    if (is_numeric_vec(x) && is_numeric_atom(y)) {
        K res = ktn(KF_VEC, x->n);
        double sy = get_scalar_float(y);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = (sy == 0.0) ? INFINITY : (get_vec_float(x, i) / sy);
        }
        return res;
    }
    if (is_numeric_vec(x) && is_numeric_vec(y)) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            double d = get_vec_float(y, i);
            kF(res)[i] = (d == 0.0) ? INFINITY : (get_vec_float(x, i) / d);
        }
        return res;
    }
    return kerror("type");
}

K v_mod(K x, K y) {
    if (!x || !y) return kerror("type");
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        int64_t m = get_scalar_int(x);
        int64_t n = get_scalar_int(y);
        if (m == 0) return kj(0);
        return kj(n % m);
    }
    return kerror("type");
}

K v_min(K x, K y) {
    if (!x || !y) return kerror("type");
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        double fx = get_scalar_float(x);
        double fy = get_scalar_float(y);
        return kf(fx < fy ? fx : fy);
    }
    if (is_numeric_vec(x) && is_numeric_vec(y)) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            double a = get_vec_float(x, i);
            double b = get_vec_float(y, i);
            kF(res)[i] = (a < b) ? a : b;
        }
        return res;
    }
    return kerror("type");
}

K v_max(K x, K y) {
    if (!x || !y) return kerror("type");
    if (is_numeric_atom(x) && is_numeric_atom(y)) {
        double fx = get_scalar_float(x);
        double fy = get_scalar_float(y);
        return kf(fx > fy ? fx : fy);
    }
    if (is_numeric_vec(x) && is_numeric_vec(y)) {
        if (x->n != y->n) return kerror("length");
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            double a = get_vec_float(x, i);
            double b = get_vec_float(y, i);
            kF(res)[i] = (a > b) ? a : b;
        }
        return res;
    }
    return kerror("type");
}

K v_neg(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) {
        if (x->t == KF_ATOM) return kf(-x->f);
        return kj(-get_scalar_int(x));
    }
    if (is_numeric_vec(x)) {
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = -get_vec_float(x, i);
        }
        return res;
    }
    return kerror("type");
}

K v_abs(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) {
        double v = fabs(get_scalar_float(x));
        return (x->t == KF_ATOM) ? kf(v) : kj((int64_t)v);
    }
    if (is_numeric_vec(x)) {
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = fabs(get_vec_float(x, i));
        }
        return res;
    }
    return kerror("type");
}

K v_sqrt(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) {
        return kf(sqrt(get_scalar_float(x)));
    }
    if (is_numeric_vec(x)) {
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = sqrt(get_vec_float(x, i));
        }
        return res;
    }
    return kerror("type");
}

K v_exp(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) {
        return kf(exp(get_scalar_float(x)));
    }
    if (is_numeric_vec(x)) {
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = exp(get_vec_float(x, i));
        }
        return res;
    }
    return kerror("type");
}

K v_log(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) {
        return kf(log(get_scalar_float(x)));
    }
    if (is_numeric_vec(x)) {
        K res = ktn(KF_VEC, x->n);
        for (int64_t i = 0; i < x->n; i++) {
            kF(res)[i] = log(get_vec_float(x, i));
        }
        return res;
    }
    return kerror("type");
}

K v_til(K x) {
    if (!x || !is_numeric_atom(x)) return kerror("type");
    int64_t n = get_scalar_int(x);
    if (n < 0) n = 0;
    K res = ktn(KJ_VEC, n);
    for (int64_t i = 0; i < n; i++) {
        kJ(res)[i] = i;
    }
    return res;
}

K v_sum(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) return r1(x);
    if (is_numeric_vec(x)) {
        double acc = 0.0;
        for (int64_t i = 0; i < x->n; i++) {
            acc += get_vec_float(x, i);
        }
        return kf(acc);
    }
    return kerror("type");
}

K v_avg(K x) {
    if (!x) return kerror("type");
    if (is_numeric_atom(x)) return kf(get_scalar_float(x));
    if (is_numeric_vec(x)) {
        if (x->n == 0) return kf(0.0);
        double acc = 0.0;
        for (int64_t i = 0; i < x->n; i++) {
            acc += get_vec_float(x, i);
        }
        return kf(acc / (double)x->n);
    }
    return kerror("type");
}

K v_var(K x) {
    if (!x || !is_numeric_vec(x) || x->n < 2) return kf(0.0);
    double m = 0.0;
    for (int64_t i = 0; i < x->n; i++) m += get_vec_float(x, i);
    m /= (double)x->n;
    double ss = 0.0;
    for (int64_t i = 0; i < x->n; i++) {
        double d = get_vec_float(x, i) - m;
        ss += d * d;
    }
    return kf(ss / (double)(x->n - 1));
}

K v_dev(K x) {
    K v = v_var(x);
    double f = v->f;
    r0(v);
    return kf(sqrt(f));
}

static int compare_doubles(const void* a, const void* b) {
    double fa = *(const double*)a;
    double fb = *(const double*)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

K v_med(K x) {
    if (!x || !is_numeric_vec(x) || x->n == 0) return kf(0.0);
    double* buf = (double*)malloc(sizeof(double) * (size_t)x->n);
    for (int64_t i = 0; i < x->n; i++) buf[i] = get_vec_float(x, i);
    qsort(buf, (size_t)x->n, sizeof(double), compare_doubles);
    double res;
    if (x->n % 2 == 1) {
        res = buf[x->n / 2];
    } else {
        res = 0.5 * (buf[(x->n / 2) - 1] + buf[x->n / 2]);
    }
    free(buf);
    return kf(res);
}
