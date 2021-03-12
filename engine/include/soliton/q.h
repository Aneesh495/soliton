#ifndef SOLITON_Q_H
#define SOLITON_Q_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Type codes matching classic kdb+/q conventions.
 * Negative codes denote scalar atoms; positive codes denote homogeneous vectors.
 * Zero denotes a general mixed list of K objects.
 */
#define KB_ATOM   (-1)
#define KG_ATOM   (-4)
#define KH_ATOM   (-5)
#define KI_ATOM   (-6)
#define KJ_ATOM   (-7)
#define KE_ATOM   (-8)
#define KF_ATOM   (-9)
#define KC_ATOM   (-10)
#define KS_ATOM   (-11)
#define KP_ATOM   (-12)  /* timestamp */
#define KD_ATOM   (-14)  /* date */
#define KT_ATOM   (-15)  /* time */
#define KN_ATOM   (-16)  /* timespan */

#define KM_LIST     0   /* mixed list */
#define KB_VEC      1   /* boolean vector */
#define KG_VEC      4   /* byte vector */
#define KH_VEC      5   /* int16 vector */
#define KI_VEC      6   /* int32 vector */
#define KJ_VEC      7   /* int64 vector */
#define KE_VEC      8   /* float32 vector */
#define KF_VEC      9   /* float64 vector */
#define KC_VEC     10   /* char vector (string) */
#define KS_VEC     11   /* symbol vector */
#define KP_VEC     12   /* timestamp vector */
#define KD_VEC     14   /* date vector */
#define KT_VEC     15   /* time vector */
#define KN_VEC     16   /* timespan vector */

#define XT_TABLE   98   /* table (flip of dictionary) */
#define XD_DICT    99   /* dictionary or keyed table */
#define XF_LAMBDA 100   /* lambda closure */
#define XV_PRIM   101   /* primitive operator */
#define XP_PROJ   102   /* projection */
#define XE_ERR   (-128) /* error object */

typedef char* S;

typedef struct k0 {
    int8_t m;   /* memory control / flags */
    int8_t t;   /* type code */
    int16_t u;  /* attribute flags (1: sorted, 2: unique, 3: parted, 4: grouped) */
    int32_t r;  /* reference count */
    union {
        /* Atoms */
        uint8_t g;       /* byte / bool */
        int16_t h;       /* short */
        int32_t i;       /* int */
        int64_t j;       /* long */
        float e;         /* real */
        double f;        /* float */
        char c;          /* char */
        S s;             /* symbol pointer */
        
        /* Vectors and Complex Structures */
        struct {
            int64_t n;   /* element count */
            union {
                uint8_t G0[1];
                int16_t H0[1];
                int32_t I0[1];
                int64_t J0[1];
                float E0[1];
                double F0[1];
                char C0[1];
                S S0[1];
                struct k0* K0[1];
            };
        };
        
        /* Dictionary / Table representation */
        struct {
            struct k0* k; /* keys (symbol list or key table) */
            struct k0* v; /* values (value list or value table) */
        } dict;
        
        /* Lambda / Closure representation */
        struct {
            char* body;
            struct k0* params;
            struct k0* env;
        } lambda;
        
        /* Projection representation */
        struct {
            struct k0* fn;
            struct k0* args;
        } proj;
        
        /* Primitive function */
        struct {
            int id;
            const char* name;
            int arity;
        } prim;
    };
} *K;

/* Vector accessors */
#define kG(x) ((uint8_t*)((x)->G0))
#define kH(x) ((int16_t*)((x)->H0))
#define kI(x) ((int32_t*)((x)->I0))
#define kJ(x) ((int64_t*)((x)->J0))
#define kE(x) ((float*)((x)->E0))
#define kF(x) ((double*)((x)->F0))
#define kC(x) ((char*)((x)->C0))
#define kS(x) ((S*)((x)->S0))
#define kK(x) ((K*)((x)->K0))

/* Reference counting helpers */
static inline K r1(K x) {
    if (x) {
        x->r++;
    }
    return x;
}

void r0(K x);

/* Constructors for scalar atoms */
K kb(bool b);
K kg(uint8_t b);
K kh(int16_t h);
K ki(int32_t i);
K kj(int64_t j);
K ke(float e);
K kf(double f);
K kc(char c);
K ks(const char* s);
K kerror(const char* msg);

/* Vector and collection constructors */
K ktn(int type, int64_t length);
K kp(const char* s);       /* string (char vector) */
K kpn(const char* s, int64_t n);
K xD(K keys, K vals);      /* dictionary */
K xT(K dict);              /* table (flip of dictionary) */
K k_mixed_list(int64_t n);

/* Symbol interning table */
void sym_init(void);
void sym_cleanup(void);
S sn(const char* s, size_t n);
S ss(const char* s);

/* Memory subsystem */
void k_mem_init(void);
void k_mem_cleanup(void);
size_t k_allocated_bytes(void);

#ifdef __cplusplus
}
#endif

#endif /* SOLITON_Q_H */
