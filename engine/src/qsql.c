#include <stdlib.h>
#include <string.h>
#include "soliton/q.h"
#include "soliton/eval.h"

K table_get_column(K tbl, S col_name);
K table_filter_rows(K tbl, K bool_mask);
int64_t table_row_count(K tbl);

K qsql_select(Env* env, K tbl, K cols, K by_cols, K where_clause) {
    (void)env;
    (void)by_cols;
    if (!tbl) return kerror("type");
    K src_tbl = tbl;
    
    /* Apply where clause filter if present */
    if (where_clause && where_clause->t == KB_VEC) {
        src_tbl = table_filter_rows(tbl, where_clause);
    } else {
        src_tbl = r1(tbl);
    }
    
    if (!cols || cols->t != KS_VEC) {
        return src_tbl; /* select * from tbl */
    }
    
    /* Project selected columns */
    K proj_keys = ktn(KS_VEC, cols->n);
    K proj_vals = k_mixed_list(cols->n);
    for (int64_t i = 0; i < cols->n; i++) {
        S col_name = kS(cols)[i];
        kS(proj_keys)[i] = col_name;
        K c = table_get_column(src_tbl, col_name);
        kK(proj_vals)[i] = c;
    }
    r0(src_tbl);
    
    return xT(xD(proj_keys, proj_vals));
}

K qsql_update(Env* env, K tbl, K cols, K by_cols, K where_clause) {
    (void)env;
    (void)by_cols;
    (void)where_clause;
    if (!tbl || !cols || cols->t != XD_DICT) return kerror("type");
    K dict = (tbl->t == XT_TABLE) ? tbl->dict.k : tbl;
    if (dict->t != XD_DICT) return kerror("type");
    
    K keys = dict->dict.k;
    K vals = dict->dict.v;
    K upd_keys = cols->dict.k;
    K upd_vals = cols->dict.v;
    
    K new_vals = k_mixed_list(vals->n);
    for (int64_t i = 0; i < keys->n; i++) {
        S sym = kS(keys)[i];
        K replacement = NULL;
        for (int64_t u = 0; u < upd_keys->n; u++) {
            if (kS(upd_keys)[u] == sym) {
                replacement = kK(upd_vals)[u];
                break;
            }
        }
        kK(new_vals)[i] = replacement ? r1(replacement) : r1(kK(vals)[i]);
    }
    
    return xT(xD(r1(keys), new_vals));
}

K qsql_delete(Env* env, K tbl, K cols, K where_clause) {
    (void)env;
    if (!tbl) return kerror("type");
    if (where_clause && where_clause->t == KB_VEC) {
        /* Invert boolean mask to drop rows */
        K inv_mask = ktn(KB_VEC, where_clause->n);
        for (int64_t i = 0; i < where_clause->n; i++) {
            kG(inv_mask)[i] = kG(where_clause)[i] ? 0 : 1;
        }
        K res = table_filter_rows(tbl, inv_mask);
        r0(inv_mask);
        return res;
    }
    if (cols && cols->t == KS_VEC) {
        /* Drop columns */
        K dict = (tbl->t == XT_TABLE) ? tbl->dict.k : tbl;
        K keys = dict->dict.k;
        K vals = dict->dict.v;
        
        int64_t keep = 0;
        for (int64_t i = 0; i < keys->n; i++) {
            bool drop = false;
            for (int64_t d = 0; d < cols->n; d++) {
                if (kS(keys)[i] == kS(cols)[d]) { drop = true; break; }
            }
            if (!drop) keep++;
        }
        K new_keys = ktn(KS_VEC, keep);
        K new_vals = k_mixed_list(keep);
        int64_t dst = 0;
        for (int64_t i = 0; i < keys->n; i++) {
            bool drop = false;
            for (int64_t d = 0; d < cols->n; d++) {
                if (kS(keys)[i] == kS(cols)[d]) { drop = true; break; }
            }
            if (!drop) {
                kS(new_keys)[dst] = kS(keys)[i];
                kK(new_vals)[dst] = r1(kK(vals)[i]);
                dst++;
            }
        }
        return xT(xD(new_keys, new_vals));
    }
    return r1(tbl);
}

K qsql_exec(Env* env, K tbl, K col, K by_cols, K where_clause) {
    (void)env;
    (void)by_cols;
    if (!tbl || !col || col->t != KS_ATOM) return kerror("type");
    K filtered = (where_clause && where_clause->t == KB_VEC) ? table_filter_rows(tbl, where_clause) : r1(tbl);
    K res = table_get_column(filtered, col->s);
    r0(filtered);
    return res;
}
