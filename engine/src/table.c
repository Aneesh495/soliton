#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "soliton/q.h"
#include "soliton/eval.h"

K table_get_column(K tbl, S col_name) {
    if (!tbl) return kerror("type");
    K dict = (tbl->t == XT_TABLE) ? tbl->dict.k : tbl;
    if (dict->t != XD_DICT) return kerror("type");
    
    K keys = dict->dict.k;
    K vals = dict->dict.v;
    if (keys->t != KS_VEC || vals->t != KM_LIST) return kerror("type");
    
    for (int64_t i = 0; i < keys->n; i++) {
        if (kS(keys)[i] == col_name) {
            return r1(kK(vals)[i]);
        }
    }
    return kerror("column");
}

int64_t table_row_count(K tbl) {
    if (!tbl) return 0;
    K dict = (tbl->t == XT_TABLE) ? tbl->dict.k : tbl;
    if (dict->t != XD_DICT) return 0;
    K vals = dict->dict.v;
    if (vals && vals->n > 0 && kK(vals)[0]) {
        return kK(vals)[0]->n;
    }
    return 0;
}

K table_filter_rows(K tbl, K bool_mask) {
    if (!tbl || !bool_mask) return kerror("type");
    if (bool_mask->t != KB_VEC) return kerror("type");
    
    K dict = (tbl->t == XT_TABLE) ? tbl->dict.k : tbl;
    if (dict->t != XD_DICT) return kerror("type");
    K keys = dict->dict.k;
    K vals = dict->dict.v;
    
    int64_t keep_cnt = 0;
    for (int64_t i = 0; i < bool_mask->n; i++) {
        if (kG(bool_mask)[i]) keep_cnt++;
    }
    
    K new_vals = k_mixed_list(vals->n);
    for (int64_t c = 0; c < vals->n; c++) {
        K col = kK(vals)[c];
        K new_col = ktn(col->t, keep_cnt);
        int64_t dst = 0;
        for (int64_t r = 0; r < bool_mask->n; r++) {
            if (kG(bool_mask)[r]) {
                if (col->t == KF_VEC) kF(new_col)[dst] = kF(col)[r];
                else if (col->t == KJ_VEC) kJ(new_col)[dst] = kJ(col)[r];
                else if (col->t == KI_VEC) kI(new_col)[dst] = kI(col)[r];
                else if (col->t == KS_VEC) kS(new_col)[dst] = kS(col)[r];
                else if (col->t == KB_VEC) kG(new_col)[dst] = kG(col)[r];
                else if (col->t == KM_LIST) kK(new_col)[dst] = r1(kK(col)[r]);
                dst++;
            }
        }
        kK(new_vals)[c] = new_col;
    }
    
    K new_dict = xD(r1(keys), new_vals);
    return xT(new_dict);
}

void table_print(K tbl, int64_t max_rows) {
    if (!tbl || tbl->t != XT_TABLE) {
        printf("[Not a table]\n");
        return;
    }
    K dict = tbl->dict.k;
    K keys = dict->dict.k;
    K vals = dict->dict.v;
    int64_t rows = table_row_count(tbl);
    if (max_rows > 0 && rows > max_rows) rows = max_rows;
    
    /* Print column headers */
    for (int64_t c = 0; c < keys->n; c++) {
        printf("%-16s ", kS(keys)[c]);
    }
    printf("\n");
    for (int64_t c = 0; c < keys->n; c++) {
        printf("---------------- ");
    }
    printf("\n");
    
    /* Print rows */
    for (int64_t r = 0; r < rows; r++) {
        for (int64_t c = 0; c < keys->n; c++) {
            K col = kK(vals)[c];
            if (col->t == KF_VEC) printf("%-16.6f ", kF(col)[r]);
            else if (col->t == KJ_VEC) printf("%-16lld ", (long long)kJ(col)[r]);
            else if (col->t == KI_VEC) printf("%-16d ", kI(col)[r]);
            else if (col->t == KS_VEC) printf("%-16s ", kS(col)[r]);
            else if (col->t == KB_VEC) printf("%-16s ", kG(col)[r] ? "1b" : "0b");
            else printf("%-16s ", "...");
        }
        printf("\n");
    }
    if (rows < table_row_count(tbl)) {
        printf("... (%lld rows total)\n", (long long)table_row_count(tbl));
    }
}
