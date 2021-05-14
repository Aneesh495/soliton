#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include "soliton/q.h"
#include "soliton/eval.h"
#include "soliton/ipc.h"

void table_print(K tbl, int64_t max_rows);

static void print_k(K x) {
    if (!x) { printf("::\n"); return; }
    if (x->t == XE_ERR) {
        printf("'%s\n", x->s);
        return;
    }
    if (x->t == XT_TABLE) {
        table_print(x, 25);
        return;
    }
    if (x->t == XD_DICT) {
        printf("[dict with %lld keys]\n", (long long)x->dict.k->n);
        return;
    }
    if (x->t < 0) {
        switch (x->t) {
            case KB_ATOM: printf("%sb\n", x->g ? "1" : "0"); break;
            case KG_ATOM: printf("0x%02x\n", x->g); break;
            case KH_ATOM: printf("%dh\n", x->h); break;
            case KI_ATOM: printf("%di\n", x->i); break;
            case KJ_ATOM: printf("%lld\n", (long long)x->j); break;
            case KF_ATOM: printf("%.6f\n", x->f); break;
            case KC_ATOM: printf("\"%c\"\n", x->c); break;
            case KS_ATOM: printf("`%s\n", x->s); break;
            default: printf("(atom %d)\n", x->t); break;
        }
        return;
    }
    if (x->t == KC_VEC) {
        printf("\"%.*s\"\n", (int)x->n, kC(x));
        return;
    }
    if (x->t == KS_VEC) {
        for (int64_t i = 0; i < x->n; i++) {
            printf("`%s", kS(x)[i]);
            if (i < x->n - 1) printf(" ");
        }
        printf("\n");
        return;
    }
    if (x->t == KJ_VEC) {
        for (int64_t i = 0; i < x->n; i++) {
            printf("%lld", (long long)kJ(x)[i]);
            if (i < x->n - 1) printf(" ");
        }
        printf("\n");
        return;
    }
    if (x->t == KF_VEC) {
        for (int64_t i = 0; i < x->n; i++) {
            printf("%.6f", kF(x)[i]);
            if (i < x->n - 1) printf(" ");
        }
        printf("\n");
        return;
    }
    if (x->t == KM_LIST) {
        printf("(");
        for (int64_t i = 0; i < x->n; i++) {
            print_k(kK(x)[i]);
            if (i < x->n - 1) printf("; ");
        }
        printf(")\n");
        return;
    }
    printf("(vector t=%d n=%lld)\n", x->t, (long long)x->n);
}

int main(int argc, char** argv) {
    k_mem_init();
    Env* global_env = env_new(NULL);
    
    int port = 0;
    const char* script_file = NULL;
    const char* eval_expr = NULL;
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            port = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-e") == 0 && i + 1 < argc) {
            eval_expr = argv[++i];
        } else if (argv[i][0] != '-') {
            script_file = argv[i];
        }
    }
    
    if (eval_expr) {
        K res = q_eval_string(global_env, eval_expr);
        print_k(res);
        r0(res);
        env_free(global_env);
        k_mem_cleanup();
        return 0;
    }
    
    if (script_file) {
        K res = q_load_file(global_env, script_file);
        if (res->t == XE_ERR) {
            fprintf(stderr, "Error loading %s: '%s\n", script_file, res->s);
            r0(res);
            env_free(global_env);
            k_mem_cleanup();
            return 1;
        }
        r0(res);
    }
    
    if (port > 0) {
        printf("Soliton IPC server listening on port %d...\n", port);
        int sfd = ipc_listen(port);
        if (sfd < 0) {
            fprintf(stderr, "Failed to bind port %d\n", port);
            return 1;
        }
        while (1) {
            int cfd = accept(sfd, NULL, NULL);
            if (cfd < 0) continue;
            uint8_t msg_type;
            K msg = ipc_recv(cfd, &msg_type);
            if (msg && msg->t != XE_ERR) {
                K resp = NULL;
                if (msg->t == KC_VEC) {
                    char* code = (char*)malloc((size_t)msg->n + 1);
                    memcpy(code, kC(msg), (size_t)msg->n);
                    code[msg->n] = '\0';
                    resp = q_eval_string(global_env, code);
                    free(code);
                } else {
                    resp = r1(msg);
                }
                if (msg_type == IPC_SYNC) {
                    size_t out_len;
                    uint8_t* out_buf = ipc_serialize(resp, IPC_RESPONSE, &out_len);
                    write(cfd, out_buf, out_len);
                    free(out_buf);
                }
                r0(resp);
                r0(msg);
            }
            close(cfd);
        }
    }
    
    if (!script_file && !eval_expr && port == 0) {
        printf("Soliton Array Engine (q REPL)\n");
        printf("Type expressions or \\\\ to exit.\n\n");
        char line[4096];
        while (1) {
            printf("q) ");
            fflush(stdout);
            if (!fgets(line, sizeof(line), stdin)) break;
            if (strncmp(line, "\\\\", 2) == 0) break;
            size_t len = strlen(line);
            while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
                line[--len] = '\0';
            }
            if (len == 0) continue;
            K res = q_eval_string(global_env, line);
            print_k(res);
            r0(res);
        }
    }
    
    env_free(global_env);
    k_mem_cleanup();
    return 0;
}
