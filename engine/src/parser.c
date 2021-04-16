#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "soliton/q.h"
#include "soliton/eval.h"

void lexer_init(Lexer* lex, const char* src) {
    lex->src = src;
    lex->pos = 0;
    lex->len = strlen(src);
    lexer_next(lex);
}

static char peek_char(Lexer* lex) {
    if (lex->pos >= lex->len) return '\0';
    return lex->src[lex->pos];
}

static char get_char(Lexer* lex) {
    if (lex->pos >= lex->len) return '\0';
    return lex->src[lex->pos++];
}

static void skip_whitespace(Lexer* lex) {
    while (lex->pos < lex->len) {
        char c = lex->src[lex->pos];
        if (c == ' ' || c == '\t' || c == '\r') {
            lex->pos++;
        } else if (c == '/' && (lex->pos == 0 || lex->src[lex->pos - 1] == '\n' || lex->src[lex->pos - 1] == ' ')) {
            /* Line comment starting with / */
            while (lex->pos < lex->len && lex->src[lex->pos] != '\n') {
                lex->pos++;
            }
        } else {
            break;
        }
    }
}

Token lexer_next(Lexer* lex) {
    skip_whitespace(lex);
    Token tok;
    memset(&tok, 0, sizeof(tok));
    
    if (lex->pos >= lex->len) {
        tok.type = TOK_EOF;
        lex->curr = tok;
        return tok;
    }
    
    tok.text = &lex->src[lex->pos];
    char c = get_char(lex);
    
    if (c == '\n') {
        tok.type = TOK_SEMI;
        tok.len = 1;
        lex->curr = tok;
        return tok;
    }
    if (c == ';') {
        tok.type = TOK_SEMI;
        tok.len = 1;
        lex->curr = tok;
        return tok;
    }
    if (c == ':') {
        tok.type = TOK_COLON;
        tok.len = 1;
        lex->curr = tok;
        return tok;
    }
    if (c == ',') {
        tok.type = TOK_COMMA;
        tok.len = 1;
        lex->curr = tok;
        return tok;
    }
    if (c == '(') { tok.type = TOK_LPAREN; tok.len = 1; lex->curr = tok; return tok; }
    if (c == ')') { tok.type = TOK_RPAREN; tok.len = 1; lex->curr = tok; return tok; }
    if (c == '[') { tok.type = TOK_LBRACKET; tok.len = 1; lex->curr = tok; return tok; }
    if (c == ']') { tok.type = TOK_RBRACKET; tok.len = 1; lex->curr = tok; return tok; }
    if (c == '{') { tok.type = TOK_LBRACE; tok.len = 1; lex->curr = tok; return tok; }
    if (c == '}') { tok.type = TOK_RBRACE; tok.len = 1; lex->curr = tok; return tok; }
    
    /* Symbol literal: `sym */
    if (c == '`') {
        size_t start = lex->pos;
        while (lex->pos < lex->len && (isalnum((unsigned char)lex->src[lex->pos]) || lex->src[lex->pos] == '_' || lex->src[lex->pos] == '.')) {
            lex->pos++;
        }
        size_t slen = lex->pos - start;
        tok.type = TOK_SYMBOL;
        tok.len = lex->pos - (start - 1);
        tok.sval = sn(&lex->src[start], slen);
        lex->curr = tok;
        return tok;
    }
    
    /* String literal: "hello" */
    if (c == '"') {
        size_t start = lex->pos;
        while (lex->pos < lex->len && lex->src[lex->pos] != '"') {
            if (lex->src[lex->pos] == '\\' && lex->pos + 1 < lex->len) {
                lex->pos += 2;
            } else {
                lex->pos++;
            }
        }
        size_t slen = lex->pos - start;
        if (lex->pos < lex->len) lex->pos++; /* skip closing quote */
        tok.type = TOK_STRING;
        tok.len = slen;
        tok.text = &lex->src[start];
        lex->curr = tok;
        return tok;
    }
    
    /* Numbers: int or float */
    if (isdigit((unsigned char)c) || (c == '-' && isdigit((unsigned char)peek_char(lex)))) {
        size_t start = lex->pos - 1;
        bool is_float = false;
        while (lex->pos < lex->len && (isdigit((unsigned char)lex->src[lex->pos]) || lex->src[lex->pos] == '.' || lex->src[lex->pos] == 'e' || lex->src[lex->pos] == 'E')) {
            if (lex->src[lex->pos] == '.') is_float = true;
            lex->pos++;
        }
        tok.len = lex->pos - start;
        char num_buf[64];
        if (tok.len < sizeof(num_buf)) {
            memcpy(num_buf, &lex->src[start], tok.len);
            num_buf[tok.len] = '\0';
            if (is_float) {
                tok.type = TOK_FLOAT;
                tok.fval = atof(num_buf);
            } else {
                tok.type = TOK_INT;
                tok.ival = atoll(num_buf);
            }
        }
        lex->curr = tok;
        return tok;
    }
    
    /* Identifiers: variable or function name */
    if (isalpha((unsigned char)c) || c == '.') {
        size_t start = lex->pos - 1;
        while (lex->pos < lex->len && (isalnum((unsigned char)lex->src[lex->pos]) || lex->src[lex->pos] == '_' || lex->src[lex->pos] == '.')) {
            lex->pos++;
        }
        tok.len = lex->pos - start;
        tok.type = TOK_IDENT;
        tok.sval = sn(&lex->src[start], tok.len);
        lex->curr = tok;
        return tok;
    }
    
    /* Single-character primitive verb */
    tok.type = TOK_VERB;
    tok.len = 1;
    tok.cval = c;
    lex->curr = tok;
    return tok;
}

static K parse_expr(Lexer* lex, Env* env);

static K parse_atom_or_vec(Lexer* lex, Env* env) {
    if (lex->curr.type == TOK_INT) {
        int64_t val = lex->curr.ival;
        lexer_next(lex);
        /* Check for consecutive ints to form vector: e.g. 1 2 3 */
        if (lex->curr.type == TOK_INT || lex->curr.type == TOK_FLOAT) {
            int64_t cap = 16;
            int64_t count = 0;
            double* buf = (double*)malloc(sizeof(double) * (size_t)cap);
            buf[count++] = (double)val;
            while (lex->curr.type == TOK_INT || lex->curr.type == TOK_FLOAT) {
                if (count >= cap) {
                    cap *= 2;
                    buf = (double*)realloc(buf, sizeof(double) * (size_t)cap);
                }
                if (lex->curr.type == TOK_INT) buf[count++] = (double)lex->curr.ival;
                else buf[count++] = lex->curr.fval;
                lexer_next(lex);
            }
            K vec = ktn(KF_VEC, count);
            for (int64_t i = 0; i < count; i++) kF(vec)[i] = buf[i];
            free(buf);
            return vec;
        }
        return kj(val);
    }
    
    if (lex->curr.type == TOK_FLOAT) {
        double val = lex->curr.fval;
        lexer_next(lex);
        if (lex->curr.type == TOK_INT || lex->curr.type == TOK_FLOAT) {
            int64_t cap = 16;
            int64_t count = 0;
            double* buf = (double*)malloc(sizeof(double) * (size_t)cap);
            buf[count++] = val;
            while (lex->curr.type == TOK_INT || lex->curr.type == TOK_FLOAT) {
                if (count >= cap) {
                    cap *= 2;
                    buf = (double*)realloc(buf, sizeof(double) * (size_t)cap);
                }
                if (lex->curr.type == TOK_INT) buf[count++] = (double)lex->curr.ival;
                else buf[count++] = lex->curr.fval;
                lexer_next(lex);
            }
            K vec = ktn(KF_VEC, count);
            for (int64_t i = 0; i < count; i++) kF(vec)[i] = buf[i];
            free(buf);
            return vec;
        }
        return kf(val);
    }
    
    if (lex->curr.type == TOK_SYMBOL) {
        S s = lex->curr.sval;
        lexer_next(lex);
        if (lex->curr.type == TOK_SYMBOL) {
            int64_t cap = 16;
            int64_t count = 0;
            S* buf = (S*)malloc(sizeof(S) * (size_t)cap);
            buf[count++] = s;
            while (lex->curr.type == TOK_SYMBOL) {
                if (count >= cap) {
                    cap *= 2;
                    buf = (S*)realloc(buf, sizeof(S) * (size_t)cap);
                }
                buf[count++] = lex->curr.sval;
                lexer_next(lex);
            }
            K svec = ktn(KS_VEC, count);
            for (int64_t i = 0; i < count; i++) kS(svec)[i] = buf[i];
            free(buf);
            return svec;
        }
        return ks(s);
    }
    
    if (lex->curr.type == TOK_STRING) {
        K str = kpn(lex->curr.text, (int64_t)lex->curr.len);
        lexer_next(lex);
        return str;
    }
    
    if (lex->curr.type == TOK_IDENT) {
        S sym = lex->curr.sval;
        lexer_next(lex);
        if (lex->curr.type == TOK_COLON) {
            /* Assignment: var: expr */
            lexer_next(lex);
            K val = parse_expr(lex, env);
            env_set(env, sym, r1(val));
            return val;
        }
        /* Variable lookup */
        K val = env_get(env, sym);
        if (val) return r1(val);
        return kerror("value");
    }
    
    if (lex->curr.type == TOK_LPAREN) {
        lexer_next(lex);
        K res = parse_expr(lex, env);
        if (lex->curr.type == TOK_RPAREN) lexer_next(lex);
        return res;
    }
    
    if (lex->curr.type == TOK_LBRACE) {
        /* Lambda definition: { ... } */
        size_t start = lex->pos;
        int depth = 1;
        while (lex->pos < lex->len && depth > 0) {
            char c = lex->src[lex->pos++];
            if (c == '{') depth++;
            else if (c == '}') depth--;
        }
        size_t body_len = (lex->pos > start) ? (lex->pos - start - 1) : 0;
        char* body = (char*)malloc(body_len + 1);
        memcpy(body, &lex->src[start], body_len);
        body[body_len] = '\0';
        
        K lam = (K)calloc(1, sizeof(struct k0));
        lam->t = XF_LAMBDA;
        lam->r = 1;
        lam->lambda.body = body;
        lam->lambda.params = NULL;
        lam->lambda.env = NULL;
        lexer_next(lex);
        return lam;
    }
    
    return kerror("syntax");
}

static K parse_expr(Lexer* lex, Env* env) {
    if (lex->curr.type == TOK_EOF || lex->curr.type == TOK_SEMI) {
        return ki(0);
    }
    
    K left = parse_atom_or_vec(lex, env);
    if (left->t == XE_ERR) return left;
    
    if (lex->curr.type == TOK_VERB) {
        char verb = lex->curr.cval;
        lexer_next(lex);
        K right = parse_expr(lex, env);
        K res;
        switch (verb) {
            case '+': res = v_add(left, right); break;
            case '-': res = v_sub(left, right); break;
            case '*': res = v_mul(left, right); break;
            case '%': res = v_div(left, right); break;
            case '&': res = v_and(left, right); break;
            case '|': res = v_or(left, right); break;
            case '=': res = v_eq(left, right); break;
            case '<': res = v_lt(left, right); break;
            case '>': res = v_gt(left, right); break;
            case ',': res = v_join(left, right); break;
            case '#': res = v_take(left, right); break;
            case '_': res = v_drop(left, right); break;
            case '?': res = v_find(left, right); break;
            case '@': res = v_at(left, right); break;
            case '.': res = v_dot(left, right); break;
            case '!': res = xD(left, right); break;
            default: res = kerror("verb"); break;
        }
        r0(left);
        r0(right);
        return res;
    }
    
    return left;
}

K q_eval_string(Env* env, const char* code) {
    Lexer lex;
    lexer_init(&lex, code);
    K last = ki(0);
    while (lex.curr.type != TOK_EOF) {
        if (lex.curr.type == TOK_SEMI) {
            lexer_next(&lex);
            continue;
        }
        r0(last);
        last = parse_expr(&lex, env);
        if (last->t == XE_ERR) break;
    }
    return last;
}
