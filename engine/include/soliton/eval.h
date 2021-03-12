#ifndef SOLITON_EVAL_H
#define SOLITON_EVAL_H

#include "soliton/q.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Execution environment and namespace mappings.
 * Variables and global tables are stored in a symbol-keyed dictionary.
 */
typedef struct Env {
    K dict;             /* symbol -> K value mapping */
    struct Env* parent; /* enclosing lexical frame */
} Env;

Env* env_new(Env* parent);
void env_free(Env* env);
K env_get(Env* env, S sym);
void env_set(Env* env, S sym, K val);

/* Lexer token types */
typedef enum {
    TOK_EOF = 0,
    TOK_INT,
    TOK_FLOAT,
    TOK_CHAR,
    TOK_STRING,
    TOK_SYMBOL,
    TOK_IDENT,
    TOK_VERB,
    TOK_ADVERB,
    TOK_LPAREN,
    TOK_RPAREN,
    TOK_LBRACKET,
    TOK_RBRACKET,
    TOK_LBRACE,
    TOK_RBRACE,
    TOK_SEMI,
    TOK_COLON,
    TOK_COMMA
} TokenType;

typedef struct {
    TokenType type;
    const char* text;
    size_t len;
    int64_t ival;
    double fval;
    char cval;
    S sval;
} Token;

/* Tokenizer interface */
typedef struct {
    const char* src;
    size_t pos;
    size_t len;
    Token curr;
} Lexer;

void lexer_init(Lexer* lex, const char* src);
Token lexer_next(Lexer* lex);

/* Parser and Evaluator */
K q_eval_string(Env* env, const char* code);
K q_eval(Env* env, K ast);
K q_load_file(Env* env, const char* filepath);

/* Builtin Verbs and Adverbs */
K v_add(K x, K y);
K v_sub(K x, K y);
K v_mul(K x, K y);
K v_div(K x, K y);
K v_mod(K x, K y);
K v_min(K x, K y);
K v_max(K x, K y);
K v_eq(K x, K y);
K v_neq(K x, K y);
K v_lt(K x, K y);
K v_gt(K x, K y);
K v_lte(K x, K y);
K v_gte(K x, K y);
K v_match(K x, K y);
K v_and(K x, K y);
K v_or(K x, K y);
K v_not(K x);
K v_neg(K x);
K v_abs(K x);
K v_sqrt(K x);
K v_exp(K x);
K v_log(K x);
K v_count(K x);
K v_first(K x);
K v_last(K x);
K v_til(K x);
K v_reverse(K x);
K v_take(K x, K y);
K v_drop(K x, K y);
K v_join(K x, K y);
K v_find(K x, K y);
K v_at(K x, K y);
K v_dot(K x, K y);
K v_sum(K x);
K v_avg(K x);
K v_dev(K x);
K v_var(K x);
K v_med(K x);
K v_enlist(K x);
K v_flip(K x);
K v_keys(K x);
K v_values(K x);

/* Adverbs */
K adv_each(Env* env, K fn, K list);
K adv_over(Env* env, K fn, K list);
K adv_scan(Env* env, K fn, K list);
K adv_prior(Env* env, K fn, K list);
K adv_each_right(Env* env, K fn, K x, K y);
K adv_each_left(Env* env, K fn, K x, K y);

/* qSQL primitives */
K qsql_select(Env* env, K tbl, K cols, K by_cols, K where_clause);
K qsql_update(Env* env, K tbl, K cols, K by_cols, K where_clause);
K qsql_delete(Env* env, K tbl, K cols, K where_clause);
K qsql_exec(Env* env, K tbl, K col, K by_cols, K where_clause);

#ifdef __cplusplus
}
#endif

#endif /* SOLITON_EVAL_H */
