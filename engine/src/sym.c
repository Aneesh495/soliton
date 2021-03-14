#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "soliton/q.h"

#define SYM_BUCKETS 65536

typedef struct SymNode {
    char* str;
    uint32_t hash;
    struct SymNode* next;
} SymNode;

static SymNode* g_buckets[SYM_BUCKETS];
static bool g_sym_initialized = false;

static uint32_t hash_str(const char* s, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; i++) {
        h ^= (uint8_t)s[i];
        h *= 16777619u;
    }
    return h;
}

void sym_init(void) {
    if (g_sym_initialized) return;
    memset(g_buckets, 0, sizeof(g_buckets));
    g_sym_initialized = true;
    /* Pre-intern empty symbol ` */
    ss("");
}

void sym_cleanup(void) {
    if (!g_sym_initialized) return;
    for (size_t i = 0; i < SYM_BUCKETS; i++) {
        SymNode* node = g_buckets[i];
        while (node) {
            SymNode* next = node->next;
            free(node->str);
            free(node);
            node = next;
        }
        g_buckets[i] = NULL;
    }
    g_sym_initialized = false;
}

S sn(const char* s, size_t n) {
    if (!g_sym_initialized) sym_init();
    if (!s) return ss("");
    
    uint32_t h = hash_str(s, n);
    size_t idx = h % SYM_BUCKETS;
    
    SymNode* node = g_buckets[idx];
    while (node) {
        if (node->hash == h && strlen(node->str) == n && memcmp(node->str, s, n) == 0) {
            return node->str;
        }
        node = node->next;
    }
    
    SymNode* new_node = (SymNode*)malloc(sizeof(SymNode));
    new_node->str = (char*)malloc(n + 1);
    memcpy(new_node->str, s, n);
    new_node->str[n] = '\0';
    new_node->hash = h;
    new_node->next = g_buckets[idx];
    g_buckets[idx] = new_node;
    
    return new_node->str;
}

S ss(const char* s) {
    if (!s) return sn("", 0);
    return sn(s, strlen(s));
}
