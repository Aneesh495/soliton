#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdio.h>
#include "soliton/q.h"
#include "soliton/ipc.h"

typedef struct {
    uint8_t* data;
    size_t size;
    size_t cap;
} ByteBuffer;

static void buf_init(ByteBuffer* b, size_t initial_cap) {
    b->cap = initial_cap ? initial_cap : 256;
    b->data = (uint8_t*)malloc(b->cap);
    b->size = 0;
}

static void buf_write(ByteBuffer* b, const void* src, size_t n) {
    while (b->size + n > b->cap) {
        b->cap *= 2;
        b->data = (uint8_t*)realloc(b->data, b->cap);
    }
    memcpy(b->data + b->size, src, n);
    b->size += n;
}

static void buf_write_byte(ByteBuffer* b, uint8_t byte) {
    buf_write(b, &byte, 1);
}

static void buf_write_int(ByteBuffer* b, int32_t val) {
    buf_write(b, &val, 4);
}

static void buf_write_long(ByteBuffer* b, int64_t val) {
    buf_write(b, &val, 8);
}

static void buf_write_double(ByteBuffer* b, double val) {
    buf_write(b, &val, 8);
}

static void serialize_k(ByteBuffer* b, K x) {
    if (!x) {
        buf_write_byte(b, (uint8_t)KB_ATOM);
        buf_write_byte(b, 0);
        return;
    }
    
    buf_write_byte(b, (uint8_t)x->t);
    
    if (x->t < 0) {
        /* Atom */
        switch (x->t) {
            case KB_ATOM:
            case KG_ATOM: buf_write_byte(b, x->g); break;
            case KH_ATOM: buf_write(b, &x->h, 2); break;
            case KI_ATOM: buf_write_int(b, x->i); break;
            case KJ_ATOM:
            case KP_ATOM:
            case KN_ATOM: buf_write_long(b, x->j); break;
            case KE_ATOM: buf_write(b, &x->e, 4); break;
            case KF_ATOM: buf_write_double(b, x->f); break;
            case KC_ATOM: buf_write_byte(b, (uint8_t)x->c); break;
            case KS_ATOM: {
                size_t slen = strlen(x->s);
                buf_write(b, x->s, slen + 1);
                break;
            }
            default: buf_write_long(b, 0); break;
        }
        return;
    }
    
    if (x->t == XD_DICT) {
        serialize_k(b, x->dict.k);
        serialize_k(b, x->dict.v);
        return;
    }
    
    if (x->t == XT_TABLE) {
        buf_write_byte(b, 0); /* flip attribute */
        serialize_k(b, x->dict.k);
        return;
    }
    
    /* Vectors */
    buf_write_byte(b, (uint8_t)x->u); /* attribute */
    buf_write_int(b, (int32_t)x->n);
    
    switch (x->t) {
        case KM_LIST: {
            for (int64_t i = 0; i < x->n; i++) serialize_k(b, kK(x)[i]);
            break;
        }
        case KB_VEC:
        case KG_VEC: buf_write(b, kG(x), (size_t)x->n); break;
        case KH_VEC: buf_write(b, kH(x), (size_t)x->n * 2); break;
        case KI_VEC: buf_write(b, kI(x), (size_t)x->n * 4); break;
        case KJ_VEC:
        case KP_VEC:
        case KN_VEC: buf_write(b, kJ(x), (size_t)x->n * 8); break;
        case KE_VEC: buf_write(b, kE(x), (size_t)x->n * 4); break;
        case KF_VEC: buf_write(b, kF(x), (size_t)x->n * 8); break;
        case KC_VEC: buf_write(b, kC(x), (size_t)x->n); break;
        case KS_VEC: {
            for (int64_t i = 0; i < x->n; i++) {
                size_t slen = strlen(kS(x)[i]);
                buf_write(b, kS(x)[i], slen + 1);
            }
            break;
        }
        default: break;
    }
}

uint8_t* ipc_serialize(K x, uint8_t msg_type, size_t* out_len) {
    ByteBuffer b;
    buf_init(&b, 512);
    
    /* 8-byte kdb+ header placeholder */
    uint8_t hdr[8] = { 1, msg_type, 0, 0, 0, 0, 0, 0 };
    buf_write(&b, hdr, 8);
    
    serialize_k(&b, x);
    
    uint32_t total_len = (uint32_t)b.size;
    memcpy(b.data + 4, &total_len, 4);
    
    *out_len = b.size;
    return b.data;
}

typedef struct {
    const uint8_t* data;
    size_t pos;
    size_t len;
} ByteReader;

static K deserialize_k(ByteReader* r) {
    if (r->pos >= r->len) return kerror("eof");
    int8_t type = (int8_t)r->data[r->pos++];
    
    if (type < 0) {
        switch (type) {
            case KB_ATOM: return kb(r->data[r->pos++] != 0);
            case KG_ATOM: return kg(r->data[r->pos++]);
            case KH_ATOM: {
                int16_t v; memcpy(&v, &r->data[r->pos], 2); r->pos += 2;
                return kh(v);
            }
            case KI_ATOM: {
                int32_t v; memcpy(&v, &r->data[r->pos], 4); r->pos += 4;
                return ki(v);
            }
            case KJ_ATOM: {
                int64_t v; memcpy(&v, &r->data[r->pos], 8); r->pos += 8;
                return kj(v);
            }
            case KF_ATOM: {
                double v; memcpy(&v, &r->data[r->pos], 8); r->pos += 8;
                return kf(v);
            }
            case KC_ATOM: return kc((char)r->data[r->pos++]);
            case KS_ATOM: {
                const char* s = (const char*)&r->data[r->pos];
                size_t slen = strlen(s);
                r->pos += slen + 1;
                return ks(s);
            }
            default: return kerror("type");
        }
    }
    
    if (type == XD_DICT) {
        K keys = deserialize_k(r);
        K vals = deserialize_k(r);
        return xD(keys, vals);
    }
    
    if (type == XT_TABLE) {
        r->pos++; /* skip flip attr */
        K dict = deserialize_k(r);
        return xT(dict);
    }
    
    r->pos++; /* skip attribute */
    int32_t n;
    memcpy(&n, &r->data[r->pos], 4);
    r->pos += 4;
    
    if (type == KM_LIST) {
        K res = k_mixed_list(n);
        for (int32_t i = 0; i < n; i++) {
            kK(res)[i] = deserialize_k(r);
        }
        return res;
    }
    
    K res = ktn(type, n);
    if (type == KF_VEC) {
        memcpy(kF(res), &r->data[r->pos], (size_t)n * 8);
        r->pos += (size_t)n * 8;
    } else if (type == KJ_VEC) {
        memcpy(kJ(res), &r->data[r->pos], (size_t)n * 8);
        r->pos += (size_t)n * 8;
    } else if (type == KI_VEC) {
        memcpy(kI(res), &r->data[r->pos], (size_t)n * 4);
        r->pos += (size_t)n * 4;
    } else if (type == KC_VEC) {
        memcpy(kC(res), &r->data[r->pos], (size_t)n);
        r->pos += (size_t)n;
    } else if (type == KS_VEC) {
        for (int32_t i = 0; i < n; i++) {
            const char* s = (const char*)&r->data[r->pos];
            size_t slen = strlen(s);
            kS(res)[i] = ss(s);
            r->pos += slen + 1;
        }
    }
    return res;
}

K ipc_deserialize(const uint8_t* buf, size_t len, uint8_t* out_msg_type) {
    if (len < 8) return kerror("short");
    if (out_msg_type) *out_msg_type = buf[1];
    ByteReader r = { buf, 8, len };
    return deserialize_k(&r);
}

int ipc_listen(int port) {
    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) return -1;
    int opt = 1;
    setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    
    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons((uint16_t)port);
    
    if (bind(fd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        close(fd);
        return -1;
    }
    if (listen(fd, 16) < 0) {
        close(fd);
        return -1;
    }
    return fd;
}

int ipc_connect(const char* host, int port) {
    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) return -1;
    
    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons((uint16_t)port);
    inet_pton(AF_INET, host, &addr.sin_addr);
    
    if (connect(fd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        close(fd);
        return -1;
    }
    return fd;
}

int ipc_send_async(int fd, K x) {
    size_t len = 0;
    uint8_t* buf = ipc_serialize(x, IPC_ASYNC, &len);
    ssize_t w = write(fd, buf, len);
    free(buf);
    return (w == (ssize_t)len) ? 0 : -1;
}

K ipc_send_sync(int fd, K x) {
    size_t len = 0;
    uint8_t* buf = ipc_serialize(x, IPC_SYNC, &len);
    ssize_t w = write(fd, buf, len);
    free(buf);
    if (w != (ssize_t)len) return kerror("send");
    uint8_t resp_type;
    return ipc_recv(fd, &resp_type);
}

K ipc_recv(int fd, uint8_t* out_msg_type) {
    uint8_t hdr[8];
    ssize_t r = read(fd, hdr, 8);
    if (r < 8) return kerror("recv");
    
    uint32_t total_len;
    memcpy(&total_len, hdr + 4, 4);
    if (total_len < 8 || total_len > 64 * 1024 * 1024) return kerror("length");
    
    uint8_t* buf = (uint8_t*)malloc(total_len);
    memcpy(buf, hdr, 8);
    size_t read_bytes = 8;
    while (read_bytes < total_len) {
        ssize_t chunk = read(fd, buf + read_bytes, total_len - read_bytes);
        if (chunk <= 0) { free(buf); return kerror("recv"); }
        read_bytes += (size_t)chunk;
    }
    
    K res = ipc_deserialize(buf, total_len, out_msg_type);
    free(buf);
    return res;
}
