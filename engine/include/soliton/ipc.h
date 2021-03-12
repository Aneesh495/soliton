#ifndef SOLITON_IPC_H
#define SOLITON_IPC_H

#include "soliton/q.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * kdb+ IPC wire format specification:
 * Header:
 *   byte 0: endianness (1: little endian, 0: big endian)
 *   byte 1: message type (0: async, 1: sync, 2: response)
 *   byte 2: compression (0: uncompressed, 1: compressed)
 *   byte 3: reserved (0)
 *   bytes 4-7: message length (uint32 including 8-byte header)
 */
#define IPC_ASYNC    0
#define IPC_SYNC     1
#define IPC_RESPONSE 2

typedef struct {
    uint8_t endian;
    uint8_t msg_type;
    uint8_t compressed;
    uint8_t reserved;
    uint32_t length;
} IpcHeader;

/* Serialization and deserialization */
uint8_t* ipc_serialize(K x, uint8_t msg_type, size_t* out_len);
K ipc_deserialize(const uint8_t* buf, size_t len, uint8_t* out_msg_type);

/* IPC Server and Client */
int ipc_listen(int port);
int ipc_connect(const char* host, int port);
K ipc_send_sync(int fd, K x);
int ipc_send_async(int fd, K x);
K ipc_recv(int fd, uint8_t* out_msg_type);

#ifdef __cplusplus
}
#endif

#endif /* SOLITON_IPC_H */
