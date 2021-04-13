#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "soliton/q.h"

/* Days between Unix epoch (1970.01.01) and kdb epoch (2000.01.01) */
#define KDB_EPOCH_DAYS 10957
#define NANOS_PER_SEC 1000000000LL
#define NANOS_PER_DAY (86400LL * NANOS_PER_SEC)

int64_t q_now_timestamp(void) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    int64_t unix_nanos = (int64_t)ts.tv_sec * NANOS_PER_SEC + ts.tv_nsec;
    int64_t kdb_nanos = unix_nanos - ((int64_t)KDB_EPOCH_DAYS * NANOS_PER_DAY);
    return kdb_nanos;
}

int32_t q_now_date(void) {
    time_t t = time(NULL);
    int64_t days_unix = t / 86400;
    return (int32_t)(days_unix - KDB_EPOCH_DAYS);
}

int32_t q_now_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    time_t sec = ts.tv_sec % 86400;
    return (int32_t)(sec * 1000 + ts.tv_nsec / 1000000);
}

void q_format_timestamp(int64_t nanos, char* buf, size_t buf_sz) {
    int64_t unix_nanos = nanos + ((int64_t)KDB_EPOCH_DAYS * NANOS_PER_DAY);
    time_t sec = (time_t)(unix_nanos / NANOS_PER_SEC);
    int64_t rem_nanos = unix_nanos % NANOS_PER_SEC;
    if (rem_nanos < 0) {
        rem_nanos += NANOS_PER_SEC;
        sec -= 1;
    }
    struct tm tm_buf;
    gmtime_r(&sec, &tm_buf);
    snprintf(buf, buf_sz, "%04d.%02d.%02dD%02d:%02d:%02d.%09lld",
             tm_buf.tm_year + 1900, tm_buf.tm_mon + 1, tm_buf.tm_mday,
             tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec,
             (long long)rem_nanos);
}

int64_t q_parse_timestamp(const char* s) {
    int year, month, day, hour, min, sec;
    long long nanos = 0;
    if (sscanf(s, "%04d.%02d.%02dD%02d:%02d:%02d.%09lld",
               &year, &month, &day, &hour, &min, &sec, &nanos) >= 6) {
        struct tm tm_val;
        memset(&tm_val, 0, sizeof(tm_val));
        tm_val.tm_year = year - 1900;
        tm_val.tm_mon = month - 1;
        tm_val.tm_mday = day;
        tm_val.tm_hour = hour;
        tm_val.tm_min = min;
        tm_val.tm_sec = sec;
        time_t t = timegm(&tm_val);
        int64_t unix_nanos = (int64_t)t * NANOS_PER_SEC + nanos;
        return unix_nanos - ((int64_t)KDB_EPOCH_DAYS * NANOS_PER_DAY);
    }
    return 0;
}
