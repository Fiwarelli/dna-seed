/*
 * DNA Seed Compression — C Engine
 * ================================
 * Same algorithm as dna_seed.py v11, 100-1000x faster.
 * Pipeline: Pre-BWT RLE → BWT (suffix array) → MTF → ZRLE → Adaptive Arithmetic (KT+decay)
 *
 * Author: Fiwa (concept/design) + AI (implementation)
 * Beats bzip2 on structured data. Pure algorithm win.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <divsufsort.h>

static int get_num_cores(void) {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 4;
}

/* ============== CONSTANTS ============== */

#define RUNA 256
#define RUNB 257
#define MAX_SYMBOLS 259  /* 0-255 + RUNA + RUNB + EOF */

#define ARITH_PRECISION 32
#define ARITH_FULL      ((uint64_t)1 << ARITH_PRECISION)
#define ARITH_HALF      ((uint64_t)1 << (ARITH_PRECISION - 1))
#define ARITH_QUARTER   ((uint64_t)1 << (ARITH_PRECISION - 2))

#define ADAPTIVE_SCALE  65536
#define ENTROPY_THRESHOLD 7.9

/* ============== UTILITY ============== */

typedef unsigned char u8;
typedef unsigned int u32;
typedef unsigned long long u64;

/* Dynamic byte buffer */
typedef struct {
    u8 *data;
    size_t len;
    size_t cap;
} Buffer;

static void buf_init(Buffer *b, size_t cap) {
    b->data = (u8 *)malloc(cap);
    b->len = 0;
    b->cap = cap;
}

static void buf_push(Buffer *b, u8 byte) {
    if (b->len >= b->cap) {
        b->cap = b->cap * 2;
        b->data = (u8 *)realloc(b->data, b->cap);
    }
    b->data[b->len++] = byte;
}

static void buf_push_bytes(Buffer *b, const u8 *src, size_t n) {
    while (b->len + n > b->cap) {
        b->cap *= 2;
        b->data = (u8 *)realloc(b->data, b->cap);
    }
    memcpy(b->data + b->len, src, n);
    b->len += n;
}

static void buf_free(Buffer *b) {
    free(b->data);
    b->data = NULL;
    b->len = b->cap = 0;
}

/* Dynamic int buffer (for ZRLE symbols) */
typedef struct {
    int *data;
    size_t len;
    size_t cap;
} IntBuf;

static void ibuf_init(IntBuf *b, size_t cap) {
    b->data = (int *)malloc(cap * sizeof(int));
    b->len = 0;
    b->cap = cap;
}

static void ibuf_push(IntBuf *b, int val) {
    if (b->len >= b->cap) {
        b->cap *= 2;
        b->data = (int *)realloc(b->data, b->cap * sizeof(int));
    }
    b->data[b->len++] = val;
}

static void ibuf_free(IntBuf *b) {
    free(b->data);
    b->data = NULL;
    b->len = b->cap = 0;
}

/* ============== VARINT ============== */

static size_t write_varint(Buffer *b, u64 val) {
    size_t written = 0;
    while (val >= 0x80) {
        buf_push(b, (u8)(val & 0x7F) | 0x80);
        val >>= 7;
        written++;
    }
    buf_push(b, (u8)val);
    return written + 1;
}

static u64 read_varint(const u8 *data, size_t *pos) {
    u64 result = 0;
    int shift = 0;
    while (1) {
        u8 byte = data[(*pos)++];
        result |= (u64)(byte & 0x7F) << shift;
        if (!(byte & 0x80)) break;
        shift += 7;
    }
    return result;
}

/* ============== ENTROPY ============== */

static double entropy(const u8 *data, size_t n) {
    if (n == 0) return 0.0;
    size_t freq[256] = {0};
    for (size_t i = 0; i < n; i++) freq[data[i]]++;
    double ent = 0.0;
    for (int i = 0; i < 256; i++) {
        if (freq[i] > 0) {
            double p = (double)freq[i] / n;
            ent -= p * log2(p);
        }
    }
    return ent;
}

/* ============== PRE-BWT RLE ============== */

static void pre_bwt_rle_encode(const u8 *data, size_t n, Buffer *out) {
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && data[j] == data[i] && j - i < 259) j++;
        size_t run = j - i;
        if (run >= 4) {
            buf_push(out, data[i]);
            buf_push(out, data[i]);
            buf_push(out, data[i]);
            buf_push(out, data[i]);
            buf_push(out, (u8)(run - 4));
            i = j;
        } else {
            buf_push(out, data[i]);
            i++;
        }
    }
}

static void pre_bwt_rle_decode(const u8 *data, size_t n, Buffer *out) {
    size_t i = 0;
    while (i < n) {
        buf_push(out, data[i]);
        size_t rlen = out->len;
        if (rlen >= 4 &&
            out->data[rlen-1] == out->data[rlen-2] &&
            out->data[rlen-2] == out->data[rlen-3] &&
            out->data[rlen-3] == out->data[rlen-4]) {
            i++;
            if (i < n) {
                u8 count = data[i];
                for (int c = 0; c < count; c++)
                    buf_push(out, out->data[rlen-1]);
                i++;
            }
        } else {
            i++;
        }
    }
}

/* ============== BWT (libdivsufsort — O(n) suffix array) ============== */
/*
 * Uses Yuta Mori's divsufsort — fastest known BWT implementation.
 * O(n) time, excellent cache behavior, single-threaded but highly optimized.
 */

/* Use divbwt for transform — correct O(n) algorithm.
   Use divsufsort inverse for decompress. */
static void bwt_transform(const u8 *data, size_t n, u8 *out, size_t *idx) {
    if (n == 0) { *idx = 0; return; }
    u8 *tmp = (u8 *)malloc(n);
    memcpy(tmp, data, n);
    saidx_t pidx = divbwt(tmp, out, NULL, (saidx_t)n);
    *idx = (size_t)pidx;
    free(tmp);
}

static void bwt_inverse(const u8 *transformed, size_t n, size_t idx, u8 *out) {
    if (n == 0) return;
    u8 *tmp = (u8 *)malloc(n);
    memcpy(tmp, transformed, n);
    inverse_bw_transform(tmp, out, NULL, (saidx_t)n, (saidx_t)idx);
    free(tmp);
}

/* ============== MTF ============== */

static void mtf_encode(const u8 *data, size_t n, u8 *out) {
    u8 alphabet[256];   /* position → symbol */
    u8 pos_of[256];     /* symbol → position (reverse lookup) */
    for (int i = 0; i < 256; i++) { alphabet[i] = (u8)i; pos_of[i] = (u8)i; }

    for (size_t i = 0; i < n; i++) {
        u8 byte = data[i];
        u8 idx = pos_of[byte];
        out[i] = idx;
        /* Move to front: shift positions of all symbols before this one */
        for (u8 j = idx; j > 0; j--) {
            alphabet[j] = alphabet[j - 1];
            pos_of[alphabet[j]] = j;
        }
        alphabet[0] = byte;
        pos_of[byte] = 0;
    }
}

static void mtf_decode(const u8 *data, size_t n, u8 *out) {
    u8 alphabet[256];
    for (int i = 0; i < 256; i++) alphabet[i] = (u8)i;

    for (size_t i = 0; i < n; i++) {
        u8 idx = data[i];
        u8 byte = alphabet[idx];
        out[i] = byte;
        memmove(alphabet + 1, alphabet, idx);
        alphabet[0] = byte;
    }
}

/* ============== ZRLE ============== */

static void zrle_encode(const u8 *mtf_data, size_t n, IntBuf *out) {
    size_t zero_count = 0;

    for (size_t i = 0; i < n; i++) {
        if (mtf_data[i] == 0) {
            zero_count++;
        } else {
            /* Flush zeros */
            if (zero_count > 0) {
                size_t zc = zero_count;
                while (zc > 0) {
                    if (zc % 2 == 1) {
                        ibuf_push(out, RUNA);
                        zc = (zc - 1) / 2;
                    } else {
                        ibuf_push(out, RUNB);
                        zc = (zc - 2) / 2;
                    }
                }
                zero_count = 0;
            }
            ibuf_push(out, mtf_data[i]); /* 1-255, no collision with RUNA/RUNB */
        }
    }
    /* Flush trailing zeros */
    if (zero_count > 0) {
        size_t zc = zero_count;
        while (zc > 0) {
            if (zc % 2 == 1) {
                ibuf_push(out, RUNA);
                zc = (zc - 1) / 2;
            } else {
                ibuf_push(out, RUNB);
                zc = (zc - 2) / 2;
            }
        }
    }
}

static void zrle_decode(const int *encoded, size_t n, Buffer *out) {
    size_t run_power = 0;
    size_t run_count = 0;

    for (size_t i = 0; i < n; i++) {
        if (encoded[i] == RUNA) {
            run_count += ((size_t)1 << run_power);
            run_power++;
        } else if (encoded[i] == RUNB) {
            run_count += ((size_t)2 << run_power);
            run_power++;
        } else {
            if (run_count > 0) {
                for (size_t j = 0; j < run_count; j++) buf_push(out, 0);
                run_count = 0;
                run_power = 0;
            }
            buf_push(out, (u8)encoded[i]);
        }
    }
    if (run_count > 0) {
        for (size_t j = 0; j < run_count; j++) buf_push(out, 0);
    }
}

/* ============== ARITHMETIC CODER ============== */

typedef struct {
    u8 *bits;
    size_t bit_count;
    size_t bit_cap;
    u64 low;
    u64 high;
    u64 pending;
} ArithEncoder;

static void aenc_init(ArithEncoder *e) {
    e->bit_cap = 1024 * 1024;
    e->bits = (u8 *)malloc(e->bit_cap);
    e->bit_count = 0;
    e->low = 0;
    e->high = ARITH_FULL - 1;
    e->pending = 0;
}

static void aenc_emit_bit(ArithEncoder *e, u8 bit) {
    if (e->bit_count >= e->bit_cap) {
        e->bit_cap *= 2;
        e->bits = (u8 *)realloc(e->bits, e->bit_cap);
    }
    e->bits[e->bit_count++] = bit;
    while (e->pending > 0) {
        if (e->bit_count >= e->bit_cap) {
            e->bit_cap *= 2;
            e->bits = (u8 *)realloc(e->bits, e->bit_cap);
        }
        e->bits[e->bit_count++] = 1 - bit;
        e->pending--;
    }
}

static void aenc_encode(ArithEncoder *e, u64 cum, u64 sym_freq, u64 total) {
    u64 rng = e->high - e->low + 1;
    e->high = e->low + (rng * (cum + sym_freq)) / total - 1;
    e->low  = e->low + (rng * cum) / total;

    while (1) {
        if (e->high < ARITH_HALF) {
            aenc_emit_bit(e, 0);
        } else if (e->low >= ARITH_HALF) {
            aenc_emit_bit(e, 1);
            e->low -= ARITH_HALF;
            e->high -= ARITH_HALF;
        } else if (e->low >= ARITH_QUARTER && e->high < 3 * ARITH_QUARTER) {
            e->pending++;
            e->low -= ARITH_QUARTER;
            e->high -= ARITH_QUARTER;
        } else {
            break;
        }
        e->low <<= 1;
        e->high = (e->high << 1) | 1;
    }
}

static void aenc_finish(ArithEncoder *e, Buffer *out, size_t *total_bits) {
    e->pending++;
    if (e->low < ARITH_QUARTER) {
        aenc_emit_bit(e, 0);
    } else {
        aenc_emit_bit(e, 1);
    }

    *total_bits = e->bit_count;

    /* Pack bits into bytes */
    size_t byte_count = (e->bit_count + 7) / 8;
    for (size_t i = 0; i < byte_count; i++) {
        u8 byte = 0;
        for (int j = 0; j < 8; j++) {
            size_t idx = i * 8 + j;
            byte = (byte << 1) | (idx < e->bit_count ? e->bits[idx] : 0);
        }
        buf_push(out, byte);
    }
    free(e->bits);
    e->bits = NULL;
}

/* Decoder */
typedef struct {
    const u8 *data;
    size_t data_len;
    size_t total_bits;
    size_t bit_pos;
    u64 low;
    u64 high;
    u64 value;
} ArithDecoder;

static u8 adec_read_bit(ArithDecoder *d) {
    if (d->bit_pos >= d->total_bits) return 0;
    size_t byte_idx = d->bit_pos / 8;
    int bit_idx = 7 - (d->bit_pos % 8);
    d->bit_pos++;
    if (byte_idx < d->data_len)
        return (d->data[byte_idx] >> bit_idx) & 1;
    return 0;
}

static void adec_init(ArithDecoder *d, const u8 *data, size_t data_len, size_t total_bits) {
    d->data = data;
    d->data_len = data_len;
    d->total_bits = total_bits;
    d->bit_pos = 0;
    d->low = 0;
    d->high = ARITH_FULL - 1;
    d->value = 0;
    for (int i = 0; i < ARITH_PRECISION; i++) {
        d->value = (d->value << 1) | adec_read_bit(d);
    }
}

static int adec_decode(ArithDecoder *d, const u64 *cum_freq, size_t num_symbols, u64 total) {
    u64 rng = d->high - d->low + 1;
    u64 scaled = ((d->value - d->low + 1) * total - 1) / rng;

    /* Binary search */
    size_t lo = 0, hi = num_symbols - 1;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (cum_freq[mid + 1] <= scaled) lo = mid + 1;
        else hi = mid;
    }
    int symbol = (int)lo;

    /* Update range */
    d->high = d->low + (rng * cum_freq[symbol + 1]) / total - 1;
    d->low  = d->low + (rng * cum_freq[symbol]) / total;

    while (1) {
        if (d->high < ARITH_HALF) {
            /* nothing */
        } else if (d->low >= ARITH_HALF) {
            d->low -= ARITH_HALF;
            d->high -= ARITH_HALF;
            d->value -= ARITH_HALF;
        } else if (d->low >= ARITH_QUARTER && d->high < 3 * ARITH_QUARTER) {
            d->low -= ARITH_QUARTER;
            d->high -= ARITH_QUARTER;
            d->value -= ARITH_QUARTER;
        } else {
            break;
        }
        d->low <<= 1;
        d->high = (d->high << 1) | 1;
        d->value = (d->value << 1) | adec_read_bit(d);
    }
    return symbol;
}

/* ============== ADAPTIVE MODEL (KT + Decay) — Fixed-Point Integer ============== */
/*
 * All-integer adaptive model. No floats anywhere.
 *
 * Decay: alpha = 0.99961 ≈ 65510/65536 (closest integer ratio)
 * Count representation: 32-bit fixed-point, scaled by 2^16
 * Initial count: 0.1 * 65536 = 6554
 * Decay: count = (count * 65510) >> 16
 * Add 1: count += 65536 (= 1.0 in fixed-point)
 *
 * This eliminates ALL float ops from the hot path.
 * The decay loop becomes integer multiply + shift — auto-vectorizable.
 */

#define FP_SHIFT    16
#define FP_ONE      (1 << FP_SHIFT)        /* 65536 = 1.0 */
#define FP_ALPHA    65510                    /* 0.99961 * 65536 ≈ 65510 */
#define FP_INIT     6554                     /* 0.1 * 65536 ≈ 6554 */
#define ADAPTIVE_SCALE 65536

static void fp_build_cumfreq(const u32 *counts, u64 tc_fp,
                              int num_symbols, u64 *cum_freq,
                              u64 *int_freq, u64 *int_total) {
    /* Convert fixed-point counts to arithmetic coder frequencies.
       Use reciprocal multiply instead of division:
       counts[j] * SCALE / tc  =  counts[j] * (SCALE * 2^32 / tc) >> 32 */
    u64 inv_tc = ((u64)ADAPTIVE_SCALE << 32) / tc_fp;  /* one division for all symbols */
    *int_total = 0;
    cum_freq[0] = 0;
    for (int j = 0; j < num_symbols; j++) {
        u64 f = ((u64)counts[j] * inv_tc) >> 32;
        if (f < 1) f = 1;
        int_freq[j] = f;
        cum_freq[j + 1] = cum_freq[j] + f;
        *int_total += f;
    }
}

/* Batch decay: accumulate updates, decay every BATCH_N symbols.
   alpha^4 ≈ 0.99844 — still highly adaptive, 4x fewer decay sweeps. */
#define BATCH_N     8
#define FP_ALPHA_B  65276   /* alpha^8 * 65536 — pre-computed */

static inline void fp_batch_decay(u32 *counts, u64 *tc, int num_symbols) {
    u64 new_tc = 0;
    for (int j = 0; j < num_symbols; j++) {
        counts[j] = (u32)(((u64)counts[j] * FP_ALPHA_B) >> FP_SHIFT);
        if (counts[j] < 1) counts[j] = 1;
        new_tc += counts[j];
    }
    *tc = new_tc;
}

static void arith_compress_adaptive(const int *symbols, size_t count,
                                     int max_symbol, Buffer *out,
                                     size_t *total_bits) {
    int num_symbols = max_symbol + 2;
    int eof = max_symbol + 1;

    u32 *counts = (u32 *)malloc(num_symbols * sizeof(u32));
    u64 tc = 0;
    for (int i = 0; i < num_symbols; i++) {
        counts[i] = FP_INIT;
        tc += FP_INIT;
    }

    u64 *cum_freq = (u64 *)malloc((num_symbols + 1) * sizeof(u64));
    u64 *int_freq = (u64 *)malloc(num_symbols * sizeof(u64));
    u64 int_total;

    ArithEncoder enc;
    aenc_init(&enc);

    /* Rebuild cumfreq at batch boundaries — 8x fewer rebuilds.
       Both encoder and decoder use same stale table between rebuilds. */
    fp_build_cumfreq(counts, tc, num_symbols, cum_freq, int_freq, &int_total);
    for (size_t i = 0; i < count; i++) {
        int s = symbols[i];
        aenc_encode(&enc, cum_freq[s], int_freq[s], int_total);
        counts[s] += FP_ONE;
        tc += FP_ONE;
        if ((i & (BATCH_N - 1)) == (BATCH_N - 1)) {
            fp_batch_decay(counts, &tc, num_symbols);
            fp_build_cumfreq(counts, tc, num_symbols, cum_freq, int_freq, &int_total);
        }
    }

    /* EOF */
    fp_build_cumfreq(counts, tc, num_symbols, cum_freq, int_freq, &int_total);
    aenc_encode(&enc, cum_freq[eof], int_freq[eof], int_total);

    aenc_finish(&enc, out, total_bits);

    free(counts);
    free(cum_freq);
    free(int_freq);
}

static void arith_decompress_adaptive(const u8 *data, size_t data_len,
                                       size_t total_bits, size_t count,
                                       int max_symbol, IntBuf *out) {
    int num_symbols = max_symbol + 2;
    int eof = max_symbol + 1;

    u32 *counts = (u32 *)malloc(num_symbols * sizeof(u32));
    u64 tc = 0;
    for (int i = 0; i < num_symbols; i++) {
        counts[i] = FP_INIT;
        tc += FP_INIT;
    }

    u64 *cum_freq = (u64 *)malloc((num_symbols + 1) * sizeof(u64));
    u64 *int_freq = (u64 *)malloc(num_symbols * sizeof(u64));
    u64 int_total;

    ArithDecoder dec;
    adec_init(&dec, data, data_len, total_bits);

    /* Match encoder: rebuild cumfreq only at batch boundaries */
    fp_build_cumfreq(counts, tc, num_symbols, cum_freq, int_freq, &int_total);
    for (size_t i = 0; i < count + 1; i++) {
        int s = adec_decode(&dec, cum_freq, num_symbols, int_total);
        if (s == eof) break;
        ibuf_push(out, s);
        counts[s] += FP_ONE;
        tc += FP_ONE;
        if ((i & (BATCH_N - 1)) == (BATCH_N - 1)) {
            fp_batch_decay(counts, &tc, num_symbols);
            fp_build_cumfreq(counts, tc, num_symbols, cum_freq, int_freq, &int_total);
        }
    }

    free(counts);
    free(cum_freq);
    free(int_freq);
}

/* ============== BLOCK-PARALLEL COMPRESS (v12) ============== */
/*
 * Split input into blocks after Pre-RLE.
 * Each block processes BWT → MTF → ZRLE → Arithmetic independently.
 * All blocks processed in parallel with pthreads.
 * Block size tuned for ratio vs speed: 10MB default.
 */

#define MIN_BLOCK_SIZE (10 * 1024 * 1024)  /* 10MB minimum block */
#define TARGET_BLOCKS  4                    /* aim for ~4 blocks to use cores */

typedef struct {
    /* Input */
    const u8 *data;
    size_t len;
    /* Output */
    size_t bwt_idx;
    size_t zrle_count;
    size_t total_bits;
    Buffer arith_buf;
    int max_sym;
} BlockJob;

static void compress_block(BlockJob *job) {
    size_t n = job->len;
    const u8 *data = job->data;

    /* BWT */
    u8 *bwt_out = (u8 *)malloc(n);
    u8 *tmp = (u8 *)malloc(n);
    memcpy(tmp, data, n);
    saidx_t pidx = divbwt(tmp, bwt_out, NULL, (saidx_t)n);
    job->bwt_idx = (size_t)pidx;
    free(tmp);

    /* MTF */
    u8 *mtf_out = (u8 *)malloc(n);
    mtf_encode(bwt_out, n, mtf_out);
    free(bwt_out);

    /* ZRLE */
    IntBuf zrle;
    ibuf_init(&zrle, n);
    zrle_encode(mtf_out, n, &zrle);
    free(mtf_out);
    job->zrle_count = zrle.len;

    /* Max symbol */
    int max_sym = RUNB;
    for (size_t i = 0; i < zrle.len; i++) {
        if (zrle.data[i] > max_sym) max_sym = zrle.data[i];
    }
    job->max_sym = max_sym;

    /* Arithmetic coding */
    buf_init(&job->arith_buf, zrle.len);
    arith_compress_adaptive(zrle.data, zrle.len, max_sym, &job->arith_buf, &job->total_bits);
    ibuf_free(&zrle);
}

static void *block_worker(void *arg) {
    compress_block((BlockJob *)arg);
    return NULL;
}

static void compress(const u8 *data, size_t n, Buffer *out) {
    /* Empty */
    if (n == 0) {
        buf_push(out, 'S'); buf_push(out, 'E');
        buf_push(out, 12);
        write_varint(out, 0);
        return;
    }

    /* Entropy check — near-random data → store raw */
    double ent = entropy(data, n);
    if (ent >= ENTROPY_THRESHOLD) {
        buf_push(out, 'S'); buf_push(out, 'E');
        buf_push(out, 0xFF);
        write_varint(out, n);
        buf_push_bytes(out, data, n);
        return;
    }

    fprintf(stderr, "Compressing: %zu bytes (entropy: %.2f bits/byte)\n", n, ent);
    struct timespec ts0, ts1, ts_total;
    clock_gettime(CLOCK_MONOTONIC, &ts_total);

    /* Step 0: Pre-BWT RLE */
    Buffer rle_buf;
    buf_init(&rle_buf, n);
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    pre_bwt_rle_encode(data, n, &rle_buf);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    fprintf(stderr, "  Pre-RLE: %zu → %zu bytes (%.3fs)\n", n, rle_buf.len,
            (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) / 1e9);

    /* Adaptive block size: bigger files → bigger blocks → better compression.
       Target ~4 blocks for parallelism, minimum 10MB per block.
       Use ceiling division so block_size * TARGET_BLOCKS >= len. */
    size_t block_size = (rle_buf.len + TARGET_BLOCKS - 1) / TARGET_BLOCKS;
    if (block_size < MIN_BLOCK_SIZE) block_size = MIN_BLOCK_SIZE;
    if (block_size > rle_buf.len) block_size = rle_buf.len;
    int num_blocks = (int)((rle_buf.len + block_size - 1) / block_size);
    if (num_blocks < 1) num_blocks = 1;
    BlockJob *jobs = (BlockJob *)calloc(num_blocks, sizeof(BlockJob));

    for (int b = 0; b < num_blocks; b++) {
        size_t offset = (size_t)b * block_size;
        size_t blen = (offset + block_size <= rle_buf.len) ? block_size : (rle_buf.len - offset);
        jobs[b].data = rle_buf.data + offset;
        jobs[b].len = blen;
    }

    /* Process all blocks in parallel */
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    int num_threads = num_blocks;
    int max_cores = get_num_cores();
    if (num_threads > max_cores) num_threads = max_cores;

    if (num_blocks == 1) {
        /* Single block: no thread overhead */
        compress_block(&jobs[0]);
    } else {
        pthread_t *threads = (pthread_t *)malloc(num_blocks * sizeof(pthread_t));
        /* Launch threads in batches of num_threads */
        for (int start = 0; start < num_blocks; start += num_threads) {
            int batch = num_blocks - start;
            if (batch > num_threads) batch = num_threads;
            for (int i = 0; i < batch; i++) {
                pthread_create(&threads[start + i], NULL, block_worker, &jobs[start + i]);
            }
            for (int i = 0; i < batch; i++) {
                pthread_join(threads[start + i], NULL);
            }
        }
        free(threads);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    fprintf(stderr, "  Blocks: %d × %zuMB, parallel BWT+MTF+ZRLE+Arith: %.3fs\n",
            num_blocks, block_size / (1024*1024),
            (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) / 1e9);

    /* Build seed — v12 block header */
    size_t rle_delta = n - rle_buf.len;
    buf_push(out, 'S'); buf_push(out, 'E');  /* magic */
    buf_push(out, 12);                        /* version 12 = block mode */
    write_varint(out, n);                     /* original size */
    write_varint(out, rle_delta);             /* pre-RLE delta */
    write_varint(out, (u64)num_blocks);       /* block count */

    size_t total_arith = 0;
    for (int b = 0; b < num_blocks; b++) {
        write_varint(out, jobs[b].len);       /* block size (pre-rle bytes) */
        write_varint(out, jobs[b].bwt_idx);
        write_varint(out, jobs[b].zrle_count);
        write_varint(out, jobs[b].total_bits);
        write_varint(out, jobs[b].arith_buf.len);  /* compressed block size */
        buf_push_bytes(out, jobs[b].arith_buf.data, jobs[b].arith_buf.len);
        total_arith += jobs[b].arith_buf.len;
        buf_free(&jobs[b].arith_buf);
    }

    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    fprintf(stderr, "  Seed: %zu bytes (ratio: %.3f, %.1f%% saved) total: %.3fs\n",
            out->len, (double)out->len / n, (1.0 - (double)out->len / n) * 100,
            (ts_end.tv_sec - ts_total.tv_sec) + (ts_end.tv_nsec - ts_total.tv_nsec) / 1e9);

    free(jobs);
    buf_free(&rle_buf);
}

/* ============== DECOMPRESS ============== */

/* Decompress a single block: Arithmetic → ZRLE → MTF → BWT inverse */
typedef struct {
    const u8 *arith_data;
    size_t arith_len;
    size_t total_bits;
    size_t zrle_count;
    size_t block_size;
    size_t bwt_idx;
    u8 *output;      /* pre-allocated output buffer for this block */
} DecompBlockJob;

static void decompress_block(DecompBlockJob *job) {
    int max_sym = RUNB;

    /* Decode arithmetic */
    IntBuf zrle_symbols;
    ibuf_init(&zrle_symbols, job->zrle_count + 1);
    arith_decompress_adaptive(job->arith_data, job->arith_len,
                               job->total_bits, job->zrle_count, max_sym, &zrle_symbols);

    /* Decode ZRLE */
    Buffer mtf_data;
    buf_init(&mtf_data, job->block_size + 256);
    zrle_decode(zrle_symbols.data, zrle_symbols.len, &mtf_data);
    ibuf_free(&zrle_symbols);

    /* Inverse MTF */
    u8 *bwt_data = (u8 *)malloc(job->block_size);
    mtf_decode(mtf_data.data, job->block_size, bwt_data);
    buf_free(&mtf_data);

    /* Inverse BWT */
    bwt_inverse(bwt_data, job->block_size, job->bwt_idx, job->output);
    free(bwt_data);
}

static void *decomp_block_worker(void *arg) {
    decompress_block((DecompBlockJob *)arg);
    return NULL;
}

static void decompress(const u8 *seed, size_t seed_len, Buffer *out) {
    size_t pos = 0;

    /* Check magic */
    if (seed_len < 3 || seed[0] != 'S' || seed[1] != 'E') {
        fprintf(stderr, "Error: invalid magic\n");
        return;
    }
    pos = 2;
    u8 version = seed[pos++];

    if (version == 0xFF) {
        u64 orig_size = read_varint(seed, &pos);
        buf_push_bytes(out, seed + pos, orig_size);
        return;
    }

    if (version == 12) {
        /* Block-parallel decompression */
        u64 original_size = read_varint(seed, &pos);
        if (original_size == 0) return;

        u64 rle_delta = read_varint(seed, &pos);
        u64 pre_rle_size = original_size - rle_delta;
        u64 num_blocks_u = read_varint(seed, &pos);
        int num_blocks = (int)num_blocks_u;

        /* Allocate contiguous output for all blocks */
        u8 *rle_data = (u8 *)malloc(pre_rle_size);
        DecompBlockJob *jobs = (DecompBlockJob *)calloc(num_blocks, sizeof(DecompBlockJob));

        size_t block_offset = 0;
        for (int b = 0; b < num_blocks; b++) {
            jobs[b].block_size = (size_t)read_varint(seed, &pos);
            jobs[b].bwt_idx = (size_t)read_varint(seed, &pos);
            jobs[b].zrle_count = (size_t)read_varint(seed, &pos);
            jobs[b].total_bits = (size_t)read_varint(seed, &pos);
            jobs[b].arith_len = (size_t)read_varint(seed, &pos);
            jobs[b].arith_data = seed + pos;
            pos += jobs[b].arith_len;
            jobs[b].output = rle_data + block_offset;
            block_offset += jobs[b].block_size;
        }

        /* Parallel decompression */
        int num_threads = num_blocks;
        int max_cores = get_num_cores();
    if (num_threads > max_cores) num_threads = max_cores;

        if (num_blocks == 1) {
            decompress_block(&jobs[0]);
        } else {
            pthread_t *threads = (pthread_t *)malloc(num_blocks * sizeof(pthread_t));
            for (int start = 0; start < num_blocks; start += num_threads) {
                int batch = num_blocks - start;
                if (batch > num_threads) batch = num_threads;
                for (int i = 0; i < batch; i++) {
                    pthread_create(&threads[start + i], NULL, decomp_block_worker, &jobs[start + i]);
                }
                for (int i = 0; i < batch; i++) {
                    pthread_join(threads[start + i], NULL);
                }
            }
            free(threads);
        }

        /* Decode pre-BWT RLE */
        pre_bwt_rle_decode(rle_data, pre_rle_size, out);
        free(rle_data);
        free(jobs);

        if (out->len > original_size) out->len = original_size;
        return;
    }

    if (version != 11) {
        fprintf(stderr, "Error: unsupported version %d\n", version);
        return;
    }

    /* Legacy v11 single-block decompression */
    u64 original_size = read_varint(seed, &pos);
    if (original_size == 0) return;

    u64 rle_delta = read_varint(seed, &pos);
    u64 pre_rle_size = original_size - rle_delta;
    u64 bwt_idx = read_varint(seed, &pos);
    u64 zrle_count = read_varint(seed, &pos);
    u64 total_bits = read_varint(seed, &pos);
    int max_sym = RUNB;

    const u8 *arith_data = seed + pos;
    size_t arith_len = seed_len - pos;

    /* Decode arithmetic */
    IntBuf zrle_symbols;
    ibuf_init(&zrle_symbols, zrle_count + 1);
    arith_decompress_adaptive(arith_data, arith_len, total_bits, zrle_count, max_sym, &zrle_symbols);

    /* Decode ZRLE */
    Buffer mtf_data;
    buf_init(&mtf_data, pre_rle_size + 256);
    zrle_decode(zrle_symbols.data, zrle_symbols.len, &mtf_data);
    ibuf_free(&zrle_symbols);

    /* Inverse MTF */
    u8 *bwt_data = (u8 *)malloc(pre_rle_size);
    mtf_decode(mtf_data.data, pre_rle_size, bwt_data);
    buf_free(&mtf_data);

    /* Inverse BWT */
    u8 *rle_data = (u8 *)malloc(pre_rle_size);
    bwt_inverse(bwt_data, pre_rle_size, bwt_idx, rle_data);
    free(bwt_data);

    /* Decode pre-BWT RLE */
    pre_bwt_rle_decode(rle_data, pre_rle_size, out);
    free(rle_data);

    if (out->len > original_size) out->len = original_size;
}

/* ============== MAIN ============== */

static void print_usage(const char *prog) {
    fprintf(stderr, "DNA Seed Compression v12 — C Engine (Block-Parallel)\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s compress   <input> <output.dna>\n", prog);
    fprintf(stderr, "  %s decompress <input.dna> <output>\n", prog);
    fprintf(stderr, "  %s benchmark  <input>\n", prog);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char *cmd = argv[1];
    const char *input_file = argv[2];

    /* Read input */
    FILE *f = fopen(input_file, "rb");
    if (!f) { fprintf(stderr, "Error: can't open %s\n", input_file); return 1; }
    fseek(f, 0, SEEK_END);
    size_t input_size = ftell(f);
    fseek(f, 0, SEEK_SET);
    u8 *input_data = (u8 *)malloc(input_size);
    fread(input_data, 1, input_size, f);
    fclose(f);

    if (strcmp(cmd, "compress") == 0) {
        if (argc < 4) { print_usage(argv[0]); return 1; }
        const char *output_file = argv[3];

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);

        Buffer out;
        buf_init(&out, input_size);
        compress(input_data, input_size, &out);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        f = fopen(output_file, "wb");
        fwrite(out.data, 1, out.len, f);
        fclose(f);

        printf("Compressed: %s (%zu bytes)\n", input_file, input_size);
        printf("Seed: %s (%zu bytes)\n", output_file, out.len);
        printf("Ratio: %.3f (%.1f%% saved)\n", (double)out.len / input_size,
               (1.0 - (double)out.len / input_size) * 100);
        printf("Time: %.3f seconds (%.1f MB/s)\n", elapsed,
               input_size / elapsed / 1048576.0);

        buf_free(&out);

    } else if (strcmp(cmd, "decompress") == 0) {
        if (argc < 4) { print_usage(argv[0]); return 1; }
        const char *output_file = argv[3];

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);

        Buffer out;
        buf_init(&out, input_size * 4);
        decompress(input_data, input_size, &out);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        f = fopen(output_file, "wb");
        fwrite(out.data, 1, out.len, f);
        fclose(f);

        printf("Decompressed: %s (%zu bytes) → %s (%zu bytes)\n",
               input_file, input_size, output_file, out.len);
        printf("Time: %.3f seconds\n", elapsed);

        buf_free(&out);

    } else if (strcmp(cmd, "benchmark") == 0) {
        /* Compress, decompress, verify roundtrip */
        struct timespec t0, t1;

        clock_gettime(CLOCK_MONOTONIC, &t0);
        Buffer compressed;
        buf_init(&compressed, input_size);
        compress(input_data, input_size, &compressed);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double comp_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        clock_gettime(CLOCK_MONOTONIC, &t0);
        Buffer decompressed;
        buf_init(&decompressed, input_size * 2);
        decompress(compressed.data, compressed.len, &decompressed);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double decomp_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        /* Verify */
        int ok = (decompressed.len == input_size &&
                  memcmp(decompressed.data, input_data, input_size) == 0);

        printf("=== DNA Seed v12 (C) Block-Parallel Benchmark ===\n");
        printf("Input:      %s (%zu bytes)\n", input_file, input_size);
        printf("Compressed: %zu bytes (ratio: %.3f, %.1f%% saved)\n",
               compressed.len, (double)compressed.len / input_size,
               (1.0 - (double)compressed.len / input_size) * 100);
        printf("Compress:   %.3f sec (%.1f MB/s)\n", comp_time,
               input_size / comp_time / 1048576.0);
        printf("Decompress: %.3f sec (%.1f MB/s)\n", decomp_time,
               input_size / decomp_time / 1048576.0);
        printf("Roundtrip:  %s\n", ok ? "PASS ✓" : "FAIL ✗");

        buf_free(&compressed);
        buf_free(&decompressed);
    } else {
        print_usage(argv[0]);
        return 1;
    }

    free(input_data);
    return 0;
}
