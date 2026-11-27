/*
 * pmscan.c -- which principal minors of the Fourier matrix F_M vanish in the
 * finite field K = F_p[x]/(h), where h is an irreducible factor of the M-th
 * cyclotomic polynomial modulo p and zeta = x mod h?
 *
 * Every subset S of {0,...,M-1} is visited exactly once by a Schur-complement
 * recursion: index t is either excluded (delete its row and column) or
 * included.  An included index with nonzero pivot is eliminated; an included
 * index with zero pivot is carried along ("deferred") and the determinant of
 * the deferred block is taken at the leaf.  At a leaf with eliminated set E and
 * deferred set D the current matrix is the Schur complement of A[E+D] with
 * respect to E, and det A[E] != 0, so det A[E+D] = 0 iff the deferred block is
 * singular.
 *
 * Build:  cc -O2 -o pmscan pmscan.c            (log/Zech tables, small fields)
 *         cc -O2 -DPOLY -o pmscan_poly pmscan.c (polynomial arithmetic, any field)
 * Run:    pmscan M p f h_0 ... h_{f-1} [list] [first] [lowrank R]
 *         with h(x) = x^f + h_{f-1} x^{f-1} + ... + h_0.
 *         `lowrank R` replaces the full scan by a direct test (one Gaussian
 *         elimination each) of every subset of size at most R that contains 0;
 *         vanishing ones are printed as `Z mask`.  `check MASK` (hexadecimal)
 *         tests the single subset MASK by direct elimination.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXM 60
#define MAXF 32

static int M, P, F;
static int Hc[MAXF];
static uint32_t RED[MAXF][MAXF]; /* RED[k][j]: coefficient of x^j in x^(F+k) mod h */

typedef struct { uint32_t c[MAXF]; } pe;

static void pe_zero(pe *a) { memset(a, 0, sizeof(*a)); }
static void pe_one(pe *a) { pe_zero(a); a->c[0] = 1 % P; }
static int pe_iszero(const pe *a) { for (int i = 0; i < F; i++) if (a->c[i]) return 0; return 1; }
static int pe_isone(const pe *a) { if (a->c[0] != 1) return 0; for (int i = 1; i < F; i++) if (a->c[i]) return 0; return 1; }

static void pe_mul(const pe *a, const pe *b, pe *r) {
    uint64_t t[2 * MAXF];
    int n = 2 * F - 1;
    for (int i = 0; i < n; i++) t[i] = 0;
    for (int i = 0; i < F; i++) {
        uint64_t ai = a->c[i];
        if (!ai) continue;
        for (int j = 0; j < F; j++) t[i + j] += ai * b->c[j];
    }
    for (int i = 0; i < n; i++) t[i] %= (uint64_t)P;
    for (int k = 0; k < F - 1; k++) {
        uint64_t tk = t[F + k];
        if (!tk) continue;
        for (int j = 0; j < F; j++) t[j] += tk * RED[k][j];
    }
    for (int j = 0; j < F; j++) r->c[j] = (uint32_t)(t[j] % (uint64_t)P);
    for (int j = F; j < MAXF; j++) r->c[j] = 0;
}

static void pe_pow(const pe *a, uint64_t e, pe *r) {
    pe base = *a, acc, t;
    pe_one(&acc);
    while (e) {
        if (e & 1) { pe_mul(&acc, &base, &t); acc = t; }
        pe_mul(&base, &base, &t); base = t;
        e >>= 1;
    }
    *r = acc;
}

/* a^(p^f - 2) = a^(p-2) * prod_{k=1}^{f-1} (a^(p-1))^(p^k) */
static void pe_inv(const pe *a, pe *r) {
    pe acc, t, u;
    pe_pow(a, (uint64_t)(P - 2), &acc);
    pe_mul(&acc, a, &t);
    for (int k = 1; k < F; k++) {
        pe_pow(&t, (uint64_t)P, &u); t = u;
        pe_mul(&acc, &t, &u); acc = u;
    }
    *r = acc;
}

static void setup_reduction(void) {
    for (int j = 0; j < F; j++) RED[0][j] = (uint32_t)((P - Hc[j] % P) % P);
    for (int k = 0; k + 1 < F - 1 + 1 && k + 1 < MAXF; k++) {
        uint64_t top = RED[k][F - 1];
        for (int j = F - 1; j >= 1; j--) RED[k + 1][j] = (uint32_t)((RED[k][j - 1] + top * RED[0][j]) % (uint64_t)P);
        RED[k + 1][0] = (uint32_t)((top * RED[0][0]) % (uint64_t)P);
    }
}

#ifdef POLY
typedef pe elt;
#define EISZ(a) pe_iszero(&(a))
static inline void e_mul(elt *r, const elt *a, const elt *b) { pe t; pe_mul(a, b, &t); *r = t; }
static inline void e_inv(elt *r, const elt *a) { pe t; pe_inv(a, &t); *r = t; }
static inline void e_mulsub(elt *r, const elt *a, const elt *fac, const elt *b) {
    pe t; pe_mul(fac, b, &t);
    for (int i = 0; i < F; i++) { uint32_t x = a->c[i] + (uint32_t)P - t.c[i]; if (x >= (uint32_t)P) x -= (uint32_t)P; r->c[i] = x; }
    for (int i = F; i < MAXF; i++) r->c[i] = 0;
}
#else
typedef int32_t elt;
static int32_t Q1, ZS, NEGL;
static int32_t *ZECH;
#define EISZ(a) ((a) == ZS)
static inline int32_t t_mul(int32_t a, int32_t b) { if (a == ZS || b == ZS) return ZS; int32_t s = a + b; if (s >= Q1) s -= Q1; return s; }
static inline int32_t t_add(int32_t a, int32_t b) {
    if (a == ZS) return b;
    if (b == ZS) return a;
    int32_t d = b - a; if (d < 0) d += Q1;
    int32_t z = ZECH[d]; if (z == ZS) return ZS;
    int32_t s = a + z; if (s >= Q1) s -= Q1; return s;
}
static inline void e_mul(elt *r, const elt *a, const elt *b) { *r = t_mul(*a, *b); }
static inline void e_inv(elt *r, const elt *a) { *r = (*a == 0) ? 0 : Q1 - *a; }
static inline void e_mulsub(elt *r, const elt *a, const elt *fac, const elt *b) {
    int32_t t = t_mul(*fac, *b);
    if (t != ZS) { t += NEGL; if (t >= Q1) t -= Q1; }
    *r = t_add(*a, t);
}
#endif

static elt BUF[MAXM + 2][MAXM * MAXM];
static uint64_t zero_by_rank[MAXM + 1], total_by_rank[MAXM + 1], nzero, nleaf;
static int list_mode, stop_first;

static int block_is_singular(const elt *C, int n) {
    elt W[16 * 16];
    if (n > 16) { fprintf(stderr, "deferred block too large (%d)\n", n); exit(3); }
    memcpy(W, C, sizeof(elt) * (size_t)n * (size_t)n);
    for (int col = 0; col < n; col++) {
        int r = col;
        while (r < n && EISZ(W[r * n + col])) r++;
        if (r == n) return 1;
        if (r != col) for (int j = 0; j < n; j++) { elt t = W[r * n + j]; W[r * n + j] = W[col * n + j]; W[col * n + j] = t; }
        elt inv; e_inv(&inv, &W[col * n + col]);
        for (int i = col + 1; i < n; i++) {
            if (EISZ(W[i * n + col])) continue;
            elt fac; e_mul(&fac, &W[i * n + col], &inv);
            for (int j = col + 1; j < n; j++) e_mulsub(&W[i * n + j], &W[i * n + j], &fac, &W[col * n + j]);
        }
    }
    return 0;
}

static void rec(int depth, int n, int nund, int ndef, uint64_t mask) {
    elt *C = BUF[depth];
    if (nund == 0) {
        int rk = __builtin_popcountll(mask);
        total_by_rank[rk]++; nleaf++;
        int z = (ndef == 0) ? 0 : (ndef == 1 ? EISZ(C[0]) : block_is_singular(C, n));
        if (z) {
            zero_by_rank[rk]++; nzero++;
            if (list_mode) printf("Z %llx\n", (unsigned long long)mask);
            if (stop_first) { printf("FIRSTZERO %llx rank %d\n", (unsigned long long)mask, rk); exit(0); }
        }
        return;
    }
    elt *D = BUF[depth + 1];
    int m = n - 1;
    /* exclude index `depth` */
    for (int i = 1; i < n; i++) memcpy(&D[(i - 1) * m], &C[i * n + 1], sizeof(elt) * (size_t)m);
    rec(depth + 1, m, nund - 1, ndef, mask);
    /* include index `depth` */
    if (!EISZ(C[0])) {
        elt inv; e_inv(&inv, &C[0]);
        for (int i = 1; i < n; i++) {
            elt fac; e_mul(&fac, &C[i * n], &inv);
            for (int j = 1; j < n; j++) e_mulsub(&D[(i - 1) * m + (j - 1)], &C[i * n + j], &fac, &C[j]);
        }
        rec(depth + 1, m, nund - 1, ndef, mask | (1ULL << depth));
    } else {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) D[i * n + j] = C[((i + 1) % n) * n + ((j + 1) % n)];
        rec(depth + 1, n, nund - 1, ndef + 1, mask | (1ULL << depth));
    }
}

/* Direct test of one subset given as a bit mask. */
static int subset_is_singular(uint64_t mask) {
    int idx[MAXM], n = 0;
    for (int i = 0; i < M; i++) if (mask >> i & 1) idx[n++] = i;
    elt *W = malloc(sizeof(elt) * (size_t)(n ? n * n : 1));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) W[i * n + j] = BUF[0][idx[i] * M + idx[j]];
    int singular = 0;
    for (int col = 0; col < n && !singular; col++) {
        int r = col;
        while (r < n && EISZ(W[r * n + col])) r++;
        if (r == n) { singular = 1; break; }
        if (r != col) for (int j = 0; j < n; j++) { elt t = W[r * n + j]; W[r * n + j] = W[col * n + j]; W[col * n + j] = t; }
        elt inv; e_inv(&inv, &W[col * n + col]);
        for (int i = col + 1; i < n; i++) {
            if (EISZ(W[i * n + col])) continue;
            elt fac; e_mul(&fac, &W[i * n + col], &inv);
            for (int j = col + 1; j < n; j++) e_mulsub(&W[i * n + j], &W[i * n + j], &fac, &W[col * n + j]);
        }
    }
    free(W);
    return singular;
}

/* Direct test of all subsets of size <= R containing index 0. */
static void lowrank_scan(int R) {
    int idx[16];
    elt W[16 * 16];
    if (R > 16) R = 16;
    for (int r = 1; r <= R && r <= M; r++) {
        idx[0] = 0;
        for (int i = 1; i < r; i++) idx[i] = i;
        for (;;) {
            for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) W[i * r + j] = BUF[0][idx[i] * M + idx[j]];
            total_by_rank[r]++; nleaf++;
            if (block_is_singular(W, r)) {
                uint64_t mask = 0;
                for (int i = 0; i < r; i++) mask |= 1ULL << idx[i];
                zero_by_rank[r]++; nzero++;
                if (list_mode) printf("Z %llx\n", (unsigned long long)mask);
            }
            int k = r - 1;
            while (k >= 1 && idx[k] == M - r + k) k--;
            if (k < 1) break;
            idx[k]++;
            for (int i = k + 1; i < r; i++) idx[i] = idx[i - 1] + 1;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s M p f h_0 ... h_{f-1} [list] [first]\n", argv[0]); return 2; }
    M = atoi(argv[1]); P = atoi(argv[2]); F = atoi(argv[3]);
    if (M < 1 || M > MAXM || F < 1 || F > MAXF || P < 2 || P > 65521 || argc < 4 + F) { fprintf(stderr, "bad arguments\n"); return 2; }
    for (int i = 0; i < F; i++) Hc[i] = ((atoi(argv[4 + i]) % P) + P) % P;
    int lowrank = 0, check_mode = 0;
    uint64_t check_mask = 0;
    for (int i = 4 + F; i < argc; i++) {
        if (!strcmp(argv[i], "list")) list_mode = 1;
        if (!strcmp(argv[i], "first")) stop_first = 1;
        if (!strcmp(argv[i], "lowrank") && i + 1 < argc) lowrank = atoi(argv[++i]);
        if (!strcmp(argv[i], "check") && i + 1 < argc) { check_mode = 1; check_mask = strtoull(argv[++i], NULL, 16); }
    }
    setup_reduction();

    /* zeta = x mod h and its powers */
    pe zeta; pe_zero(&zeta);
    if (F == 1) zeta.c[0] = (uint32_t)((P - Hc[0]) % P); else zeta.c[1] = 1;
    static pe zp[MAXM + 1];
    pe_one(&zp[0]);
    for (int k = 1; k <= M; k++) pe_mul(&zp[k - 1], &zeta, &zp[k]);
    if (!pe_isone(&zp[M])) { fprintf(stderr, "zeta^M != 1: h is not a factor of Phi_M mod p\n"); return 4; }
    for (int d = 1; d < M; d++) if (M % d == 0 && pe_isone(&zp[d])) { fprintf(stderr, "zeta has order %d < M\n", d); return 4; }

#ifdef POLY
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) BUF[0][i * M + j] = zp[(i * j) % M];
#else
    uint64_t q = 1;
    for (int i = 0; i < F; i++) { q *= (uint64_t)P; if (q > (1ULL << 24)) { fprintf(stderr, "field too large for tables\n"); return 5; } }
    Q1 = (int32_t)(q - 1); ZS = Q1; NEGL = (P == 2) ? 0 : Q1 / 2;
    uint32_t *EXPI = malloc(sizeof(uint32_t) * q);
    int32_t *LOGI = malloc(sizeof(int32_t) * q);
    ZECH = malloc(sizeof(int32_t) * q);
    /* find a generator of K^* */
    uint64_t fac[64]; int nf = 0; { uint64_t r = q - 1; for (uint64_t d = 2; d * d <= r; d++) if (r % d == 0) { fac[nf++] = d; while (r % d == 0) r /= d; } if (r > 1) fac[nf++] = r; }
    pe g; int found = 0;
    for (uint64_t cand = 1; cand < q && !found; cand++) {
        uint64_t v = cand; pe_zero(&g);
        for (int i = 0; i < F; i++) { g.c[i] = (uint32_t)(v % (uint64_t)P); v /= (uint64_t)P; }
        found = 1;
        for (int i = 0; i < nf; i++) { pe t; pe_pow(&g, (q - 1) / fac[i], &t); if (pe_isone(&t)) { found = 0; break; } }
        if (q == 2) found = 1;
    }
    if (!found) { fprintf(stderr, "no generator found\n"); return 6; }
    pe e; pe_one(&e);
    for (uint64_t i = 0; i < q; i++) LOGI[i] = ZS;
    for (uint64_t k = 0; k + 1 < q; k++) {
        uint64_t idx = 0, pw = 1;
        for (int i = 0; i < F; i++) { idx += e.c[i] * pw; pw *= (uint64_t)P; }
        EXPI[k] = (uint32_t)idx;
        if (LOGI[idx] != ZS) { fprintf(stderr, "generator table collision\n"); return 6; }
        LOGI[idx] = (int32_t)k;
        pe t; pe_mul(&e, &g, &t); e = t;
    }
    for (uint64_t k = 0; k + 1 < q; k++) {
        uint64_t idx = EXPI[k];
        uint64_t plus1 = (idx % (uint64_t)P == (uint64_t)(P - 1)) ? idx - (uint64_t)(P - 1) : idx + 1;
        ZECH[k] = LOGI[plus1];
    }
    uint64_t zidx = 0, pw = 1;
    for (int i = 0; i < F; i++) { zidx += zeta.c[i] * pw; pw *= (uint64_t)P; }
    int64_t zlog = LOGI[zidx];
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) BUF[0][i * M + j] = (int32_t)((zlog * (int64_t)((i * j) % M)) % (int64_t)Q1);
#endif

    if (check_mode) {
        printf("CHECK %llx %s\n", (unsigned long long)check_mask, subset_is_singular(check_mask) ? "singular" : "nonsingular");
        return 0;
    }
    if (lowrank > 0) lowrank_scan(lowrank); else rec(0, M, M, 0, 0ULL);

    printf("%s M=%d p=%d f=%d h=", lowrank > 0 ? "LOWRANK" : "RESULT", M, P, F);
    for (int i = 0; i < F; i++) printf("%d%s", Hc[i], i + 1 < F ? "," : "");
    printf(" subsets=%llu zeros=%llu byrank=", (unsigned long long)nleaf, (unsigned long long)nzero);
    for (int r = 0; r <= M; r++) printf("%llu%s", (unsigned long long)zero_by_rank[r], r < M ? "," : "");
    printf("\n");
    return 0;
}
