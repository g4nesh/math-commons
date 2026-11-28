/*
 * window.c -- the window polynomial of a principal minor of F_(qM), exactly.
 *
 * M is a positive integer, q a prime not dividing M, h an irreducible factor of
 * the M-th cyclotomic polynomial modulo q of degree f (a "chart"), and
 * FF = F_q[x]/(h) with eta = x.  An index set K in Z_M x Z_q is given by its
 * fibres K_a.  With E(z) = sum_{m=0}^{q-1} z^m / m! in F_q[z], the program
 * computes the determinant of
 *
 *     H[(a,k), (a',b')] = eta^(a a') * sum_{i>=0} tau^i b'^(k+i) h_i(b_0..b_k) / (k+i)!
 *
 * (rows: a in Z_M and 0 <= k < |K_a|, with b_0 < b_1 < ... the elements of K_a
 * and h_i the complete homogeneous symmetric polynomial; columns: (a', b') in K;
 * terms with k + i > q - 1 are absent) as a power series in tau modulo tau^L.
 * H is obtained from the matrix ( eta^(a a') E(tau b b') ) by Newton divided
 * differences in b along every fibre, followed by division of row (a,k) by
 * tau^k, so
 *
 *     det( eta^(a a') E(tau b b') )_K = prod_a V(K_a) * tau^e * det H,
 *
 * e = sum_a C(|K_a|, 2), V = product of differences.  See README.md, Theorem A.
 *
 * The determinant is computed by elimination over the discrete valuation ring
 * FF[[tau]]: the pivot of a column is an entry of least valuation, so every
 * multiplier is integral; a pivot of valuation v costs v digits of absolute
 * precision in the remaining matrix, and the product of the pivots is still
 * known modulo tau^L.
 *
 * Build:  cc -O2 -o window window.c
 * Run:    window M q L f h_0 ... h_{f-1}  < sets
 *         h(x) = x^f + h_{f-1} x^{f-1} + ... + h_0.  Every input line holds M
 *         integers, the fibres K_0 ... K_(M-1) as bit masks over Z_q.
 * Output: one line per set: L integers, the coefficients of tau^0 ... tau^(L-1)
 *         of det H.  A field element sum_i c_i x^i is printed as sum_i c_i q^i.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXF 16
#define MAXL 64

static int M, Q, L, F;
static int Hc[MAXF];
static int64_t QF;            /* q^f */
static int32_t Q1, ZS;        /* q^f - 1; ZS = Q1 encodes zero in the log representation */
static int32_t *EXPT, *LOGT, *ZECH;

static inline int32_t t_mul(int32_t a, int32_t b) {
    if (a == ZS || b == ZS) return ZS;
    int32_t s = a + b; if (s >= Q1) s -= Q1; return s;
}
static inline int32_t t_add(int32_t a, int32_t b) {
    if (a == ZS) return b;
    if (b == ZS) return a;
    int32_t d = b - a; if (d < 0) d += Q1;
    int32_t z = ZECH[d]; if (z == ZS) return ZS;
    int32_t s = a + z; if (s >= Q1) s -= Q1; return s;
}
static int32_t NEGL;
static inline int32_t t_neg(int32_t a) { if (a == ZS) return ZS; int32_t s = a + NEGL; if (s >= Q1) s -= Q1; return s; }
static inline int32_t t_inv(int32_t a) { return a == 0 ? 0 : Q1 - a; }

/* multiplication of codes (digit vectors base q) modulo h, only for building the tables */
static int64_t code_mul(int64_t a, int64_t b) {
    int da[MAXF], db[MAXF], t[2 * MAXF];
    for (int i = 0; i < F; i++) { da[i] = (int)(a % Q); a /= Q; db[i] = (int)(b % Q); b /= Q; }
    for (int i = 0; i < 2 * F; i++) t[i] = 0;
    for (int i = 0; i < F; i++) for (int j = 0; j < F; j++) t[i + j] = (t[i + j] + da[i] * db[j]) % Q;
    for (int k = 2 * F - 2; k >= F; k--) {
        int c = t[k]; t[k] = 0;
        if (c) for (int j = 0; j < F; j++) t[k - F + j] = ((t[k - F + j] - c * Hc[j]) % Q + Q) % Q;
    }
    int64_t r = 0;
    for (int i = F - 1; i >= 0; i--) r = r * Q + t[i];
    return r;
}
static int64_t code_pow(int64_t a, int64_t e) {
    int64_t r = 1;
    while (e) { if (e & 1) r = code_mul(r, a); a = code_mul(a, a); e >>= 1; }
    return r;
}

static void build_tables(void) {
    QF = 1; for (int i = 0; i < F; i++) QF *= Q;
    if (QF > (1 << 24)) { fprintf(stderr, "field too large for tables\n"); exit(2); }
    Q1 = (int32_t)(QF - 1); ZS = Q1;
    int64_t pr[64]; int np = 0; int64_t r = Q1;
    for (int64_t d = 2; d * d <= r; d++) if (r % d == 0) { pr[np++] = d; while (r % d == 0) r /= d; }
    if (r > 1) pr[np++] = r;
    int64_t g = 0;
    for (int64_t cand = 1; cand < QF; cand++) {
        int ok = 1;
        for (int i = 0; i < np && ok; i++) if (code_pow(cand, Q1 / pr[i]) == 1) ok = 0;
        if (ok) { g = cand; break; }
    }
    if (!g) { fprintf(stderr, "no generator\n"); exit(2); }
    EXPT = malloc(sizeof(int32_t) * (size_t)QF); LOGT = malloc(sizeof(int32_t) * (size_t)QF); ZECH = malloc(sizeof(int32_t) * (size_t)QF);
    int64_t x = 1;
    for (int32_t k = 0; k < Q1; k++) { EXPT[k] = (int32_t)x; LOGT[x] = k; x = code_mul(x, g); }
    if (x != 1) { fprintf(stderr, "h is not irreducible\n"); exit(2); }
    LOGT[0] = ZS;
    for (int32_t k = 0; k < Q1; k++) {
        int64_t c = EXPT[k];
        int64_t c1 = (c % Q == Q - 1) ? c - (Q - 1) : c + 1;
        ZECH[k] = LOGT[c1];
    }
    NEGL = (Q == 2) ? 0 : Q1 / 2;
}

static inline int s_val(const int32_t *a, int cur) { for (int i = 0; i < cur; i++) if (a[i] != ZS) return i; return cur; }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: window M q L f h_0 ... h_{f-1} < sets\n"); return 2; }
    M = atoi(argv[1]); Q = atoi(argv[2]); L = atoi(argv[3]); F = atoi(argv[4]);
    if (argc != 5 + F || F < 1 || F > MAXF || L < 1 || L > MAXL || Q > 31) { fprintf(stderr, "bad arguments\n"); return 2; }
    for (int i = 0; i < F; i++) Hc[i] = ((atoi(argv[5 + i]) % Q) + Q) % Q;
    build_tables();
    int32_t ilog[32];
    for (int n = 0; n < Q; n++) ilog[n] = LOGT[n];
    int32_t eta = (F == 1) ? LOGT[(Q - Hc[0]) % Q] : LOGT[Q];
    int32_t *ep = malloc(sizeof(int32_t) * (size_t)M);
    for (int d = 0; d < M; d++) ep[d] = (int32_t)(((int64_t)eta * d) % Q1);
    if (M > 1 && (((int64_t)eta * M) % Q1 != 0)) { fprintf(stderr, "x is not an M-th root of unity\n"); return 2; }
    int invfact[32]; { int fct = 1; for (int n = 0; n < Q; n++) { if (n) fct = fct * n % Q; int inv = 1; for (int e = 0; e < Q - 2; e++) inv = inv * fct % Q; invfact[n] = (Q == 2) ? 1 : inv; } }
    int maxn = M * Q;
    int32_t *A = malloc(sizeof(int32_t) * (size_t)maxn * (size_t)maxn * (size_t)L);
    int *fib = malloc(sizeof(int) * (size_t)M);
    int *ra = malloc(sizeof(int) * (size_t)maxn), *rk = malloc(sizeof(int) * (size_t)maxn), *cb = malloc(sizeof(int) * (size_t)maxn);
    int32_t U[MAXL], u[MAXL], ui[MAXL], mu[MAXL], tmp[MAXL];
    for (;;) {
        int ok = 1;
        for (int a = 0; a < M; a++) if (scanf("%d", &fib[a]) != 1) { ok = 0; break; }
        if (!ok) break;
        int n = 0;
        for (int a = 0; a < M; a++) { int k = 0; for (int b = 0; b < Q; b++) if (fib[a] >> b & 1) { ra[n] = a; rk[n] = k++; cb[n] = b; n++; } }
        /* build H */
        for (int r = 0; r < n; r++) {
            int a = ra[r], k = rk[r];
            int nodes[32], nn = 0;
            for (int b = 0; b < Q && nn <= k; b++) if (fib[a] >> b & 1) nodes[nn++] = b;
            int hs[MAXL];                      /* h_i(nodes[0..k]) mod q */
            for (int i = 0; i < L; i++) hs[i] = (i == 0);
            for (int j = 0; j <= k; j++) for (int i = 1; i < L; i++) hs[i] = (hs[i] + nodes[j] * hs[i - 1]) % Q;
            for (int c = 0; c < n; c++) {
                int a2 = ra[c], b2 = cb[c];
                int32_t *e = A + ((size_t)r * (size_t)n + (size_t)c) * (size_t)L;
                int32_t root = ep[(a * a2) % M];
                int pw = 1;
                for (int j = 0; j < k; j++) pw = pw * b2 % Q;   /* b2^k, with 0^0 = 1 */
                for (int i = 0; i < L; i++) {
                    int v = (k + i <= Q - 1) ? pw * hs[i] % Q * invfact[k + i] % Q : 0;
                    e[i] = v ? t_mul(ilog[v], root) : ZS;
                    pw = pw * b2 % Q;
                }
            }
        }
        /* elimination over FF[[tau]] */
        int cur = L, V = 0, neg = 0, zero = 0;
        U[0] = 0; for (int i = 1; i < L; i++) U[i] = ZS;
        for (int c = 0; c < n && !zero; c++) {
            int best = -1, bv = cur;
            for (int r = c; r < n; r++) { int v = s_val(A + ((size_t)r * n + c) * L, cur); if (v < bv) { bv = v; best = r; } }
            if (best < 0) { zero = 1; break; }
            if (best != c) {
                for (int j = c; j < n; j++) {
                    int32_t *x = A + ((size_t)c * n + j) * L, *y = A + ((size_t)best * n + j) * L;
                    for (int i = 0; i < cur; i++) { int32_t t = x[i]; x[i] = y[i]; y[i] = t; }
                }
                neg ^= 1;
            }
            int nc = cur - bv;
            const int32_t *piv = A + ((size_t)c * n + c) * L;
            for (int i = 0; i < nc; i++) u[i] = piv[i + bv];
            /* ui = 1/u mod tau^nc */
            ui[0] = t_inv(u[0]);
            for (int i = 1; i < nc; i++) {
                int32_t s = ZS;
                for (int j = 1; j <= i; j++) s = t_add(s, t_mul(u[j], ui[i - j]));
                ui[i] = t_neg(t_mul(s, ui[0]));
            }
            for (int i = 0; i < nc; i++) { int32_t s = ZS; for (int j = 0; j <= i; j++) s = t_add(s, t_mul(U[j], u[i - j])); tmp[i] = s; }
            for (int i = 0; i < nc; i++) U[i] = tmp[i];
            V += bv;
            for (int r = c + 1; r < n; r++) {
                const int32_t *x = A + ((size_t)r * n + c) * L;
                if (s_val(x, cur) >= cur) continue;
                int allz = 1;
                for (int i = 0; i < nc; i++) { int32_t s = ZS; for (int j = 0; j <= i; j++) s = t_add(s, t_mul(x[j + bv], ui[i - j])); mu[i] = s; if (s != ZS) allz = 0; }
                if (allz) continue;
                for (int j = c + 1; j < n; j++) {
                    const int32_t *p = A + ((size_t)c * n + j) * L;
                    int32_t *y = A + ((size_t)r * n + j) * L;
                    for (int i = 0; i < nc; i++) {
                        int32_t s = ZS;
                        for (int t = 0; t <= i; t++) if (mu[t] != ZS && p[i - t] != ZS) s = t_add(s, t_mul(mu[t], p[i - t]));
                        if (s != ZS) y[i] = t_add(y[i], t_neg(s));
                    }
                }
            }
            cur = nc;
        }
        for (int j = 0; j < L; j++) {
            int32_t v = ZS;
            if (!zero && j >= V && j - V < cur) v = neg ? t_neg(U[j - V]) : U[j - V];
            printf(j ? " %d" : "%d", v == ZS ? 0 : EXPT[v]);
        }
        printf("\n");
    }
    return 0;
}
