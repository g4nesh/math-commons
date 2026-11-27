"""Checker B: principal minors of (zeta^(ij)), zeta of order 14 in F_1331, via Zech/log tables.
Elements are stored as discrete logs base g (ZERO = -1). Determinants: Bareiss (fraction-free)
elimination; rechecks by bitmask-DP Laplace (cofactor) expansion. Standard library only."""
import itertools, json, os, random, time
from math import comb

P, Q = 11, 1330  # characteristic, order of multiplicative group of F_1331
ZERO = -1
T0 = time.perf_counter()


# ---------- field construction ----------
def xtimes(v, c):  # multiply (c0,c1,c2) by x modulo x^3 + a x^2 + b x + d
    a, b, d = c
    c0, c1, c2 = v
    return ((-c2 * d) % P, (c0 - c2 * b) % P, (c1 - c2 * a) % P)


def x_order(c):
    v, k = (0, 1, 0), 1
    while v != (1, 0, 0):
        v = xtimes(v, c)
        k += 1
        if k > Q:
            return None
    return k


def find_modulus():
    for a, b, d in itertools.product(range(P - 1, -1, -1), repeat=3):  # descending lex
        if any((r**3 + a * r * r + b * r + d) % P == 0 for r in range(P)):
            continue  # has a root -> reducible
        if x_order((a, b, d)) == Q:
            return (a, b, d)


C = find_modulus()
print(f"modulus: x^3 + {C[0]} x^2 + {C[1]} x + {C[2]} over F_11 (first irreducible primitive cubic, descending lex)")
vec = [(1, 0, 0)]
for _ in range(Q - 1):
    vec.append(xtimes(vec[-1], C))
enc = lambda v: v[0] + P * v[1] + P * P * v[2]
EXP = [enc(v) for v in vec]  # exp[k] = x^k as integer code
LOG = {e: k for k, e in enumerate(EXP)}
assert len(LOG) == Q and sorted(EXP) == list(range(1, P**3)), "log table not a bijection onto nonzero elements"
ZECH = []  # ZECH[d] = log(1 + x^d), ZERO if 1+x^d = 0
for d in range(Q):
    e = EXP[d]
    e1 = (e % P + 1) % P + (e - e % P)
    ZECH.append(LOG[e1] if e1 else ZERO)

HALF = Q // 2  # log(-1)


def mul(a, b): return ZERO if a < 0 or b < 0 else (a + b) % Q
def neg(a): return a if a < 0 else (a + HALF) % Q


def add(a, b):
    if a < 0: return b
    if b < 0: return a
    z = ZECH[(b - a) % Q]
    return ZERO if z < 0 else (a + z) % Q


def sub(a, b): return add(a, neg(b))
def fdiv(a, b): return ZERO if a < 0 else (a - b) % Q  # b != ZERO


def bareiss(M):
    """Fraction-free (Bareiss) determinant of a log-matrix; returns log or ZERO. Exact."""
    n = len(M)
    if n == 0:
        return 0  # empty determinant = 1 = g^0
    M = [r[:] for r in M]
    flip = False
    prev = 0
    for k in range(n - 1):
        if M[k][k] < 0:
            p = next((i for i in range(k + 1, n) if M[i][k] >= 0), None)
            if p is None:
                return ZERO
            M[k], M[p] = M[p], M[k]
            flip = not flip
        pk, rk = M[k][k], M[k]
        for i in range(k + 1, n):
            ri = M[i]
            ik = ri[k]
            for j in range(k + 1, n):
                ri[j] = fdiv(sub(mul(ri[j], pk), mul(ik, rk[j])), prev)
        prev = pk
    d = M[n - 1][n - 1]
    return neg(d) if flip else d


def cofactor_det(M):
    """Laplace expansion by rows, DP over column subsets (independent of Bareiss)."""
    n = len(M)
    f = {0: 0}
    for r in range(n):
        g = {}
        for cols in itertools.combinations(range(n), r + 1):
            mask = sum(1 << c for c in cols)
            tot = ZERO
            for pos, j in enumerate(cols):
                t = mul(M[r][j], f[mask ^ (1 << j)])
                tot = add(tot, neg(t) if (r + pos) & 1 else t)  # cofactor sign (-1)^(row+col)
            g[mask] = tot
        f = g
    return f[(1 << n) - 1] if n else 0


def build(Nn, zl, S):  # log-matrix of zeta^(ij), zeta = g^zl of order Nn (Nn | 1330)
    return [[(zl * ((i * j) % Nn)) % Q for j in S] for i in S]


# ---------- self-tests ----------
def pmul(u, v):  # independent polynomial-basis multiplication mod c
    a, b, d = C
    p = [0] * 5
    for i in range(3):
        for j in range(3):
            p[i + j] = (p[i + j] + u[i] * v[j]) % P
    for k in (4, 3):
        t = p[k]
        p[k] = 0
        p[k - 3] = (p[k - 3] - t * d) % P
        p[k - 2] = (p[k - 2] - t * b) % P
        p[k - 1] = (p[k - 1] - t * a) % P
    return tuple(p[:3])


def ppow(y, e):
    r = (1, 0, 0)
    while e:
        if e & 1:
            r = pmul(r, y)
        y = pmul(y, y)
        e >>= 1
    return r


results = {}


def check(name, ok, extra=""):
    results[name] = bool(ok)
    print(f"SELFTEST {'PASS' if ok else 'FAIL'}: {name} {extra}")
    assert ok, name


rs = random.Random(1)
allv = [(a, b, c) for a in range(P) for b in range(P) for c in range(P) if (a, b, c) != (0, 0, 0)]
check("T1a log table is a bijection onto the 1330 nonzero elements", all(LOG[enc(vec[k])] == k for k in range(Q)))
check("T1b Frobenius y^1331 == y on all 1331 elements (polynomial-basis pow)",
      all(ppow(y, P**3) == y for y in allv) and ppow((0, 0, 0), P**3) == (0, 0, 0))
check("T1c exp table agrees with polynomial multiplication (500 random pairs)",
      all(vec[(k + m) % Q] == pmul(vec[k], vec[m]) for k, m in [(rs.randrange(Q), rs.randrange(Q)) for _ in range(500)]))


def vadd_log(a, b):
    s = tuple((x + y) % P for x, y in zip(vec[a], vec[b]))
    return LOG[enc(s)] if any(s) else ZERO


check("T1d Zech addition agrees with vector addition (2000 random pairs)",
      all(vadd_log(a, b) == add(a, b) for a, b in [(rs.randrange(Q), rs.randrange(Q)) for _ in range(2000)]))

ZL = 95  # zeta = g^95
zp = lambda k: ppow(vec[ZL], k)
check("zeta = g^95: zeta^14 = 1, zeta^7 != 1, zeta^2 != 1", zp(14) == (1, 0, 0) and zp(7) != (1, 0, 0) and zp(2) != (1, 0, 0))
ROOTS = {"zeta": ZL, "zeta^-1": (-ZL) % Q}
N = 14

check("T2a empty set -> 1", all(bareiss(build(N, z, [])) == 0 and cofactor_det([]) == 0 for z in ROOTS.values()))
check("T2b singleton {i} -> zeta^(i^2) (nonzero; =1 for i=0)", all(bareiss(build(N, z, [i])) == (z * i * i) % Q for z in ROOTS.values() for i in range(N)))
ok = True
for _ in range(300):
    S = sorted(rs.sample(range(N), rs.randint(2, 9)))
    for z in ROOTS.values():
        M = build(N, z, S)
        ok &= bareiss(M) == cofactor_det(M)
check("T2c Bareiss == cofactor expansion on 600 random minors (size 2..9)", ok)
Mr = [[rs.randrange(Q) for _ in range(6)] for _ in range(6)]
Mr[2][3] = ZERO
Mr[0][0] = ZERO
Mr[1][0] = ZERO
check("T2d Bareiss == cofactor on random 6x6 with zero entries (pivot swaps)", bareiss(Mr) == cofactor_det(Mr))

full = list(range(N))
Mz, Mi = build(N, ZL, full), build(N, (-ZL) % Q, full)
dz, di = bareiss(Mz), bareiss(Mi)
log3 = LOG[3]  # 14 = 3 in F_11
check("T3a det M(zeta)*det M(zeta^-1) == 14^14 (nonzero)", dz >= 0 and di >= 0 and mul(dz, di) == (14 * log3) % Q, f"(logs {dz},{di})")
prod_ok = True
for i in range(N):
    for j in range(N):
        s = ZERO
        for k in range(N):
            s = add(s, mul(Mz[i][k], Mi[k][j]))
        prod_ok &= (s == log3) if i == j else (s == ZERO)
check("T3b matrix product M(zeta) M(zeta^-1) == 14 I, entrywise", prod_ok)
check("T3c full-set Bareiss equals cofactor DP (14x14) for both roots", cofactor_det(Mz) == dz and cofactor_det(Mi) == di)
R = [[rs.randrange(Q) for _ in range(5)] for _ in range(4)]
R.insert(3, R[0][:])
check("T4a singular by construction (repeated row, 5x5) -> 0 in both algorithms", bareiss(R) == ZERO and cofactor_det(R) == ZERO)
check("T4b [[0,1],[1,0]] needs a pivot swap, det = -1", bareiss([[ZERO, 0], [0, ZERO]]) == neg(0))
z2 = (2 * ZL) % Q  # zeta^2 has order 7; rows for i=0,7 coincide on columns {0,7}
check("T4c non-primitive root zeta^2 with S={0,7}: det = 0", bareiss([[(z2 * ((i * j) % 14)) % Q for j in (0, 7)] for i in (0, 7)]) == ZERO)


def analog(Nn, zl, label):
    zs = []
    for m in range(1 << Nn):
        S = [i for i in range(Nn) if m >> i & 1]
        M = build(Nn, zl, S)
        d = bareiss(M)
        if d == ZERO:
            assert cofactor_det(M) == ZERO
            zs.append(S)
    print(f"  analog {label}: {1 << Nn} subsets, {len(zs)} zero principal minors" + (f", e.g. {zs[:3]}" if zs else ""))
    return zs


print("T5 analogues (same pipeline; N=11 skipped: no primitive 11th root in characteristic 11):")
a7 = [analog(7, 190 * u % Q, f"N=7 in F_1331, zeta=g^(190*{u})") for u in range(1, 7)]
a10 = analog(10, 133, "N=10, roots in F_11, zeta=g^133")
results["T5 analogues ran (N=7 six roots, N=10)"] = True
results["T5 N=7 zero counts per root"] = [len(z) for z in a7]
results["T5 N=10 zero count"] = len(a10)
print(f"  (self-tests done at {time.perf_counter() - T0:.1f}s)")

# ---------- main run ----------
UNITS = [1, 3, 5, 9, 11, 13]


def image(m, a, t):
    r = 0
    for i in range(N):
        if m >> i & 1:
            r |= 1 << ((a * i + t) % N)
    return r


NM = 1 << N
par = list(range(NM))  # union-find on generators: translation by 1, dilation by 3 (3 generates the units)


def find(x):
    while par[x] != x:
        par[x] = par[par[x]]
        x = par[x]
    return x


for m in range(NM):
    for (a, t) in ((1, 1), (3, 0)):
        u, v = find(m), find(image(m, a, t))
        if u != v:
            par[max(u, v)] = min(u, v)
rep = [find(m) for m in range(NM)]
reps = sorted(set(rep))
orb_size = {r: 0 for r in reps}
for m in range(NM):
    orb_size[rep[m]] += 1
fix = sum(1 for a in UNITS for t in range(N) for m in range(NM) if image(m, a, t) == m)  # Burnside, 84 maps
check("orbit count: union-find == Burnside over all 84 affine maps", fix % 84 == 0 and fix // 84 == len(reps), f"({len(reps)} orbits)")
pop = lambda m: bin(m).count("1")
orb_by_size = [0] * (N + 1)
for r in reps:
    orb_by_size[pop(r)] += 1
print("orbits by subset size:", orb_by_size)


def status(m, z):  # True if nonzero
    return bareiss(build(N, z, [i for i in range(N) if m >> i & 1])) != ZERO


t = time.perf_counter()
q_stat = {name: {r: status(r, z) for r in reps} for name, z in ROOTS.items()}
tq = time.perf_counter() - t
print(f"quotient run: {len(reps)} reps x 2 roots = {2 * len(reps)} determinants in {tq:.2f}s")
t = time.perf_counter()
f_stat = {name: [status(m, z) for m in range(NM)] for name, z in ROOTS.items()}
tf = time.perf_counter() - t
print(f"full run: {NM} subsets x 2 roots = {2 * NM} determinants in {tf:.2f}s")

bad = [(name, m) for name in ROOTS for m in range(NM) if f_stat[name][m] != q_stat[name][rep[m]]]
check("symmetry claim: every orbit member has same zero/nonzero status as its rep, per root", not bad, f"({len(bad)} mismatches over {2 * NM} cases)")
bad2 = sum(f_stat[name][m] != f_stat[name][image(m, a, tt)] for name in ROOTS for m in range(0, NM, 7) for a in UNITS for tt in (0, 1, 5, 13))
check("symmetry claim spot-check over explicit (a,t) images (all 6 units)", bad2 == 0)

counts, zero_cases = {}, []
for name in ROOTS:
    nz, zr, qnz, qzr = ([0] * (N + 1) for _ in range(4))
    for m in range(NM):
        if f_stat[name][m]:
            nz[pop(m)] += 1
        else:
            zr[pop(m)] += 1
            zero_cases.append((name, m))
    for r in reps:
        (qnz if q_stat[name][r] else qzr)[pop(r)] += orb_size[r]
    assert nz == qnz and zr == qzr
    counts[name] = {"nonzero_by_size": nz, "zero_by_size": zr, "total_nonzero": sum(nz), "total_zero": sum(zr)}
    print(f"\ncounts by size, root {name} (full run; quotient weighted by orbit size gives identical rows)")
    print("size    " + " ".join(f"{s:5d}" for s in range(N + 1)))
    print("C(14,s) " + " ".join(f"{comb(N, s):5d}" for s in range(N + 1)))
    print("nonzero " + " ".join(f"{x:5d}" for x in nz))
    print("zero    " + " ".join(f"{x:5d}" for x in zr))
    print(f"total nonzero {sum(nz)}, zero {sum(zr)}")
check("quotient and full runs agree (per size, per root)", True)

for name, m in zero_cases:
    S = [i for i in range(N) if m >> i & 1]
    M = build(N, ROOTS[name], S)
    print(f"ZERO MINOR: root {name}, S={S}, Bareiss={bareiss(M)}, cofactor={cofactor_det(M)} (-1 = field element 0)")

verdict = not zero_cases
print("\nVERDICT: " + ("Hyp(14,11) HOLDS: all 2^14 principal minors nonzero in F_1331 for both zeta and zeta^-1 (32768 determinants)"
                       if verdict else "Hyp(14,11) FAILS: see ZERO MINOR lines"))
total = time.perf_counter() - T0
print(f"total runtime {total:.1f}s")
json.dump({"modulus_x3_ax2_bx_d": list(C), "zeta": "g^95", "orbits": len(reps), "orbits_by_size": orb_by_size,
           "counts": counts, "zero_cases": [[n, [i for i in range(N) if m >> i & 1]] for n, m in zero_cases],
           "verdict_all_nonzero": verdict, "quotient_full_agree": True, "selftests": results,
           "runtime_seconds": {"quotient": round(tq, 2), "full": round(tf, 2), "total": round(total, 2)}},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "counts.json"), "w"), indent=1)
