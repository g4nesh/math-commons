"""Generalization map n = 4 -> n = 5 for the Zheng MOC(4) proof (audit label n5).

Own code only. Reads nothing from the supplement. Deterministic seeds.
Run one part at a time:  python -I -B generalization_n5.py PART
PART in {alg, burnside, facets, uo}.

Levels: functions marked EXACT use integer / rational arithmetic (python-flint) or
arithmetic modulo the prime p = 2^31 - 1 at random points (stated as such);
functions marked NUMERIC use floating point.
"""
import itertools, math, sys, time, os
from collections import Counter
from fractions import Fraction

import numpy as np
import flint

PRIME = 2147483647


# ----------------------------------------------------------------------------- basics
def perms(n):
    return list(itertools.permutations(range(n)))  # lexicographic one-line order


def sgn(p):
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True; j = p[j]; L += 1
            if L % 2 == 0:
                s = -s
    return s


def subsets(n, k):
    return list(itertools.combinations(range(n), k))


def indicator_rows(n, ks):
    """chi_{I,J}(sigma) = [sigma(I) = J] for |I| = |J| = k in ks; rows over perms(n)."""
    P = perms(n)
    rows, keys = [], []
    for k in ks:
        for I in subsets(n, k):
            for J in subsets(n, k):
                Js = set(J)
                rows.append([1 if {s[i] for i in I} == Js else 0 for s in P])
                keys.append((I, J))
    return rows, keys


def zrank(rows):
    if not rows:
        return 0
    return flint.fmpz_mat(rows).rank()


# ----------------------------------------------------------------------------- brackets
def bracket(i, j):
    """[ij] = x_i y_j - x_j y_i as {mask: coeff}; bit l set <=> x_l used."""
    return {1 << i: 1, 1 << j: -1}


def pmul(p, q):
    out = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = m1 | m2
            out[m] = out.get(m, 0) + c1 * c2
    return {m: c for m, c in out.items() if c}


def matching_poly(pairs):
    p = {0: 1}
    for i, j in pairs:
        p = pmul(p, bracket(min(i, j), max(i, j)))
    return p


def all_matchings(labels):
    if not labels:
        yield []
        return
    a = labels[0]
    for t in range(1, len(labels)):
        rest = labels[1:t] + labels[t + 1:]
        for m in all_matchings(rest):
            yield [(a, labels[t])] + m


def noncrossing(labels):
    if not labels:
        yield []
        return
    a = labels[0]
    for t in range(1, len(labels), 2):
        for m1 in noncrossing(labels[1:t]):
            for m2 in noncrossing(labels[t + 1:]):
                yield [(a, labels[t])] + m1 + m2


def poly_rows(polys, n):
    cols = [m for m in range(1 << (2 * n)) if bin(m).count("1") == n]
    idx = {m: i for i, m in enumerate(cols)}
    rows = []
    for p in polys:
        r = [0] * len(cols)
        for m, c in p.items():
            r[idx[m]] = c
        rows.append(r)
    return rows


def catalan(n):
    return math.comb(2 * n, n) // (n + 1)


def bip_pairs(n, s):
    return [(i, n + s[i]) for i in range(n)]


# ----------------------------------------------------------------------------- PART alg
def exact_P(n):
    """EXACT: P_n with p_sigma = sum_a P[sigma,a] B_a (B = noncrossing matchings)."""
    NC = list(noncrossing(list(range(2 * n))))
    Nrows = poly_rows([matching_poly(m) for m in NC], n)
    Brows = poly_rows([matching_poly(bip_pairs(n, s)) for s in perms(n)], n)
    N = flint.fmpq_mat(flint.fmpz_mat(Nrows))
    B = flint.fmpq_mat(flint.fmpz_mat(Brows))
    NNt = N * N.transpose()
    Pq = B * N.transpose() * NNt.inv()
    assert Pq * N == B, "p_sigma not in the span of the noncrossing basis"
    P = [[int(Pq[i, j]) for j in range(Pq.ncols())] for i in range(Pq.nrows())]
    assert all(Pq[i, j] == P[i][j] for i in range(Pq.nrows()) for j in range(Pq.ncols()))
    return NC, P


def solve_rows_q(P, target):
    """EXACT least-squares-free solve P c = target (checked)."""
    A = flint.fmpq_mat(flint.fmpz_mat(P))
    At = A.transpose()
    tq = flint.fmpq_mat(len(target), 1, [flint.fmpq(Fraction(t).numerator, Fraction(t).denominator) for t in target])
    c = (At * A).inv() * At * tq
    assert A * c == tq, "target not in image of P"
    return [c[i, 0] for i in range(c.nrows())]


def part_alg():
    print("=" * 78)
    print("PART alg: matching cone X_n, Catalan ranks, P_n, kappa, entry subspace, dim X_n")
    print("=" * 78)
    for n in (3, 4, 5, 6):
        t0 = time.time()
        NC = list(noncrossing(list(range(2 * n))))
        nc_rank = zrank(poly_rows([matching_poly(m) for m in NC], n))
        bip = [matching_poly(bip_pairs(n, s)) for s in perms(n)]
        bip_rank = zrank(poly_rows(bip, n))
        both = zrank(poly_rows([matching_poly(m) for m in NC] + bip, n))
        line = (f"n={n}: Catalan C_n={catalan(n)}; #noncrossing={len(NC)} rank={nc_rank}; "
                f"#bipartite p_sigma={len(bip)} rank={bip_rank}; rank(NC u bip)={both}")
        if n <= 5:
            allm = list(all_matchings(list(range(2 * n))))
            line += f"; all {len(allm)} perfect matchings rank={zrank(poly_rows([matching_poly(m) for m in allm], n))}"
        print(line + f"  [EXACT, {time.time() - t0:.1f}s]")
        # indicator span on S_n
        if n <= 5:
            chi_all, _ = indicator_rows(n, range(0, n + 1))
            chi1, _ = indicator_rows(n, [0, 1])
            chi12, _ = indicator_rows(n, [0, 1, 2])
            print(f"      indicator functions chi_IJ on S_{n}: rank(all k)={zrank(chi_all)}, "
                  f"rank(k<=2)={zrank(chi12)}, rank(entry k<=1)={zrank(chi1)} "
                  f"(expected (n-1)^2+1={(n - 1) ** 2 + 1})  [EXACT]")

    for n in (4, 5):
        print("-" * 78)
        t0 = time.time()
        NC, P = exact_P(n)
        print(f"n={n}: P_n is an integer {len(P)}x{len(P[0])} matrix, entries in "
              f"{sorted(set(x for r in P for x in r))}, rank {zrank(P)}  [EXACT, {time.time() - t0:.1f}s]")
        chi_all, _ = indicator_rows(n, range(0, n + 1))
        colsP = [list(r) for r in zip(*P)]
        print(f"      rank[columns of P ; all chi_IJ] = {zrank(colsP + chi_all)}  "
              f"(= rank P means im P = span of the chi_IJ)  [EXACT]")
        kappa = solve_rows_q(P, [1] * len(P))
        print(f"      kappa (P kappa = 1) in noncrossing coordinates: "
              f"{[str(x) for x in kappa]}  [EXACT]")
        # dim X_n via Jacobian at an integer point
        rng = np.random.default_rng(100 + n)
        pt = [(int(rng.integers(-9, 10)), int(rng.integers(1, 10))) for _ in range(2 * n)]
        polys = [matching_poly(m) for m in NC]

        def evalmono(mask, skip=None):
            v = 1
            for l in range(2 * n):
                if l == skip:
                    continue
                v *= pt[l][0] if (mask >> l) & 1 else pt[l][1]
            return v
        J = []
        for p in polys:
            row = []
            for l in range(2 * n):
                dx = sum(c * evalmono(m, l) for m, c in p.items() if (m >> l) & 1)
                dy = sum(c * evalmono(m, l) for m, c in p.items() if not (m >> l) & 1)
                row += [dx, dy]
            J.append(row)
        rJ = zrank(J)
        # kernel fields: sl2 (3) and torus with prod t = 1 (2n-1)
        K = []
        for X in ([[0, 1], [0, 0]], [[0, 0], [1, 0]], [[1, 0], [0, -1]]):
            vec = []
            for (x, y) in pt:
                vec += [X[0][0] * x + X[0][1] * y, X[1][0] * x + X[1][1] * y]
            K.append(vec)
        for l in range(2 * n - 1):
            vec = [0] * (4 * n)
            vec[2 * l], vec[2 * l + 1] = pt[l]
            vec[2 * l + 2], vec[2 * l + 3] = -pt[l + 1][0], -pt[l + 1][1]
            K.append(vec)
        JK = flint.fmpz_mat(J) * flint.fmpz_mat(K).transpose()
        print(f"      dim X_{n}: Jacobian rank at integer point {pt} is {rJ}; the {len(K)} "
              f"invariance fields are independent (rank {zrank(K)}) and J*K = 0: {JK.is_zero()}; "
              f"so generic rank = {4 * n} - {len(K)} = {4 * n - len(K)} = 2n-2  [EXACT]")
        quadric_analysis(n, NC, P, kappa)


def nc_values_mod(NC, z):
    out = []
    for m in NC:
        v = 1
        for i, j in m:
            v = v * (z[i] - z[j]) % PRIME
        out.append(v)
    return out


def quadric_analysis(n, NC, P, kappa):
    """Quadrics through X_n (mod p, random points) and K = {x : grad f(kappa).x = 0}."""
    t0 = time.time()
    D = len(NC)
    rng = np.random.default_rng(7 + n)
    monos = [(a, b) for a in range(D) for b in range(a, D)]
    npts = len(monos) + 60

    def randz():
        return [int(rng.integers(1, PRIME)) for _ in range(2 * n)]
    rows = []
    for _ in range(npts):
        m = nc_values_mod(NC, randz())
        rows.append([m[a] * m[b] % PRIME for a, b in monos])
    M = flint.nmod_mat(rows, PRIME)
    rk = M.rank()
    X, nul = M.nullspace()
    quads = [[int(X[r, c]) for r in range(len(monos))] for c in range(nul)]
    # fresh-point check of every quadric
    bad = 0
    for _ in range(40):
        m = nc_values_mod(NC, randz())
        mv = [m[a] * m[b] % PRIME for a, b in monos]
        for q in quads:
            if sum(x * y for x, y in zip(q, mv)) % PRIME:
                bad += 1
    print(f"      quadrics through X_{n}: {len(monos)} quadratic monomials in {D} coordinates; "
          f"evaluation rank {rk} at {npts} random points => dim I_2 = {nul}, "
          f"dim R_2 = {rk}; {bad} failures at 40 fresh points  "
          f"[EXACT mod p={PRIME}, randomized; {time.time() - t0:.1f}s]")
    # gradients at kappa
    kap = [int(x.p) * pow(int(x.q), -1, PRIME) % PRIME for x in kappa]
    G = []
    for q in quads:
        g = [0] * D
        for (a, b), c in zip(monos, q):
            if c:
                g[a] = (g[a] + c * kap[b]) % PRIME
                g[b] = (g[b] + c * kap[a]) % PRIME
        G.append(g)
    Gm = flint.nmod_mat(G, PRIME)
    rG = Gm.rank()
    KX, knul = Gm.nullspace()
    Kbasis = [[int(KX[r, c]) for r in range(D)] for c in range(knul)]
    # entry subspace E_n in noncrossing coordinates
    chi1, keys1 = indicator_rows(n, [1])
    E = []
    for row in chi1:
        e = solve_rows_q(P, row)
        E.append([int(x.p) * pow(int(x.q), -1, PRIME) % PRIME for x in e])
    rE = flint.nmod_mat(E, PRIME).rank()
    rKE = flint.nmod_mat(Kbasis + E, PRIME).rank()
    kap_in_K = flint.nmod_mat(Kbasis + [kap], PRIME).rank() == knul
    print(f"      K = common kernel of grad f(kappa), f in I_2: rank of gradients {rG}, dim K = {knul}; "
          f"dim E_{n} = {rE}; rank[K;E] = {rKE} -> "
          f"{'K = E' if (rKE == knul == rE) else 'K != E'}; kappa in K: {kap_in_K}  [EXACT mod p]")
    if n == 5 and not (rKE == knul == rE):
        print("      (K strictly larger than E: see analysis below)")
    tangent_tests(n, NC, P, Kbasis, E)


def nc_jac_mod(NC, z):
    """d m / d z_l for v_l = (z_l, 1): [ij] = z_i - z_j."""
    L = len(z)
    J = []
    for m in NC:
        row = [0] * L
        facs = [(i, j, (z[i] - z[j]) % PRIME) for i, j in m]
        for t, (i, j, _) in enumerate(facs):
            rest = 1
            for u, (_, _, f) in enumerate(facs):
                if u != t:
                    rest = rest * f % PRIME
            row[i] = (row[i] + rest) % PRIME
            row[j] = (row[j] - rest) % PRIME
        J.append(row)
    return J


def tangent_tests(n, NC, P, Kbasis, E):
    rng = np.random.default_rng(31 + n)
    D = len(NC)
    kap_vals = None
    # kappa in noncrossing coords (mod p) is the last-solved vector; recompute
    kap = solve_rows_q(P, [1] * len(P))
    kap = [int(x.p) * pow(int(x.q), -1, PRIME) % PRIME for x in kap]
    res = []
    for trial in range(4):
        z = [int(rng.integers(1, PRIME)) for _ in range(2 * n)]
        if trial == 3:
            # first spectrum with n-1 equal entries: z_1 = .. = z_{n-1}? use a_2..a_n equal
            for l in range(2, n):
                z[l] = z[1]
        m = nc_values_mod(NC, z)
        J = nc_jac_mod(NC, z)
        T = [list(c) for c in zip(*J)] + [m]  # tangent space spanned by columns + Euler
        rT = flint.nmod_mat(T, PRIME).rank()
        rTk = flint.nmod_mat(T + [kap], PRIME).rank()
        inK = flint.nmod_mat(Kbasis + [m], PRIME).rank() == len(Kbasis)
        inE = flint.nmod_mat(E + [m], PRIME).rank() == flint.nmod_mat(E, PRIME).rank()
        tag = "a_2=..=a_n (n-1 equal)" if trial == 3 else "random point"
        res.append(f"        {tag}: dim T_xX = {rT}, kappa in T_xX: {rTk == rT}, x in K: {inK}, x in E: {inE}")
    print(f"      tangent test on X_{n} (mod p, random points):")
    print("\n".join(res))


# ----------------------------------------------------------------------------- PART burnside
def group_action(n):
    P = perms(n)
    idx = {p: i for i, p in enumerate(P)}
    inv = [tuple(np.argsort(p)) for p in P]
    acts = []
    for r in P:
        for r2 in P:
            for eps in (0, 1):
                a = []
                for s in P:
                    t = inv[idx[s]] if eps else s
                    a.append(idx[tuple(r[t[r2[i]]] for i in range(n))])
                acts.append(a)
    return np.array(acts, dtype=np.int32)


def cycle_type(a):
    seen = np.zeros(len(a), bool)
    ct = []
    for i in range(len(a)):
        if not seen[i]:
            L, j = 0, i
            while not seen[j]:
                seen[j] = True; j = a[j]; L += 1
            ct.append(L)
    return tuple(sorted(ct))


def orbit_counts(n, kmax):
    G = group_action(n)
    assert len({tuple(g) for g in G}) == len(G), "action not faithful"
    types = Counter(cycle_type(g) for g in G)
    total = [0] * (kmax + 1)
    for ct, mult in types.items():
        poly = [1] + [0] * kmax
        for L in ct:
            new = poly[:]
            for d in range(kmax + 1 - L):
                new[d + L] += poly[d]
            poly = new
        for k in range(kmax + 1):
            total[k] += mult * poly[k]
    return len(G), [t // len(G) for t in total], [t % len(G) for t in total]


def part_burnside():
    print("=" * 78)
    print("PART burnside: orbits of k-subsets of S_n under G_n = (S_n x S_n) x| C_2  [EXACT]")
    print("=" * 78)
    for n, kmax in ((4, 12), (5, 12)):
        t0 = time.time()
        order, orb, rem = orbit_counts(n, kmax)
        N = math.factorial(n)
        assert all(r == 0 for r in rem)
        print(f"n={n}: |G| = {order} (faithful on S_{n}); {time.time() - t0:.1f}s")
        for k in range(1, kmax + 1):
            print(f"   k={k:2d}: C({N},{k}) = {math.comb(N, k):>20d}   orbits = {orb[k]}")


# ----------------------------------------------------------------------------- PART facets
def vertex_coords(n):
    """EXACT: vertices v_sigma of Q_n in 42 (resp. 14) independent indicator coordinates."""
    chi, keys = indicator_rows(n, [0, 1, 2])
    # greedy independent rows
    chosen, cur = [], []
    r = 0
    for i, row in enumerate(chi):
        if zrank(cur + [row]) > r:
            cur.append(row); chosen.append(i); r += 1
    V = [list(c) for c in zip(*cur)]  # n! x r
    return V, [keys[i] for i in chosen]


def facet_sampling(n, nrays, seed):
    from scipy.optimize import linprog
    V, _ = vertex_coords(n)
    Vf = np.array(V, float)
    N, D = Vf.shape
    b = Vf.mean(0)
    U_, S_, Wt = np.linalg.svd(Vf - b)
    B = Wt[: D - 1]  # affine directions (D-1)
    U = (Vf - b) @ B.T
    rng = np.random.default_rng(seed)
    G = group_action(n)
    found, classes = {}, {}
    nondeg = 0
    scale = np.linalg.norm(U, axis=1).mean()
    for ray in range(nrays):
        if ray % 2 == 0:  # uniform random direction
            d = rng.standard_normal(D - 1)
        else:  # towards the centroid of a random small vertex set, slightly perturbed
            T = rng.choice(N, size=int(rng.integers(2, 9)), replace=False)
            d = U[T].mean(0) + 0.02 * scale * rng.standard_normal(D - 1)
        A_eq = np.vstack([np.hstack([U.T, -d[:, None]]), np.hstack([np.ones(N), [0.0]])])
        b_eq = np.zeros(D); b_eq[-1] = 1.0
        c = np.zeros(N + 1); c[-1] = -1.0
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * N + [(0, None)], method="highs")
        if res.status != 0:
            continue
        y = res.eqlin.marginals
        vals = U @ y[:-1] + y[-1]
        Z = [i for i in range(N) if abs(vals[i]) < 1e-7]
        rows = [V[i] for i in Z]
        if zrank(rows) != D - 1:
            continue
        nondeg += 1
        ns, nul = flint.fmpz_mat(rows).nullspace()
        cvec = [int(ns[i, 0]) for i in range(D)]
        val = [sum(a * c_ for a, c_ in zip(V[s], cvec)) for s in range(N)]
        if min(val) < 0:
            val = [-x for x in val]
        assert min(val) == 0 and all(x >= 0 for x in val)
        g = math.gcd(*val)
        val = tuple(x // g for x in val)
        Zt = tuple(i for i in range(N) if val[i] == 0)
        found[Zt] = val
        if Zt not in classes:
            imgs = np.sort(G[:, list(Zt)], axis=1)
            canon = tuple(min(map(tuple, imgs)))
            classes[Zt] = canon
    # group facets by class
    byclass = {}
    for Zt, canon in classes.items():
        byclass.setdefault(canon, []).append(Zt)
    return V, found, byclass, nondeg


def is_minor_type(val, n):
    chi, keys = indicator_rows(n, range(1, n))
    # minor type: zero set = {sigma(I) != J}, i.e. val proportional to chi_IJ
    sv = set(val)
    if sv == {0, 1}:
        for row in chi:
            if tuple(row) == tuple(val):
                return True
    return False


def part_facets():
    print("=" * 78)
    print("PART facets: dimension of Q_n, minor facets, facet discovery by random rays")
    print("=" * 78)
    for n, nrays in ((4, 400), (5, 600)):
        t0 = time.time()
        V, keys = vertex_coords(n)
        D = len(V[0])
        print(f"n={n}: {len(V)} vertices, linear rank {zrank(V)} => dim Q_{n} = {zrank(V) - 1}  [EXACT]")
        # minor vectors as facets
        chi, ks = indicator_rows(n, [1, 2])
        seen, nfac, nnot = set(), 0, 0
        for row, (I, J) in zip(chi, ks):
            if tuple(row) in seen:
                continue
            seen.add(tuple(row))
            Z = [V[s] for s in range(len(V)) if row[s] == 0]
            if zrank(Z) == D - 1:
                nfac += 1
            else:
                nnot += 1
        print(f"      distinct minor vectors n_IJ (|I|=1,2, complements identified): {len(seen)}; "
              f"facets: {nfac}; not facets: {nnot}  [EXACT]")
        V, found, byclass, nondeg = facet_sampling(n, nrays, seed=11 + n)
        print(f"      random-ray facet discovery: {nrays} rays, {nondeg} nondegenerate exits, "
              f"{len(found)} distinct facets in {len(byclass)} G-classes (lower bounds)  "
              f"[LP NUMERIC proposes, facet property EXACT; {time.time() - t0:.1f}s]")
        for canon, Zs in sorted(byclass.items(), key=lambda kv: -len(kv[0])):
            val = found[Zs[0]]
            prof = sorted(Counter(val).items())
            print(f"        class: |zero set| = {len(Zs[0])}, value profile {prof}, "
                  f"minor type: {is_minor_type(val, n)}, distinct facets hit: {len(Zs)}")
    return None


def nonminor_facet_functionals(n=5, nrays=600, seed=16):
    V, found, byclass, _ = facet_sampling(n, nrays, seed)
    out = []
    for canon, Zs in sorted(byclass.items(), key=lambda kv: -len(kv[0])):
        val = found[Zs[0]]
        if not is_minor_type(val, n):
            out.append(np.array(val, float))
    return out


# ----------------------------------------------------------------------------- PART uo
class MinorFunctional:
    """l(U) = w0 + sum w1[i,j] |U_ij|^2 + sum W2[R,C] |det U[R,C]|^2 (R,C 2-subsets).
    k >= 3 minors are folded onto complements (|det U[I,J]| = |det U[I^c,J^c]| for
    unitary U), so every functional on r_n(U) can be written this way for n <= 5."""

    def __init__(self, n, w0, w1, W2):
        self.n, self.w0, self.w1, self.W2 = n, w0, np.array(w1, float), np.array(W2, float)
        pairs = subsets(n, 2)
        self.pairs = pairs
        self.r1 = np.array([p[0] for p in pairs]); self.r2 = np.array([p[1] for p in pairs])

    def __call__(self, U):
        a, b = self.r1, self.r2
        M2 = (U[a][:, a] * U[b][:, b] - U[a][:, b] * U[b][:, a])
        return self.w0 + np.sum(self.w1 * np.abs(U) ** 2) + np.sum(self.W2 * np.abs(M2) ** 2)

    def values(self):
        P = perms(self.n)
        out = []
        for s in P:
            v = self.w0 + sum(self.w1[i, s[i]] for i in range(self.n))
            for ri, R in enumerate(self.pairs):
                img = tuple(sorted((s[R[0]], s[R[1]])))
                ci = self.pairs.index(img)
                v += self.W2[ri, ci]
            out.append(v)
        return np.array(out)


def functional_from_terms(n, terms):
    """terms: list of (coef, I, J) with 0-based index tuples of any size."""
    pairs = subsets(n, 2)
    w0, w1, W2 = 0.0, np.zeros((n, n)), np.zeros((len(pairs), len(pairs)))
    full = set(range(n))
    for c, I, J in terms:
        I, J = tuple(sorted(I)), tuple(sorted(J))
        k = len(I)
        if k in (0, n):
            w0 += c; continue
        if k > n - k:
            I, J = tuple(sorted(full - set(I))), tuple(sorted(full - set(J)))
            k = n - k
        if k == 1:
            w1[I[0], J[0]] += c
        elif k == 2:
            W2[pairs.index(I), pairs.index(J)] += c
        else:
            raise ValueError("k > 2 after folding")
    return MinorFunctional(n, w0, w1, W2)


def h_terms(pi, eps, rows=(0, 1, 2, 3), cols=(0, 1, 2, 3), extra=()):
    """h_{pi,eps} of Sect. 4.1 on the 4x4 block rows x cols; extra = index pair appended to
    every I and J (complement-style embedding)."""
    R = lambda S: tuple(rows[i] for i in S)
    Cc = lambda S: tuple(cols[pi[i]] for i in S)
    t = [(0.5, R((i,)), Cc((i,))) for i in range(4)]
    d = [(+0.5 * eps, R((0, 1)), Cc((0, 2))), (+0.5 * eps, R((0, 1)), Cc((1, 3))),
         (-0.5 * eps, R((0, 2)), Cc((0, 1))), (-0.5 * eps, R((0, 2)), Cc((2, 3)))]
    out = []
    for c, I, J in t + d:
        if extra:
            I, J = I + (extra[0],), J + (extra[1],)
        out.append((c, I, J))
    return out


def haar_unitary(n, rng):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def haar_so(n, rng):
    Z = rng.standard_normal((n, n))
    Q, R = np.linalg.qr(Z)
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def exp_herm(p, n):
    H = np.zeros((n, n), complex)
    H[np.diag_indices(n)] = p[:n]
    iu = np.triu_indices(n, 1)
    m = len(iu[0])
    H[iu] = p[n:n + m] + 1j * p[n + m:]
    H = H + np.triu(H, 1).conj().T
    lam, W = np.linalg.eigh(H)
    return (W * np.exp(1j * lam)) @ W.conj().T


def exp_skew(p, n):
    S = np.zeros((n, n))
    iu = np.triu_indices(n, 1)
    S[iu] = p
    S = S - S.T
    lam, W = np.linalg.eigh(1j * S)
    return ((W * np.exp(-1j * lam)) @ W.conj().T).real


def minimize_over(f, n, group, starts, rng):
    """NUMERIC multistart L-BFGS-B in the chart X0 * exp(.) around Haar starts; the
    best point is re-centred and polished. Returns (best value, best matrix)."""
    from scipy.optimize import minimize
    npar = n * n if group == "U" else n * (n - 1) // 2
    ex = (lambda p: exp_herm(p, n)) if group == "U" else (lambda p: exp_skew(p, n))
    best, bestX = np.inf, None
    for _ in range(starts):
        X0 = haar_unitary(n, rng) if group == "U" else haar_so(n, rng)
        g = lambda p: float(f(X0 @ ex(p)))
        r = minimize(g, np.zeros(npar), method="L-BFGS-B",
                     options={"ftol": 1e-15, "gtol": 1e-10, "maxiter": 2000})
        if r.fun < best:
            best, bestX = r.fun, X0 @ ex(r.x)
    for _ in range(3):  # polish
        X0 = bestX
        g = lambda p: float(f(X0 @ ex(p)))
        r = minimize(g, np.zeros(npar), method="L-BFGS-B",
                     options={"ftol": 1e-16, "gtol": 1e-12, "maxiter": 4000})
        if r.fun <= best:
            best, bestX = r.fun, X0 @ ex(r.x)
    return best, bestX


def functional_from_values(n, v):
    """Minor functional (k<=2 coordinates) with value vector v at permutation matrices."""
    chi, keys = indicator_rows(n, [0, 1, 2])
    A = np.array(chi, float).T  # n! x m
    w, *_ = np.linalg.lstsq(A, v, rcond=None)
    assert np.max(np.abs(A @ w - v)) < 1e-9, "value vector not in the indicator span"
    terms = [(w[i], I, J) for i, (I, J) in enumerate(keys)]
    return functional_from_terms(n, terms)


def leibniz(v, O, P, S):
    n = O.shape[0]
    prod = np.ones(len(P))
    for i in range(n):
        prod = prod * O[i, P[:, i]]
    return float(np.sum(S * v * prod))


def compare(name, n, f, starts_O, starts_U, seed, vals=None):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    mO, XO = minimize_over(f, n, "O", starts_O, rng)
    mU, XU = minimize_over(f, n, "U", starts_U, rng)
    mperm = float(np.min(vals)) if vals is not None else float("nan")
    flag = "U BELOW O" if mU < mO - 1e-7 else "equal within 1e-7"
    print(f"   {name}: min over perms = {mperm:.10f}; min SO({n}) = {mO:.10f}; "
          f"min U({n}) = {mU:.10f}; diff U-O = {mU - mO:.2e} -> {flag}  "
          f"[NUMERIC, {starts_O}/{starts_U} starts, seed {seed}, {time.time() - t0:.0f}s]")
    return mO, mU, XO, XU


def r_vector(U, n):
    """all k = 1..n-1 minors squared (no folding), for rank tests."""
    out = []
    for k in range(1, n):
        for I in subsets(n, k):
            for J in subsets(n, k):
                out.append(abs(np.linalg.det(U[np.ix_(I, J)])) ** 2)
    return np.array(out)


def membership(n, U, rng, nsample=3000, iters=40, starts=12):
    """NUMERIC column generation: is r(U) in conv{r(O)}? Uses the k<=2 coordinates."""
    from scipy.optimize import linprog
    chi, keys = indicator_rows(n, [1, 2])

    def y(X):
        a = np.abs(X) ** 2
        out = [a[I[0], J[0]] if len(I) == 1 else abs(X[I[0], J[0]] * X[I[1], J[1]] - X[I[0], J[1]] * X[I[1], J[0]]) ** 2
               for (I, J) in keys]
        return np.array(out)
    Ys = [y(haar_so(n, rng)) for _ in range(nsample)]
    for s in perms(n):
        Pm = np.zeros((n, n)); Pm[range(n), s] = 1
        Ys.append(y(Pm))
    Ys = np.array(Ys)
    mean = Ys.mean(0)
    _, sv, Wt = np.linalg.svd(Ys - mean, full_matrices=False)
    r = int(np.sum(sv > 1e-8 * sv[0]))
    B = Wt[:r]
    proj = lambda Y: (Y - mean) @ B.T
    yU = proj(y(U))
    Z = proj(Ys)
    hist = []
    for it in range(iters):
        # separation LP: max t s.t. c.(z_k - yU) >= t, -1 <= c <= 1
        A_ub = np.hstack([-(Z - yU), np.ones((len(Z), 1))])
        res = linprog(np.r_[np.zeros(r), -1.0], A_ub=A_ub, b_ub=np.zeros(len(Z)),
                      bounds=[(-1, 1)] * r + [(None, None)], method="highs")
        c, t = res.x[:r], res.x[-1]
        if t < 1e-9:
            hist.append((it, t, None))
            return True, r, it, t, hist
        # minimise c.y(O) over SO(n)
        w = B.T @ c
        f = lambda X: float(np.dot(w, y(X)))
        found = []
        for _ in range(starts):
            m, X = minimize_over(f, n, "O", 1, rng)
            found.append((m, X))
        gap = min(m for m, _ in found) - float(np.dot(w, y(U)))
        hist.append((it, t, gap))
        if gap > 1e-7:
            return False, r, it, t, hist + [("separating functional (k<=2 coords)", w)]
        Z = np.vstack([Z] + [proj(y(X))[None, :] for _, X in found])
    return None, r, iters, t, hist


def part_uo(which="1"):
    print("=" * 78)
    print(f"PART uo{which}: is r_n(U) in conv r_n(O(n))?  support functions and membership LP")
    print("=" * 78)
    if which == "2":
        return part_uo_realizable()
    if which == "3":
        return part_uo_membership()
    # span checks
    for n in (3, 4, 5):
        rng = np.random.default_rng(500 + n)
        RU = np.array([r_vector(haar_unitary(n, rng), n) for _ in range(300)])
        RO = np.array([r_vector(haar_so(n, rng), n) for _ in range(300)])
        rk = lambda M: int(np.sum(np.linalg.svd(np.c_[M, np.ones(len(M))], compute_uv=False) > 1e-9))
        print(f"n={n}: numerical rank of [1, r_n(.)] over 300 Haar samples: U({n}) {rk(RU)}, SO({n}) {rk(RO)}, "
              f"joint {rk(np.r_[RU, RO])}; Catalan {catalan(n)}  [NUMERIC, tol 1e-9 relative]")
    # Leibniz cross-check of folding and (4.2) for n = 5
    n = 5
    P5 = np.array(perms(5)); S5 = np.array([sgn(p) for p in perms(5)], float)
    rng = np.random.default_rng(77)
    errs = []
    for _ in range(20):
        terms = [(rng.standard_normal(), I, J) for k in (1, 2, 3, 4) for I in subsets(5, k) for J in subsets(5, k)
                 if rng.random() < 0.05]
        f = functional_from_terms(5, terms)
        O = haar_so(5, rng)
        direct = sum(c * abs(np.linalg.det(O[np.ix_(I, J)])) ** 2 for c, I, J in terms)
        errs.append(max(abs(f(O) - direct), abs(leibniz(f.values(), O, P5, S5) - direct)))
        U = haar_unitary(5, rng)
        directU = sum(c * abs(np.linalg.det(U[np.ix_(I, J)])) ** 2 for c, I, J in terms)
        errs.append(abs(f(U) - directU))
    print(f"n=5: folding k>=3 onto complements and Leibniz sum = sum w det(O[I,J])^2 on SO(5): "
          f"max error {max(errs):.1e} over 20 random functionals  [NUMERIC]")
    # functional vanishing at all permutation matrices vanishes on U(5)
    chi, keys = indicator_rows(5, [0, 1, 2])
    A = np.array(chi, float).T
    _, sv, Wt = np.linalg.svd(A)
    ker = Wt[np.sum(sv > 1e-9):]
    worst = 0.0
    for _ in range(10):
        w = ker.T @ rng.standard_normal(ker.shape[0])
        f = functional_from_terms(5, [(w[i], I, J) for i, (I, J) in enumerate(keys)])
        worst = max(worst, abs(f(haar_unitary(5, rng))))
    print(f"n=5: {ker.shape[0]}-dim space of k<=2 functionals vanishing at all 120 permutation matrices; "
          f"max |value| at random U(5): {worst:.1e}  [NUMERIC]")

    print("-- controls n = 3, 4")
    for j in range(2):
        rng = np.random.default_rng(900 + j)
        f = MinorFunctional(3, 0.0, rng.standard_normal((3, 3)), np.zeros((3, 3)))
        compare(f"n=3 random functional #{j}", 3, f, 20, 20, 900 + j, f.values())
    for eps in (-1, +1):
        f = functional_from_terms(4, h_terms((0, 1, 2, 3), eps))
        vals = f.values()
        compare(f"n=4 h_(id,{eps:+d})", 4, f, 40, 40, 950 + eps, vals)
    # cross-check h against the shared library's value vector (P * h_vec)
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
        import moc4lib as L
        libv = [float(x) for x in L.values(L.h_vec(tuple(range(4)), -1))]
        mine = functional_from_terms(4, h_terms((0, 1, 2, 3), -1)).values()
        print(f"   n=4 h_(id,-1): value vector equals moc4lib values(h_vec(id,-1)): {np.allclose(libv, mine)}")
    except Exception as e:
        print(f"   moc4lib cross-check skipped: {e!r}")

    print("-- n = 5 test functionals")
    tests = []
    P5l = perms(5)

    def lifted(name, terms, R, C):
        """terms + c (1 - r_{R,C}), c the least value making all permutation values >= 0;
        the correction only changes permutations with sigma(R) != C."""
        v = functional_from_terms(5, terms).values()
        off = [v[i] for i, s in enumerate(P5l) if s[R] != C]
        c = max(0.0, -min(off))
        f = functional_from_terms(5, terms + [(c, (), ()), (-c, (R,), (C,))])
        return (f"{name} + {c:g}(1 - r_{R + 1}{C + 1})", f)
    tests.append(lifted("E1 h_(id,-1) on block rows1-4 x cols1-4", h_terms((0, 1, 2, 3), -1), 4, 4))
    tests.append(lifted("E1' h_(id,+1) on block rows1-4 x cols1-4", h_terms((0, 1, 2, 3), +1), 4, 4))
    tests.append(lifted("E1'' h_(2143,-1) on block rows1-4 x cols1-4", h_terms((1, 0, 3, 2), -1), 4, 4))
    tests.append(lifted("E2 h_(id,-1) with 5 appended to every I,J", h_terms((0, 1, 2, 3), -1, extra=(4, 4)), 4, 4))
    tests.append(lifted("E3 h_(id,-1) on rows1-4 x cols2-5", h_terms((0, 1, 2, 3), -1, cols=(1, 2, 3, 4)), 4, 0))
    tests.append(lifted("E4 E1 + E2", h_terms((0, 1, 2, 3), -1) + h_terms((0, 1, 2, 3), -1, extra=(4, 4)), 4, 4))
    tests.append(lifted("E5 E1 + E3", h_terms((0, 1, 2, 3), -1) + h_terms((0, 1, 2, 3), -1, cols=(1, 2, 3, 4)), 4, 4))
    V5, _ = vertex_coords(5)
    results = []
    for k, (name, f) in enumerate(tests):
        v = f.values()
        vi = np.rint(2 * v).astype(int)
        assert np.allclose(2 * v, vi)
        Z = [V5[i] for i in range(120) if vi[i] == 0]
        print(f"   {name}: values x2 profile {sorted(Counter(vi.tolist()).items())}; zero set {len(Z)} "
              f"vertices spanning a face of dimension {zrank(Z) - 1} of Q_5 (facet iff 40)  [EXACT]")
        mO, mU, XO, XU = compare(name, 5, f, 40, 40, 1000 + k, f.values())
        results.append((name, mO, mU))


def part_uo_realizable():
    print("-- n = 5 positive combinations of G-images of the non-minor facet E4 = E1 + E2")
    G = group_action(5)
    v4 = functional_from_terms(5, h_terms((0, 1, 2, 3), -1) + h_terms((0, 1, 2, 3), -1, extra=(4, 4))).values()
    rng = np.random.default_rng(4321)
    for j in range(4):
        v = np.zeros(120)
        for g in rng.choice(len(G), size=3, replace=False):
            w = np.zeros(120); w[G[g]] = v4
            v += rng.uniform(0.2, 1.0) * w
        f = functional_from_values(5, v)
        compare(f"facet-image combo #{j}", 5, f, 40, 40, 1500 + j, f.values())
    print("-- n = 5 realizable functionals Re(theta det(D_a + U D_b U*))")
    rng = np.random.default_rng(1234)
    for j in range(6):
        a = rng.standard_normal(5) + 1j * rng.standard_normal(5)
        b = rng.standard_normal(5) + 1j * rng.standard_normal(5)
        th = np.exp(1j * rng.uniform(0, 2 * np.pi))
        Da, Db = np.diag(a), np.diag(b)
        f = lambda U, Da=Da, Db=Db, th=th: float((th * np.linalg.det(Da + U @ Db @ U.conj().T)).real)
        zs = np.array([(th * np.prod([a[i] + b[s[i]] for i in range(5)])).real for s in perms(5)])
        compare(f"realizable #{j}", 5, f, 30, 30, 2000 + j, zs)


def combo_functional(j_target):
    G = group_action(5)
    v4 = functional_from_terms(5, h_terms((0, 1, 2, 3), -1) + h_terms((0, 1, 2, 3), -1, extra=(4, 4))).values()
    rng = np.random.default_rng(4321)
    for j in range(j_target + 1):
        v, gs = np.zeros(120), []
        for g in rng.choice(len(G), size=3, replace=False):
            w = np.zeros(120); w[G[g]] = v4
            c = rng.uniform(0.2, 1.0)
            v += c * w; gs.append((int(g), c))
    return v, gs


def part_uo_combo3():
    """Follow-up on facet-image combo #3 (U(5) minimum below SO(5) minimum)."""
    from scipy.linalg import polar
    v, gs = combo_functional(3)
    f = functional_from_values(5, v)
    P5 = perms(5)
    print(f"combo #3 = sum of c_g * (g . E4) over (g index, c_g) = {[(g, round(c, 6)) for g, c in gs]}")
    print(f"   value vector at permutations: min {v.min():.3g}, #zeros {int(np.sum(np.abs(v) < 1e-12))}")
    rng = np.random.default_rng(77003)
    t0 = time.time()
    mO, XO = minimize_over(f, 5, "O", 600, rng)
    print(f"   SO(5): 600 Haar multistarts: min = {mO:.3e}  [{time.time() - t0:.0f}s]")
    t0 = time.time()
    mU, XU = minimize_over(f, 5, "U", 120, rng)
    XU, _ = polar(XU)  # re-unitarize
    val = f(XU)
    print(f"   U(5): 120 Haar multistarts: min = {mU:.6e}; value at re-unitarized minimizer "
          f"{val:.6e}; ||U*U^H - I|| = {np.linalg.norm(XU @ XU.conj().T - np.eye(5)):.1e}  [{time.time() - t0:.0f}s]")
    A = np.abs(XU)
    dist = min(np.linalg.norm(A - np.eye(5)[list(s)]) for s in P5)
    print(f"   |U*| entries (rounded):\n{np.round(A, 4)}")
    print(f"   distance of |U*| to nearest permutation pattern: {dist:.3f}")
    # orthogonal starts derived from U*: real part, polar-projected, and its sign variants
    best = np.inf
    from scipy.optimize import minimize
    for M in (XU.real, XU.imag, (XU * np.exp(-1j * np.angle(XU[0, 0]))).real):
        try:
            O0, _ = polar(M)
        except Exception:
            continue
        if np.linalg.det(O0) < 0:
            O0[:, 0] *= -1
        g = lambda p: float(f(O0 @ exp_skew(p, 5)))
        r = minimize(g, np.zeros(10), method="L-BFGS-B", options={"ftol": 1e-16, "gtol": 1e-12})
        best = min(best, r.fun)
    print(f"   SO(5) local searches started from projections of U*: min = {best:.3e}")
    # the minimum over the phase torus at fixed |U*|-pattern is not separately explored
    print(f"   verdict: min_U - min_O <= {val - mO:.3e}  [NUMERIC]")


def part_uo_membership():
    print("-- membership of r_n(U) for Haar U in conv r_n(O(n)) (column generation)")
    for n, cnt in ((3, 1), (4, 2), (5, 3)):
        for j in range(cnt):
            rng = np.random.default_rng(3000 + 10 * n + j)
            U = haar_unitary(n, rng)
            t0 = time.time()
            inside, r, it, t, hist = membership(n, U, rng)
            last = [h for h in hist if isinstance(h[0], int)]
            print(f"   n={n} Haar U #{j}: affine coords {r}; verdict "
                  f"{'INSIDE (sample hull contains r(U))' if inside else ('OUTSIDE: separating functional found' if inside is False else 'UNDECIDED')} "
                  f"after {it} rounds; last LP margin {t:.2e}; gaps {[None if h[2] is None else round(h[2], 9) for h in last[-4:]]}  "
                  f"[NUMERIC, {time.time() - t0:.0f}s]")


if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "alg"
    {"alg": part_alg, "burnside": part_burnside, "facets": part_facets,
     "uo1": lambda: part_uo("1"), "uo2": lambda: part_uo("2"), "uo3": lambda: part_uo("3"), "uo4": part_uo_combo3}[part]()
