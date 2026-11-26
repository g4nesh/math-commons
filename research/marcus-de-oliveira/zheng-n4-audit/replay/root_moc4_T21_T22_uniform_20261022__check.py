"""Audit replay for root_moc4_T21_T22_uniform_20261022
(paper cases T_22 = {2,3,4,8,10,12} and T_23 = {2,3,5,9,11,12}).

Own code only; uses ../lib/moc4lib.py.  Nothing from the supplement is read or executed.
Prints one line per check; exits nonzero if an EXACT check fails.
"""
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


def vec(cnt):
    out = [0] * 14
    for j, c in cnt.items():
        for r in range(14):
            out[r] = out[r] + c * L.F_VECS[j][r]
    return out


def eqv(lhs, rhs):
    return vec(Counter(lhs)) == vec(Counter(rhs))


T22 = frozenset({2, 3, 4, 8, 10, 12})
T23 = frozenset({2, 3, 5, 9, 11, 12})
check("contact_sets_equal_paper_T22_T23", L.CONTACT_SETS[22] == T22 and L.CONTACT_SETS[23] == T23)
S21 = [0, 1, 2, 6, 10, 53, 54, 55, 56, 67, 69, 70, 72]
S22 = [0, 1, 2, 7, 9, 12, 39, 40, 41, 67, 68, 72]
check("generators_T22", L.generators(T22) == S21, str(L.generators(T22)))
check("generators_T23", L.generators(T23) == S22, str(L.generators(T23)))
NM21 = [j for j in S21 if not L.F_IS_MINOR[j]]
NM22 = [j for j in S22 if not L.F_IS_MINOR[j]]
check("nonminor_T22", NM21 == [54, 56, 69, 70], str(NM21))
check("nonminor_T23", NM22 == [40, 68], str(NM22))
print(f"INFO rank generators T22 = {L.rank([L.F_VECS[j] for j in S21])}, T23 = {L.rank([L.F_VECS[j] for j in S22])}")

PAIRS = [([54, 56], [0, 2, 10, 53, 55]), ([54, 69], [2, 6, 10, 53, 67]),
         ([54, 70], [0, 6, 55, 67]), ([56, 69], [1, 2, 2, 10, 53, 72]),
         ([56, 70], [0, 1, 2, 55, 72]), ([69, 70], [1, 2, 6, 67, 72]),
         ([40, 68], [0, 7, 41, 67])]
for l, r in PAIRS:
    check(f"pair_{l}={r}", eqv(l, r) and all(L.F_IS_MINOR[j] for j in r))
# every pair of nonminor types in each family has a replacement
check("all_nonminor_pairs_listed",
      {tuple(sorted(l)) for l, _ in PAIRS} ==
      {(a, b) for fam in (NM21, NM22) for a in fam for b in fam if a < b})

ABS = [([54, 1], [0, 6, 53]), ([54, 72], [10, 55, 67]), ([56, 6], [1, 2, 10, 55]),
       ([56, 67], [0, 2, 53, 72]), ([69, 0], [1, 2, 10, 67]), ([69, 55], [2, 6, 53, 72]),
       ([70, 10], [0, 6, 72]), ([70, 53], [1, 55, 67]), ([40, 1], [0, 7, 39]),
       ([40, 9], [2, 7, 41]), ([40, 72], [12, 41, 67]), ([68, 2], [0, 9, 67]),
       ([68, 12], [0, 7, 72]), ([68, 39], [1, 41, 67])]
res = {}
for l, r in ABS:
    check(f"absorb_{l}={r}", eqv(l, r) and all(L.F_IS_MINOR[j] for j in r))
    res.setdefault(l[0], set()).add(l[1])
for N, R in res.items():
    rhs = set().union(*[set(r) for l, r in ABS if l[0] == N])
    check(f"absorb_no_regeneration_{N}", not (rhs & R), f"resources={sorted(R)}")

NF = [(21, 54, -1, (0, 2, 6, 10, 53, 55, 67, 54), "P24"),
      (21, 56, 1, (1, 2, 10, 0, 55, 53, 72, 56), "P24"),
      (21, 69, 1, (1, 2, 10, 6, 67, 53, 72, 69), "P24"),
      (21, 70, -1, (0, 2, 6, 1, 72, 55, 67, 70), "P24"),
      (22, 40, 1, (12, 7, 2, 0, 67, 41, 40, 39), "P33"),
      (22, 68, 1, (1, 0, 9, 7, 41, 67, 68, 72), "P33")]
for fam, N, sg, lam, Pn in NF:
    S = S21 if fam == 21 else S22
    NM = NM21 if fam == 21 else NM22
    left = set(S) - (set(NM) - {N}) - res[N]
    check(f"normal_form_support_{fam}_{N}", left == set(lam), f"{sorted(left)}")

r = sp.symbols("r0:8", nonnegative=True)
d = r[7]
P24 = (d**3 + d**2 * (r[0] + r[2] + r[3] + r[4] + r[5] + r[6])
       + d * (r[5]*r[2] + r[5]*r[3] + r[5]*r[6] + r[5]*r[0] + r[2]*r[4] + r[2]*r[6] + r[2]*r[0]
              + r[3]*r[4] + r[3]*r[6] + r[4]*r[0] + r[6]*r[0])
       + r[5]*r[2]*r[6] + r[5]*r[2]*r[0] + r[5]*r[3]*r[6] + r[5]*r[6]*r[0] + r[2]*r[4]*r[0]
       + r[2]*r[6]*r[0])
N_ = (r[2] + r[6]) * (r[0] + r[1] + r[6]) + r[3] * (r[1] + r[6])
P33 = (r[5] + r[6] + r[7]) * N_ + r[4] * ((r[0] + r[1] + r[3] + r[6]) * (r[5] + r[6]) + r[1] * r[3])
POLYS = {"P24": sp.expand(P24), "P33": sp.expand(P33)}
for nm, poly in POLYS.items():
    pp = sp.Poly(poly, *r)
    check(f"{nm}_all_coeffs_positive", all(c > 0 for c in pp.coeffs()))
dvar = {"P24": r[7], "P33": r[6]}
for nm, poly in POLYS.items():
    pp = sp.Poly(poly, *r)
    mon = [0] * 8
    mon[r.index(dvar[nm])] = 3
    check(f"{nm}_has_d^3_coeff_1", pp.coeff_monomial(tuple(mon)) == 1)

for fam, N, sg, lam, Pn in NF:
    m = [0] * 14
    for i, j in enumerate(lam):
        for k in range(14):
            m[k] = m[k] + r[i] * L.F_VECS[j][k]
    C = sp.expand(L.cubic_C(m))
    check(f"cubic_identity_{fam}_{N}_sigma={sg}_{Pn}", sp.expand(sg * C - POLYS[Pn]) == 0,
          "" if sp.expand(sg * C - POLYS[Pn]) == 0 else f"ratio={sp.simplify(C / POLYS[Pn])}")
    check(f"nonminor_position_{fam}_{N}", lam[7 if Pn == 'P24' else 6] == N)

# rectangles and swaps (c = -b convention: F_sigma = prod (A_i - C_sigma(i)))
A = sp.symbols("A0:4")
Cc = sp.symbols("C0:4")


def Fs(k):
    p = L.PERMS[k]
    return sp.Mul(*[A[i] - Cc[p[i]] for i in range(4)])


check("rectangle_F2F10=F4F8", sp.expand(Fs(2) * Fs(10) - Fs(4) * Fs(8)) == 0)
check("rectangle_F3F11=F5F9", sp.expand(Fs(3) * Fs(11) - Fs(5) * Fs(9)) == 0)
for i, j in [(2, 4), (2, 8), (3, 5), (3, 9)]:
    f = sp.factor_list(sp.expand(Fs(i) - Fs(j)))
    lin = [fa for fa, e in f[1]]
    ok = len(lin) == 4 and all(sp.Poly(fa, *A, *Cc).total_degree() == 1 and len(fa.free_symbols) == 2 for fa in lin)
    check(f"swap_F{i}-F{j}_is_product_of_4_label_differences", ok, str(sp.factor(Fs(i) - Fs(j))))

print(f"SUMMARY exact checks failed: {len(FAIL)} {FAIL}")

if "--numeric" in sys.argv:
    import numpy as np
    for T, ps, tag in [(T22, 2, "T22"), (T23, 2, "T23")]:
        worst, found, nd = np.inf, 0, 0
        for seed in range(12):
            pt = L.random_point_on_contact_stratum(T, seed=seed, tries=10, pi_star=ps)
            if pt is None:
                continue
            found += 1
            labels = list(pt["a"]) + [-x for x in pt["b"]]
            sep = min(abs(labels[i] - labels[j]) for i in range(8) for j in range(i + 1, 8))
            nd += sep > 1e-3
            c = np.linalg.lstsq(L.P_numpy(), pt["gaps"], rcond=None)[0]
            val, _ = L.min_phi_over_SO4(c, restarts=8, seed=seed)
            worst = min(worst, val)
            print(f"NUMERIC {tag} seed={seed} minLabelSep={sep:.2e} contact={pt['contact']} minPhi={val:.3e}")
        print(f"NUMERIC {tag} points={found}/12 eight-distinct(sep>1e-3)={nd} worst minPhi={worst:.3e}")

sys.exit(1 if FAIL else 0)
