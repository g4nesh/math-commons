"""Audit replay for root_moc4_T20_uniform_20261022 (paper case T_21 = {0,3,5,9,11,12}).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Run:  python -I -B root_moc4_T20_uniform_20261022__check.py
"""
import os
import sys
import random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


def vec_comb(coeffs):
    """sum c_j f_j as a 14-vector (sympy or Fraction coefficients)."""
    out = [0] * 14
    for j, c in coeffs.items():
        for r in range(14):
            out[r] = out[r] + c * L.F_VECS[j][r]
    return out


def gval(vec, s):
    """(P vec)_s for a 14-vector with arbitrary-ring entries."""
    return sum((L.P[s][r] * vec[r] for r in range(14)), 0)


T = frozenset({0, 3, 5, 9, 11, 12})

# 1. contact set = paper T_21
check("contact_set_equals_paper_T21", L.CONTACT_SETS[21] == T)

# 2. complete source list
gens = L.generators(T)
check("generators_T21", gens == [0, 1, 2, 7, 9, 12, 29, 39, 40, 41], str(gens))
nonminor = [j for j in gens if not L.F_IS_MINOR[j]]
check("only_f40_nonminor", nonminor == [40], str(nonminor))
rk = L.rank([L.F_VECS[j] for j in gens])
print(f"INFO rank of the 10 generators = {rk}")

# 3. detector delta = g8 - g2 - g6 - g16 on generators
dv = {j: L.F_VALUES[j][8] - L.F_VALUES[j][2] - L.F_VALUES[j][6] - L.F_VALUES[j][16] for j in gens}
want = {j: (1 if j == 40 else -1 if j in (1, 9) else 0) for j in gens}
check("detector_delta_values_on_generators", dv == want, str(dv))

# 4. replacement identities
def eqv(lhs, rhs):
    return vec_comb({j: lhs.count(j) for j in set(lhs)}) == vec_comb({j: rhs.count(j) for j in set(rhs)})


check("f40+f1=f0+f7+f39", eqv([40, 1], [0, 7, 39]))
check("f40+f9=f2+f7+f41", eqv([40, 9], [2, 7, 41]))
check("replacement_rhs_all_minor", all(L.F_IS_MINOR[j] for j in [0, 7, 39, 2, 41]))

# 5. product (rectangle) relation and swap factorizations, symbolic spectra
A = sp.symbols("A0:4")
Cc = sp.symbols("C0:4")


def Fs(k):
    p = L.PERMS[k]
    return sp.Mul(*[A[i] - Cc[p[i]] for i in range(4)])


check("rectangle_F3F11=F5F9", sp.expand(Fs(3) * Fs(11) - Fs(5) * Fs(9)) == 0)
f35 = sp.factor(Fs(3) - Fs(5))
f39 = sp.factor(Fs(3) - Fs(9))
t35 = (A[0] - Cc[0]) * (A[1] - A[2]) * (A[3] - Cc[1]) * (Cc[2] - Cc[3])
t39 = (A[0] - A[3]) * (A[1] - Cc[2]) * (A[2] - Cc[3]) * (Cc[0] - Cc[1])
check("F3-F5_factorization", sp.expand(Fs(3) - Fs(5) - t35) == 0, str(f35))
check("F3-F9_factorization", sp.expand(Fs(3) - Fs(9) - t39) == 0, str(f39))
check("F3-F5_up_to_sign", sp.expand(Fs(3) - Fs(5) + t35) == 0 or sp.expand(Fs(3) - Fs(5) - t35) == 0)

# 6. H = 0 argument: (H+iq3)(H+iq11) = (H+iq5)(H+iq9) with H != 0 forces q3 in {q5,q9}
H, q3, q5, q9 = sp.symbols("H q3 q5 q9", real=True)
q11 = q5 + q9 - q3   # from imaginary part / H
re_part = sp.expand(sp.re(sp.expand((H + sp.I * q3) * (H + sp.I * q11) - (H + sp.I * q5) * (H + sp.I * q9))))
check("real_part_forces_(q3-q5)(q3-q9)=0", sp.expand(re_part + (q3 - q5) * (q3 - q9)) == 0
      or sp.expand(re_part - (q3 - q5) * (q3 - q9)) == 0, str(sp.factor(re_part)))

# 7. lemma C(Re x) = 0 on X: exact Gaussian-rational matching vectors
rng = random.Random(20261125)
ok = True
for _ in range(12):
    a = [L.QI(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
    b = [L.QI(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
    al = L.QI(rng.randint(-3, 3), rng.randint(1, 3))
    m = L.matching_vector(a, b, al)
    y = [L.re_part(t) for t in m]
    ok = ok and L.cubic_C(y) == 0
check("C(Re x)=0_on_12_gaussian_rational_matching_vectors", ok)

# 8. cubic identity (2)
q, a_, b_, c_, e_, f_, d_, h_ = sp.symbols("q a b c e f d h", nonnegative=True)
m = vec_comb({0: q, 2: a_, 7: b_, 12: c_, 29: e_, 39: f_, 40: d_, 41: h_})
Cval = sp.expand(L.cubic_C(m))
M = (b_ + d_) * (q + a_ + d_) + c_ * (a_ + d_)
claim = (f_ + d_ + h_) * M + e_ * ((q + a_ + c_ + d_) * (f_ + d_) + a_ * c_)
ratio = sp.simplify(Cval / sp.expand(claim)) if Cval != 0 else 0
check("cubic_identity_(2)_up_to_nonzero_constant", ratio != 0 and ratio.is_number, f"C(m)/claim = {ratio}")

# also the vector m is the zero-height lift: (P m) on T vanishes, and P m = g.
check("lift_vanishes_on_T", all(sp.expand(gval(m, s)) == 0 for s in T))

# 9. exhaustiveness of the delta split: in any nonnegative representation, delta = w40 - w1 - w9;
# replacements preserve the vector; checked symbolically by the detector values above.

print(f"SUMMARY exact checks failed: {len(FAIL)} {FAIL}")

# NUMERIC falsification on the contact stratum (best effort)
if "--numeric" in sys.argv:
    import numpy as np
    worst = np.inf
    worstdelta = -np.inf
    found = 0
    for seed in range(40):
        pt = L.random_point_on_contact_stratum(T, seed=seed, tries=20, pi_star=0)
        if pt is None:
            continue
        found += 1
        g = pt["gaps"]
        # the real lift: c = Re(alpha*m) - Re(F_0) kappa ; check Re F_0 ~ 0 relative scale
        Fv = L.F_values_np(pt["a"], pt["b"], pt["alpha"])
        H0 = Fv[0].real / max(1.0, np.max(np.abs(Fv)))
        delta = g[8] - g[2] - g[6] - g[16]
        worstdelta = max(worstdelta, delta / max(1.0, np.max(np.abs(g))))
        c = np.linalg.lstsq(L.P_numpy(), g, rcond=None)[0]
        val, _ = L.min_phi_over_SO4(c, restarts=10, seed=seed)
        worst = min(worst, val)
        labels = list(pt["a"]) + [-x for x in pt["b"]]
        sep = min(abs(labels[i] - labels[j]) for i in range(8) for j in range(i + 1, 8))
        print(f"NUMERIC seed={seed} minLabelSep={sep:.2e} contact={pt['contact']} relReF0={H0:.2e} "
              f"rel_delta={delta / max(1.0, np.max(np.abs(g))):.3e} minPhi={val:.3e}")
    print(f"NUMERIC stratum points found={found}/40 max rel delta={worstdelta:.3e} worst minPhi={worst:.3e}")

sys.exit(1 if FAIL else 0)
