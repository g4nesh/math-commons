"""Audit replay for independent_moc4_T23_g17_uniform_20261022 (seven-contact set T_24 + {17}).

Own code only; uses ../lib/moc4lib.py.  The only supplement input is the JSON data file
independent_moc4_T23_g17_full_radial_20261022.json, loaded with json.load as DATA (the
degree-8 polynomial G); its role is checked independently below.  Location from MOC4_SUPP
(the supplement root) or argv[1].  Prints one line per check; exits nonzero if an EXACT
check fails.  Options: --numeric (floating-point tests), --bernstein (the 42-table replay).
"""
import os
import sys
import json
import random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info="", level="EXACT"):
    print(f"{'PASS' if ok else 'FAIL'} {level} {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


args = [a for a in sys.argv[1:] if not a.startswith("--")]
SUPP = args[0] if args else os.environ.get("MOC4_SUPP")

T = frozenset({2, 3, 4, 6, 7, 8})
T7 = T | {17}
gens = L.generators(T7)
check("generators_T24+17", gens == [0, 2, 3, 5, 53, 54, 67], str(gens))
check("rank_7_independent", L.rank([L.F_VECS[j] for j in gens]) == 7)
check("only_f54_nonminor", [j for j in gens if not L.F_IS_MINOR[j]] == [54])
hits = [k for k, v in L.AUX7.items() if L.find_g(v, T7) is not None]
print(f"INFO T_24+17 is in the orbit of Table-5 seven-set(s): {hits}")

# coefficient formulas: a=g1 (f67), b=g22 (f0), e=g12 (f3), p=g13 (f2), k=g16 (f53), j=g10 (f5), d=g9 (f54)
CO = {67: 1, 0: 22, 3: 12, 2: 13, 53: 16, 5: 10, 54: 9}
ok = True
okd = True
for jj in gens:
    vv = L.F_VALUES[jj]
    rec = [sum(vv[CO[k]] * L.F_VALUES[k][s] for k in CO) for s in range(24)]
    ok = ok and rec == list(vv)
    okd = okd and vv[9] == vv[5] - vv[16] == vv[0] - vv[1] - vv[16]
check("canonical_coefficients_are_gaps", ok)
check("d=g9=g5-g16=g0-g1-g16", okd)

a, b, e, p, k, d, j, h = sp.symbols("a b e p k d j h")
coeff = {67: a, 0: b, 3: e, 2: p, 53: k, 54: d, 5: j}
y = [sum(coeff[f] * L.F_VECS[f][r] for f in coeff) + h * L.KAPPA[r] for r in range(14)]
Cy = sp.expand(L.cubic_C(y))
U = (a + d) * (b + d) + d * e
W = U + k * (b + d + e)
Q = (e - p) * j - (e + d) * (k + d)
K4 = a * b * (e - d) - d**2 * (a + b + e + d)
check("identity_(2)_C(g+h kappa)=-dW+hQ", sp.expand(Cy - (-d * W + h * Q)) == 0)

# Hessian of C (linear in m)
ms = sp.symbols("m0:14")
Csym = sp.expand(L.cubic_C(ms))
HESS = [[sp.diff(Csym, ms[r], ms[s]) for s in range(14)] for r in range(14)]
HESS_f = [[sp.lambdify(ms, HESS[r][s]) for s in range(14)] for r in range(14)]


def hess_at(vals):
    return [[Fraction(HESS_f[r][s](*vals)) if HESS[r][s].free_symbols else Fraction(int(HESS[r][s]))
             for s in range(14)] for r in range(14)]


# lemma: det Hess C(Re x) = 0 and C(Re x) = 0 on X (Gaussian-rational matching vectors)
rng = random.Random(17)
okh = True
for _ in range(8):
    aa = [L.QI(rng.randint(-4, 4), rng.randint(-4, 4)) for _ in range(4)]
    bb = [L.QI(rng.randint(-4, 4), rng.randint(-4, 4)) for _ in range(4)]
    m = L.matching_vector(aa, bb, L.QI(rng.randint(-2, 2), 1))
    yr = [L.re_part(t) for t in m]
    okh = okh and L.cubic_C(yr) == 0 and L.det(hess_at(yr)) == 0
check("C(Re x)=0_and_detHess(Re x)=0_on_8_gaussian_rational_points", okh)

# G from the supplement JSON (DATA ONLY); variable order there is (a,b,c,d,e,f,v,z)
G = None
if SUPP:
    fn = os.path.join(SUPP, "proofs", "manin_transfer_20261021", "independent_moc4_T23_g17_full_radial_20261022.json")
    with open(fn) as fh:
        D = json.load(fh)
    # identify the variables from the stored cubic/Q terms (independently re-derived above)
    jv = sp.symbols("ja jb jc jd je jf jv jz")

    def poly_of(terms):
        return sp.Add(*[sp.Rational(c) * sp.Mul(*[v**ex for v, ex in zip(jv, mono)]) for mono, c in terms])
    MAP = {jv[0]: b, jv[1]: p, jv[2]: e, jv[3]: j, jv[4]: k, jv[5]: a, jv[6]: d, jv[7]: h}
    Cj = sp.expand(poly_of(D["C_terms"]).subs(MAP, simultaneous=True))
    Qj = sp.expand(poly_of(D["Q_terms"]).subs(MAP, simultaneous=True))
    check("json_variable_map_matches_-dW_and_Q", sp.expand(Cj + d * W) == 0 and sp.expand(Qj - Q) == 0)
    facs = D["radial_factors"]
    G = sp.expand(poly_of([t for f in facs if f["degree"] == 8 for t in f["terms"]]).subs(MAP, simultaneous=True))
    print(f"INFO G: {len(sp.Poly(G, a, b, e, p, k, d, j, h).terms())} terms, total degree "
          f"{sp.Poly(G, a, b, e, p, k, d, j, h).total_degree()}, involves h: {G.has(h)}, degree in j: {sp.degree(G, j)}")
    # (3) at random rational points, exactly
    Gf = sp.lambdify((a, b, e, p, k, d, j), G)
    ok3 = True
    rr = random.Random(3)
    for _ in range(25):
        A_, B_, E_, P_, K_, D_, J_ = [Fraction(rr.randint(1, 40), rr.randint(1, 9)) for _ in range(7)]
        Wv = (A_ + D_) * (B_ + D_) + D_ * E_ + K_ * (B_ + D_ + E_)
        Qv = (E_ - P_) * J_ - (E_ + D_) * (K_ + D_)
        if Qv == 0:
            continue
        hv = D_ * Wv / Qv
        sub = {a: A_, b: B_, e: E_, p: P_, k: K_, d: D_, j: J_, h: hv}
        yv = [Fraction(sp.Rational(sp.sympify(t).subs(sub))) if not isinstance(t, (int,)) else Fraction(t) for t in y]
        dH = L.det(hess_at(yv))
        Gv = Fraction(sp.Rational(G.subs(sub)))
        ok3 = ok3 and (Qv**4 * dH == -4 * D_**2 * Qv**2 * Gv**2)
    check("identity_(3)_Q^4 detHess = -4 d^2 Q^2 G^2 at 25 random rational points", ok3,
          level="EXACT-AT-POINTS")
    k0 = sp.expand(G.subs(j, 0))
    k0claim = d * (b + d) * W * ((a + d + k) * (-K4) + (a + d) * (p * W + a * e * k))
    check("identity_(7)_k0", sp.expand(k0 - k0claim) == 0)
    Gpe0 = sp.Poly(sp.expand(G.subs({p: 0, e: 0})), a, b, k, d, j)
    check("p=e=0: all coefficients of G positive", all(c > 0 for c in Gpe0.coeffs()),
          f"d^8 coeff={Gpe0.coeff_monomial(d**8)}")
else:
    print("SKIP G-dependent checks: no supplement path (set MOC4_SUPP)")

# capacity transfer witness: some g maps the facets (59,53,55,67,54,2) to (5,3,67,0,54,2)
src = [59, 53, 55, 67, 54, 2]
dst = [5, 3, 67, 0, 54, 2]
T25 = L.CONTACT_SETS[25]
print(f"INFO generators of paper T_25 = {L.generators(T25)}")
found = []
for g in L.G_ELEMENTS:
    t = L.g_index_map(g)
    ok = True
    for s_, d_ in zip(src, dst):
        img = [None] * 24
        for sg in range(24):
            img[t[sg]] = L.F_VALUES[s_][sg]
        if tuple(img) != L.F_VALUES[d_]:
            ok = False
            break
    if ok:
        found.append(g)
check("capacity_transfer_group_element_exists", bool(found), f"{len(found)} elements, e.g. {found[:1]}")
if found:
    g = found[0]
    print(f"INFO that element maps T_25 to {sorted(L.g_act_set(g, T25))}; T_24+17 = {sorted(T7)}")

print(f"SUMMARY exact checks failed: {len(FAIL)} {FAIL}")

if "--bernstein" in sys.argv and G is not None:
    import flint
    names = ("A", "T", "V", "U", "k", "j")
    ctx = flint.fmpq_mpoly_ctx.get(names, "lex")
    A_, T_, V_, U_, k_, j_ = ctx.gens()
    one = ctx.from_dict({(0,) * 6: 1})
    Gd1 = sp.Poly(sp.expand(G.subs(d, 1)), a, b, e, p, k, j)
    deg = {v: sp.degree(Gd1.as_expr(), v) for v in (a, b, e, p)}

    def chart(name):
        """returns dict var -> (num, den) as fmpq_mpoly, and J* = (N, D)."""
        if name == "I-":
            av = (A_, one); bv = (T_, A_ * (1 + T_)); ev = (V_, one); pv = (V_ + U_, one)
            Js = (pv[1], pv[0])
        elif name in ("I+", "II"):
            av = (A_, one); bv = (1 + T_, A_)
            # L = (a+1)(b+1)/(ab-1) = (A+1)(1+T+A)/(A T)
            Ln, Ld = (A_ + 1) * (1 + T_ + A_), A_ * T_
            ev = (Ln * V_, Ld * (1 + V_))
            if name == "I+":
                pv = (ev[0] + U_ * ev[1], ev[1]); Js = (pv[1], pv[0])
            else:
                pv = (ev[0] * U_, ev[1] * (1 + U_))
                # 1/e - 1/L = Ld/(Ln V)
                Js = (Ld, Ln * V_)
        elif name in ("IIIa-", "IIIa+"):
            bv = (T_, 1 + T_); av = (T_ * V_, (1 + T_) * (1 + V_)); ev = (A_, one)
            if name == "IIIa-":
                pv = (A_ * U_, (1 + A_) * (1 + U_)); Js = (1 + A_, A_)
            else:
                pv = (A_ * (1 + U_) + A_ * A_ * U_, (1 + A_) * (1 + U_)); Js = (pv[1], pv[0])
        elif name in ("IIIb+", "IIIb-"):
            bv = (1 + T_, one); ev = (A_, one)
            if name == "IIIb+":
                av = (V_, (1 + T_) * (1 + V_))
                pv = (A_ * (1 + T_) * (1 + U_) + A_ * A_ * U_, (A_ + 1 + T_) * (1 + U_)); Js = (pv[1], pv[0])
            else:
                # w = V/(1+V), a = w^2/b
                av = (V_ * V_, (1 + V_) * (1 + V_) * (1 + T_))
                pv = (A_ * (1 + T_) * U_, (A_ + 1 + T_) * (1 + U_))
                # J* = 1/e + b(1-w)^2/(b-w)^2 = 1/A + (1+T)/((1+T)(1+V)-V)^2
                bw = (1 + T_) * (1 + V_) - V_
                Js = (bw * bw + A_ * (1 + T_), A_ * bw * bw)
        return {a: av, b: bv, e: ev, p: pv}, Js

    def subst(poly_sp, cmap, swap):
        """poly in a,b,e,p,k (d=1) -> numerator after multiplying by prod den^deg."""
        out = ctx.from_dict({})
        P = sp.Poly(poly_sp, a, b, e, p, k)
        dg = {v: max(P.degree(v), 0) for v in (a, b, e, p)}
        vmap = dict(cmap)
        if swap:
            vmap = {a: cmap[b], b: cmap[a], e: cmap[e], p: cmap[p]}
            dg = {a: dg[a], b: dg[b], e: dg[e], p: dg[p]}
        for mono, c in P.terms():
            t = ctx.from_dict({(0,) * 6: int(c)}) if c.q == 1 else ctx.from_dict({(0,) * 6: flint.fmpq(int(c.p), int(c.q))})
            for v, ex in zip((a, b, e, p), mono[:4]):
                nu, de = vmap[v]
                t = t * nu**ex * de**(dg[v] - ex)
            t = t * k_**mono[4]
            out = out + t
        return out, dg

    def all_pos(poly):
        cs = poly.coeffs()
        return len(cs) > 0 and all(c > 0 for c in cs), len(cs)

    k0s = Gd1.as_expr().coeff(j, 0)
    k1s = Gd1.as_expr().coeff(j, 1)
    k2s = Gd1.as_expr().coeff(j, 2)
    # use a common degree vector for k0,k1,k2 so that numerators share the denominator
    Pall = sp.Poly(Gd1.as_expr(), a, b, e, p, k, j)
    for reg in ("I-", "I+", "II", "IIIa-", "IIIa+", "IIIb+", "IIIb-"):
        cmap, (N, Dn) = chart(reg)
        for swap in (False, True):
            # common denominators: multiply every coefficient by the same prod den^maxdeg
            dgmax = {v: Pall.degree(v) for v in (a, b, e, p)}

            def sub_common(expr):
                P = sp.Poly(expr, a, b, e, p, k)
                out = ctx.from_dict({})
                vmap = cmap if not swap else {a: cmap[b], b: cmap[a], e: cmap[e], p: cmap[p]}
                for mono, c in P.terms():
                    t = ctx.from_dict({(0,) * 6: flint.fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q))})
                    for v, ex in zip((a, b, e, p), mono[:4]):
                        nu, de = vmap[v]
                        t = t * nu**ex * de**(dgmax[v] - ex)
                    out = out + t * k_**mono[4]
                return out
            # G(b,a) on the chart == G with a <- chart(b), b <- chart(a) (done inside sub_common)
            K0 = sub_common(k0s)
            K1 = sub_common(k1s)
            K2 = sub_common(k2s)
            mid = 2 * Dn * K0 + N * K1
            end = Dn * Dn * K0 + N * Dn * K1 + N * N * K2
            pm, nm = all_pos(mid)
            pe, ne = all_pos(end)
            tag = "G(b,a)" if swap else "G(a,b)"
            check(f"bernstein_{reg}_{tag}_middle_coeffs_positive", pm, f"terms={nm}")
            check(f"bernstein_{reg}_{tag}_endpoint_coeffs_positive", pe, f"terms={ne}")
    print(f"SUMMARY after Bernstein: exact checks failed: {len(FAIL)} {FAIL}")

if "--numeric" in sys.argv:
    import numpy as np
    Pn = L.P_numpy()
    V = {f: Pn @ np.array([float(t) for t in L.F_VECS[f]]) for f in coeff}

    def minphi(cf, seed):
        vals = sum(cf[f] * V[f] for f in cf)
        cvec = np.linalg.lstsq(Pn, vals, rcond=None)[0]
        return L.min_phi_over_SO4(cvec, restarts=8, seed=seed)[0]
    rng = np.random.default_rng(5)
    worst = {"SOC": np.inf, "oldH": np.inf, "cap5": np.inf, "cap6": np.inf}
    for it in range(12):
        av, bv, ev, pv = rng.uniform(0.05, 3, 4)
        base = {67: av, 0: bv, 3: ev, 2: pv, 53: 0.0, 54: 1.0}
        worst["SOC"] = min(worst["SOC"], minphi({**base, 5: 1 / pv}, it))
        worst["oldH"] = min(worst["oldH"], minphi({**base, 5: 1 + 1 / ev}, it))
        if av * bv > 1:
            Lv = (av + 1) * (bv + 1) / (av * bv - 1)
            worst["cap5"] = min(worst["cap5"], minphi({**base, 2: 0.0, 5: max(0.0, 1 / ev - 1 / Lv)}, it))
        else:
            A2, B2 = sorted([av, bv])
            wv = np.sqrt(A2 * B2)
            worst["cap6"] = min(worst["cap6"], minphi({**base, 67: A2, 0: B2, 2: 0.0,
                                                      5: 1 / ev + B2 * (1 - wv)**2 / (B2 - wv)**2}, it))
    print(f"NUMERIC capacities at their thresholds (d=1, k=0, 12 draws, seed 5, 8 restarts): "
          + ", ".join(f"{kk} worst minPhi={vv:.3e}" for kk, vv in worst.items()))

sys.exit(1 if FAIL else 0)
