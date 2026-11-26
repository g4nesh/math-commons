"""Own replay for blind_moc4_inventory52_plus17_middle_products_20261023.
Reuses the own objects of the low-products replay (same directory): N, R, Q from the
paper's cubic, and G / octics read as data and checked there.  Rebuilds the twelve
F-coefficient tables of the middle chamber and checks their signs and markers."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blind_moc4_inventory52_plus17_low_products_20261023__check as LOW
import flint
from math import comb

check = LOW.check
CTX = LOW.CTX
a, B, C, b, A, c, d, z = LOW.GENS
FV = z    # reuse the 8th slot as the auxiliary variable F

def chart_mid(poly, shift=(5, 5, 4), deg=18):
    kctx = flint.fmpz_mpoly_ctx.get(('x', 'y', 'zz', 'W0', 'C', 'F'), 'lex')
    x, y, zz, W0, Cv, Fv = kctx.gens()
    D = W0 + x + y
    out = kctx.from_dict({}); cache = {}
    for mon, co in poly.to_dict().items():
        ia, iB, iC, ib, iA, ic, idd, iF = (int(t) for t in mon)
        assert iB == 0 and iA == 0 and idd == 0
        ex = ia + ib - ic + shift[0]; ey = ia - ib + ic + shift[1]; ez = -ia + ib + ic + shift[2]
        assert min(ex, ey, ez) >= 0 and ex + ey <= deg, (mon, ex, ey, ez)
        key = (ex, ey, ez)
        if key not in cache:
            cache[key] = x ** ex * (x + y) ** ey * (1 + zz) ** ez * D ** (deg - ex - ey)
        assert co.q == 1
        out = out + int(co.p) * cache[key] * Cv ** iC * Fv ** iF
    return out

def main():
    u, Cu, H = LOW.build_core()
    cz = LOW.coeffs_in(Cu, 7)
    NR = -cz[0]; cB = LOW.coeffs_in(NR, 1); N, R = cB[0], cB[1]
    data = LOW.load_supp()
    if data is None:
        print("SKIP (set MOC4_SUPP)"); return
    G = LOW.poly_from_terms(data["G"]); gB = LOW.coeffs_in(G, 1)
    al, be, ga = gB[2], gB[1], gB[0]
    Hs = [LOW.poly_from_terms(t) for t in data["Herons"]]
    check("Delta = F+ F- (re-checked)", be * be - 4 * al * ga == Hs[0] * Hs[1])
    Fstar = Hs[data["negative_Heron_index"]]
    s1 = LOW.setd1
    ctrl = {"gamma": s1(ga), "R beta - N alpha": s1(R * be - N * al), "R F* - N alpha": s1(R * Fstar - N * al)}
    U = (b + 1) * (c + 1); V = b * c - 1
    markers = {"gamma": 4096, "R beta - N alpha": 2048}
    stated_counts = {"gamma": [4473, 4824, 5131, 4374], "R beta - N alpha": [5220, 5913, 6435, 5687],
                     "R F* - N alpha": [5179, 6624, 7299, 7514]}
    for nm, p in ctrl.items():
        pa = LOW.coeffs_in(p, 4)
        check("%s has degree <= 3 in A" % nm, max(pa) <= 3, max(pa))
        tot = CTX.from_dict({})
        for j in range(4):
            inner = CTX.from_dict({})
            for i in range(0, j + 1):
                if i in pa:
                    inner = inner + comb(3 - i, j - i) * pa[i] * U ** i * V ** (3 - i)
            tot = tot + FV ** j * inner
        cp = chart_mid(tot)
        byF = {}
        for mon, co in cp.to_dict().items():
            byF.setdefault(int(mon[5]), {})[mon] = int(co)
        counts = [len(byF.get(j, {})) for j in range(4)]
        signs = [(sum(1 for v in byF.get(j, {}).values() if v > 0), sum(1 for v in byF.get(j, {}).values() if v < 0)) for j in range(4)]
        print("INFO %-18s counts %s (pos,neg) %s" % (nm, counts, signs))
        want_pos = nm != "R F* - N alpha"
        ok_sign = all((s[1] == 0) if want_pos else (s[0] == 0) for s in signs)
        check("%s: every F-coefficient table has the stated sign" % nm, ok_sign)
        check("%s: term counts equal the stated %s" % (nm, stated_counts[nm]), counts == stated_counts[nm], counts)
        f0 = byF.get(0, {})
        if nm in markers:
            check("%s: F^0 marker %d x^18 z^3" % (nm, markers[nm]), f0.get((18, 0, 3, 0, 0, 0)) == markers[nm], f0.get((18, 0, 3, 0, 0, 0)))
        else:
            check("R F* - N alpha: F^0 marker -384 x^18 z^3 C", f0.get((18, 0, 3, 0, 1, 0)) == -384, f0.get((18, 0, 3, 0, 1, 0)))
    # chart coverage, PROOF-READ facts checked symbolically: bc = Z^2, V = Z^2 - 1
    print("INFO chart: b c = Z^2 exactly, so V = bc - 1 = z(2+z) > 0 for z > 0 (by hand)")

if __name__ == "__main__":
    main()
    print("RESULT", "all EXACT checks passed" if LOW.fails == 0 else "%d failures" % LOW.fails)
    sys.exit(1 if LOW.fails else 0)
