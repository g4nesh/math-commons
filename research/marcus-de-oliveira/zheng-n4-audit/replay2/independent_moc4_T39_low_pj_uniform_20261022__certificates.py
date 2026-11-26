"""Certificate replay (auditor cert-T39) for independent_moc4_T39_low_pj_uniform_20261022.

Imports the shared CORE (own code) from the high_p replay in this directory, which rebuilds
the Hessian identity, G, the octics, P*, pseudo-division, trace-norm identity and bordered
minors from the paper's definitions.  Then replays the four low-(p,j) coefficient
certificates under p = P/(1+P), j = J/(1+J).  Supplement JSON is read as DATA only.

usage: MOC4_SUPP=<supplement root> python -I -B this_file.py
"""
import importlib.util
import os
import sys
import time

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hp", os.path.join(HERE, "independent_moc4_T39_high_p_uniform_20261022__certificates.py"))
hp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hp)
check, substitute, all_negative, all_nonneg = hp.check, hp.substitute, hp.all_negative, hp.all_nonneg
q, p, j, c, e, b, a = hp.q, hp.p, hp.j, hp.c, hp.e, hp.b, hp.a


def main():
    R = hp.core(full=True)
    A_, B_, N_, Nj, Q, k0, k1, k2, Fa, Fb = (R[k] for k in
                                            ["A", "B", "N", "Nj", "Q", "gamma", "beta", "alpha", "Fa", "Fb"])
    # sign guards on the region 0 <= p, j < 1
    check("N - ea has nonnegative coefficients and constant 1 (so N_j = (N - ea) + (1-j)ea > 0)",
          all_nonneg(N_ - e * a) and (N_ - e * a).to_dict().get((0,) * 8) == 1)
    rest = -Q - (1 - p) * c - (1 - j) * (e + a) - p * j - 1
    check("-Q = (1-p)c + (1-j)(e+a) + pj + 1 + (polynomial with nonnegative coefficients)",
          all_nonneg(rest))
    disp = (q + e + a + b + (1 - p)) * c + q * (e + a + 1) + b * (e + a + 1) + (1 - j) * (e + a) + 1 + p * j
    check("-Q equals the displayed expression of p.1", -Q == disp)
    octx = flint.fmpz_mpoly_ctx.get(["q", "P", "J", "e", "b", "a"])
    oq, oP, oJ, oe, ob, oa = octx.gens()
    one = octx.constant(1)
    images = [oq, oP, oJ, None, oe, ob, oa, None]
    dens = [None, one + oP, one + oJ, None, None, None, None, None]
    data = hp.load_supp_json("independent_moc4_T39_low_pj_uniform_exact_20261022.json")
    stored = {cc["name"]: cc["coefficients"] for cc in data["strict_certificates"]} if data else {}
    if data:
        print("INFO stored certificate names: " + ", ".join(stored))
    certs = [("k0", k0, 696),
             ("AB_beta_minus_Nj_alpha", A_ * B_ * k1 - Nj * k2, 1985),
             ("Fa_strict_negative_from_alpha", A_ * B_ * Fa + N_ * k2, 1457),
             ("Fb_strict_negative_from_alpha", A_ * B_ * Fb + N_ * k2, 1457)]
    outs = {}
    for name, pol, want in certs:
        out = substitute(pol, octx, images, dens)
        cst = int(out.to_dict().get((0,) * 6, 0))
        check("low-p,j certificate %s: %d terms (claimed %d), every coefficient < 0, constant %d < 0"
              % (name, len(out), want, cst),
              len(out) == want and all_negative(out) and cst < 0)
        outs[name] = out
    # data comparison (F+/F- labels may be swapped between my factor order and theirs)
    if data:
        for name in ["k0", "AB_beta_minus_Nj_alpha"]:
            hp.compare_with_data(name, outs[name], None, stored.get(name))
        hp.compare_factor_pair("F+-certificates", [outs["Fa_strict_negative_from_alpha"],
                                                   outs["Fb_strict_negative_from_alpha"]],
                               [stored.get("Fa_strict_negative_from_alpha"),
                                stored.get("Fb_strict_negative_from_alpha")])
        hp.compare_factor_pair("octics F+, F-", [Fa, Fb], [data.get("Fa_terms"), data.get("Fb_terms")])
        hp.compare_factor_pair("factors L, M of r1", [R["L"], R["M"]],
                               [data.get("U16_terms"), data.get("V16_terms")])
    print("INFO chart inverse: P = p/(1-p), J = j/(1-j) are >= 0 and finite exactly for 0 <= p, j < 1; "
          "every P, J >= 0 gives 0 <= p, j < 1; the multiplier (1+P)^deg_p (1+J)^deg_j is positive.")
    print("SUMMARY %d failures; elapsed %.1f s" % (len(hp.FAILS), time.monotonic() - hp.T0))
    return 1 if hp.FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
