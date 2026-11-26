"""High-precision check of lifted block-edge contacts (n = 5). Usage: python -I -B hp_lift.py DPS
For representative records (edge_s501.json best[0], and the two adv_s508.json block-edge results with
the largest lambda): rebuild U0 = (R_ij(phi) P_s) (+) 1 exactly unitary at DPS digits (float spectra
and c taken as exact binary rationals), compute d, all 120 z, the edge normal of the lifted edge,
the signed distance of d to that line, the largest projection of any other z beyond d, and
lambda_max / lambda_min of Re(e^{-i theta} W) (W = second variation in the off-block directions).
Writes out/hp_lift.json."""
import sys, os, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
from lift_lib import W_matrix_mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 50
mp.mp.dps = DPS
P5 = list(itertools.permutations(range(5)))


def build(a, b, s, i, j, c):
    a = [mp.mpc(*t) for t in a]; b = [mp.mpc(*t) for t in b]
    c = mp.mpf(c)
    U = mp.zeros(5, 5)
    for r in range(5):
        U[r, s[r]] = 1
    R = mp.eye(5)
    co, si = mp.sqrt(c), mp.sqrt(1 - c)
    R[i, i], R[i, j], R[j, i], R[j, j] = co, si, -si, co
    return a, b, R * U


def analyse(a, b, U, s, i, j):
    M = mp.diag(a) + U * mp.diag(b) * U.transpose_conj()
    d = mp.det(M)
    z = [mp.fprod([a[k] + b[p[k]] for k in range(5)]) for p in P5]
    s5 = tuple(s) + (4,)
    t = list(s5); t[i], t[j] = t[j], t[i]; t = tuple(t)
    zs, zt = z[P5.index(s5)], z[P5.index(t)]
    e = zt - zs
    nrm = -1j * e / abs(e)
    proj = [mp.re((w - zs) * mp.conj(nrm)) for w in z]
    if max(proj) > -min(proj):          # orient outward: all z must lie on the nonpositive side
        nrm = -nrm; proj = [-p for p in proj]
    diam = max(abs(z[p] - z[q]) for p in range(120) for q in range(p))
    dd = mp.re((d - zs) * mp.conj(nrm)) / diam
    others = sorted((proj[k] / diam, k) for k in range(120) if k not in (P5.index(s5), P5.index(t)))
    unitary_defect = mp.mnorm(U * U.transpose_conj() - mp.eye(5), 1)
    th = mp.arg(nrm)
    W = W_matrix_mp(a, b, U, dps=DPS)
    H = mp.matrix(8, 8)
    for p in range(8):
        for q in range(8):
            H[p, q] = mp.re(mp.exp(-1j * th) * W[p, q])
    ev = mp.eigsy(H, eigvals_only=True)
    ev = sorted([ev[k] / diam for k in range(8)])
    return dict(signed_dist_d_to_edge_line_over_diam=mp.nstr(dd, 8),
                max_other_projection_over_diam=mp.nstr(others[-1][0], 12), other_perm=list(P5[others[-1][1]]),
                n_other_within_1e_30=sum(1 for o in others if o[0] > -mp.mpf(10) ** -30),
                unitary_defect=mp.nstr(unitary_defect, 5), theta=mp.nstr(th, 20),
                lam_max_over_diam=mp.nstr(ev[-1], 25), lam_min_over_diam=mp.nstr(ev[0], 15),
                eigs_over_diam=[mp.nstr(v, 15) for v in ev])


def main():
    out = []
    e = json.load(open(os.path.join(HERE, "out", "edge_s501.json")))["best"][0]
    a, b, U = build(e["a"], e["b"], e["edge"]["s"] + [4], e["edge"]["i"], e["edge"]["j"], e["edge"]["c"])
    out.append(dict(src="edge_s501 best[0]", edge=e["edge"], **analyse(a, b, U, e["edge"]["s"], e["edge"]["i"], e["edge"]["j"])))
    adv = [r for r in json.load(open(os.path.join(HERE, "out", "adv_s508.json")))["results"]
           if r["kind"] == "block-edge" and r["final"] is not None and r["lam_final"] > -1e-3]
    adv.sort(key=lambda r: -r["lam_final"])
    for r in adv[:3]:
        a, b, U = build(r["a"], r["b"], r["edge"]["s"] + [4], r["edge"]["i"], r["edge"]["j"], r["c"])
        out.append(dict(src="adv_s508 block-edge (lam float %.3e)" % r["lam_final"], edge=r["edge"], c=r["c"],
                        **analyse(a, b, U, r["edge"]["s"], r["edge"]["i"], r["edge"]["j"])))
    for o in out:
        print(json.dumps({k: o[k] for k in ("src", "signed_dist_d_to_edge_line_over_diam", "max_other_projection_over_diam",
                                             "n_other_within_1e_30", "unitary_defect", "lam_max_over_diam", "lam_min_over_diam")}))
    json.dump(dict(dps=DPS, results=out), open(os.path.join(HERE, "out", "hp_lift.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
