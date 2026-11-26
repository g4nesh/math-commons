"""Calibration runs for the toolkit. Usage: python calibrate.py STAGE [out_dir]
Stages: n2, n3, n4, psi, fakes, verifier, hp. Each writes calib_<stage>.json."""
import json
import os
import subprocess
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from moc_core import margin, margin_points, z_all, det_d, perms, hull_indices          # noqa: E402
from moc_param import (cayley, skew_from_real, random_unitary, random_spectra,  # noqa: E402
                       perm_matrix)
from moc_opt import run_restarts                                          # noqa: E402
import moc_exact as ex                                                    # noqa: E402
from moc_hp import hp_margin                                              # noqa: E402

O_STAR = [[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]


def psi_np(V):
    """Zheng's Example 1.3 functional (as transcribed in the repo's literature lane)."""
    def m(I, J):
        return abs(np.linalg.det(V[np.ix_([i - 1 for i in I], [j - 1 for j in J])])) ** 2
    return 0.5 * (sum(abs(V[i, i]) ** 2 for i in range(4)) - m([1, 2], [1, 3]) - m([1, 2], [2, 4])
                  + m([1, 3], [1, 2]) + m([1, 3], [3, 4]))


def psi_exact_ostar():
    O = [[Fr(v, 7) for v in r] for r in O_STAR]
    def m(I, J):
        r = [i - 1 for i in I]
        c = [j - 1 for j in J]
        return (O[r[0]][c[0]] * O[r[1]][c[1]] - O[r[0]][c[1]] * O[r[1]][c[0]]) ** 2
    return Fr(1, 2) * (sum(O[i][i] ** 2 for i in range(4)) - m([1, 2], [1, 3]) - m([1, 2], [2, 4])
                       + m([1, 3], [1, 2]) + m([1, 3], [3, 4]))


def stage_restarts(n, R, seed, out, **opts):
    return run_restarts(n, R, seed, out_json=out, time_limit=540, **opts)["summary"]


def stage_psi(out, restarts=200, seed=2026):
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    vals = []
    for r in range(restarts):
        U0 = random_unitary(4, rng)
        for _ in range(3):   # recentred stages: U = U0 cayley(K(X))
            f = lambda x: psi_np(U0 @ cayley(skew_from_real(x.reshape(4, 4))))
            res = minimize(f, np.zeros(16), method="BFGS", options=dict(gtol=1e-11))
            U0 = U0 @ cayley(skew_from_real(res.x.reshape(4, 4)))
        v = psi_np(U0)
        vals.append(v)
        if v < best[0]:
            best = (v, U0)
    Pv = sorted(set(round(psi_np(perm_matrix(s)), 12) for s in perms(4)))
    O = np.array(O_STAR) / 7.0
    res = dict(psi_Ostar_exact=str(psi_exact_ostar()), psi_Ostar_float=psi_np(O),
               psi_on_permutation_matrices=Pv, restarts=restarts, seed=seed,
               min_psi_over_U4=best[0], abs_U_sq_at_min=(np.abs(best[1]) ** 2).round(6).tolist(),
               n_within_1e6_of_min=int(np.sum(np.array(vals) < best[0] + 1e-6)))
    json.dump(res, open(out, "w"), indent=1)
    return res


def stage_fakes(out, scratch):
    rng = np.random.default_rng(5)
    res = {}
    # (1) non-unitary "U": find one that pushes d outside H
    for t in range(200):
        a, b = random_spectra(4, rng)
        Ufake = random_unitary(4, rng) * 1.3
        r = margin(a, b, Ufake)
        if r["m"] > 1e-3:
            break
    res["nonunitary_float"] = dict(m=r["m"], unitarity_defect=r["unitarity_defect"], ok_flag=r["ok"])
    Ue = [[ex.rat(v, 10 ** 6) for v in row] for row in Ufake]
    cert = ex.build_certificate(a, b, None, r["theta"], U_exact=Ue,
                                out_json=os.path.join(scratch, "fake_nonunitary_cert.json"))
    chk = subprocess.run([sys.executable, "-I", os.path.join(HERE, "check_cert.py"),
                          os.path.join(scratch, "fake_nonunitary_cert.json")], capture_output=True, text=True)
    res["nonunitary_exact"] = dict(builder_status=cert["status"], builder_gap_positive=Fr(cert["gap"]) > 0,
                                   checker_exit=chk.returncode, checker_msg=chk.stdout.strip())
    # (2) wrong hull: honest inside config, drop the active permutations
    a, b = random_spectra(4, rng)
    U = random_unitary(4, rng)
    r = margin(a, b, U)
    z = z_all(a, b)
    d = det_d(a, b, U)
    keep = np.ones(len(z), bool)
    proj = (np.exp(-1j * r["theta"]) * z).real
    order = np.argsort(-proj)
    keep[order[:6]] = False
    try:
        margin_points(d, z[keep], 4)
        guard = "NOT RAISED"
    except ValueError as e:
        guard = "raised: " + str(e)
    rw = margin_points(d, z[keep], 4, allow_partial=True)
    res["wrong_hull_float"] = dict(true_m=r["m"], m_with_6_perms_dropped=rw["m"], completeness_guard=guard)
    cert = ex.build_certificate(a, b, U, rw["theta"], out_json=None)
    # forge: certificate whose gap was computed against the partial list only
    Lz = [Fr(cert["direction"][0]) * z_.x + Fr(cert["direction"][1]) * z_.y
          for k, (s, z_) in enumerate(ex.z_exact([ex.rat(v, 10 ** 12) for v in a],
                                                  [ex.rat(v, 10 ** 12) for v in b])) if keep[k]]
    dd = ex.G(Fr(cert["d"][0]), Fr(cert["d"][1]))
    cert["gap"] = str(Fr(cert["direction"][0]) * dd.x + Fr(cert["direction"][1]) * dd.y - max(Lz))
    p = os.path.join(scratch, "fake_wronghull_cert.json")
    json.dump(cert, open(p, "w"))
    chk = subprocess.run([sys.executable, "-I", os.path.join(HERE, "check_cert.py"), p],
                         capture_output=True, text=True)
    res["wrong_hull_exact"] = dict(forged_gap=cert["gap"][:40], forged_gap_positive=Fr(cert["gap"]) > 0,
                                   checker_exit=chk.returncode, checker_msg=chk.stdout.strip()[:200])
    # (3) trivial contact noise: permutation matrix perturbed by 1e-9 (unitary)
    a, b = random_spectra(5, rng)
    zz = z_all(a, b)
    vtx = perms(5)[hull_indices(zz)[0]]          # a permutation whose z is a hull vertex
    P = perm_matrix(vtx).astype(complex)
    E = 1e-9 * (rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5)))
    Up = P @ cayley(0.5 * (E - E.conj().T))
    rf = margin(a, b, Up)
    rh = hp_margin(a, b, Up, dps=60)
    res["near_permutation_noise"] = dict(float_m=rf["m"], hp_m=rh["m_str"])
    json.dump(res, open(out, "w"), indent=1, default=str)
    return res


def stage_verifier(out, scratch, repo_verify):
    rng = np.random.default_rng(44)
    a, b = random_spectra(4, rng)
    U = random_unitary(4, rng)
    r = margin(a, b, U)
    K, Ue, R = ex.exact_unitary_from_float(U, 12)    # small denominators: verifier caps 64 digits
    ae = [ex.rat(v, 12) for v in a]
    be = [ex.rat(v, 12) for v in b]
    d = ex.d_exact(ae, be, Ue)
    memb = ex.membership(d, ex.z_exact(ae, be))
    w = ex.witness_json(ae, be, Ue, "toolkit calibration n=4 inside instance")
    wp = os.path.join(scratch, "verifier_witness_n4.json")
    rp = os.path.join(scratch, "verifier_result_n4.json")
    json.dump(w, open(wp, "w"))
    chk = subprocess.run([sys.executable, "-I", repo_verify, wp, "--output", rp], capture_output=True, text=True)
    rv = json.load(open(rp)) if os.path.exists(rp) else {}
    Uf = np.array([[complex(v) for v in row] for row in Ue])
    rr = margin([complex(v) for v in ae], [complex(v) for v in be], Uf)
    res = dict(float_m_of_rationalized=rr["m"], toolkit_d=d.enc(), toolkit_membership=memb, unitary=ex.is_unitary(Ue),
               verifier_exit=chk.returncode, verifier_stderr=chk.stderr[-300:],
               verifier_d=rv.get("determinant"), verifier_membership=rv.get("membership"),
               d_agree=(rv.get("determinant") == d.enc()))
    json.dump(res, open(out, "w"), indent=1)
    return res


def load_config(rec):
    a = np.array([complex(*v) for v in rec["a"]])
    b = np.array([complex(*v) for v in rec["b"]])
    U = np.array(rec["U_re"]) + 1j * np.array(rec["U_im"])
    return a, b, U


def stage_hp(out, outd, top=5):
    res = {}
    for st in ("n2", "n3", "n4", "n3c", "n4c"):
        p = os.path.join(outd, "calib_%s.json" % st)
        if not os.path.exists(p):
            continue
        rows = []
        for rec in json.load(open(p))["best"][:top]:
            a, b, U = load_config(rec)
            h = hp_margin(a, b, U, dps=60)
            rows.append(dict(restart=rec["restart"], float_m=rec["m"], hp_m=h["m_str"],
                             dist_genperm=rec["dist_genperm"], dist_blockdiag=rec["dist_blockdiag"]))
        res[st] = rows
    json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    stage = sys.argv[1]
    outd = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "calibration")
    scratch = sys.argv[3] if len(sys.argv) > 3 else outd
    os.makedirs(outd, exist_ok=True)
    o = os.path.join(outd, "calib_%s.json" % stage)
    t = time.time()
    if stage == "n2":
        r = stage_restarts(2, 2000, 101, o)
    elif stage == "n3":
        r = stage_restarts(3, 2000, 102, o)
    elif stage == "n3c":
        r = stage_restarts(3, 500, 103, o, dperm=0.15, dblock=0.15)
    elif stage == "n4":
        r = stage_restarts(4, 1000, 104, o)
    elif stage == "n4c":
        r = stage_restarts(4, 300, 105, o, dperm=0.15, dblock=0.15)
    elif stage == "psi":
        r = stage_psi(o)
    elif stage == "fakes":
        r = stage_fakes(o, scratch)
    elif stage == "hp":
        r = stage_hp(o, outd)
    elif stage == "verifier":
        r = stage_verifier(o, scratch, sys.argv[4])
    print(json.dumps(r, indent=1, default=str)[:3000])
    print("seconds %.1f" % (time.time() - t))
