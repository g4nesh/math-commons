"""Lifted-contact analysis, n = 4 -> 5. Usage: python -I -B run_lift.py STAGE [args]
Stages: check, edge SEED R, vertex SEED, adv SEED R, direct SEED R
Writes JSON into out/ (next to this file)."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from lift_lib import (K_off, w_quad, w_fd, W_matrix, lift_curvature, normals_at)
from moc_core import z_all, det_d, margin, perms, hull_indices, dist_genperm, dist_blockdiag
from moc_param import random_unitary, random_spectra, embed, perm_matrix, canonical_spectra
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
T0 = time.time()
TLIM = 540.0


def cj(z):
    return [[float(np.real(v)), float(np.imag(v))] for v in np.atleast_1d(z)]


def save(name, obj):
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(obj, f, indent=1)


def stage_check():
    rng = np.random.default_rng(7)
    res = []
    for t in range(5):
        a, b = random_spectra(5, rng, "gaussian")
        U4 = random_unitary(4, rng)
        U0 = np.eye(5, dtype=complex); U0[:4, :4] = U4
        x = rng.normal(size=4) + 1j * rng.normal(size=4)
        wa = w_quad(a, b, U0, x)
        wf = w_fd(a, b, U0, x)
        # evenness: first derivative
        K = K_off(x, 5); h = 1e-6
        d1 = (det_d(a, b, U0 @ expm(h * K)) - det_d(a, b, U0 @ expm(-h * K))) / (2 * h)
        W = W_matrix(a, b, U0)
        v = np.concatenate([x.real, x.imag])
        res.append(dict(w_analytic=cj(wa), w_fd=cj(wf), rel_err=float(abs(wa - wf) / abs(wa)),
                        first_deriv=float(abs(d1)), quadform_err=float(abs(v @ W @ v - wa) / abs(wa))))
    # trivial lift U0 = I: w(x) = sum |x_i|^2 (z_(i,4) - z_id)
    a, b = random_spectra(5, rng, "gaussian")
    x = rng.normal(size=4) + 1j * rng.normal(size=4)
    P = [tuple(int(t) for t in p) for p in perms(5)]
    z = z_all(a, b)
    pred = 0
    for i in range(4):
        s = list(range(5)); s[i], s[4] = 4, i
        pred += abs(x[i]) ** 2 * (z[P.index(tuple(s))] - z[0])
    wa = w_quad(a, b, np.eye(5, dtype=complex), x)
    res.append(dict(identity_lift_err=float(abs(wa - pred) / abs(pred))))
    save("check.json", res)
    print(json.dumps(res, indent=1))


def fifth_pairs(rng, scale, k):
    """Grid/random sample of (a5, b5): mixture of small, comparable and large magnitudes."""
    out = []
    for j in range(k):
        r = scale * [0.3, 1.0, 3.0, 10.0][j % 4]
        a5 = r * (rng.normal() + 1j * rng.normal())
        b5 = scale * (rng.normal() + 1j * rng.normal()) * [0.1, 1.0][j % 2]
        out.append((a5, b5))
    return out


def block_edge_contacts(a, b):
    """n=4 edge contacts of block type: hull edges [z_s, z_t] with t = s o (i j)."""
    P = [tuple(int(t) for t in p) for p in perms(4)]
    z = z_all(a, b)
    idx = hull_indices(z)
    out = []
    for k in range(len(idx)):
        s, t = P[idx[k]], P[idx[(k + 1) % len(idx)]]
        diff = [i for i in range(4) if s[i] != t[i]]
        if len(diff) == 2:
            out.append((s, diff[0], diff[1]))
    return out


def U_block_edge(s, i, j, c):
    """U = R_ij(phi) P_s with cos^2 phi = c: d = c z_s + (1 - c) z_t."""
    Ps = perm_matrix(s).astype(complex)
    R = np.eye(4, dtype=complex)
    co, si = np.sqrt(c), np.sqrt(1 - c)
    R[i, i], R[i, j], R[j, i], R[j, j] = co, si, -si, co
    return R @ Ps


def analyze_lift(a4, b4, U4, a5, b5, tol=1e-9):
    a, b, U0 = embed(a4, b4, U4, a5, b5)
    kind, info, out, W = lift_curvature(a, b, U0, tol)
    if kind == "not-boundary":
        return None
    lam = max(o["lam_max"] for o in out)
    return dict(kind=kind, lam_max=lam, per_normal=[dict(theta=o["theta"], lam_max=o["lam_max"],
                lam_min=o["lam_min"]) for o in out], info={k: (v if not isinstance(v, tuple) else list(v))
                for k, v in info.items()}, a=cj(a), b=cj(b),
                U_re=U0.real.tolist(), U_im=U0.imag.tolist())


def stage_edge(seed, R):
    rng = np.random.default_rng(seed)
    recs, lams, nb, ntot = [], [], 0, 0
    for r in range(R):
        if time.time() - T0 > TLIM:
            break
        kind = ["gaussian", "circle", "mixed"][r % 3]
        a4, b4 = random_spectra(4, rng, kind)
        a4, b4 = canonical_spectra(a4, b4)
        for (s, i, j) in block_edge_contacts(a4, b4):
            for c in (0.2, 0.5, 0.8):
                U4 = U_block_edge(s, i, j, c)
                for (a5, b5) in fifth_pairs(rng, 1.0, 8):
                    ntot += 1
                    rec = analyze_lift(a4, b4, U4, a5, b5)
                    if rec is None:
                        continue
                    nb += 1
                    rec.update(restart=r, spectra_kind=kind, edge=dict(s=list(s), i=i, j=j, c=c))
                    lams.append(rec["lam_max"])
                    recs.append(rec)
                    recs.sort(key=lambda t: -t["lam_max"])
                    del recs[20:]
    lams = np.array(lams)
    summ = dict(stage="edge", seed=seed, restarts=r + 1, lifts_tried=ntot, lifts_on_boundary=nb,
                lam_max_overall=float(lams.max()) if len(lams) else None,
                n_lam_pos=int(np.sum(lams > 1e-12)), quantiles=np.quantile(lams, [0, .5, .9, .99, 1]).tolist() if len(lams) else [],
                seconds=time.time() - T0)
    save("edge_s%d.json" % seed, dict(summary=summ, best=recs))
    print(json.dumps(summ, indent=1))


def load_vertex_contacts():
    j = json.load(open(os.path.join(HERE, "..", "lib", "contact_n4.json")))
    cs = []
    for c in j["numerical_contacts"]:
        a = np.array([complex(*v) for v in c["a"]]); b = np.array([complex(*v) for v in c["b"]])
        U = np.array(c["U_re"]) + 1j * np.array(c["U_im"])
        cs.append(dict(src=c["source"], a=a, b=b, U=U))
    extra = os.path.join(OUT, "n4_contacts.json")
    if os.path.exists(extra):
        for c in json.load(open(extra))["contacts"]:
            a = np.array([complex(*v) for v in c["a"]]); b = np.array([complex(*v) for v in c["b"]])
            U = np.array(c["U_re"]) + 1j * np.array(c["U_im"])
            cs.append(dict(src=c["src"], a=a, b=b, U=U))
    return cs


def stage_vertex(seed, K=400):
    rng = np.random.default_rng(seed)
    cs = load_vertex_contacts()
    allres = []
    for ci, c in enumerate(cs):
        r4 = margin(c["a"], c["b"], c["U"])
        k4, th4, inf4 = normals_at(det_d(c["a"], c["b"], c["U"]), z_all(c["a"], c["b"]), tol=1e-7)
        lams, best, nb = [], None, 0
        for (a5, b5) in fifth_pairs(rng, 1.0, K):
            rec = analyze_lift(c["a"], c["b"], c["U"], a5, b5, tol=1e-6)
            if rec is None:
                continue
            nb += 1
            lams.append(rec["lam_max"])
            if best is None or rec["lam_max"] > best["lam_max"]:
                best = rec
        lams = np.array(lams)
        allres.append(dict(contact=ci, src=c["src"], m4=r4["m"], kind4=k4, info4={k: (list(v) if isinstance(v, tuple) else v) for k, v in inf4.items()},
                           dist_genperm=float(dist_genperm(c["U"])[0]), dist_blockdiag=float(dist_blockdiag(c["U"])[0]),
                           lifts=K, on_boundary=nb,
                           lam_max_overall=float(lams.max()) if nb else None,
                           lam_quantiles=np.quantile(lams, [0, .5, .9, 1]).tolist() if nb else [], best=best))
        print(ci, c["src"], k4, nb, allres[-1]["lam_max_overall"])
    save("vertex_s%d.json" % seed, dict(seed=seed, results=allres))


if __name__ == "__main__":
    st = sys.argv[1]
    if st == "check":
        stage_check()
    elif st == "edge":
        stage_edge(int(sys.argv[2]), int(sys.argv[3]))
    elif st == "vertex":
        stage_vertex(int(sys.argv[2]))


def stage_direct(seed, R, sources=("edge_s501.json", "vertex_s502.json")):
    """Direct n = 5 ascent from lifted contacts: U = U0 expm(eps K) with K random (off-block
    and in-block parts), then toolkit smooth_maximize (spectra move too) and polish.
    Every result is re-scored with the robust hull (lift_lib.robust_margin)."""
    from moc_opt import smooth_maximize, polish
    from lift_lib import robust_margin
    rng = np.random.default_rng(seed)
    seeds = []
    for src in sources:
        p = os.path.join(OUT, src)
        if not os.path.exists(p):
            continue
        j = json.load(open(p))
        if "best" in j:
            lst = j["best"]
        else:
            lst = [r.get("best") or r.get("final") for r in j["results"]]
            lst = [t for t in lst if t]
        for rec in lst:
            a = np.array([complex(*v) for v in rec["a"]]); b = np.array([complex(*v) for v in rec["b"]])
            U0 = np.array(rec["U_re"]) + 1j * np.array(rec["U_im"])
            seeds.append((src, rec["lam_max"], a, b, U0))
    recs, ms, nev = [], [], 0
    for r in range(R):
        if time.time() - T0 > TLIM:
            break
        src, lam, a, b, U0 = seeds[r % len(seeds)]
        eps = [0.02, 0.1, 0.3][(r // len(seeds)) % 3]
        X = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
        K = 0.5 * (X - X.conj().T)
        K[:4, :4] *= 0.3            # mostly off-block
        U = U0 @ expm(eps * K / np.linalg.norm(K))
        m_start = robust_margin(det_d(a, b, U), z_all(a, b))
        try:
            a2, b2, U2, th, info = smooth_maximize(a, b, U)
            nev += info["nev"]
            try:
                a3, b3, U3, th3, _ = polish(a2, b2, U2, th)
                if robust_margin(det_d(a3, b3, U3), z_all(a3, b3)) > robust_margin(det_d(a2, b2, U2), z_all(a2, b2)):
                    a2, b2, U2 = a3, b3, U3
            except Exception:
                pass
        except Exception as e:
            continue
        m_rob = robust_margin(det_d(a2, b2, U2), z_all(a2, b2))
        m_tk = margin(a2, b2, U2)["m"]
        S = np.abs(a2[:, None] + b2[None, :]).min()
        ms.append(m_rob)
        recs.append(dict(restart=r, src=src, seed_lam=lam, eps=eps, m_start=m_start, m_robust=m_rob, m_toolkit=m_tk,
                         dist_genperm=float(dist_genperm(U2)[0]), dist_blockdiag=float(dist_blockdiag(U2)[0]),
                         min_abs_a_plus_b=float(S), a=cj(a2), b=cj(b2), U_re=U2.real.tolist(), U_im=U2.imag.tolist()))
        recs.sort(key=lambda t: -t["m_robust"])
        del recs[20:]
    ms = np.array(ms)
    summ = dict(stage="direct", seed=seed, runs=len(ms), seeds_available=len(seeds), nev=nev,
                max_m_robust=float(ms.max()) if len(ms) else None, n_gt_1e7=int(np.sum(ms > 1e-7)),
                quantiles=np.quantile(ms, [0, .5, .9, 1]).tolist() if len(ms) else [], seconds=time.time() - T0)
    save("direct_s%d.json" % seed, dict(summary=summ, best=recs))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__" and sys.argv[1] == "direct":
    stage_direct(int(sys.argv[2]), int(sys.argv[3]),
                 sources=tuple(sys.argv[4].split(",")) if len(sys.argv) > 4 else ("edge_s501.json", "vertex_s502.json"))


def _sep_ok(a, b, tol=0.05):
    s = np.sqrt(np.mean(np.abs(np.concatenate([a - a.mean(), b - b.mean()])) ** 2))
    sab = np.abs(a[:, None] + b[None, :]).min()
    sa = min(abs(a[i] - a[j]) for i in range(len(a)) for j in range(i))
    sb = min(abs(b[i] - b[j]) for i in range(len(b)) for j in range(i))
    return min(sab, sa, sb) / s >= tol


def stage_adv(seed, K=12, maxfev=1500, tol=1e-12):
    """Adversarial: maximize lambda_max over the spectra (all 5 pairs) and the mixing c, for the
    block-type edge contacts with the largest lambda_max in edge_s501.json (edge (s,i,j) fixed),
    and over (a5, b5) only for the lifted degenerate-spectrum contact from vertex_s502.json.
    Nondegeneracy guard: all |a_i - a_j|, |b_i - b_j|, |a_i + b_j| >= 0.05 x rms spread."""
    from lift_lib import robust_margin
    rng = np.random.default_rng(seed)
    j = json.load(open(os.path.join(OUT, "edge_s501.json")))
    res = []
    seen = set()
    for rec in j["best"]:
        key = (tuple(rec["edge"]["s"]), rec["edge"]["i"], rec["edge"]["j"], rec["restart"])
        if key in seen:
            continue
        seen.add(key)
        if len(seen) > K or time.time() - T0 > TLIM:
            break
        s, i, jj = tuple(rec["edge"]["s"]), rec["edge"]["i"], rec["edge"]["j"]
        a0 = np.array([complex(*v) for v in rec["a"]]); b0 = np.array([complex(*v) for v in rec["b"]])
        c0 = rec["edge"]["c"]
        p0 = np.concatenate([a0.real, a0.imag, b0.real, b0.imag, [np.arcsin(np.sqrt(c0))]])

        def unpack(p):
            a = p[0:5] + 1j * p[5:10]; b = p[10:15] + 1j * p[15:20]
            c = float(np.sin(p[20]) ** 2)
            return a, b, c

        def obj(p):
            a, b, c = unpack(p)
            if not _sep_ok(a, b) or not (0.02 < c < 0.98):
                return 1.0
            U4 = U_block_edge(s, i, jj, c)
            r = analyze_lift(a[:4], b[:4], U4, a[4], b[4], tol=tol)
            if r is None:
                aa, bb, U0 = embed(a[:4], b[:4], U4, a[4], b[4])
                return 1.0 - robust_margin(det_d(aa, bb, U0), z_all(aa, bb))
            return -r["lam_max"]
        f0 = obj(p0)
        r1 = minimize(obj, p0, method="Nelder-Mead", options=dict(maxfev=maxfev, xatol=1e-10, fatol=1e-14))
        a, b, c = unpack(r1.x)
        fin = analyze_lift(a[:4], b[:4], U_block_edge(s, i, jj, c), a[4], b[4], tol=tol)
        res.append(dict(kind="block-edge", edge=dict(s=list(s), i=i, j=jj), lam_start=-f0, lam_final=-r1.fun,
                        nfev=int(r1.nfev), c=c, a=cj(a), b=cj(b), final=fin))
        print("block-edge", s, i, jj, "lam %.3e -> %.3e" % (-f0, -r1.fun), flush=True)
    # degenerate-spectrum edge contact: vary (a5, b5) only
    v = json.load(open(os.path.join(OUT, "vertex_s502.json")))["results"][0]
    best = v["best"]
    a = np.array([complex(*t) for t in best["a"]]); b = np.array([complex(*t) for t in best["b"]])
    U0 = np.array(best["U_re"]) + 1j * np.array(best["U_im"]); U4 = U0[:4, :4]
    for trial in range(4):
        q0 = np.array([a[4].real, a[4].imag, b[4].real, b[4].imag]) + (0 if trial == 0 else 0.3 * rng.normal(size=4))

        def objq(q):
            r = analyze_lift(a[:4], b[:4], U4, q[0] + 1j * q[1], q[2] + 1j * q[3], tol=tol)
            if r is None:
                return 1.0
            return -r["lam_max"]
        f0 = objq(q0)
        if f0 >= 1.0:
            continue
        r1 = minimize(objq, q0, method="Nelder-Mead", options=dict(maxfev=maxfev, xatol=1e-10, fatol=1e-14))
        fin = analyze_lift(a[:4], b[:4], U4, r1.x[0] + 1j * r1.x[1], r1.x[2] + 1j * r1.x[3], tol=tol)
        res.append(dict(kind="collision-edge (a0+b3=0)", lam_start=-f0, lam_final=-r1.fun, nfev=int(r1.nfev),
                        a5=[r1.x[0], r1.x[1]], b5=[r1.x[2], r1.x[3]], final=fin))
        print("collision", trial, "lam %.3e -> %.3e" % (-f0, -r1.fun), flush=True)
    save("adv_s%d.json" % seed, dict(seed=seed, tol=tol, results=res, seconds=time.time() - T0))


if __name__ == "__main__" and sys.argv[1] == "adv":
    stage_adv(int(sys.argv[2]), tol=float(sys.argv[3]) if len(sys.argv) > 3 else 1e-12)
