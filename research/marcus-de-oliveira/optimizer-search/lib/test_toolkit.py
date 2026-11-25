"""Unit tests: python -I -B test_toolkit.py  (exit 0 iff all pass)."""
import itertools
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from moc_core import (perms, z_all, z_chunks, det_d, margin, margin_points, hull_lean,  # noqa: E402
                      dist_genperm, dist_blockdiag, hull_indices, diameter, signed_distance)
from moc_param import (cayley, inverse_cayley, skew_from_real, real_from_skew, expm_unitary,  # noqa: E402
                       log_unitary, cayley_coords, random_unitary, random_spectra, canonical_spectra,
                       embed, perm_matrix)
from moc_opt import Problem, pack, smooth_maximize, polish  # noqa: E402
import moc_exact as ex  # noqa: E402
from moc_hp import hp_margin  # noqa: E402
from calibrate import psi_exact_ostar, psi_np, O_STAR  # noqa: E402

FAILS = []


def check(name, cond, info=""):
    print(("ok   " if cond else "FAIL ") + name + ("  " + str(info) if info else ""))
    if not cond:
        FAILS.append(name)


rng = np.random.default_rng(0)

# permutations and z
for n in range(1, 8):
    P = perms(n)
    check("perms(%d) complete and distinct" % n, len({tuple(p) for p in P}) == len(P) == np.prod(range(1, n + 1)))
a, b = random_spectra(5, rng)
z = z_all(a, b)
zb = [np.prod([a[i] + b[s[i]] for i in range(5)]) for s in itertools.permutations(range(5))]
check("z_all matches brute force (n=5)", np.allclose(z, zb))
a7, b7 = random_spectra(7, rng)
zc = np.concatenate([zz for _, zz in z_chunks(a7, b7)])
check("z_chunks == z_all in order (n=7)", np.allclose(zc, z_all(a7, b7)))
V, S = hull_lean(a7, b7)
Vf = z_all(a7, b7)[hull_indices(z_all(a7, b7))]
check("hull_lean == full hull (n=7)", len(V) == len(Vf) and np.allclose(sorted(V, key=lambda v: (v.real, v.imag)), sorted(Vf, key=lambda v: (v.real, v.imag))))

# d at permutation matrices
for s in [(0, 1, 2, 3, 4), (2, 0, 4, 1, 3)]:
    Pm = perm_matrix(s)
    check("d(P_sigma) = z_sigma %s" % (s,), np.isclose(det_d(a, b, Pm), np.prod([a[i] + b[s[i]] for i in range(5)])))
r = margin(a, b, perm_matrix((2, 0, 4, 1, 3)))
check("margin at a permutation is <= 0 (z_sigma may be interior)", r["m"] <= 1e-12, r["m"])
vtx = tuple(int(v) for v in perms(5)[hull_indices(z)[0]])
r = margin(a, b, perm_matrix(vtx))
check("margin at a hull-vertex permutation is ~0", abs(r["m"]) < 1e-12, r["m"])

# signed distance on a square
Vsq = np.array([0, 1, 1 + 1j, 1j])
check("signed distance inside", np.isclose(signed_distance(0.5 + 0.25j, Vsq)[0], -0.25))
dd, th = signed_distance(2 + 2j, Vsq)
check("signed distance outside (vertex region)", np.isclose(dd, np.sqrt(2)) and np.isclose(th, np.pi / 4))
check("signed distance outside (edge region)", np.isclose(signed_distance(0.5 - 1j, Vsq)[0], 1.0))
try:
    margin_points(0, z[:-1], 5)
    check("incomplete permutation list rejected", False)
except ValueError:
    check("incomplete permutation list rejected", True)

# Goyal example (Kovacec v2 Sect. 4 as transcribed in the repo): d = 2/3 + 2i/9, inside
ag = np.array([0, 1, 1j])
bg = np.array([-1 - 1j, 1, -1 + 1j])
Ug = np.array([[-1, 2, 2], [2, -1, 2], [2, 2, -1]]) / 3
rg = margin(ag, bg, Ug)
check("Goyal d", np.isclose(rg["d"], 2 / 3 + 2j / 9), rg["d"])
check("Goyal inside (m<0)", rg["m"] < 0, rg["m"])

# parametrizations
K = skew_from_real(rng.standard_normal((4, 4)))
check("skew_from_real skew-Hermitian", np.allclose(K, -K.conj().T))
check("real_from_skew inverse", np.allclose(skew_from_real(real_from_skew(K)), K))
U = cayley(K)
check("cayley unitary", np.allclose(U @ U.conj().T, np.eye(4)))
check("inverse_cayley", np.allclose(inverse_cayley(U), K))
check("expm/log unitary", np.allclose(expm_unitary(log_unitary(U)), U))
Kc, R, gap = cayley_coords(-np.eye(3) + 0j)
check("cayley_coords handles eigenvalue -1", np.allclose(cayley(Kc), -np.eye(3) * R[None, :]) and gap > 0.1)
a4, b4 = random_spectra(4, rng)
U4 = random_unitary(4, rng)
ca, cb = canonical_spectra(a4, b4)
check("canonical_spectra keeps m", np.isclose(margin(a4, b4, U4)["m"], margin(ca, cb, U4)["m"]))
check("canonical_spectra mean b = 0", abs(cb.mean()) < 1e-12)
L = np.diag(np.exp(1j * rng.random(4)))
Rr = np.diag(np.exp(1j * rng.random(4)))
check("d invariant under U -> L U R", np.isclose(det_d(a4, b4, L @ U4 @ Rr), det_d(a4, b4, U4)))
ae, be, Ue = embed(ag, bg, Ug, [2 + 1j], [0.5j])
check("embed: d factorizes", np.isclose(det_d(ae, be, Ue), det_d(ag, bg, Ug) * (2 + 1.5j)))
check("dist_genperm(P)=0, dist_blockdiag(embed)=0", dist_genperm(perm_matrix((1, 0, 2)))[0] < 1e-15 and dist_blockdiag(Ue)[0] < 1e-15)
check("dist_genperm(Haar)>0", dist_genperm(U4)[0] > 1e-3)

# gradient
for n, opts in [(3, {}), (5, dict(sep=1.0, dperm=0.9, dblock=0.9))]:
    an, bn = random_spectra(n, rng)
    pr = Problem(n, U0=random_unitary(n, rng), tau=0.05, **opts)
    x = pack(an, bn, 0.3 * rng.standard_normal((n, n)), 0.7)
    f, g = pr.fg(x)
    fd = np.array([(pr.fg(x + h)[0] - pr.fg(x - h)[0]) / 2e-6 for h in 1e-6 * np.eye(len(x))])
    check("analytic gradient n=%d %s" % (n, opts), np.abs(fd - g).max() < 1e-5 * max(1, np.abs(g).max()), np.abs(fd - g).max())

# optimizer on n = 3 lands on m <= ~0
an, bn = random_spectra(3, rng)
a2, b2, U2, th, _ = smooth_maximize(an, bn, random_unitary(3, rng))
a2, b2, U2, th, _ = polish(a2, b2, U2, th)
check("n=3 optimizer max m <= 1e-12", margin(a2, b2, U2)["m"] <= 1e-12, margin(a2, b2, U2)["m"])

# exact: psi(O*) and O* orthogonal
check("psi(O*) = -158/2401 exact", psi_exact_ostar() == Fr(-158, 2401))
check("psi float agrees", np.isclose(psi_np(np.array(O_STAR) / 7.0), -158 / 2401))

# exact: Goyal through the exact pipeline
gA = [ex.G(0), ex.G(1), ex.G(0, 1)]
gB = [ex.G(-1, -1), ex.G(1), ex.G(-1, 1)]
gU = [[ex.G(Fr(v, 3)) for v in r] for r in [[-1, 2, 2], [2, -1, 2], [2, 2, -1]]]
dG = ex.d_exact(gA, gB, gU)
check("exact Goyal d = 2/3 + 2i/9", dG == ex.G(Fr(2, 3), Fr(2, 9)))
check("exact Goyal inside", ex.membership(dG, ex.z_exact(gA, gB)) == "inside")
K3, U3, R3 = ex.exact_unitary_from_float(random_unitary(3, rng), 10 ** 6)
check("exact Cayley unitary", ex.is_unitary(U3))

# exact checker: a genuine outside instance is impossible to plant with unitary U if MOC holds
# for n = 3, so test the checker on (i) a non-unitary matrix and (ii) a valid inside instance.
with tempfile.TemporaryDirectory() as td:
    cert = ex.build_certificate(ag, bg, Ug, rg["theta"], out_json=os.path.join(td, "c.json"))
    out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "check_cert.py"), os.path.join(td, "c.json")],
                         capture_output=True, text=True)
    check("checker rejects inside instance", out.returncode == 1 and cert["status"] == "not-certified", out.stdout.strip())
    Ubad = [[ex.G(Fr(13, 10) if i == j else 0) for j in range(3)] for i in range(3)]
    cert = ex.build_certificate(ag, bg, None, 0.0, U_exact=Ubad, out_json=os.path.join(td, "c2.json"))
    out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "check_cert.py"), os.path.join(td, "c2.json")],
                         capture_output=True, text=True)
    check("checker rejects non-unitary U", out.returncode == 1 and "not unitary" in out.stdout, out.stdout.strip())
    # checker accepts a correct separation when the 'points' are honest: check the
    # arithmetic path with a 1x1 case where d = a + b equals z (gap 0 -> reject)
    c1 = dict(n=1, a=[["1", "0"]], b=[["0", "1"]], U=[[["1", "0"]]], direction=["1", "0"])
    json.dump(c1, open(os.path.join(td, "c3.json"), "w"))
    out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "check_cert.py"), os.path.join(td, "c3.json")],
                         capture_output=True, text=True)
    check("checker rejects gap 0 (n=1)", out.returncode == 1, out.stdout.strip())

# high precision
rh = hp_margin(ag, bg, Ug, dps=50)
check("hp Goyal agrees with float", abs(float(rh["m"]) - rg["m"]) < 1e-12, rh["m_str"])
check("hp U unitary to 1e-45", float(rh["unitarity_defect"]) < 1e-45)

print("\n%d failures" % len(FAILS) if FAILS else "\nall passed")
sys.exit(1 if FAILS else 0)
