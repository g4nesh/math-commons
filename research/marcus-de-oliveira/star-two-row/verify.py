"""Exact regression for the analytic two-row compression; not a global proof."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import platform
import time

HERE = Path(__file__).resolve().parent
HELPER = HERE.parent / 'star-small-ratio' / 'verify.py'
HELPER_SHA = '710cd066cd4f0045089bbd845df9280113a31f9e335e1f157bf3bb12a0efa0bc'
if hashlib.sha256(HELPER.read_bytes()).hexdigest() != HELPER_SHA:
    raise RuntimeError('The pinned exact-arithmetic helper has changed; inspect before replay.')
spec = importlib.util.spec_from_file_location('small_ratio_exact', HELPER)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
TIMES = (F(1, 100), F(1, 5), F(1), F(5))


def scalars(U):
    if len(U) != 4 or any(len(row) != 4 for row in U) or not r.unitary(U):
        raise ValueError('Expected an exact 4-by-4 unitary.')
    if any(im != 0 or re < 0 for re, im in U[0]):
        raise ValueError('Normalize the first row to real nonnegative entries.')
    x, a, _, _ = U
    u, v, w, z = (entry[0] for entry in x)
    A, C, g0 = r.quantities(U)
    m14 = r.minor(U, (0, 1), (0, 3))
    m34 = r.minor(U, (0, 1), (2, 3))
    k = A + u*u + w*w + z*z - r.norm(m34) - r.norm(m14)
    h = 1 - w*w - r.norm(a[2]) - w*w - z*z + r.norm(m34)
    d = u*u + z*z - r.norm(m14)
    e = r.norm(r.add(r.scale(1-w*w, a[1]), r.scale(v*w, a[2]))) - z*z*r.norm(a[1])
    coeff = [A*A, A*(2*k+h), k*k+k*h+A*(d-w*w), k*d+e]
    return {'A': A, 'C': C, 'g0': g0, 'k': k, 'h': h, 'd': d, 'e': e,
            'coefficients': coeff}


def direct_compression(U):
    """Compare scalar formulas with the actual Hermitian matrix in the b,c basis."""
    data = scalars(U)
    x, _, b, c = U
    u, _, w, z = (entry[0] for entry in x)
    E = [r.z(), r.z(), r.z(1), r.z()]
    W = [r.z(-z), r.z(), r.z(), r.z(u)]
    V = [r.z(), r.z(), r.z(-z), r.z(w)]
    for t in TIMES:
        H = []
        for row in (b, c):
            entries = []
            for col in (b, c):
                def rank_one(vector):
                    return r.mul(r.dot(row, vector), r.conj(r.dot(col, vector)))
                entries.append(r.sub(r.add(rank_one(E), r.scale(t, rank_one(W))), rank_one(V)))
            H.append(entries)
        assert H[0][0][1] == H[1][1][1] == 0
        assert H[1][0] == r.conj(H[0][1])
        trace = H[0][0][0] + H[1][1][0]
        determinant = r.mul(H[0][0], H[1][1])[0] - r.norm(H[0][1])
        assert trace == data['h'] + t*data['d']
        assert determinant == t*data['e'] - w*w*data['A']
        K = data['A']/t + data['k']
        assert K + H[0][0][0] == data['g0'] + (1+1/t)*data['A'] + t*data['C']
        polynomial = sum(co*t**i for i, co in enumerate(data['coefficients']))
        assert polynomial == t*t*(K*K+K*trace+determinant)
        assert t*(2*K+trace) == 2*data['A']+(2*data['k']+data['h'])*t+data['d']*t*t
    return data


def read_fixture(name):
    payload = json.loads((HERE / name).read_text())
    if payload.get('schema') != 'exact-gaussian-rational-matrix-v1':
        raise ValueError('Missing exact matrix schema.')
    U = [[r.z(F(re), F(im)) for re, im in row] for row in payload['U']]
    data = direct_compression(U)
    assert list(map(str, data['coefficients'])) == payload['coefficients']
    for key in ('k', 'h', 'd', 'e'):
        assert str(data[key]) == payload[key]
    assert data['coefficients'][2] < 0
    if name.startswith('u11-'):
        assert U[0][0] == r.z()
    return U, data


def run():
    start = time.process_time()
    pairs = [(F(3, 5), F(4, 5)), (F(5, 13), F(12, 13))]
    mix = [(F(0), F(1)), (F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(1), F(0))]
    phases = [r.z(1), r.z(0, 1), r.z(F(3, 5), F(4, 5))]
    general = zero_u11 = zero_u14 = phase_cases = 0
    for ip, it, jt, ia, ib, ic in product(range(2), range(2), range(2), range(4), range(3), range(4)):
        p, q = pairs[ip]
        co, si = pairs[it]
        ce, te = pairs[jt]
        beta = r.scale(mix[ia][0], phases[ia % 3])
        gamma = r.scale(mix[ia][1], phases[(ia+1) % 3])
        delta = r.scale(mix[ib][0], phases[(ib+2) % 3])
        eta = r.z(mix[ib][1])
        alpha = r.scale(mix[ic][0], phases[ic % 3])
        rr = mix[ic][1]
        U = r.boundary(p, q, co, si, ce, te, beta, gamma, delta, eta, alpha, rr)
        direct_compression(U)
        general += 1

        # Entire-stratum theorem has no sampled assumptions; these replay its identity.
        U0 = r.boundary(p, q, co, si, F(1), F(0), beta, gamma, delta, eta, alpha, rr)
        data0 = direct_compression(U0)
        x, a, b, c = U0
        u, _, w, z = (entry[0] for entry in x)
        assert z == 0
        assert data0['g0'] == u*u*(1-r.norm(a[3]))+r.norm(b[2])+w*w*r.norm(c[3]) >= 0
        zero_u14 += 1

        U1 = r.boundary(p, q, F(0), F(1), ce, te, beta, gamma, delta, eta, alpha, rr)
        data1 = direct_compression(U1)
        v, w, z = p, q*ce, q*te
        s, rr_gamma = r.norm(alpha), rr*rr*r.norm(gamma)
        beta_actual, gamma_actual = r.scale(rr, beta), r.scale(rr, gamma)
        R = r.norm(r.add(r.scale(v*w, beta_actual), r.scale(z, gamma_actual)))/(q*q)
        A, k, h, d, e = (data1[key] for key in ('A', 'k', 'h', 'd', 'e'))
        assert A == v*v*s
        assert k == (1-z*z)*s+q*q*rr_gamma
        assert k+h == 2*A+R
        assert d == z*z*(1-s)
        assert e == z*z*(R-q*q*rr_gamma)
        # The magnitudes in this exact frame are known rational numbers.
        alpha_abs = mix[ic][0]
        beta_abs, gamma_abs = rr*mix[ia][0], rr*mix[ia][1]
        Rmin = (v*w*beta_abs-z*gamma_abs)**2/(q*q)
        assert R >= Rmin
        Umin = r.boundary(p, q, F(0), F(1), ce, te,
                          r.z(mix[ia][0]), r.z(-mix[ia][1]), delta, eta,
                          r.z(alpha_abs), rr)
        data_min = direct_compression(Umin)
        assert all(entry[1] == 0 for row in Umin[:2] for entry in row)
        for t in TIMES:
            original = sum(co*t**i for i, co in enumerate(data1['coefficients']))
            minimum = sum(co*t**i for i, co in enumerate(data_min['coefficients']))
            assert original-minimum == (R-Rmin)*(A*t+k*t*t+z*z*t**3) >= 0
            assert 2*A+(k+2*A+R)*t+d*t*t >= 0
        zero_u11 += 1
        phase_cases += 1

    fixtures = [read_fixture(name)[1] for name in ('negative-c2.json', 'u11-negative-c2.json')]
    # Include permutation degeneracies and the q=0 case, without any division by q.
    from itertools import permutations
    for permutation in permutations(range(4)):
        U = [[r.z(j == permutation[i]) for j in range(4)] for i in range(4)]
        direct_compression(U)
    malformed = [[r.z(i == j) for j in range(4)] for i in range(4)]
    malformed[0][0] = r.z(2)
    try:
        scalars(malformed)
    except ValueError:
        pass
    else:
        raise AssertionError('Nonunitary input was accepted.')
    return {'general_exact_unitaries': general, 'u14_zero_exact_unitaries': zero_u14,
            'u11_zero_exact_unitaries': zero_u11, 'opposite_real_phase_comparisons': phase_cases,
            'permutation_degeneracies': 24, 'saved_negative_coefficient_fixtures': 2,
            't_values': list(map(str, TIMES)), 'all_exact_identities_pass': True,
            'negative_c2_exact': [str(item['coefficients'][2]) for item in fixtures],
            'nonunitary_input_rejected': True, 'helper_sha256': HELPER_SHA,
            'python': platform.python_version(), 'cpu_seconds': time.process_time()-start,
            'scope': 'Exact finite regression, not a universal positivity proof.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-reports', action='store_true')
    arguments = parser.parse_args()
    report = run()
    if arguments.write_reports:
        (HERE / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
