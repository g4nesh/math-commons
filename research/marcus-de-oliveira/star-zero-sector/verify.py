"""Exact symbolic replay of the U11=0 proof, using only the standard library."""
from fractions import Fraction
from itertools import product
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import platform
import time

N = 7
ZERO = (0,)*N


class Poly(dict):
    """Sparse polynomials over Z in V,W,Z,t,s,K0,R0; exact equality is dict equality."""
    @staticmethod
    def lift(value):
        if isinstance(value, Poly):
            return value
        if isinstance(value, int):
            return Poly({ZERO: value}) if value else Poly()
        raise TypeError('Only integer polynomial coefficients are supported.')

    def __add__(self, other):
        result = Poly(self)
        for monomial, value in self.lift(other).items():
            result[monomial] = result.get(monomial, 0)+value
        return Poly({m: c for m, c in result.items() if c})

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.items()})

    def __sub__(self, other):
        return self + -self.lift(other)

    def __rsub__(self, other):
        return self.lift(other) + -self

    def __mul__(self, other):
        result = {}
        for left, a in self.items():
            for right, b in self.lift(other).items():
                powers = tuple(x+y for x, y in zip(left, right))
                result[powers] = result.get(powers, 0)+a*b
        return Poly({m: c for m, c in result.items() if c})

    __rmul__ = __mul__

    def __pow__(self, degree):
        if not isinstance(degree, int) or degree < 0:
            raise ValueError('Polynomial power must be a nonnegative integer.')
        result = self.lift(1)
        for _ in range(degree):
            result = result*self
        return result

    def substitute(self, index, replacement):
        result = Poly()
        for powers, coefficient in self.items():
            residual = list(powers)
            residual[index] = 0
            result += Poly({tuple(residual): coefficient})*replacement**powers[index]
        return result


def variable(index):
    powers = [0]*N
    powers[index] = 1
    return Poly({tuple(powers): 1})


def square_sum(V, W, Z, t):
    return (V*W*(V-Z*t)**2
            +(2*V**3*W+V**3*Z+V**2*W*Z+2*V**2*Z**2)*t
            +(V**3*W+2*V**3*Z+6*V**2*Z**2+2*V*W*Z**2)*t**2
            +(V**3*Z+2*V**2*W*Z+6*V**2*Z**2+V*W**2*Z+5*V*W*Z**2+W**2*Z**2)*t**3
            +(2*V**2*Z**2+V*W*Z**2)*t**4)


def matrix_entries(V, W, Z, t):
    Q, K1 = W+Z, V+W
    return (V**2*W+V*(V*W+Z**2)*t+Z*(V*W+Q*K1)*t**2,
            V*(Z+Q**2)+(Z*K1+Q*V*(W+3*Z))*t+V*Z*(Q+Z)*t**2,
            V+K1*t+Z*t**2)


def symbolic_checks():
    V, W, Z, t, s, K0, R0 = (variable(i) for i in range(N))
    Q, K1, m = W+Z, V+W, 1-s
    A, k, R, d = V*s, s*K1+m*K0, m*R0, Z*m
    F = (A*A+A*(k+2*A+R)*t+(k*(2*A+R)+A*(d-W))*t*t
         +(k*d+Z*(R-m*K0))*t**3)
    F0 = t*t*R0*(K0+Z*t)
    F1 = V*V+V*(3*V+W)*t+V*(2*V+W)*t*t
    Fc = (V*(K0+R0)*t+(2*V*K0+K1*R0+V*(Z-W))*t*t
          +Z*(R0+K1-K0)*t**3)
    assert F == m*m*F0+s*m*Fc+s*s*F1, 'Bernstein identity failed.'

    m11, m22, L = matrix_entries(V, W, Z, t)
    # Assemble coefficients of |B|², |G|² and 2 vwz Re(conj(B)G)
    # independently from Fc/t, using |B|²+|G|²=1.
    raw11 = V*W*L+Q*(V*(Z-W)*t+Z*K1*t*t)
    raw22 = Z*L+Q*Q*(V+2*V*t-Z*t*t)+Q*(V*(Z-W)*t+Z*K1*t*t)
    assert raw11 == m11
    assert not (raw22-m22).substitute(0, 1-W-Z)
    assert V+K1*t+Z*t*t == L
    determinant = m11*m22-V*W*Z*L*L
    residual = determinant-Q*Q*square_sum(V, W, Z, t)
    assert not residual.substitute(0, 1-W-Z), 'Determinant square identity failed.'

    # The positivity certificate consists of one weighted square plus polynomials
    # with nonnegative integer coefficients, rather than a sign sample.
    remainder = square_sum(V, W, Z, t)-V*W*(V-Z*t)**2
    assert all(coefficient >= 0 for coefficient in remainder.values())
    assert all(coefficient >= 0 for coefficient in m11.values())
    assert all(coefficient >= 0 for coefficient in m22.values())
    return {'bernstein_identity': True, 'quadratic_form_assembly': True,
            'determinant_square_identity_modulo_V_plus_W_plus_Z_equals_1': True,
            'nonnegative_certificate_coefficients': True,
            'arithmetic': 'Sparse integer polynomials; no sampling or tolerance.'}


HERE = Path(__file__).resolve().parent
TWO_ROW = HERE.parent / 'star-two-row' / 'verify.py'
TWO_ROW_SHA = '9b50005b4bbd40b0b985c54af5c7b0a0ce91e7bd3860d57015c2718a9a3ee6f5'


def rational_sqrt(value):
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    if numerator*numerator != value.numerator or denominator*denominator != value.denominator:
        raise ValueError('Regression frame must have an exact rational square root.')
    return Fraction(numerator, denominator)


def matrix_replay():
    if hashlib.sha256(TWO_ROW.read_bytes()).hexdigest() != TWO_ROW_SHA:
        raise RuntimeError('Pinned two-row exact helper changed; inspect before replay.')
    spec = importlib.util.spec_from_file_location('two_row_exact', TWO_ROW)
    previous = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(previous)
    r = previous.r
    F = Fraction
    saved, _ = previous.read_fixture('u11-negative-c2.json')
    frames = [saved]
    first_rows = [(F(0), F(1), F(3, 5), F(4, 5)),
                  (F(1), F(0), F(3, 5), F(4, 5)),
                  (F(3, 5), F(4, 5), F(3, 5), F(4, 5)),
                  (F(5, 13), F(12, 13), F(4, 5), F(3, 5)),
                  (F(3, 5), F(4, 5), F(0), F(1)),
                  (F(3, 5), F(4, 5), F(1), F(0))]
    for (p, q, co, si), (alpha, rr) in product(first_rows,
            [(F(0), F(1)), (F(3, 5), F(4, 5)), (F(1), F(0))]):
        frames.append(r.boundary(p, q, F(0), F(1), co, si,
                                 r.z(0, F(3, 5)), r.z(F(4, 5)),
                                 r.z(F(3, 5)), r.z(0, F(4, 5)),
                                 r.z(0, alpha), rr))
    q_zero = s_one = middle = 0
    for U in frames:
        data = previous.direct_compression(U)
        x, a, _, _ = U
        assert x[0] == r.z()
        v, w, z = (entry[0] for entry in x[1:])
        V, W, Z = v*v, w*w, z*z
        Q, K1, s = W+Z, V+W, r.norm(a[0])
        if Q == 0:
            q_zero += 1
            assert data['g0'] >= 0
            continue
        for t in previous.TIMES:
            direct = sum(co*t**i for i, co in enumerate(data['coefficients']))
            F1 = V*V+V*(3*V+W)*t+V*(2*V+W)*t*t
            if s == 1:
                assert direct == F1 >= 0
                continue
            q, radius = rational_sqrt(Q), rational_sqrt(1-s)
            beta = r.scale(1/q, r.sub(r.scale(z, a[2]), r.scale(w, a[3])))
            gamma = r.scale(1/q, a[1])
            B, G = r.scale(1/radius, beta), r.scale(1/radius, gamma)
            assert r.norm(B)+r.norm(G) == 1
            K0 = Q*r.norm(G)
            R0 = r.norm(r.add(r.scale(v*w, B), r.scale(z, G)))/Q
            F0 = t*t*R0*(K0+Z*t)
            Fc = (V*(K0+R0)*t+(2*V*K0+K1*R0+V*(Z-W))*t*t
                  +Z*(R0+K1-K0)*t**3)
            m11, m22, L = matrix_entries(V, W, Z, t)
            off = v*w*z*L
            form = m11*r.norm(B)+m22*r.norm(G)+2*off*r.mul(r.conj(B), G)[0]
            assert Q*Fc/t == form
            assert m11*m22-off*off == Q*Q*square_sum(V, W, Z, t) >= 0
            assert m11 >= 0 and m22 >= 0 and form >= 0
            assert direct == (1-s)**2*F0+s*(1-s)*Fc+s*s*F1 >= 0
            assert data['g0']+(1+1/t)*data['A']+t*data['C'] >= 0
        if s == 1:
            s_one += 1
        else:
            middle += 1
    return {'exact_unitary_frames': len(frames), 'q_zero_frames': q_zero,
            's_one_nonzero_q_frames': s_one, 'middle_normalized_frames': middle,
            't_values': list(map(str, previous.TIMES)),
            'earlier_negative_coefficient_witness_replayed': True,
            'helper_sha256': TWO_ROW_SHA,
            'scope': 'Finite regression of the universal symbolic proof, not a search.'}


def run():
    start = time.process_time()
    return {'symbolic': symbolic_checks(), 'regression': matrix_replay(),
            'python': platform.python_version(), 'cpu_seconds': time.process_time()-start,
            'claim': 'U11=0 only; global star and MDO remain open.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true')
    arguments = parser.parse_args()
    report = run()
    if arguments.write_report:
        (HERE / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
