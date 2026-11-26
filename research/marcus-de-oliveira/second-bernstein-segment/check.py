"""Replay and evaluate the entire fourth-plane certificate for a2 in [0,1]."""
import argparse
import json
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'four-cycle-spectral'))
import spectral as sp

PERMS = tuple(permutations(range(4)))


@lru_cache(None)
def terms():
    data = json.loads((HERE / 'certificate.json').read_text())
    if data['degree'] != 3 or len(data['terms']) != 252:
        raise ValueError('Expected 252 cubic Bernstein control terms')
    output = []
    seen = set()
    floor = [Fraction(0) for _ in range(4)]
    for control, p, vector, q in data['terms']:
        if type(control) is not int or control not in range(4) or type(p) is not int or p not in range(24):
            raise ValueError('Invalid control or permutation index')
        if len(vector) != 3 or any(type(c) is not int or c not in (-1, 0, 1) for c in vector):
            raise ValueError('Invalid square')
        form = tuple(vector)
        mass = Fraction(q)
        key = (control, p, form)
        if not any(form) or key in seen or mass <= 0:
            raise ValueError('Expected distinct, positive nonzero square terms')
        seen.add(key)
        output.append((control, p, form, mass))
        if form == (1, 0, 0):
            floor[control] += mass
    if floor != [1, 1, 1, 1]:
        raise ValueError('Denominator floor does not replay')
    return tuple(output)


def weights_at(eta, real_fourth, imag_fourth):
    """Exact rationals in, 24 rational probabilities in column-to-row order out."""
    values = (eta, real_fourth, imag_fourth)
    if any(isinstance(v, bool) or not isinstance(v, (int, Fraction)) for v in values):
        raise TypeError('Supply exact integer or Fraction coordinates')
    eta, s, t = map(Fraction, values)
    if not 0 <= eta <= 1:
        raise ValueError('The second eigenvalue must lie on [0,1]')
    basis = [comb(3, k)*eta**k*(1-eta)**(3-k) for k in range(4)]
    nums = [Fraction(0)]*24
    for k, p, (a, b, c), mass in terms():
        nums[p] += basis[k]*mass*(a+b*s+c*t)**2
    den = sum(nums)
    if den < 1 or min(nums) < 0:
        raise ArithmeticError('Invalid numerator or denominator')
    return [w/den for w in nums]


def check():
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    count = 0
    for j in range(5):
        eta = Fraction(j, 4)
        for s, t in [(s, t) for s in range(4) for t in range(4-s)]:
            weights = weights_at(eta, s, t)
            assert sum(weights) == 1
            for x, y in product((0, 1), repeat=2):
                a = [sp.Q(x), sp.Q(eta), sp.Q(y), sp.Q(s, t)]
                target = sp.leibniz(sp.matrix_sum(a))
                actual = sum((sp.Q(w)*sp.vertex(a, p) for w, p in zip(weights, PERMS)), sp.Z)
                assert target == actual, (eta, s, t, x, y)
                count += 1
    assert count == 200
    print('252 positive coefficients; four denominator floors; 200 exact original-determinant '
          'identities on a unisolvent grid.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coordinates', nargs='*', type=Fraction, metavar='Q')
    tokens = sys.argv[1:]
    if tokens in (['-h'], ['--help']):
        parser.print_help()
        sys.exit(0)
    if tokens[:1] == ['--']:
        tokens = tokens[1:]
    args = parser.parse_args(['--', *tokens])
    if not args.coordinates:
        check()
    elif len(args.coordinates) == 3:
        print(json.dumps([str(q) for q in weights_at(*args.coordinates)]))
    else:
        parser.error('Supply eta, real(a4), imag(a4), or no arguments for the exact replay.')
