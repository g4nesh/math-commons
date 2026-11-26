"""Exact common-parent certificate for every orientation of one row plane.

Python 3.9+ standard library; no optimization package is used in replay.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import prod
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('orientation_split_base', HERE.parent / 'split-one-eigenvalue-pair' / 'check.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
Q, exact = base.Q, base.exact
require, dot = base.require, base.dot


def load():
    return json.loads((HERE / 'certificate.json').read_text())


def rotate(r, u):
    result = [row[:] for row in u]
    for i in range(2):
        result[4 + i] = [sum((r[i][j] * u[4 + j][k] for j in range(2)), Q(0)) for k in range(6)]
    return result


def canonical_frame(a):
    """Divide complex grouped Leibniz sums by det(a), then take real part."""
    det = exact.determinant(a)
    require(det != Q(0), 'nonsingular frame')
    s = base.structure()
    q = [Q(0)] * 33
    for p in base.PERMS:
        q[s['lookup'][base.table(p)]] += base.parity(p) * prod(a[i][j] for i, j in enumerate(p)) / det
    return [x.real for x in q]


@lru_cache(None)
def response():
    """Algebraic constant,x,y,z coefficients of ALL 33 canonical masses."""
    u = base.fixture(base.load())
    base.check_unitary(u)
    det = exact.determinant(u)
    s = base.structure()
    out = [[Q(0)] * 4 for _ in range(33)]
    for p in base.PERMS:
        factor = base.parity(p) * prod(u[i][p[i]] for i in range(4)) / det
        a, b, c, d = u[4][p[4]], u[5][p[4]], u[4][p[5]], u[5][p[5]]
        # For R=[[u,v],[-conj(v),conj(u)]], R* sigma_z R=n.sigma.
        coefficients = ((a*d-b*c)/2, (-a*c+b*d)/2,
                        Q(0,-1)*(a*c+b*d)/2, (a*d+b*c)/2)
        row = out[s['lookup'][base.table(p)]]
        for j, value in enumerate(coefficients):
            row[j] += factor * value
    return tuple(tuple(x.real for x in row) for row in out)


def canonical_at(n):
    return [row[0] + dot(row[1:], n) for row in response()]


def rational_direction(n):
    require(isinstance(n, (list, tuple)) and len(n) == 3, 'three exact coordinates')
    require(all(type(x) in (int, F) for x in n), 'rational coordinates, not floats or booleans')
    n = tuple(map(F, n))
    require(dot(n, n) == 1, 'unit Bloch direction')
    return n


def direction_frame(n):
    n = rational_direction(n)
    x, y, z = n
    if z == -1:
        return [[Q(0), Q(1)], [Q(-1), Q(0)]]
    return [[Q(1+z), Q(x,-y)], [Q(-x,-y), Q(1+z)]]


def inverse(a):
    n = len(a)
    x = [[F(v) for v in row] + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if x[i][j]), None)
        if pivot is None:
            return None
        x[j], x[pivot] = x[pivot], x[j]
        c = x[j][j]
        x[j] = [v/c for v in x[j]]
        for i in range(n):
            if i != j:
                c = x[i][j]
                x[i] = [u-c*v for u, v in zip(x[i], x[j])]
    return [row[n:] for row in x]


def imbalances_at(n, data):
    """Bounded exact vertex search: at most C(12,6)*64 candidates."""
    s = base.structure()
    masses = list(map(F, data['parent_masses']))
    floor = F(data['child_floor'])
    cap = [masses[p]-2*floor for p, _ in s['doubles']]
    target = [dot(list(map(F, row)), n) for row in data['moment_response']]
    attempts = 0
    for free in combinations(range(12), 6):
        inv = inverse([[row[j] for j in free] for row in s['b']])
        if inv is None:
            continue
        fixed = [j for j in range(12) if j not in free]
        for signs in product((-1, 1), repeat=6):
            attempts += 1
            d = [F(0)] * 12
            for j, sign in zip(fixed, signs):
                d[j] = sign*cap[j]
            rhs = [target[i]-sum(s['b'][i][j]*d[j] for j in fixed) for i in range(6)]
            values = [dot(row, rhs) for row in inv]
            if all(abs(x) <= cap[j] for j, x in zip(free, values)):
                for j, x in zip(free, values):
                    d[j] = x
                return d, attempts
    raise ValueError('no lift found within the finite vertex enumeration')


def weights_at(n, data=None):
    """Return 33 exact masses for a rational unit Bloch vector."""
    n = rational_direction(n)
    data = load() if data is None else data
    verify(data, fixtures=False)
    d, _ = imbalances_at(n, data)
    result = base.lift(canonical_at(n), list(map(F, data['repairs'])), d)
    require(min(result) >= F(data['child_floor']) and sum(result) == 1, 'uniform child floor')
    require(all(dot(row, [a-b for a,b in zip(result, canonical_at(n))]) == 0
                for row in base.structure()['matrix']), 'all 192 polynomial coefficients')
    return result


def nearby_weights(u, rotation, data=None):
    """Exact Gaussian-rational U and U(2) rotation in the uniform open family."""
    data = load() if data is None else data
    base.check_unitary(u)
    require(len(rotation)==2 and all(len(row)==2 and all(isinstance(x,Q) for x in row) for row in rotation), 'exact 2 by 2 rotation')
    require(exact.multiply(exact.adjoint(rotation),rotation)==exact.identity(2), 'unitary rotation')
    old = base.fixture(base.load())
    radius = F(1,53760000)
    require(all((x-y).real**2+(x-y).imag**2 <= radius**2
                for row, before in zip(u,old) for x,y in zip(row,before)), 'outside uniform neighborhood')
    h = exact.multiply(exact.adjoint(rotation),exact.multiply([[Q(1),Q(0)],[Q(0),Q(-1)]],rotation))
    n = (h[0][1].real,-h[0][1].imag,h[0][0].real)
    initial = weights_at(n,data)
    q = base.canonical(rotate(rotation,u))
    q0 = base.canonical(rotate(rotation,old))
    result = [w+a-b for w,a,b in zip(initial,q,q0)]
    require(min(result)>=F(1,40000) and sum(result)==1, 'uniform nearby floor')
    require(all(dot(row,[a-b for a,b in zip(result,q)])==0 for row in base.structure()['matrix']), 'nearby coefficients')
    return result


def mixed_rotations():
    for parameters in ((F(1,3), F(1,4), F(1,5)), (F(2,3), F(-1,2), F(3,7)), (F(1),F(2),F(3))):
        norm = dot(parameters, parameters)
        a = (1-norm)/(1+norm)
        b, c, d = [2*x/(1+norm) for x in parameters]
        yield [[Q(a,b), Q(c,d)], [Q(-c,d), Q(a,-b)]]


def verify(data=None, fixtures=True):
    data = load() if data is None else data
    require(data['schema'] == 'mdo.split-row-orientations.v1', 'schema')
    s = base.structure()
    base.verify(minors=False)
    r = response()
    coarse = [sum(r[j][0] for j in f) for f in s['fibers']]
    require(all(sum(r[j][k] for j in f) == 0 for f in s['fibers'] for k in (1,2,3)), 'invariant parent polynomial')
    require(all(r[f[0]][0] == r[f[1]][0] for _,f in s['doubles']), 'odd constant zero')
    d = [[r[f[0]][k]-r[f[1]][k] for k in (1,2,3)] for _,f in s['doubles']]
    t = [[sum(row[j]*d[j][k] for j in range(12)) for k in range(3)] for row in s['b']]
    require(d == [list(map(F,row)) for row in data['odd_response']], '12 by 3 odd response')
    require(t == [list(map(F,row)) for row in data['moment_response']], '6 by 3 moment response')
    require(coarse == list(map(F,data['parent_canonical'])), 'fixed canonical parents')
    repairs = list(map(F,data['repairs']))
    require(len(repairs) == 6, 'six repairs')
    masses = [coarse[p]+sum(repairs[v]*s['kernel'][v][p] for v in range(6)) for p in range(21)]
    require(masses == list(map(F,data['parent_masses'])), 'parent masses')
    floor = F(data['child_floor'])
    require(floor == F(1,20000), 'claimed child floor')
    require(1344 * F(1,53760000) == floor/2 == F(1,40000), 'uniform neighborhood constants')
    require(sum(masses) == 1 and min(masses) > 0, 'strict normalized parent')
    require(all(masses[p] >= floor for p,f in enumerate(s['fibers']) if len(f)==1), 'single child floor')
    cap = [masses[p]-2*floor for p,_ in s['doubles']]
    require(min(cap) > 0, 'strict shrunken capacities')
    squares, gaps, reduced_gaps = [], [], []
    for h in s['normals']:
        magnitudes = [abs(sum(h[i]*s['b'][i][j] for i in range(6))) for j in range(12)]
        square = sum(sum(h[i]*t[i][k] for i in range(6))**2 for k in range(3))
        support = sum(masses[p]*v for v,(p,_) in zip(magnitudes,s['doubles']))
        reduced = dot(magnitudes,cap)
        require(reduced >= 0 and reduced*reduced >= square, 'full sphere support containment')
        squares.append(square); gaps.append(support*support-square); reduced_gaps.append(reduced*reduced-square)
    require(squares == list(map(F,data['norm_squares'])) and gaps == list(map(F,data['squared_slacks'])), 'archived support arithmetic')
    directions, attempts = [], []
    if fixtures:
        u = base.fixture(base.load())
        for rotation in mixed_rotations():
            require(exact.multiply(exact.adjoint(rotation),rotation) == exact.identity(2), 'SU2 fixture')
            h = exact.multiply(exact.adjoint(rotation),exact.multiply([[Q(1),Q(0)],[Q(0),Q(-1)]],rotation))
            n = rational_direction((h[0][1].real,-h[0][1].imag,h[0][0].real))
            q = base.canonical(rotate(rotation,u))
            require(q == canonical_at(n), 'exact rotated canonical coefficients')
            imbalance, count = imbalances_at(n,data)
            w = base.lift(q,repairs,imbalance)
            require(min(w) >= floor and sum(w)==1, 'mixed fixture weights')
            require(all(dot(row,[a-b for a,b in zip(w,q)])==0 for row in s['matrix']), 'mixed fixture 192 coefficients')
            require(canonical_frame(rotate(direction_frame(n),u)) == q, 'rational frame gauge')
            directions.append(list(map(str,n))); attempts.append(count)
    return {'scope':'one fixed coarse unitary, every U(2) orientation of its last row plane',
            'canonical_response_shape':[33,4], 'odd_response_shape':[12,3], 'moment_response_shape':[6,3],
            'support_normal_directions':len(s['normals']), 'uniform_child_mass_floor':str(floor),
            'all_orientation_entrywise_unitary_radius':'1/53760000', 'nearby_child_mass_floor':'1/40000',
            'minimum_parent_mass':str(min(masses)), 'minimum_reduced_squared_slack':str(min(reduced_gaps)),
            'mixed_exact_SU2_directions':directions, 'vertex_candidates_tested':attempts,
            'polynomial_coefficients_per_fixture':192}


if __name__ == '__main__':
    report = verify()
    if '--write-report' in sys.argv:
        (HERE / 'report.json').write_text(json.dumps(report,indent=2)+'\n')
    else:
        require(report == json.loads((HERE / 'report.json').read_text()), 'archived report')
    print(json.dumps(report,indent=2))
