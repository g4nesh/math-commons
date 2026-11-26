"""Strict paired law, analytic perturbation radius, and an exact fully split example."""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import prod
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
colors = runpy.run_path(str(ROOT/'three-color-coupling'/'coupling.py'))
double = runpy.run_path(str(ROOT/'double-paired-contingencies'/'check.py'))
g = runpy.run_path(str(ROOT/'double-paired-five-term'/'grassmann.py'))
Q, exact = colors['Q'], colors['e']
PAIRS = double['PAIRS']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def cross(z,w):
    return (z.conjugate()*w).imag


def interior(triangle, target):
    a,b,c = triangle
    area = cross(b-a,c-a)
    if not area:
        return False
    weights = (cross(b-target,c-target)/area, cross(c-target,a-target)/area,
               cross(a-target,b-target)/area)
    return min(weights)>0 and sum(weights)==1 and sum(w*z for w,z in zip(weights,triangle))==target


def reflector(v):
    norm = sum(colors['abs2'](z) for z in v)
    return [[Q(int(i==j))-2*x*y.conjugate()/norm for j,y in enumerate(v)] for i,x in enumerate(v)]


def fixture():
    u = exact.identity(6)
    for v in [[1,1,Q(0,1),1,Q(1,1),Q(0,-1)],
              [Q(1,1),-1,2,Q(0,1),1,1],
              [1,Q(0,2),-1,Q(2,1),1,-2]]:
        u = exact.multiply(reflector([exact.coerce(z) for z in v]),u)
    colors['check_unitary'](u)
    return u


def original_det(u, a, b):
    m = exact.multiply(exact.multiply(u,exact.diagonal(b)),exact.adjoint(u))
    for i in range(6):
        m[i][i] += a[i]
    return exact.determinant(m)


def run():
    u = fixture()
    tables = {g['table_of_permutations'](i,j):n for n,(i,j) in enumerate(PAIRS)}
    q = [F(0)]*21
    for c,v in zip(colors['COLORINGS'],colors['canon'](u)):
        table = tuple(sum(c[k]==s for k in (2*r,2*r+1)) for r in range(3) for s in range(3))
        q[tables[table]] += v
    rays = double['classify'](double['monomial_rows']())
    floor = min(F(1,42), min(double['dot'](q,y)/sum(y) for y in rays)/2)
    require(floor>0, 'generic unitary fixture')
    residual = [(v-floor)/(1-21*floor) for v in q]
    residual_law = double['construct_law'](residual)
    require(residual_law is not None, 'positive residual law')
    law = [floor+(1-21*floor)*v for v in residual_law]
    require(min(law)>=floor and sum(law)==1, 'strict law')
    a,b = (Q(6,-1),Q(1,6),Q(-6,-8)),(Q(-6,-8),Q(-1,2),Q(3,-1))
    p = [prod(a[i]+b[j] for i,j in enumerate(perm)) for perm in g['PERMS']]
    v = [p[i]*p[j] for i,j in PAIRS]
    av,bv = ([x for z in triple for x in (z,z)] for triple in (a,b))
    center = original_det(u,av,bv)
    require(sum(w*x for w,x in zip(law,v))==center, 'direct center determinant')
    # L1 edge lengths give rational upper bounds on Euclidean side lengths.
    def altitude_bound(t):
        l,m,n = t
        diameter = max(abs(z.real)+abs(z.imag) for z in (m-l,n-m,l-n))
        return abs(cross(m-l,n-l))/diameter if diameter else F(0)
    altitude = max(altitude_bound(t) for t in combinations(v,3))
    margin = floor*altitude
    # All max moduli are bounded by 10, hence ||A||+||B|| <= 20.
    require(all(colors['abs2'](z)<=100 for z in a+b), 'spectral norm bound')
    lipschitz = 12*(2+10)*22**5
    delta = F(1)
    while delta*lipschitz>=margin:
        delta /= 2
    require(delta>0 and margin>lipschitz*delta, 'strict support clearance')
    af = [z+Q(F((i+1)*delta,8),0) for i,z in enumerate(av)]
    bf = [z+Q(0,F((6-i)*delta,8)) for i,z in enumerate(bv)]
    require(len(set(af))==len(set(bf))==6, 'fully simple spectra')
    # Move U as well. The nontrivial rational unitary has
    # ||R-I||₂=2t/sqrt(1+t²) <= 2t < delta.
    t=delta/4
    rotated=exact.multiply(colors['rotation'](6,0,3,t,Q(0,1)),u)
    colors['check_unitary'](rotated)
    target=original_det(rotated,af,bf)
    # Find a triangle in the original labels and verify that the SAME three
    # actual permutation products certify the fully split, rotated example.
    reps={}
    for perm in permutations(range(6)):
        tab=tuple(sum(perm[k]//2==s for k in (2*r,2*r+1)) for r in range(3) for s in range(3))
        reps.setdefault(tables[tab],perm)
    split=[prod(af[i]+bf[j] for i,j in enumerate(reps[k])) for k in range(21)]
    triangle=next((t for t in combinations(range(21),3)
                   if interior([v[i] for i in t],center) and interior([split[i] for i in t],target)),None)
    require(triangle is not None, 'exact split-spectrum convex certificate')
    print('Minimum common mass >=',floor,'; support margin >=',margin)
    print('Certified simultaneous radius',delta,'; split triangle table indices',triangle)


if __name__=='__main__':
    run()
