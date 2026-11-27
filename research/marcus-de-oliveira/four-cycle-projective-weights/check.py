"""Exact exterior weights for the archived four-cycle. Run with Python's standard library."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
prior = runpy.run_path(str(HERE.parent/'complex-second-region'/'check.py'))
g, Q = prior['g'], prior['Q']
COLUMNS, ALLOWED = prior['COLUMNS'], prior['ALLOWED']
RADIUS = 384
GENERATORS = [v for v in product((-1,0,1),repeat=3) if next((i for i in v if i),0)>0]


def check(condition, description):
    if not condition:
        raise ValueError(description)


def chart():
    b = prior['coefficient_matrix'](Q(0))
    a = g.difference(prior['coefficient_matrix'](Q(1)),b)
    # Multiplying a complex residual by i sends real/imaginary rows to -imag/real.
    j = [[-v for v in row] for row in b[40:]]+b[:40]
    return a,b,j


def center(data):
    n = len(ALLOWED)*len(GENERATORS)
    mass = [F(1)]*n
    used = set()
    for index,value in data['extra']:
        check(type(index) is int and 0<=index<n and index not in used,'invalid extra-square index')
        weight = F(value)
        check(weight>0,'extra weights must be positive')
        mass[index] += weight
        used.add(index)
    total = g.square((0,0,0))
    for vector in GENERATORS:
        for e,c in g.square(vector).items():
            total[e] += c
    check(all(total[e]==(9 if e in ((0,0),(2,0),(0,2)) else 0) for e in total),'base floor')
    result = {key:F(0) for key in COLUMNS}
    for k,w in enumerate(mass):
        color, generator = divmod(k,len(GENERATORS))
        for exponent,c in g.square(GENERATORS[generator]).items():
            result[ALLOWED[color],exponent] += w*c
    return [result[key] for key in COLUMNS]


@lru_cache(None)
def prepare():
    data = json.loads((HERE/'center.json').read_text())
    check(data['positiveBase']=='1','incorrect base')
    a,b,j = chart()
    c = center(data)
    check(not any(g.matvec(a,c)),'wrong polynomial kernel at infinity')
    # Verify the finite-chart rank cap on which the reciprocal correction relies.
    lo,hi = prior['LOWER'],prior['UPPER']
    base = prior['coefficient_matrix'](Q(lo))
    prior['verify_relations'](prior['load']('rank.json'),base,
        g.difference(prior['coefficient_matrix'](Q(hi)),base),
        g.difference(prior['coefficient_matrix'](Q(lo,1)),base))
    rows = g.indices(data['rows'],54,80)
    cols = g.indices(data['columns'],54,108)
    select = lambda m: [[m[i][k] for k in cols] for i in rows]
    inv = g.inverse(select(a))
    x,y = g.matmul(inv,select(b)),g.matmul(inv,select(j))
    f = g.matvec(inv, g.matvec([b[i] for i in rows], c))
    h = g.matvec(inv, g.matvec([j[i] for i in rows], c))
    op = prior['norm'](x)+prior['norm'](y)
    force = max(map(abs,f))+max(map(abs,h))
    check(op<145 and force<814,'archived rational bounds failed')
    check(F(814,RADIUS-145)<F(7,2),'bound on coefficient correction')
    return a,b,j,c,cols,x,y,f,h


def weights_at(second_real,second_imag,fourth_real,fourth_imag):
    u,v,s,t = map(g.exact,(second_real,second_imag,fourth_real,fourth_imag))
    radius_sq = u*u+v*v
    check(radius_sq>=RADIUS**2,'second entry outside certified exterior')
    u,v = u/radius_sq,-v/radius_sq
    a,b,j,c,cols,x,y,f,h = prepare()
    mat = [[F(i==k)+u*x[i][k]+v*y[i][k] for k in range(54)] for i in range(54)]
    corr = g.solve(mat,[[-u*fi-v*hi] for fi,hi in zip(f,h)])
    coeff = c.copy()
    for k,(value,) in zip(cols,corr):
        coeff[k] += value
    check(all(not (p+u*q+v*r) for p,q,r in zip(g.matvec(a,coeff),g.matvec(b,coeff),g.matvec(j,coeff))),
          'full polynomial gap at evaluation')
    numerator = [F(0)]*24
    for value,(p,(k,l)) in zip(coeff,COLUMNS):
        numerator[p] += value*s**k*t**l
    check(all(numerator[p]>=2*(1+s*s+t*t) for p in ALLOWED),'numerator floor')
    total = sum(numerator)
    return tuple(w/total for w in numerator)


if __name__=='__main__':
    for u,v,s,t in ((384,0,0,1),(-384,0,7,-9),(0,384,-10**5,1),(F(192),F(-384),F(2,3),F(5,7))):
        weights = weights_at(u,v,s,t)
        g.prior.direct_check([Q(-7,3),Q(u,v),Q(11,-5),Q(s,t)],weights)
    print('Exterior |second eigenvalue| >= 384, all other three entries arbitrary complex.')
    print('Verified rank cap 54, exact reciprocal center, operator bound 145, forcing bound 814 and four original determinants.')
