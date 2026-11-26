"""Exact Grassmannian support certificates for an entire split-plane fiber."""
import argparse
import json
import runpy
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
g = runpy.run_path(str(ROOT / 'double-paired-five-term' / 'grassmann.py'))
prior = runpy.run_path(str(ROOT / 'split-one-eigenvalue-pair' / 'check.py'))
colors = runpy.run_path(str(ROOT / 'three-color-coupling' / 'coupling.py'))
Q, exact = colors['Q'], colors['e']
mm, tr, minor, parity = g['multiply'], g['transpose'], g['minor'], g['parity']
check = g['require']
PAIRS = g['PAIRS']
ROWS = g['ROW_PAIRS']
STAR = [[parity(i+j) if len(set(i+j)) == 4 else 0 for j in PAIRS] for i in PAIRS]
# Each shift comes with a coordinate fixed real by the global phase of a bivector.
SHIFTS = {8: (-4,0), 11: (8,0), 14: (-16,5), 46: (21,0), 47: (-29,5)}


def table_index(groups):
    return tuple(sum(i//2 == color for i in group) for group in groups for color in range(3))


def matrices(tables):
    """Precompute the 33 symmetric basis forms. Columns are singleton,singleton,pair,pair."""
    x, y = g['frames'](1,1,1,0,1)
    gram = mm(tr(y),y)
    middle = g['compound'](y,ROWS,PAIRS)
    tail = mm(mm(middle,STAR),g['compound'](gram,PAIRS,PAIRS))
    lookup = {tuple(t):i for i,t in enumerate(tables)}
    raw = [[[0]*6 for _ in PAIRS] for _ in tables]
    for i,j in permutations(range(6),2):
        if not x[i][0]*x[j][1]:
            continue
        remaining = sorted(set(range(6))-{i,j})
        for k in combinations(remaining,2):
            l = tuple(v for v in remaining if v not in k)
            dest = raw[lookup[table_index((k,l,(i,),(j,)))]]
            factor = 3*parity((i,j)+k+l)*x[i][0]*x[j][1]
            left, right = middle[ROWS.index(k)],tail[ROWS.index(l)]
            for r in range(6):
                for c in range(6):
                    dest[r][c] += factor*left[r]*right[c]
    return [[[a[i][j]+a[j][i] for j in range(6)] for i in range(6)] for a in raw]


def combine(basis, ray):
    return [[sum(h*a[i][j] for h,a in zip(ray,basis)) for j in range(6)] for i in range(6)]


def principal(matrix, retained):
    return [[matrix[i][j] for j in retained] for i in retained]


def strict(matrix):
    for k in range(1,len(matrix)+1):
        check(prior['determinant_integer'](principal(matrix,range(k)))>0,'positive leading minor')


def certify(ray_id, matrix):
    # Four times the shifted form minus the identity gives a 1/4 lower bound.
    shift, pinned = SHIFTS.get(ray_id,(0,0))
    for sign in (1,-1):
        a = [[4*(matrix[i][j]+sign*shift*STAR[i][j])-int(i==j)
              for j in range(6)] for i in range(6)]
        keep = [i for i in range(6) if sign==1 or i!=pinned]
        strict(principal(a,keep))


def coarse_three(system, ray):
    """Orbit 3 is pulled back after identifying original row colors 1 and 2."""
    rows = defaultdict(lambda:[0]*33)
    for exps, row in zip(system['monomials'],system['matrix']):
        key = (exps[0],exps[1]+exps[2],*exps[3:])
        for i,v in enumerate(row):
            rows[key][i] += v
    basis = list(rows.values())
    selected = prior['independent_rows'](basis)
    check(len(selected)==12,'coarse rank')
    check(len(prior['independent_rows']([basis[i] for i in selected]+[ray]))==12,'coarse pullback')
    check(min(ray)>=0,'negative coarse table score')


def compare(basis, tables, zframe, phase):
    """Evaluate the actual 720 terms on a complex flag without square-root normalizations."""
    x,y = g['frames'](1,1,1,0,1)
    gram = mm(tr(y),y)
    adj = [[entry.conjugate() for entry in row] for row in tr(zframe)]
    annihilator = mm(adj,gram)
    (a,b),(c,d) = [row[:2] for row in annihilator]
    invnumer = [[d,-b],[-c,a]]
    top = mm(invnumer,[row[2:] for row in annihilator])
    tail = [[-v/(a*d-b*c) for v in row] for row in top]+[[Q(1),Q(0)],[Q(0),Q(1)]]
    check(all(v==Q(0) for row in mm(annihilator,tail) for v in row),'orthogonal paired complement')
    for row in tail:
        row[0] *= phase
    frame = [[Q(v) for v in first]+second+third for first,second,third in zip(x,mm(y,zframe),mm(y,tail))]
    # Basis changes in either nonsingleton block cancel in these determinant ratios.
    det = exact.determinant(frame)
    check(det.imag!=0,'complex determinant phase')
    lookup = {tuple(t):i for i,t in enumerate(tables)}
    totals = [Q(0)]*33
    for p in permutations(range(6)):
        key = table_index((p[2:4],p[4:6],p[:1],p[1:2]))
        totals[lookup[key]] += parity(p)*prod(frame[p[c]][c] for c in range(6))
    z = [minor(zframe,r,(0,1)) for r in PAIRS]
    norm = minor(mm(annihilator,zframe),(0,1),(0,1))
    check(norm.imag==0 and norm.real>0,'positive metric norm')
    for a,value in zip(basis,totals):
        predicted = sum(z[i].conjugate()*a[i][j]*z[j] for i in range(6) for j in range(6))/(54*norm)
        check(predicted==Q((value/det).real),'direct grouped coefficient disagrees')
    check(sum(totals)==det,'normalization')


def main(data):
    system = json.loads((data/'system.json').read_text())
    records = json.loads((data/'orbits.json').read_text())['records']
    check([r['id'] for r in records]==list(range(62)),'expected 62 ordered orbit representatives')
    check(sum(r['invariant']['orbit_size'] for r in records)==935,'all physical rays')
    check(max(sum(r['ray']) for r in records)==22,'ray Lipschitz multiplier')
    basis = matrices(system['tables'])
    ordinary = []
    covered = []
    for r in records:
        t = combine(basis,r['ray'])
        if any(prior['determinant_integer'](principal(t,subset))<0
               for k in range(1,7) for subset in combinations(range(6),k)):
            ordinary.append(r['id'])
        if r['grouped_coefficient']:
            check(any(all(v*max(r['ray'])==w*max(row) for v,w in zip(row,r['ray']))
                      for row in system['matrix']),'grouped coefficient status')
        elif r['id']==3:
            coarse_three(system,r['ray'])
        else:
            certify(r['id'],t)
            covered.append(r['id'])
    check(ordinary==[14,46,47] and len(covered)==51,'exact scope of the repairs')
    for k in (1,2,3):
        zframe = [[Q(1),Q(0)],[Q(0),Q(1)],[Q(1,k),Q(-2,1)],[Q(2,-1),Q(k,2)]]
        compare(basis,system['tables'],zframe,Q(1,2))
    gap = Fraction(1,4*54*9)
    radius = gap/(2*22*672)
    check(radius==Fraction(1,57480192),'tube radius')
    print('935 rays in 62 orbits. Ten coefficient and one coarse orbit are globally nonnegative.')
    print('51 orbit representatives have a physical lower bound of',gap)
    print('Indefinite ordinary forms repaired on the Pluecker quadric:',ordinary)
    print('All 33 coefficients agree with the full complex determinant on three flags.')
    print('Certified entrywise tube radius:',radius)


if __name__=='__main__':
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('--ray-data',type=Path,default=ROOT/'split-dual-rays')
    main(args.parse_args().ray_data)
