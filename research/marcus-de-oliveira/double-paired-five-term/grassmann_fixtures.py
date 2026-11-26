"""Compare the chart formulas to all 720 Leibniz terms on exact complex frames."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations
from math import prod
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
g = runpy.run_path(str(HERE/'grassmann.py'))
c = runpy.run_path(str(HERE.parent/'three-color-coupling'/'coupling.py'))
double = runpy.run_path(str(HERE.parent/'double-paired-contingencies'/'check.py'))
Q, exact = c['Q'],c['e']
mm,trans,minor,require = g['multiply'],g['transpose'],g['minor'],g['require']
pairs = list(combinations_with_replacement(range(6),2))


def adjoint(a):
    return [[x.conjugate() for x in row] for row in trans(a)]


def check(parameters, zframe):
    x,y = g['frames'](*parameters)
    gram = mm(trans(y),y)
    annihilator = mm(adjoint(zframe),gram)
    (a,b),(c,d) = [row[:2] for row in annihilator]
    inverse = [[d,-b],[-c,a]]
    top = mm(inverse, [row[2:] for row in annihilator])
    tail = [[-v/(a*d-b*c) for v in row] for row in top]+[[Q(1),Q(0)],[Q(0),Q(1)]]
    require(all(entry==Q(0) for row in mm(annihilator,tail) for entry in row),
            'last two subspaces must be orthogonal')
    # A nonreal column phase catches taking a real part before complex division.
    for row in tail:
        row[0] = row[0]*Q(1,2)
    frame = [[Q(v) for v in first]+second+third for first,second,third
             in zip(x,mm(y,zframe),mm(y,tail))]
    det = exact.determinant(frame)
    require(det.imag!=0,'fixture should exercise the determinant phase')
    lookup = {g['table_of_permutations'](i,j):k for k,(i,j) in enumerate(pairs)}
    totals = [Q(0)]*21
    for perm in permutations(range(6)):
        table = tuple(sum(perm[k]//2==s for k in (2*r,2*r+1)) for r in range(3) for s in range(3))
        totals[lookup[table]] += g['parity'](perm)*prod(frame[k][p] for k,p in enumerate(perm))
    q = [(v/det).real for v in totals]
    require(sum(q)==1, 'canonical coefficients')
    z = [minor(zframe,indices,(0,1)) for indices in g['PAIRS']]
    sdet = minor(mm(trans(x),x),(0,1),(0,1))
    znorm = minor(mm(annihilator,zframe),(0,1),(0,1))
    require(znorm.imag==0 and znorm.real>0, 'norm in the induced metric')
    for ray in g['RAYS']:
        matrix = g['facet_matrix'](x,y,ray)
        value = sum(z[i].conjugate()*matrix[i][j]*z[j] for i in range(6) for j in range(6))
        predicted = value/(2*sdet*znorm)
        direct = sum(w*q[pairs.index((i,j))] for i,j,w in ray)
        require(predicted==Q(direct) and direct>=0, 'direct 720-term chart comparison')
    law = double['construct_law'](q)
    require(law is not None and min(law)>=0 and sum(law)==1, '21 nonnegative weights')
    for row in double['monomial_rows']().values():
        require(sum((u-v)*coef for u,v,coef in zip(law,q,row))==0, 'weighted polynomial')
    return len([w for w in q if w<0])


if __name__ == '__main__':
    frames = [(1,2,3,-2,1), (2,3,1,2,4), (0,1,2,0,2), (F(1,2),F(2,3),F(1,4),F(-1,3),F(4,5))]
    for t, chart in enumerate(frames,1):
        z = [[Q(1),Q(0)],[Q(0),Q(1)],[Q(1,t),Q(2,-1)],[Q(-1,2),Q(t,1)]]
        print('Fixture',t,'negative starting masses',check(chart,z))
    print('All 16 exact facet values and 4 reconstructed 141-coefficient laws agree.')
