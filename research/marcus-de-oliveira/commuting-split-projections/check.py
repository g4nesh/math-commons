"""Exact subset weights on the aligned commuting-projection family."""
import json
import runpy
from itertools import combinations, permutations
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
prev = runpy.run_path(str(ROOT/'isoclinic-split-planes'/'check.py'))
g, Q, exact, balanced = (prev[k] for k in ('g','Q','exact','balanced'))
check, mm, tr, minor = (prev[k] for k in ('check','mm','tr','minor'))
pairs, rows, star = g['PAIRS'],g['ROW_PAIRS'],balanced['STAR']
mixed,pure = prev['mixed'],prev['pure']
shifts = {8:(-2,0),10:(2,0),11:(2,3),14:(-2,5),16:(-2,2),17:(-4,2),29:(-4,2),
          34:(-2,0),36:(-2,5),46:(2,5),47:(-4,5),49:(4,5),54:(2,5),59:(2,0)}


def forms(tables, a, b, c, d, t=None):
    """Return symmetric forms with the two split columns grouped individually."""
    x,y = g['frames'](a,b,c,0,d)
    left, right = (list(v) for v in tr(x))
    if t is not None:
        r0,r1 = 1+a*a+c*c,1+b*b+d*d
        left,right = ([u+t*v for u,v in zip(left,right)],
                      [r0*v-t*r1*u for u,v in zip(left,right)])
    gram = mm(tr(y),y)
    y2 = g['compound'](y,rows,pairs)
    hodge = mm(mm(y2,star),g['compound'](gram,pairs,pairs))
    loc = {tuple(row):j for j,row in enumerate(tables)}
    result = [[[0]*6 for _ in pairs] for _ in tables]
    for i,j in permutations(range(6),2):
        if not prev['terms'](left[i]) or not prev['terms'](right[j]):
            continue
        remaining = [k for k in range(6) if k not in (i,j)]
        for k in combinations(remaining,2):
            l = tuple(q for q in remaining if q not in k)
            score = g['parity']((i,j)+k+l)*left[i]*right[j]
            dest = result[loc[balanced['table_index']((k,l,(i,),(j,)))]]
            vec, cov = y2[rows.index(k)],hodge[rows.index(l)]
            for u in range(6):
                for v in range(6):
                    dest[u][v] += score*vec[u]*cov[v]
    return [[[m[i][j]+m[j][i] for j in range(6)] for i in range(6)] for m in result]


def certify_polynomial(p, dimension, squares, strict):
    """Subtract positive multiples of (x^(u/2)-x^(v/2))^2, then check coefficients."""
    left = prev['terms'](p).copy()
    def index(e):
        check(len(e)==4 and all(type(v) is int and 0<=v<64 for v in e),'bad exponent')
        return sum(v*64**i for i,v in enumerate(e))
    for row in squares:
        low,high = row['endpoints']
        w = int(row['weight'])
        check(w>0 and all(v%2==0 for v in low+high),'invalid square')
        mid = [(x+y)//2 for x,y in zip(low,high)]
        for exp,change in ((low,-w),(mid,2*w),(high,-w)):
            key = index(exp)
            left[key] = left.get(key,0)+change
    count = 0
    for key,value in left.items():
        if not value:
            continue
        e = [(key//64**i)%64 for i in range(4)]
        check(0<=key<64**4 and sum(e)<=8*dimension,'invalid packed degree')
        check(value>0 and all(v%2==0 for v in e),'signed remainder')
        count += 1
    check(count or not strict,'strict bound needs a positive monomial')


def certify(orbit, mat, shift_unit, squares, used):
    det = prev['determinant_factory'](mat)
    check(all(not prev['terms'](mat[i][j]) for i in mixed for j in pure),'cross-sector term')
    check(not prev['terms'](mat[1][4]),'pure off-diagonal')
    if orbit==3:
        for k in range(1,7):
            for ix in combinations(range(6),k):
                certify_polynomial(det(ix),k,[],False)
        return 63
    scale,pin = shifts.get(orbit,(0,0))
    count = 0
    for sign in ((1,-1) if scale else (1,)):
        adj = [[mat[i][j]+sign*scale*star[i][j]*shift_unit for j in range(6)] for i in range(6)]
        det = prev['determinant_factory'](adj)
        sectors = [('mixed',tuple(v for v in mixed if sign==1 or v!=pin))]
        if sign==1:
            sectors.append(('pure',pure))
        for label,ids in sectors:
            for k in range(1,len(ids)+1):
                key = f'{orbit}/{label}/{k}' if sign==1 else ''
                corrections = squares.get(key,[])
                if corrections:
                    used.add(key)
                certify_polynomial(det(ids[:k]),k,corrections,True)
                count += 1
    return count


def compare(tables,a,b,c,d,t):
    basis = forms(tables,a,b,c,d,t)
    x,y = g['frames'](a,b,c,0,d)
    r0,r1 = 1+a*a+c*c,1+b*b+d*d
    flag = [[Q(1),Q(0)],[Q(0),Q(1)],[Q(1,2),Q(-2,1)],[Q(2,-1),Q(1,2)]]
    gram = mm(tr(y),y)
    ann = mm([[z.conjugate() for z in row] for row in tr(flag)],gram)
    (u,v),(w,h) = [r[:2] for r in ann]
    top = mm([[h,-v],[-w,u]],[r[2:] for r in ann])
    tail = [[-z/(u*h-v*w) for z in r] for r in top]+[[Q(1),Q(0)],[Q(0),Q(1)]]
    check(all(z==Q(0) for r in mm(ann,tail) for z in r),'complement frame')
    for row in tail:
        row[0] *= Q(2,1)
    if t is not None:
        x = [[u+t*v,r0*v-t*r1*u] for u,v in x]
    frame = [[Q(x0),Q(x1)]+mid+last for (x0,x1),mid,last in zip(x,mm(y,flag),mm(y,tail))]
    total = exact.determinant(frame)
    check(total.imag!=0,'complex phase')
    direct = [Q(0)]*33
    loc = {tuple(row):j for j,row in enumerate(tables)}
    for p in permutations(range(6)):
        j = loc[balanced['table_index']((p[2:4],p[4:],p[:1],p[1:2]))]
        direct[j] += g['parity'](p)*prod(frame[p[k]][k] for k in range(6))
    z = [minor(flag,ij,(0,1)) for ij in pairs]
    gz = minor(mm(ann,flag),(0,1),(0,1))
    denom = 2*r0*r1*gz*(1 if t is None else r0+t*t*r1)
    for m,val in zip(basis,direct):
        compressed = sum(z[i].conjugate()*m[i][j]*z[j] for i in range(6) for j in range(6))
        check(compressed/denom==Q((val/total).real),'direct 720-term comparison')


def main():
    data = ROOT/'split-dual-rays'
    tables = json.loads((data/'system.json').read_text())['tables']
    records = json.loads((data/'orbits.json').read_text())['records']
    squares = json.loads((HERE/'squares.json').read_text())
    check(len(records)==62 and sum(r['invariant']['orbit_size'] for r in records)==935,'dual classification')
    a,b,c,d = [prev['poly']({64**j:1}) for j in range(4)]
    basis = forms(tables,a,b,c,d)
    used = set()
    checks = [certify(r['id'],balanced['combine'](basis,r['ray']),a*b*c*d,squares,used)
              for r in records if not r['grouped_coefficient']]
    check(len(checks)==52 and used==set(squares),'certificate coverage')
    for vals in ((1,2,3,4,None),(0,2,3,0,None),(1,2,3,4,2)):
        compare(tables,*vals)
    print('Aligned singleton axes: 52 non-coefficient orbits and',sum(checks),'polynomial checks.')
    print(sum(map(len,squares.values())),'exact binomial squares; three direct complex-flag comparisons.')


if __name__=='__main__':
    main()
