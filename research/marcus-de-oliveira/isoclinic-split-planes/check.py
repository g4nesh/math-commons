"""Polynomial certificates for arbitrary scalar block projections of the split plane."""
import argparse
import hashlib
import json
import runpy
from functools import lru_cache
from itertools import combinations, permutations
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
balanced = runpy.run_path(str(ROOT/'split-physical-compression'/'check.py'))
g, Q, exact = balanced['g'], balanced['Q'], balanced['exact']
mm, tr, minor, check = g['multiply'], g['transpose'], g['minor'], g['require']
pairs, row_pairs, star = g['PAIRS'], g['ROW_PAIRS'], balanced['STAR']
poly = g['Polynomial']
mixed, pure = (0,2,3,5), (1,4)
# (coefficient of a^2 b^2 in the Hodge shift, index made real)
shifts = {8:(-2,0), 10:(2,0), 11:(2,0), 14:(-2,5), 16:(-2,2),
          17:(-4,2), 29:(-4,2), 36:(-2,2), 46:(4,0), 47:(-4,5), 49:(4,0)}


def form_basis(tables, a, b):
    first, complement = g['frames'](a,a,b,0,b)
    gram = mm(tr(complement),complement)
    middle = g['compound'](complement,row_pairs,pairs)
    tail = mm(mm(middle,star),g['compound'](gram,pairs,pairs))
    lookup = {tuple(t):i for i,t in enumerate(tables)}
    entries = [[[0]*6 for _ in pairs] for _ in tables]
    # Singleton 0 can use only an even row, singleton 1 only an odd row.
    for i in (0,2,4):
        for j in (1,3,5):
            remaining = [k for k in range(6) if k not in (i,j)]
            for k in combinations(remaining,2):
                last = tuple(q for q in remaining if q not in k)
                table = balanced['table_index']((k,last,(i,),(j,)))
                score = g['parity']((i,j)+k+last)*first[i][0]*first[j][1]
                left, right = middle[row_pairs.index(k)],tail[row_pairs.index(last)]
                result = entries[lookup[table]]
                for u in range(6):
                    for v in range(6):
                        result[u][v] += score*left[u]*right[v]
    return [[[m[i][j]+m[j][i] for j in range(6)] for i in range(6)] for m in entries]


def terms(p):
    if isinstance(p,int):
        return {0:p} if p else {}
    return p.terms


def determinant_factory(m):
    @lru_cache(None)
    def determinant(rows,cols):
        if not rows:
            return 1
        return sum((-1)**i*m[rows[0]][j]*determinant(rows[1:],cols[:i]+cols[i+1:])
                   for i,j in enumerate(cols))
    return lambda rows: determinant(rows,rows)


def check_polynomial(p, strict, order, digest):
    coefficients = terms(p)
    check(coefficients or not strict,'identically zero strict certificate')
    for key,v in coefficients.items():
        # Each factor has degree at most eight. At most six factors are multiplied.
        check(0<=key<64**2 and key%64+key//64<=8*order,'exponent packing')
        check(key%2==0 and (key//64)%2==0 and v>0,'expected positive even monomials')
    serial = sorted((key%64,key//64,v) for key,v in coefficients.items())
    digest.update(json.dumps(serial,separators=(',',':')).encode()+b'\n')
    return len(serial)


def certificate(orbit, m, shift_unit):
    digest = hashlib.sha256()
    counts = []
    check(all(not terms(m[i][j]) for i in mixed for j in pure),'coupling of different spin sectors')
    def add(polynomial, order, strict=True):
        counts.append(check_polynomial(polynomial,strict,order,digest))
    check(not terms(m[1][4]),'pure block must have zero unshifted off-diagonal')
    if orbit==3:
        # This representative can vanish. Use all principal minors.
        det = determinant_factory(m)
        for k in range(1,7):
            for subset in combinations(range(6),k):
                add(det(subset),k,False)
    else:
        scale, pinned = shifts.get(orbit,(0,0))
        sigma = scale*shift_unit
        for sign in ((1,-1) if scale else (1,)):
            adjusted = [[m[i][j]+sign*sigma*star[i][j] for j in range(6)] for i in range(6)]
            det = determinant_factory(adjusted)
            sub = tuple(i for i in mixed if sign==1 or i!=pinned)
            for k in range(1,len(sub)+1):
                add(det(sub[:k]),k)
            # The two signs have identical pure-block diagonals and determinant.
            if sign==1:
                for k in (1,2):
                    add(det(pure[:k]),k)
    return {'checked_polynomials':len(counts),'collected_monomials':sum(counts),
            'zero_polynomials':counts.count(0),'sha256':digest.hexdigest()}


def direct(tables, a, b, seed):
    basis = form_basis(tables,a,b)
    x,y = g['frames'](a,a,b,0,b)
    zframe = [[Q(1),Q(0)],[Q(0),Q(1)],[Q(1,seed),Q(-2,1)],[Q(2,-1),Q(seed,2)]]
    gram = mm(tr(y),y)
    adj = [[entry.conjugate() for entry in row] for row in tr(zframe)]
    ann = mm(adj,gram)
    (v,w),(k,l) = [row[:2] for row in ann]
    top = mm([[l,-w],[-k,v]],[row[2:] for row in ann])
    ztail = [[-entry/(v*l-w*k) for entry in row] for row in top]+[[Q(1),Q(0)],[Q(0),Q(1)]]
    check(all(entry==Q(0) for row in mm(ann,ztail) for entry in row),'complement frame')
    for row in ztail:
        row[0] *= Q(2,1)
    frame = [[Q(t) for t in left]+mid+right for left,mid,right in zip(x,mm(y,zframe),mm(y,ztail))]
    denom = exact.determinant(frame)
    check(denom.imag!=0,'determinant phase')
    lookup = {tuple(t):i for i,t in enumerate(tables)}
    direct_scores = [Q(0)]*33
    for p in permutations(range(6)):
        idx = lookup[balanced['table_index']((p[2:4],p[4:6],p[:1],p[1:2]))]
        direct_scores[idx] += g['parity'](p)*prod(frame[p[j]][j] for j in range(6))
    z = [minor(zframe,ij,(0,1)) for ij in pairs]
    norm = minor(mm(ann,zframe),(0,1),(0,1))
    check(norm.imag==0 and norm.real>0,'bivector metric')
    for mat,score in zip(basis,direct_scores):
        val = sum(z[i].conjugate()*mat[i][j]*z[j] for i in range(6) for j in range(6))
        check(val/(2*(1+a*a+b*b)**2*norm)==Q((score/denom).real),'physical coefficient')


def compute(data):
    tables = json.loads((data/'system.json').read_text())['tables']
    records = json.loads((data/'orbits.json').read_text())['records']
    check(sum(r['invariant']['orbit_size'] for r in records)==935 and len(records)==62,'ray scope')
    a, b = poly({1:1}), poly({64:1})
    basis = form_basis(tables,a,b)
    certificates = {}
    for rec in records:
        if not rec['grouped_coefficient']:
            mat = balanced['combine'](basis,rec['ray'])
            certificates[str(rec['id'])] = certificate(rec['id'],mat,a*a*b*b)
    check(len(certificates)==52 and len(shifts)==11,'expected certificates')
    for a,b,seed in ((1,2,1),(3,1,2),(0,2,3)):
        direct(tables,a,b,seed)
    return certificates


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ray-data',type=Path,default=ROOT/'split-dual-rays')
    parser.add_argument('--print-certificate',action='store_true')
    args = parser.parse_args()
    result = compute(args.ray_data)
    if args.print_certificate:
        print(json.dumps(result,indent=2))
    else:
        check(result==json.loads((HERE/'polynomials.json').read_text()),'certificate fingerprints')
        print('All 52 non-coefficient orbits pass for every nonnegative projection triple.')
        print(sum(r['checked_polynomials'] for r in result.values()),'exact polynomial checks across the 52 orbits.')
        print('Three nonreal determinant frames, including a boundary projection, agree on all 33 coefficients.')
