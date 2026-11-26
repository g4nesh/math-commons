"""An exact all-orientations certificate at two different projection triples."""
import json
import runpy
from fractions import Fraction as F
from itertools import permutations
from math import comb
from pathlib import Path

core = runpy.run_path(str(Path(__file__).resolve().with_name('check.py')))
prev,check,star = core['prev'],core['check'],core['star']


def interval_test(det, indices, split, lower):
    for k in range(1,len(indices)+1):
        coefficients = prev['terms'](det(indices[:k]))
        check(all(power%2==0 for power in coefficients),'unexpected odd angle term')
        if not coefficients or max(coefficients)!=2*k:
            return False
        p = [coefficients.get(2*i,0) for i in range(k+1)]
        if lower:
            # (1+y)^k P(split*y/(1+y)) for y>=0
            v = [sum(p[i]*split**i*comb(k-i,j-i) for i in range(j+1)) for j in range(k+1)]
        else:
            # P(split+y) for y>=0
            v = [sum(p[i]*split**(i-j)*comb(i,j) for i in range(j,k+1)) for j in range(k+1)]
        if min(v)<=0:
            return False
    return True


def verify(tables,records,alpha,beta):
    poly = prev['poly']
    t = poly({1:1})
    tests = 0
    for pivot,j,k in permutations(range(3)):
        a,b,c,d = [poly({0:value}) for value in (F(alpha[j],alpha[pivot]),F(beta[j],beta[pivot]),
                                                F(alpha[k],alpha[pivot]),F(beta[k],beta[pivot]))]
        raw = core['forms'](tables,a,b,c,d,t)
        r0,r1 = 1+a*a+c*c,1+b*b+d*d
        unit = a*b*c*d*(r0+r1*t*t)
        split = r0.terms[0]/r1.terms[0]
        for rec in records:
            orbit = rec['id']
            if rec['grouped_coefficient'] or orbit==3:
                continue
            m = core['balanced']['combine'](raw,rec['ray'])
            scale,preferred = core['shifts'].get(orbit,(0,0))
            for sign in ((1,-1) if scale else (1,)):
                adj = [[m[u][v]+unit*(scale*sign*star[u][v]) for v in range(6)] for u in range(6)]
                det = prev['determinant_factory'](adj)
                order = (0,2,3,5,1,4)
                for lower in (True,False):
                    pins = dict.fromkeys((preferred,*order)) if sign==-1 else (None,)
                    check(any(interval_test(det,tuple(u for u in order if u!=pin),split,lower)
                              for pin in pins),'no certified phase chart')
                    tests += 1
    return tests


if __name__=='__main__':
    data = core['ROOT']/'split-dual-rays'
    tables = json.loads((data/'system.json').read_text())['tables']
    records = json.loads((data/'orbits.json').read_text())['records']
    count = verify(tables,records,(1,1,3),(1,2,4))
    print(count,'positive half-interval matrix tests cover all singleton orientations and complement flags.')
