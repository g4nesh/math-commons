#!/usr/bin/env python3
"""Exact continuation of a fixed-support three-variable determinant certificate."""
from fractions import Fraction as F
from itertools import combinations, permutations
import importlib.util
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'four-cycle-spectral'/'spectral.py'
spec=importlib.util.spec_from_file_location('moving_source',SOURCE)
s=importlib.util.module_from_spec(spec);sys.modules[spec.name]=s;spec.loader.exec_module(s)
PARENT=json.loads((ROOT/'four-cycle-spectral'/'fiber_certificate.json').read_text())
SUPPORT=tuple(tuple(row['columnToRow']) for row in PARENT['weightsAtFixedFourthEigenvalue'])
BASE_WEIGHTS=tuple(F(row['mass']) for row in PARENT['weightsAtFixedFourthEigenvalue'])
SUBSETS=tuple(chosen for k in range(4) for chosen in combinations(range(3),k))


def coefficient_system(z,support=SUPPORT):
    base=s.matrix_sum([s.Z,s.Z,s.Z,z])
    rows=[[F(1)]*len(support)];rhs=[F(1)]
    for chosen in SUBSETS:
        remaining=[i for i in range(4) if i not in chosen]
        minor=[[base[i][j] for j in remaining] for i in remaining]
        target=s.leibniz(minor)
        assert target==s.e.determinant(minor)
        points=[s.product(s.B_SPECTRUM[j]+(z if p[j]==3 else s.Z)
                          for j in range(4) if p[j] not in chosen) for p in support]
        for part in ('real','imag'):
            rows.append([getattr(point,part) for point in points]);rhs.append(getattr(target,part))
    return rows,rhs


def solve(rows,rhs):
    n=len(rows[0]);r=[list(map(F,row))+[F(b)] for row,b in zip(rows,rhs)]
    labels=list(range(len(r)));rank=0;pivots=[];basis=[]
    for col in range(n):
        pivot=next((i for i in range(rank,len(r)) if r[i][col]),None)
        if pivot is None:continue
        r[rank],r[pivot]=r[pivot],r[rank];labels[rank],labels[pivot]=labels[pivot],labels[rank]
        factor=r[rank][col];r[rank]=[x/factor for x in r[rank]]
        for i in range(len(r)):
            if i!=rank:
                factor=r[i][col];r[i]=[x-factor*y for x,y in zip(r[i],r[rank])]
        pivots.append(col);basis.append(labels[rank]);rank+=1
    if any(all(x==0 for x in row[:-1]) and row[-1]!=0 for row in r):
        return None,rank,basis
    w=[F(0)]*n
    for i,col in enumerate(pivots):w[col]=r[i][-1]
    return tuple(w),rank,basis


if __name__=='__main__':
    for t in (F(0),F(1,100),F(-1,100),F(1),F(-1)):
        rows,rhs=coefficient_system(s.Q(1,-1)+t)
        weights,rank,basis=solve(rows,rhs)
        print(json.dumps({'t':str(t),'rank':rank,'basisRows':basis,'consistent':weights is not None,
                          'minimumWeight':str(min(weights)) if weights else None,
                          'agreesAtBase':weights==BASE_WEIGHTS if t==0 else None}))
    rows,rhs=coefficient_system(s.Q(1,-1),tuple(permutations(range(4))))
    print('full rank',solve(rows,rhs)[1])


def ptrim(p):
    p=list(map(F,p))
    while len(p)>1 and p[-1]==0:p.pop()
    return tuple(p)


def padd(p,q):
    return ptrim([(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0) for i in range(max(len(p),len(q)))])


def pmul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):r[i+j]+=x*y
    return ptrim(r)


def pdiv(p,q):
    p=list(ptrim(p));q=ptrim(q)
    if q==(F(0),):raise ZeroDivisionError('zero polynomial')
    out=[F(0)]*max(1,len(p)-len(q)+1)
    while len(p)>=len(q) and any(p):
        power=len(p)-len(q);factor=p[-1]/q[-1];out[power]=factor
        for i,c in enumerate(q):p[power+i]-=factor*c
        p=list(ptrim(p))
    return ptrim(out),ptrim(p)


def pgcd(p,q):
    while q!=(F(0),):p,q=q,pdiv(p,q)[1]
    return tuple(v/p[-1] for v in p)


def interpolate(values):
    diffs=list(map(F,values));out=(F(0),);basis=(F(1),)
    for k in range(len(values)):
        out=padd(out,tuple(diffs[0]*v for v in basis))
        diffs=[diffs[i+1]-diffs[i] for i in range(len(diffs)-1)]
        basis=tuple(v/F(k+1) for v in pmul(basis,(F(-k),F(1))))
    return ptrim(out)


def real_det(matrix):
    a=[list(map(F,row)) for row in matrix];n=len(a);out=F(1)
    for col in range(n):
        pivot=next((i for i in range(col,n) if a[i][col]),None)
        if pivot is None:return F(0)
        if pivot!=col:a[col],a[pivot]=a[pivot],a[col];out=-out
        v=a[col][col];out*=v
        for i in range(col+1,n):
            factor=a[i][col]/v
            for j in range(col+1,n):a[i][j]-=factor*a[col][j]
    return out


def derive():
    from math import gcd
    _,rank,basis=solve(*coefficient_system(s.Q(1,-1)))
    assert rank==12
    determinants=[];values=[[]for _ in SUPPORT]
    # Every chosen square-system entry and RHS is affine in t. Cramer's
    # determinants have degree <=12, hence 13 exact interpolation values suffice.
    for t in range(13):
        rows,rhs=coefficient_system(s.Q(1,-1)+t)
        determinant=real_det([rows[i]for i in basis]);assert determinant
        w,_,_=solve(rows,rhs);assert w is not None
        determinants.append(determinant)
        for col,value in enumerate(w):values[col].append(value*determinant)
    rawden=interpolate(determinants);nums=[interpolate(v)for v in values]
    common=rawden
    for p in nums:common=pgcd(common,p)
    den,rem=pdiv(rawden,common);assert rem==(F(0),)
    nums=[pdiv(p,common)[0]for p in nums]
    scale=1
    for p in [den]+nums:
        for x in p:scale=scale*x.denominator//gcd(scale,x.denominator)
    ints=[[int(x*scale)for x in p]for p in [den]+nums]
    content=0
    for p in ints:
        for x in p:content=gcd(content,abs(x))
    if ints[0][0]<0:content=-content
    ints=[[x//content for x in p]for p in ints]
    den,nums=ints[0],ints[1:]
    rows0,rhs0=coefficient_system(s.Q(1,-1));rows1,rhs1=coefficient_system(s.Q(2,-1))
    for i,row in enumerate(rows0):
        left=(F(0),)
        for j,x in enumerate(row):left=padd(left,pmul((x,rows1[i][j]-x),nums[j]))
        assert left==pmul(den,(rhs0[i],rhs1[i]-rhs0[i]))
    return {'support':[list(p)for p in SUPPORT],'numerators':nums,'denominator':den,
            'basisRows':basis,'basisDeterminant':[str(x)for x in rawden]}
