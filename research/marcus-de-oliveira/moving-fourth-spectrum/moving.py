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


def evaluate(p,t):
    out=F(0)
    for coefficient in reversed(p):out=out*t+coefficient
    return out


def bernstein(p,left,right):
    """Exact coefficients in B_i^n(x), t=left+(right-left)*x."""
    from math import comb
    n=len(p)-1
    power=[sum(F(p[k])*comb(k,j)*left**(k-j)*(right-left)**j for k in range(j,n+1))for j in range(n+1)]
    return tuple(sum(power[j]*F(comb(i,j),comb(n,j))for j in range(i+1))for i in range(n+1))


INTERVAL=(F(-3,4),F(2))
DENOMINATOR=(9765625,15625000,6640625,781250)


def weights_at(data,t):
    if isinstance(t,bool) or not isinstance(t,(int,F)):
        raise TypeError('Parameter must be an exact integer or Fraction')
    if not INTERVAL[0]<=t<=INTERVAL[1]:
        raise ValueError('Outside the certified closed interval [-3/4,2]')
    return tuple(evaluate(p,t)/evaluate(data['denominator'],t)for p in data['numerators'])


def verify(data):
    s.source_check()
    if s.e.multiply(s.U,s.e.adjoint(s.U))!=s.e.identity(4) or s.e.multiply(s.e.adjoint(s.U),s.U)!=s.e.identity(4):
        raise ArithmeticError('Archived matrix is not exactly unitary')
    if s.cross(s.B_SPECTRUM[0]-s.B_SPECTRUM[1],s.B_SPECTRUM[3]-s.B_SPECTRUM[1])!=-5:
        raise ArithmeticError('The stated noncollinear second spectrum changed')
    if data['support']!=[list(p)for p in SUPPORT] or len(data['numerators'])!=12:
        raise ValueError('Wrong permutation support')
    polynomials=[data['denominator']]+data['numerators']
    if any(not isinstance(p,list) or not p or any(type(x)is not int for x in p) for p in polynomials):
        raise ValueError('Expected exact integer polynomial coefficients')
    if tuple(data['denominator'])!=DENOMINATOR:
        raise ValueError('Unexpected certificate denominator')
    if pmul(pmul((1,1),(5,2)),(5,1))!=tuple(F(x,390625)for x in DENOMINATOR):
        raise ArithmeticError('Denominator factorization failed')
    rows0,rhs0=coefficient_system(s.Q(1,-1));rows1,rhs1=coefficient_system(s.Q(2,-1))
    for i,row in enumerate(rows0):
        left=(F(0),)
        for j,x in enumerate(row):left=padd(left,pmul((x,rows1[i][j]-x),data['numerators'][j]))
        if left!=pmul(data['denominator'],(rhs0[i],rhs1[i]-rhs0[i])):
            raise ArithmeticError('A cleared coefficient identity failed')
    controls=[bernstein(p,*INTERVAL)for p in data['numerators']]
    if any(min(c)<=0 for c in controls):
        raise ArithmeticError('Bernstein coefficients do not prove strict positivity')
    if weights_at(data,0)!=BASE_WEIGHTS:
        raise ArithmeticError('Base point does not recover archived masses')
    # Check the recorded Cramer minor polynomial using its proven degree bound.
    if len(set(data['basisRows']))!=12 or any(type(i)is not int or not 0<=i<17 for i in data['basisRows']):
        raise ValueError('Invalid basis row indices')
    determinant_polynomial=tuple(F(x)for x in data['basisDeterminant'])
    if len(determinant_polynomial)>13:raise ValueError('Basis determinant exceeds degree bound12')
    for t in range(13):
        rows,_=coefficient_system(s.Q(1,-1)+t)
        if evaluate(determinant_polynomial,t)!=real_det([rows[i]for i in data['basisRows']]):
            raise ArithmeticError('Recorded basis minor polynomial failed exact interpolation check')
    # Each direct evaluation checks a different route than coefficient matching.
    from itertools import product
    endpoint_checks=[]
    for t in (INTERVAL[0],F(0),INTERVAL[1]):
        w=weights_at(data,t);z=s.Q(1,-1)+t
        for corner in product((0,1),repeat=3):
            a=[s.Q(bit)for bit in corner]+[z]
            actual=s.leibniz(s.matrix_sum(a))
            if actual!=s.e.determinant(s.matrix_sum(a)):raise ArithmeticError('Determinants disagree')
            predicted=sum((mass*s.vertex(a,p)for mass,p in zip(w,SUPPORT)),s.Z)
            if actual!=predicted:raise ArithmeticError('Corner reconstruction failed')
        endpoint_checks.append({'t':str(t),'a3':z.encode(),'weights':list(map(str,w)),'minimumWeight':str(min(w))})
    singular=[]
    for t in (F(-1),F(-5,2),F(-5)):
        rows,rhs=coefficient_system(s.Q(1,-1)+t)
        solution,rank,_=solve(rows,rhs)
        _,augmented_rank,_=solve([row+[b]for row,b in zip(rows,rhs)],[0]*len(rows))
        singular.append({'t':str(t),'a3':(s.Q(1,-1)+t).encode(),'coefficientRank':rank,
                         'augmentedRank':augmented_rank,'specifiedSupportConsistent':solution is not None})
    return {'parameter':'a3=1-i+t, t real; a0,a1,a2 arbitrary complex',
            'closedInterval':list(map(str,INTERVAL)),'realCoefficientIdentities':17,
            'positiveMasses':12,'bernsteinCoefficients':[list(map(str,c))for c in controls],
            'minimumBernsteinCoefficients':[str(min(c))for c in controls],
            'denominatorFactorization':'390625*(t+1)*(2*t+5)*(t+5)',
            'baseMassesRecovered':True,'coefficientDegreeBound':12,'basisMinorValuesChecked':13,
            'exactEndpointAndBaseChecks':endpoint_checks,'specifiedSupportAtFormulaPoles':singular,
            'scope':'Fixed archived U,b; uniform over three complex eigenvalues on the stated fourth-eigenvalue segment. No general MDO theorem or pointwise counterexample.'}


def main():
    import argparse,resource
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--derive',action='store_true',help='Recreate rational polynomial weights by exact linear algebra; no optimizer')
    parser.add_argument('--check',type=Path,default=Path(__file__).with_name('derived.json'))
    parser.add_argument('--output',type=Path,help='Write derived coefficient data (with --derive)')
    parser.add_argument('--report',type=Path,help='Write version-independent verification details')
    args=parser.parse_args();resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    data=derive()if args.derive else json.loads(args.check.read_text())
    report=verify(data)
    if args.output:
        if not args.derive:parser.error('--output requires --derive')
        args.output.write_text(json.dumps(data,indent=2)+'\n')
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ('closedInterval','realCoefficientIdentities','positiveMasses','baseMassesRecovered','specifiedSupportAtFormulaPoles')}))

if __name__=='__main__':main()
