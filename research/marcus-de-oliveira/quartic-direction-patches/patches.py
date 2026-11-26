#!/usr/bin/env python3
"""Exact replay for nonnegative radial weights on complex angular patches."""
from fractions import Fraction as F
from itertools import permutations, product
from math import comb
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'four-cycle-spectral'))
import spectral as s
PERMS=tuple(permutations(range(4)))


def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col))for col in zip(*b)]for row in a]


def solve_many(a,b):
    n=len(a);width=len(b[0]);rows=[list(map(F,row))+list(map(F,right))for row,right in zip(a,b)]
    if len(b)!=n or any(len(row)!=n for row in a)or any(len(row)!=width for row in b):raise ValueError('Invalid solve dimensions')
    for j in range(n):
        k=next((i for i in range(j,n)if rows[i][j]),None)
        if k is None:raise ValueError('Singular exact matrix')
        rows[j],rows[k]=rows[k],rows[j];pivot=rows[j][j];rows[j]=[v/pivot for v in rows[j]]
        for i in range(n):
            if i!=j and rows[i][j]:
                ratio=rows[i][j];rows[i]=[x-ratio*y for x,y in zip(rows[i],rows[j])]
    answer=[row[n:]for row in rows]
    if matmul(a,answer)!=[[F(v)for v in row]for row in b]:raise ValueError('Exact solve residual')
    return answer


def pencil_data():
    s.source_check();rows=[]
    for x,y in product((0,1),repeat=2):
        samples=[]
        for z in (s.Z,s.O):
            diagonal=[s.Q(x),z,s.Q(y),s.Z]
            target=s.e.determinant(s.matrix_sum(diagonal))
            samples.append([s.vertex(diagonal,p)-target for p in PERMS])
        rows.append((samples[0],[b-a for a,b in zip(*samples)]))
    return rows


def primitive(centers):
    if len(centers)>2 or any(isinstance(r,bool)or not isinstance(r,int)for r in centers):raise ValueError('Expected at most two exact integer centers')
    poly=[F(1)]
    for root in centers:
        result=[F(0)]*(len(poly)+2)
        for i,x in enumerate(poly):
            for j,y in enumerate((root*root,-2*root,1)):result[i+j]+=x*y
        poly=result
    return poly+[F(0)]*(5-len(poly))


def matrices(terms,chart,center,data=None):
    if chart not in ('real','imag'):raise ValueError('Unknown angular chart')
    data=pencil_data()if data is None else data
    omega=s.Q(1,center)if chart=='real'else s.Q(center,1)
    derivative=s.Q(0,1)if chart=='real'else s.O
    columns=[[],[]]
    for term in terms:
        p=tuple(term['permutation'])
        if p not in PERMS:raise ValueError('Invalid permutation')
        pi=PERMS.index(p);poly=primitive(term['centers']);base=[];slope=[]
        for degree in range(6):
            for const,linear in data:
                value=const[pi]*s.Q(poly[degree])if degree<5 else s.Z
                deriv=s.Z
                if degree:
                    value+=linear[pi]*omega*s.Q(poly[degree-1])
                    deriv=linear[pi]*derivative*s.Q(poly[degree-1])
                base.append(value);slope.append(deriv)
        columns[0].append([v.real for v in base]+[v.imag for v in base]+[F(1)])
        columns[1].append([v.real for v in slope]+[v.imag for v in slope]+[F(0)])
    return [[list(row)for row in zip(*cols)]for cols in columns]


def bernstein(poly,left,right):
    n=len(poly)-1
    shifted=[sum(poly[j]*comb(j,i)*left**(j-i)*(right-left)**i for j in range(i,n+1))for i in range(n+1)]
    return [sum(shifted[i]*F(comb(k,i),comb(n,i))for i in range(k+1))for k in range(n+1)]


def horner(poly,x):
    value=F(0)
    for c in reversed(poly):value=value*x+c
    return value


def check_patch(patch,data=None):
    center=F(patch['center']);left,right=map(F,patch['interval'])
    if not F(1,2)<=left<right<=1:raise ValueError('Invalid closed angular interval')
    terms=patch['terms'];keys=[(tuple(t['centers']),tuple(t['permutation']))for t in terms]
    if not terms or len(set(keys))!=len(keys):raise ValueError('Missing or repeated support columns')
    denominator=list(map(F,patch['denominator']));numerators=[list(map(F,t['numerator']))for t in terms]
    if not denominator or any(not p for p in numerators):raise ValueError('Empty polynomial')
    bounds=[bernstein(p,left-center,right-center)for p in [denominator]+numerators]
    if min(min(b)for b in bounds)<=0:raise ValueError('Nonpositive Bernstein coefficient')
    m0,m1=matrices(terms,patch['chart'],center,data)
    degree=max(len(denominator),max(map(len,numerators))+1)
    for row in range(49):
        for power in range(degree):
            value=sum(m0[row][j]*(p[power]if power<len(p)else 0)+m1[row][j]*(p[power-1]if 0<power<=len(p)else 0)for j,p in enumerate(numerators))
            expected=denominator[power]if row==48 and power<len(denominator)else 0
            if value!=expected:raise ValueError('Polynomial identity failed')
    # Every radial primitive is a nonnegative product of real squares. A
    # positive constant or two squares with distinct roots certify Z>0.
    constant=sum(min(bounds[j+1])for j,t in enumerate(terms)if not t['centers'])
    squares={}
    for j,term in enumerate(terms):
        if len(term['centers'])==1:
            root=term['centers'][0];squares[root]=squares.get(root,F(0))+min(bounds[j+1])
    floor=constant
    for a,x in squares.items():
        for b,y in squares.items():
            if a<b:floor=max(floor,constant+x*y*(a-b)**2/(x+y))
    if floor<=0:raise ValueError('No certified positive radial denominator')
    return {'chart':patch['chart'],'interval':patch['interval'],'terms':len(terms),
            'denominatorDegree':len(denominator)-1,'numeratorDegree':max(map(len,numerators))-1,
            'radialDenominatorLowerBound':str(floor),'minimumBernsteinCoefficient':str(min(min(b)for b in bounds))}


def weights_at(patch,k,t):
    if any(isinstance(v,bool)or not isinstance(v,(int,F))for v in(k,t)):raise TypeError('Use exact integer or Fraction coordinates')
    k,t=F(k),F(t);left,right=map(F,patch['interval'])
    if not left<=k<=right:raise ValueError('Outside patch')
    masses=[F(0)]*24
    for term in patch['terms']:
        numerator=horner(list(map(F,term['numerator'])),k-F(patch['center']))
        masses[PERMS.index(tuple(term['permutation']))]+=numerator*horner(primitive(term['centers']),t)
    total=sum(masses)
    if min(masses)<0 or total<=0:raise ValueError('Invalid pointwise masses')
    return [v/total for v in masses]


def verify(certificate):
    if certificate.get('axis')!=1 or certificate.get('fixedFourth')!=['0','0']:raise ValueError('Wrong spectral slice')
    expected={('real',(F(1,2),F(9,16))),('real',(F(15,16),F(1)))}
    actual={(patch['chart'],tuple(map(F,patch['interval'])))for patch in certificate['patches']}
    if len(certificate['patches'])!=2 or actual!=expected:raise ValueError('The two stated angular patches must both be present')
    data=pencil_data();reports=[check_patch(p,data)for p in certificate['patches']]
    regions={}
    for patch in certificate['patches']:
        regions.setdefault(patch['chart'],[]).append(patch['interval'])
    return {'coveredRegions':regions,'patches':reports,
            'scope':'Two separate real-chart patches only; a4=0, a2=(1+i*k)t, t real, a1 and a3 arbitrary complex. No full angular-cover claim.'}



if __name__=='__main__':
    print(json.dumps(verify(json.loads((HERE/'certificate.json').read_text())),indent=2))
