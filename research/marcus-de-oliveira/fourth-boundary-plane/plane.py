#!/usr/bin/env python3
"""Exact two-dimensional circular boundary for a fixed three-variable fiber."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import importlib.util,sys,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('plane_boundary',HERE.parent/'fourth-boundary'/'boundary.py')
b=importlib.util.module_from_spec(spec);sys.modules[spec.name]=b;spec.loader.exec_module(b)
m=b.m;P=b.P
RAY=json.loads((HERE.parent/'fourth-boundary'/'positive.json').read_text())['bases'][0]
ROWS=json.loads((HERE.parent/'fourth-boundary'/'dual.json').read_text())['basisRows']
SUPPORT=tuple(tuple(p)for p in RAY['support']);IDS=tuple(P.index(p)for p in SUPPORT)
SI=(F(1,10),F(31,200));TI=(F(-1,50),F(1,50))
Q={(0,0):F(-350),(0,1):F(1079),(0,2):F(1808),(1,0):F(2187),(2,0):F(1808)}


def clean(p):return {k:F(v)for k,v in p.items()if v}
def add(p,q):
    r=dict(p)
    for k,v in q.items():r[k]=r.get(k,F(0))+v
    return clean(r)
def scale(p,a):return clean({k:v*a for k,v in p.items()})
def mul(p,q):
    r={}
    for(i,j),a in p.items():
        for(k,l),v in q.items():r[i+k,j+l]=r.get((i+k,j+l),F(0))+a*v
    return clean(r)
def evaluate(p,s,t):return sum(v*s**i*t**j for(i,j),v in p.items())
def decode(p):
    if not isinstance(p,list):raise ValueError('Polynomial must be a term list')
    out={}
    for term in p:
        if not isinstance(term,list)or len(term)!=3:raise ValueError('Wrong polynomial term')
        i,j,v=term
        if type(i)is not int or type(j)is not int or min(i,j)<0 or i+j>12 or type(v)is not str:raise ValueError('Invalid exact polynomial term or degree')
        v=F(v)
        if not v or(i,j)in out:raise ValueError('Zero or duplicate polynomial term')
        out[i,j]=v
    return out

def bernstein(p,si=SI,ti=TI):
    nx=max((i for i,j in p),default=0);ny=max((j for i,j in p),default=0);left,right=si;bottom,top=ti
    power={}
    for i in range(nx+1):
        for j in range(ny+1):
            power[i,j]=sum(v*comb(a,i)*left**(a-i)*(right-left)**i*comb(bb,j)*bottom**(bb-j)*(top-bottom)**j for(a,bb),v in p.items()if a>=i and bb>=j)
    return [sum(power[k,l]*F(comb(i,k),comb(nx,k))*F(comb(j,l),comb(ny,l))for k in range(i+1)for l in range(j+1))for i in range(nx+1)for j in range(ny+1)]

def system(s,t):return m.coefficient_system(m.s.Q(s+t,t-s),P)
def raw_values(s,t):
    A,rhs=system(s,t);C=[[A[i][j]for j in IDS]for i in ROWS];r=[rhs[i]for i in ROWS]
    d=m.real_det(C);out=[d]
    for j in range(12):out.append(m.real_det([row[:j]+[v]+row[j+1:]for row,v in zip(C,r)]))
    for j in range(24):
        # Repeating another selected column gives an identically zero determinant.
        if j in IDS:out.append(d if j==IDS[1]else F(0))
        else:out.append(m.real_det([row[:1]+[A[i][j]]+row[2:]for row,i in zip(C,ROWS)]))
    return out


def load_data():return json.loads((HERE/'certificate.json').read_text())
def parameters(s,t):
    if any(isinstance(v,bool)or not isinstance(v,(int,F))for v in(s,t)):raise TypeError('Use exact integer or Fraction coordinates')
    s,t=F(s),F(t)
    if not SI[0]<=s<=SI[1]or not TI[0]<=t<=TI[1]:raise ValueError('Outside the certified rectangle')
    return s,t

def weights_at(s,t):
    s,t=parameters(s,t)
    if evaluate(Q,s,t)<0:raise ValueError('Inside the excluded circle region')
    data=load_data();d=evaluate(decode(data['denominator']),s,t)
    if d<=0:raise ValueError('Certificate denominator is not positive')
    weights=tuple((p,evaluate(decode(n),s,t)/d)for p,n in zip(SUPPORT,data['numerators']))
    if any(w<0 for p,w in weights)or sum(w for p,w in weights)!=1:raise ValueError('Invalid convex weights')
    return weights


def boundary_bracket(t,bits=60):
    _,t=parameters(SI[0],t)
    if type(bits)is not int or not 1<=bits<=200:raise ValueError('Bisection steps must be between1 and200')
    lo,hi=SI
    for _ in range(bits):
        mid=(lo+hi)/2
        if evaluate(Q,mid,t)<0:lo=mid
        else:hi=mid
    return lo,hi


def verify(data=None):
    data=load_data()if data is None else data
    if data['support']!=[list(p)for p in SUPPORT]or data['basisRows']!=ROWS or type(data['orientation'])is not int or data['orientation']!=1:raise ValueError('Archived basis/orientation changed')
    if data['sInterval']!=list(map(str,SI))or data['tInterval']!=list(map(str,TI)):raise ValueError('Certified domain changed')
    if decode(data['boundaryPolynomial'])!=Q:raise ValueError('Circle polynomial changed')
    if len(data['numerators'])!=12 or len(data['dualColumns'])!=24:raise ValueError('Missing primal or dual columns')
    D=decode(data['denominator']);N=list(map(decode,data['numerators']));H=list(map(decode,data['dualColumns']))
    K=mul(mul({(1,0):7,(0,1):-1},{(2,0):1,(0,2):1}),mul({(1,0):1,(0,1):3,(0,0):4},{(1,0):1,(0,1):3,(0,0):4}))
    K=mul(K,{(2,0):1,(1,0):3,(0,2):1,(0,0):1})
    if N[1]!=scale(mul(K,Q),F(57344,5)):raise ValueError('Critical circle factorization failed')
    if D!=scale(mul(K,{(2,0):1,(0,2):1,(0,1):1,(0,0):-1}),-280000000):raise ValueError('Denominator factorization failed')
    if min(bernstein(K))<=0:raise ValueError('Common factor is not positive')
    A0,b0=system(F(0),F(0));As,bs=system(F(1),F(0));At,bt=system(F(0),F(1))
    for i in range(17):
        left={}
        for j,p in zip(IDS,N):left=add(left,mul(clean({(0,0):A0[i][j],(1,0):As[i][j]-A0[i][j],(0,1):At[i][j]-A0[i][j]}),p))
        right=mul(clean({(0,0):b0[i],(1,0):bs[i]-b0[i],(0,1):bt[i]-b0[i]}),D)
        if left!=right:raise ValueError('Cleared multiaffine coefficient identity failed')
    primal_controls=[bernstein(D)]+[bernstein(p)for j,p in enumerate(N)if j!=1]
    if any(min(v)<=0 for v in primal_controls):raise ValueError('A noncritical mass or denominator has no positive Bernstein certificate')
    dual_controls=[bernstein(p)for p in H]
    if any(min(v)<0 for v in dual_controls):raise ValueError('A permutation cofactor has negative Bernstein control')
    # Every raw determinant has total degree <=12; the stored polynomials are
    # checked to have that degree bound. The13x13 grid proves their identities.
    polynomials=[D]+N+H
    for t in range(13):
        for s in range(13):
            actual=raw_values(F(s),F(t))
            if actual!=[evaluate(p,F(s),F(t))for p in polynomials]:raise ValueError('Exact cofactor interpolation identity failed')
    # The reduced quadratic identifies the whole curve, not just sample roots.
    ds={(1,0):F(3616),(0,0):F(2187)};dt={(0,1):F(3616),(0,0):F(1079)}
    if min(bernstein(ds))<=0 or min(bernstein(dt))<=0:raise ValueError('Circle is not monotone in both coordinates')
    if min(bernstein(add(ds,scale(dt,-2))))<=0 or min(bernstein(add(scale(dt,3),scale(ds,-1))))<=0:raise ValueError('Boundary slope bounds failed')
    if evaluate(Q,SI[0],TI[1])>=0 or evaluate(Q,SI[1],TI[0])<=0:raise ValueError('A transverse slice is missing its unique boundary')
    if m.pmul((4,1),(-350,2187,1808))!=b.ALPHA:raise ValueError('Earlier ray cubic does not reduce to this circle')
    if m.s.e.multiply(m.s.U,m.s.e.adjoint(m.s.U))!=m.s.e.identity(4):raise ValueError('Archived U is not unitary')
    m.s.source_check()
    return {'coordinates':'a3=s*(1-i)+t*(1+i)','sInterval':list(map(str,SI)),'tInterval':list(map(str,TI)),
            'boundaryPolynomial':'1808*(s^2+t^2)+2187*s+1079*t-350',
            'originalComplexCoordinates':'904*(x^2+y^2)+1633*x-554*y-350',
            'feasibleSide':'q>=0, including the circle','excludedSide':'q<0; uniform weights only',
            'positiveMassesOnBoundary':11,'zeroMassesOnBoundary':1,'coefficientIdentities':17,
            'cofactorPolynomials':37,'cofactorTotalDegreeBound':12,'exactGridPoints':169,
            'minimumPrimalBernsteinControl':str(min(min(v)for v in primal_controls)),
            'minimumNonzeroDualBernsteinControl':str(min(min(v)for p,v in zip(H,dual_controls)if p)),
            'zeroDualColumns':sum(not p for p in H),'boundarySlopeInterval':['-1/2','-1/3'],
            'transverseEndpointRootBrackets':[{ 't':str(t),'sBracket':list(map(str,boundary_bracket(t)))}for t in(TI[0],F(0),TI[1])],
            'scope':'Fixed archived U and noncollinear b; weights uniform over three other arbitrary complex eigenvalues. No pointwise MDO counterexample or classification outside the rectangle.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--report',type=Path);parser.add_argument('--t',type=F)
    args=parser.parse_args();report=verify()
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    if args.t is not None:report['requestedBoundaryBracket']=list(map(str,boundary_bracket(args.t)))
    print(json.dumps(report,indent=2))
