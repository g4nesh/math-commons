#!/usr/bin/env python3
"""Exact, deterministic exploratory checks. Standard library; no repository imports."""
from fractions import Fraction as F
from itertools import combinations, permutations
import json

def q(x): return (F(x), F(0))
ZERO, ONE = q(0), q(1)
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))
def mul(x,y): return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def conj(x): return (x[0],-x[1])
def norm(x): return x[0]**2+x[1]**2
def div(x,y): return tuple(t/norm(y) for t in mul(x,conj(y)))
def scale(t,x): return (t*x[0],t*x[1])
def prod(xs):
    out=ONE
    for x in xs: out=mul(out,x)
    return out
def det(M):
    M=[row[:] for row in M]; out=ONE
    for j in range(len(M)):
        pivot=next((i for i in range(j,len(M)) if M[i][j]!=ZERO),None)
        if pivot is None: return ZERO
        if pivot!=j: M[j],M[pivot]=M[pivot],M[j]; out=neg(out)
        x=M[j][j]; out=mul(out,x)
        for i in range(j+1,len(M)):
            r=div(M[i][j],x)
            for k in range(j+1,len(M)): M[i][k]=sub(M[i][k],mul(r,M[j][k]))
    return out
def weights(p,c,lam):
    n=len(p); ident=tuple(range(n)); ws={ident:F(1)}
    for i,j in combinations(range(n),2):
        w=c*p[i]*p[j]*(1-lam*(1-p[i]-p[j]))
        perm=list(ident);perm[i],perm[j]=perm[j],perm[i]
        ws[tuple(perm)]=w;ws[ident]-=w
    for i,j,k in combinations(range(n),3):
        w=lam*c*p[i]*p[j]*p[k]
        for a,b,d in [(i,j,k),(i,k,j)]:
            perm=list(ident);perm[a],perm[b],perm[d]=b,d,a
            ws[tuple(perm)]=w;ws[ident]-=w
    assert sum(ws.values())==1
    return ws
def minor_formula(p,c,I,J):
    a,b=set(I)-set(J),set(J)-set(I)
    if not a: s=sum(p[i] for i in I);return 1-c*s*(1-s)
    if len(a)==1:return c*p[next(iter(a))]*p[next(iter(b))]
    return F(0)
def compound_checks(p,c,ws,U=None):
    count=0;n=len(p)
    for k in range(n+1):
        sets=list(combinations(range(n),k))
        for I in sets:
            for J in sets:
                mass=sum(w for perm,w in ws.items() if {perm[j] for j in J}==set(I))
                expected=minor_formula(p,c,I,J)
                assert mass==expected,(p,c,I,J,mass,expected)
                if U is not None: assert norm(det([[U[i][j] for j in J] for i in I]))==expected
                count+=1
    return count
def one_case(n,z,r=None):
    r=r if r is not None else [(F(i%3+1),F((2*i+1)%3-1)) for i in range(n)]
    R=sum(map(norm,r));p=[norm(x)/R for x in r];c=norm(sub(z,ONE))
    U=[[add(ONE if i==j else ZERO,scale(1/R,mul(sub(z,ONE),mul(r[i],conj(r[j]))))) for j in range(n)] for i in range(n)]
    a=[(F(i-2),F((i*i+1)%5-2)) for i in range(n)]
    b=[(F(3-i),F((3*i+2)%7-3)) for i in range(n)]
    if n%2==0:b[0]=neg(a[0])  # Exercise singular D_a+D_b without division.
    C=[[add(a[i] if i==j else ZERO,sum_complex(mul(mul(U[i][k],b[k]),conj(U[j][k])) for k in range(n))) for j in range(n)] for i in range(n)]
    P=lambda perm:prod(add(a[i],b[perm[i]]) for i in range(n))
    ws=weights(p,c,F(1));identity=tuple(range(n));value=det(C)
    assert sum_complex(scale(w,P(perm)) for perm,w in ws.items())==value
    signed=weights(p,c,F(0));assert sum_complex(scale(w,P(perm)) for perm,w in signed.items())==value
    if c<=3: assert min(ws.values())>=0
    count=compound_checks(p,c,ws,U)
    positive=sorted(x for x in p if x>0)
    if len(positive)>=3:
        lam=1/(1-positive[0]-positive[1]); improved=weights(p,c,lam)
        assert all(w>=0 for perm,w in improved.items() if perm!=identity)
        count+=compound_checks(p,c,improved)
    # Phase interpolation against the half-turn, in the exact permutation formula.
    H=weights(p,F(4),F(0)); hv=sum_complex(scale(w,P(perm)) for perm,w in H.items())
    assert value==add(scale(1-c/4,P(identity)),scale(c/4,hv))
    return count
def sum_complex(xs):
    out=ZERO
    for x in xs:out=add(out,x)
    return out

def main():
    count=0;cases=0
    for n in range(1,7):
        for z in [ONE,(F(0),F(1)),(F(-1),F(0)),(F(3,5),F(4,5)),(F(-3,5),F(4,5))]:
            count+=one_case(n,z);cases+=1
    for z in [ONE,(F(0),F(1)),(F(-1),F(0))]:
        count+=one_case(5,z,[q(0),q(1),(F(0),F(2)),q(0),q(3)]);cases+=1
    flat=[]
    for n in range(3,8):
        p=[F(1,n)]*n;c=F(3*n,n-1);lam=F(n,n-2)
        ws=weights(p,c,lam)
        assert min(ws.values())>=0
        assert ws[tuple(range(n))]==0
        assert all(w==0 for perm,w in ws.items() if sum(i!=j for i,j in enumerate(perm))==2)
        count+=compound_checks(p,c,ws)
        flat.append({'n':n,'formal_compound_threshold':str(c),'physical_threshold':str(min(F(4),c)),'half_turn_expected_moved_points':str(4*(1-sum(x*x for x in p))),'half_turn_obstructed':4*(1-sum(x*x for x in p))>3})
    # General singleton/permutation support argument: every permutation moving >=4
    # points has a subset with >=2 elements moved out. Exhaustive n<=7 sanity check.
    permutations_checked=0
    for n in range(1,8):
        for perm in permutations(range(n)):
            if sum(i!=j for i,j in enumerate(perm))<4:continue
            assert any(len(set(I)-{perm[j] for j in I})>=2 for I in combinations(range(n),2))
            permutations_checked+=1
    print(json.dumps({'status':'passed','exact_unitary_cases':cases,'compound_entries_checked':count,'permutations_with_support_at_least_four_checked':permutations_checked,'flat_cases':flat},indent=2))
if __name__=='__main__':main()
