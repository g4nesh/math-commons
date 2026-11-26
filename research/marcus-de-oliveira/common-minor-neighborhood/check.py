"""Exact edge-gap, curvature and outside-common-minor checks in dimension five."""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, prod
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
colors = runpy.run_path(str(ROOT/'three-color-coupling'/'coupling.py'))
Q, e = colors['Q'], colors['e']
A = (Q(-4,-4), Q(0,3), Q(2,-2), Q(-3,-2), Q(4,-3))
B = (Q(-1,1), Q(2,5), Q(4), Q(5,-5), Q(2,-1))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def cross(direction, value):
    return (direction.conjugate()*value).imag


def gap(left, right, z):
    return cross(right-left,z-left)


def hull(vertices):
    points = sorted(set(vertices),key=lambda z:(z.real,z.imag))
    def half(seq):
        chain = []
        for point in seq:
            while len(chain)>=2 and gap(chain[-2],chain[-1],point)<=0:
                chain.pop()
            chain.append(point)
        return chain[:-1]
    return half(points)+half(reversed(points))


def parameters(a,b):
    """Return a conservative rational threshold for the entrywise minor LP residual."""
    n=len(a)
    ps=list(permutations(range(n)))
    values=[prod(a[i]+b[j] for i,j in enumerate(p)) for p in ps]
    polygon=hull(values)
    check(len(polygon)>=3,'positive area is required')
    edges=list(zip(polygon,polygon[1:]+polygon[:1]))
    gaps=[]
    for left,right in edges:
        deficits=[gap(left,right,z) for z in values]
        check(min(deficits)==0 and deficits.count(0)==2,'supporting edge must have exactly two labels')
        length_upper=abs((right-left).real)+abs((right-left).imag)
        gaps.append(min(d for d in deficits if d>0)/length_upper)
    g=min(gaps)
    r=max(F(1),max(abs(z.real)+abs(z.imag) for z in a)+max(abs(z.real)+abs(z.imag) for z in b))
    rho=min(F(1),3*g/(16*n**3*r**n))
    eta=min(F(1,16),rho**2/(512*n))
    xi=min(eta/n,g*eta/(2*comb(2*n,n)*r**n))
    return ps,values,edges,g,xi


def add(a,b,sign=1):
    return [[x+sign*y for x,y in zip(row,col)] for row,col in zip(a,b)]


def commutator(a,b):
    return add(e.multiply(a,b),e.multiply(b,a),-1)


def derivatives(base,spectrum,x):
    """First and second coefficients of det(Da + exp(tX) base Db base* exp(-tX))."""
    n=len(x)
    conjugated=e.multiply(e.multiply(base,e.diagonal(spectrum)),e.adjoint(base))
    m=add(e.diagonal(A),conjugated)
    first=commutator(x,conjugated)
    second=commutator(x,first)
    def columns(changes):
        return [[changes.get(j,m)[i][j] for j in range(n)] for i in range(n)]
    linear=sum(e.determinant(columns({i:first})) for i in range(n))
    quadratic=sum(e.determinant(columns({i:second})) for i in range(n))/2
    quadratic+=sum(e.determinant(columns({i:first,j:first})) for i,j in combinations(range(n),2))
    return linear,quadratic


def curvature(ps,values,edges):
    n=len(A);id_perm=tuple(range(n))
    left,right=next((l,r) for l,r in edges if values[0] in (l,r)
                    and sum(x!=y for x,y in zip(ps[values.index(l)],ps[values.index(r)]))==2)
    other=ps[values.index(right if values[0]==left else left)]
    pair=[i for i in range(n) if other[i]!=id_perm[i]]
    base=colors['rotation'](n,*pair,F(2,5),Q(F(3,5),F(4,5)))
    blocks=[pair]+[[i] for i in range(n) if i not in pair]
    total=[[Q(0) for _ in A] for _ in A]
    coefficients=[]
    unscaled_gap=min(gap(left,right,z) for z in values if z not in (left,right))
    for a,b in combinations(blocks,2):
        x=[[Q(0) for _ in A] for _ in A]
        squared=0
        for i in a:
            for j in b:
                z=Q(i+j+1,i-j)
                x[i][j],x[j][i]=z,-z.conjugate()
                squared+=colors['abs2'](z)
        total=add(total,x)
        first,second=derivatives(base,B,x)
        check(cross(right-left,first)==0,'linear normal term')
        quad=cross(right-left,second)
        check(quad>=unscaled_gap*squared,'transverse lower bound')
        coefficients.append(quad)
    all_first,all_second=derivatives(base,B,total)
    check(cross(right-left,all_first)==0 and cross(right-left,all_second)==sum(coefficients),
          'mixed normal terms did not cancel')
    return len(coefficients)


def obstruction(u):
    def sq(r,c):
        return colors['abs2'](colors['minor'](u,r,c))
    return (sq((0,2),(2,3))-sq((0,),(2,))-sq((3,),(1,))-sq((0,1),(1,3))
            -sum(1-sq((i,),(i,)) for i in range(4,len(u))))


def four_cycle(n, t):
    u=e.identity(n)
    for i,j in [(0,1),(1,2),(2,3),(0,3)]:
        u=e.multiply(colors['rotation'](n,i,j,t),u)
    return u


def cycle_gap(t):
    return 16*t**4*(t**4-6*t**2+1)/(1+t*t)**6


def original(u):
    return e.determinant(add(e.diagonal(A),e.multiply(e.multiply(u,e.diagonal(B)),e.adjoint(u))))


def main():
    n=len(A)
    ps,values,edges,g,xi=parameters(A,B)
    count=curvature(ps,values,edges)
    # Rational half-angle rotations; the last gate connects the fifth site.
    t=min(F(1,128),xi/(24*n))
    # After clearing (1+t^2)^16 both polynomials have degree at most 32.
    for value in range(33):
        parameter=F(value)
        check(obstruction(four_cycle(4,parameter))==cycle_gap(parameter),'polynomial separator identity')
    u=four_cycle(n,t)
    check(obstruction(u)==cycle_gap(t) and cycle_gap(t)>0,'four-gate separator')
    u=e.multiply(colors['rotation'](n,3,4,t**5,Q(0,1)),u)
    check(e.multiply(u,e.adjoint(u))==e.identity(n),'exact unitarity')
    check(2*n*(8*t+2*t**5)<=xi,'near-identity residual threshold')
    check(obstruction(u)>0,'new unitary must lie outside the common-minor locus')
    for p in ps:
        pm=[[Q(int(p[j]==i)) for j in range(n)] for i in range(n)]
        check(obstruction(pm)<=0,'invalid permutation separator')
    # A connected nonzero support graph excludes a common reducing split for
    # these two simple spectra. Row and column vertices have separate labels.
    reached={0}
    for _ in range(2*n):
        for i,j in ((i,j) for i in range(n) for j in range(n) if u[i][j]!=Q(0)):
            if i in reached or n+j in reached:
                reached.update((i,n+j))
    check(len(reached)==2*n,'irreducible fixture')
    value=original(u)
    check(min(gap(l,r,value) for l,r in edges)>0,'direct polygon verification')
    print('Exact five-dimensional hull:',len(edges),'edges;',count,'normal block checks.')
    print('Unit-normal edge-gap lower bound:',g)
    print('Certified minor residual threshold:',xi)
    print('The irreducible five-site fixture is strictly inside every hull edge and violates the common-minor separator.')


if __name__=='__main__':
    main()
