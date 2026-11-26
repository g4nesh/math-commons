"""Exact finite replays of the transportation-table contact reduction."""
from fractions import Fraction as F
from itertools import product
from math import gcd, prod


def check(ok, why):
    if not ok:
        raise ValueError(why)


def tables(r,s):
    if not r:
        if not any(s):
            yield ()
        return
    for row in product(*(range(min(r[0],v)+1) for v in s)):
        if sum(row)==r[0]:
            for rest in tables(r[1:],tuple(x-y for x,y in zip(s,row))):
                yield (row,)+rest


def bracket(v,w):
    return v[0]*w[1]-v[1]*w[0]


def evaluate(m,left,right):
    return prod(bracket(v,w)**entry for v,row in zip(left,m) for w,entry in zip(right,row))


def combine(m,n,i=1,j=1):
    return tuple(tuple(i*x+j*y for x,y in zip(a,b)) for a,b in zip(m,n))


def chart(r,s,expected):
    family = tuple(tables(r,s))
    check(len(family)==expected,"number of distinct contact coordinates")
    p = [F(5+2*i) for i in range(len(r)-2)]
    q = [F(2+2*j) for j in range(len(s)-1)]
    left = [(0,1),(1,0)]+[(v,1) for v in p]
    right = [(1,1)]+[(w,1) for w in q]
    scale = F(17,13)
    values = {m:scale*evaluate(m,left,right) for m in family}
    def ratio(i,k,j,l):
        m = next(m for m in family if m[i][j] and m[k][l])
        n = [list(a) for a in m]
        n[i][j] -= 1
        n[k][l] -= 1
        n[i][l] += 1
        n[k][j] += 1
        return values[m]/values[tuple(tuple(row) for row in n)]
    recovered_q = [1/ratio(0,1,0,j) for j in range(1,len(s))]
    check(recovered_q==q,"right projective coordinates")
    for i in range(2,len(r)):
        h = ratio(0,i,0,1)
        check(q[0]*(h-1)/(h*q[0]-1)==p[i-2],"left projective coordinate")
    check(values[family[0]]/evaluate(family[0],left,right)==scale,"common scale")
    d = gcd(*r,*s)
    if d>1:
        small = tuple(tables(tuple(x//d for x in r),tuple(x//d for x in s)))
        zero = tuple(tuple(0 for _ in s) for _ in r)
        sums = {zero}
        for _ in range(d):
            sums = {combine(a,b) for a in sums for b in small}
        check(sums==set(family),"integer decomposition of all full tables")
        for a,b in product(small,repeat=2):
            x = evaluate(a,left,right)
            y = evaluate(b,left,right)
            ca = evaluate(combine(a,b,d,0),left,right)
            cb = evaluate(combine(a,b,d-1,1),left,right)
            cc = evaluate(combine(a,b,d-2,2),left,right)
            check(ca*cc==cb*cb and ca+cc-2*cb==x**(d-2)*(x-y)**2,"Veronese line test")
    return len(r)+len(s)-2


def signed_example():
    perms = [(0,1,2),(1,2,0),(2,0,1),(0,2,1),(1,0,2),(2,1,0)]
    def evaluate_at(t):
        a,b = (t/10,t,1+t),(0,1,1+2*t)
        p = [prod(a[i]+b[j] for i,j in enumerate(perm)) for perm in perms]
        coords = {(i,j):p[i]*p[j] for i in range(6) for j in range(i,6)}
        low = min(coords.values())
        contacts = {m for m,val in coords.items() if val==low}
        gap = (2*coords[0,1]-coords[1,1]-low)/sum(z-low for z in coords.values())
        check(gap<0,"expected signed-data violation")
        return contacts,gap,p
    start,g0,_ = evaluate_at(F(1))
    end,g1,p = evaluate_at(F(0))
    check(start=={(0,0)} and start<end and len(end)==18,"strict contact growth")
    check(p==[0,1,0,0,0,1] and g1==-F(1,3) and g1<g0,"rigid endpoint and normalized gap")
    # The first moment of 2*(P0+P1)-2*P1 is 2*P0, a transportation table.
    for r,c in product(range(3),repeat=2):
        moment = 2*(int(perms[0][r]==c)+int(perms[1][r]==c))-2*int(perms[1][r]==c)
        check(moment==2*int(r==c),"admissible first moments")
    return g0,g1


def main():
    cases = [((1,)*6,(2,2,2),90),((2,2,2),(2,2,2),21),
             ((3,3),(3,3),4),((1,1),(1,1),2)]
    bounds = [chart(*args) for args in cases]
    check(bounds==[7,4,2,2],"pattern contact bounds")
    g0,g1 = signed_example()
    print(f"Pattern contact bounds {bounds}; exact charts, layer decompositions and line tests pass.")
    print(f"Admissible signed example: 1 to 18 contact tables; normalized gap {g0} to {g1}.")


if __name__ == "__main__":
    main()
