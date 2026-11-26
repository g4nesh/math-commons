"""Exact C3-averaged six-site coloring model (standard library only)."""
from fractions import Fraction as F
from itertools import combinations, product

SITES=tuple(range(6))
SUBSETS=tuple(s for k in range(7) for s in combinations(SITES,k))
COLORS=tuple(c for c in product(range(3),repeat=6) if all(c.count(j)==2 for j in range(3)))
ROWS=tuple((s,h) for s in SUBSETS for h in product(range(3),repeat=3) if sum(h)==len(s))
ORBIT_REPS=tuple(c for c in COLORS if c[0]==0)
SIGNS=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
FACETS=tuple([(tuple(sign if k==j else 0 for k in range(4)),F(1)) for j in range(3) for sign in [-1,1]]+
             [(tuple(d)+(sign,),F(3)) for d in SIGNS for sign in [-1,1]]+
             [(tuple(-3*x for x in d)+(2*sign,),F(9)) for d in SIGNS for sign in [-1,1]])


def require(condition,message):
    if not condition: raise ValueError(message)


def parity(seq):
    return (-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))


def histogram(c,s): return tuple(sum(c[i]==j for i in s) for j in range(3))


def coefficients(s,h):
    """Coefficients of (1,d01,d02,d12,kappa), including complement rows."""
    if len(s)>3: return coefficients(tuple(i for i in SITES if i not in s),tuple(2-i for i in h))
    q=[F(0)]*5
    if len(s)==0: q[0]=F(1)
    elif len(s)==1: q[0]=F(1,3)
    elif len(s)==2:
        groups=tuple(i//2 for i in s)
        repeated=2 in h
        if groups[0]==groups[1]: q[0]=F(1,9) if repeated else F(2,9)
        else:
            sign=(-1)**sum(i%2 for i in s)
            q[0]=F(1,18) if repeated else F(5,18)
            q[1+((0,1),(0,2),(1,2)).index(groups)]=F(-sign if repeated else sign,18)
    else:
        if len({i//2 for i in s})==2: q[0]=F(1,3) if h==(1,1,1) else F(1,9)
        else:
            signs=[(-1)**(i%2) for i in s]
            if h==(1,1,1): q[0]=F(1,2);scale=F(1,6)
            else: q[0]=F(1,12);scale=-F(1,36);q[4]=F(parity(h)*signs[0]*signs[1]*signs[2],36)
            for j,(a,b) in enumerate(((0,1),(0,2),(1,2))):q[j+1]=scale*signs[a]*signs[b]
    return tuple(q)


def target(point):
    z=(F(1),)+tuple(point)
    return tuple(sum(a*b for a,b in zip(coefficients(s,h),z)) for s,h in ROWS)


def orbit_column(c):
    orbit=[tuple((i+k)%3 for i in c) for k in range(3)]
    return tuple(F(sum(histogram(o,s)==h for o in orbit),3) for s,h in ROWS)


def solve(matrix,rhs):
    """Unique exact square solve; None exactly when singular."""
    n=len(rhs);a=[[F(x) for x in row]+[F(y)] for row,y in zip(matrix,rhs)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return None
        a[j],a[pivot]=a[pivot],a[j]
        z=a[j][j];a[j]=[x/z for x in a[j]]
        for i in range(n):
            if i!=j:
                z=a[i][j]
                if z:a[i]=[x-z*y for x,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a)


def vertices():
    answer=set()
    for selected in combinations(FACETS,4):
        v=solve([a for a,b in selected],[b for a,b in selected])
        if v is not None and all(sum(x*y for x,y in zip(a,v))<=b for a,b in FACETS):answer.add(v)
    return tuple(sorted(answer))
