"""Exact paired Pfaffians, classical weights, perturbations and the order-eight limit."""
from collections import defaultdict
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location("coloring_algebra",HERE.parent/"three-color-coupling"/"coupling.py")
c = module_from_spec(spec)
spec.loader.exec_module(c)
Q, e = c.Q, c.e
S = tuple(range(6))
P = tuple(permutations(range(3)))
PAIRS = tuple(combinations_with_replacement(range(6),2))


def check(ok, reason):
    if not ok:
        raise ValueError(reason)


def symplectic_form(m):
    return [[Q(1 if j==i+1 and i%2==0 else -1 if i==j+1 and j%2==0 else 0)
             for j in range(2*m)] for i in range(2*m)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def pfaffian(a):
    n = len(a)
    if n==0:
        return Q(1)
    total = Q(0)
    for j in range(1,n):
        keep = [i for i in range(1,n) if i!=j]
        total += (-1)**(j+1) * a[0][j] * pfaffian([[a[i][k] for k in keep] for i in keep])
    return total


def block(z,w):
    return [[z,w],[-w.conjugate(),z.conjugate()]]


def gate(r,s,t,v):
    cos,sin = (1-t*t)/(1+t*t),2*t/(1+t*t)
    result = e.identity(6)
    for i in range(2):
        result[2*r+i][2*r+i] = result[2*s+i][2*s+i] = Q(cos)
        for j in range(2):
            result[2*r+i][2*s+j] = sin*v[i][j]
            result[2*s+j][2*r+i] = -sin*v[i][j].conjugate()
    return result


def fixture():
    phases = [block(Q(0,1),Q(0)),block(Q(0),Q(1)),block(Q(F(3,5)),Q(0,F(4,5)))]
    check(e.multiply(*phases[:2])!=e.multiply(*reversed(phases[:2])),"noncommuting phases")
    u = e.identity(6)
    for r,s,t,k in [(0,1,F(1,3),0),(1,2,F(2,5),1),(0,2,F(1,4),2),(0,1,F(2,7),1)]:
        u = e.multiply(gate(r,s,t,phases[k]),u)
    return u


def expanded(indices):
    return tuple(k for i in indices for k in (2*i,2*i+1))


def birkhoff(u):
    mat = [[c.minor(u,expanded([r]),expanded([s])) for s in range(3)] for r in range(3)]
    check(all(z.imag==0 for row in mat for z in row),"real block determinants")
    mat = [[z.real for z in row] for row in mat]
    check(all(sum(r)==1 for r in mat+transpose(mat)),"doubly stochastic")
    h = min(min(row) for row in mat)
    check(h>0,"strict fixture")
    residual = [[x-h for x in row] for row in mat]
    law = []
    for p in P:
        amount = min(residual[i][j] for i,j in enumerate(p))
        for i,j in enumerate(p):
            residual[i][j] -= amount
        law.append(h/2+amount)
    check(sum(law)==1 and min(law)>=h/2,"Birkhoff reserve")
    for k in range(4):
        for rows in combinations(range(3),k):
            for cols in combinations(range(3),k):
                value = sum(w for p,w in zip(P,law) if tuple(sorted(p[i] for i in rows))==cols)
                check(Q(value)==c.minor(u,expanded(rows),expanded(cols)),"paired minor")
    return h,law


def table(coloring):
    return tuple(tuple(sum(coloring[i]==j for i in (2*r,2*r+1)) for j in range(3))
                 for r in range(3))


def node(p,q):
    return tuple(y for i in range(3) for y in (p[i],q[i]))


def signed(u):
    ix = {table(node(P[i],P[j])):k for k,(i,j) in enumerate(PAIRS)}
    law = [F(0)]*21
    for coloring,w in zip(c.COLORINGS,c.canon(u)):
        law[ix[table(coloring)]] += w
    return law


def coefficients(u,law):
    classical = defaultdict(F)
    for (i,j),w in zip(PAIRS,law):
        assignment = node(P[i],P[j])
        for rows in c.SUBSETS:
            a = tuple(2-sum(k//2==g for k in rows) for g in range(3))
            b = tuple(sum(assignment[k]==g for k in rows) for g in range(3))
            classical[a,b] += w
    physical = defaultdict(F)
    for rows in c.SUBSETS:
        a = tuple(2-sum(k//2==g for k in rows) for g in range(3))
        for cols in combinations(S,len(rows)):
            b = tuple(sum(k//2==g for k in cols) for g in range(3))
            physical[a,b] += c.abs2(c.minor(u,rows,cols))
    check(classical==physical and len(physical)==141,"924-to-141 coefficient expansion")


def four_pair_limit():
    data = json.loads((HERE.parent/"single-cycle"/"four-cycle-obstruction.json").read_text())
    u = [[Q(F(x)) if i%2==j%2 else Q(0) for j,x in enumerate(
        value for value in data["U"][i//2] for _ in range(2))] for i in range(8)]
    J = symplectic_form(4)
    check(e.multiply(e.multiply(u,J),transpose(u))==J,"order-eight symplectic identity")
    cover = [(v["rows"],v["columns"]) for v in data["cover"]]
    score = sum(c.minor(u,expanded(r),expanded(t)) for r,t in cover)
    check(score==Q(F(14953,15625)),"order-eight cover gap")
    minimum = min(sum(set(p[i] for i in r)==set(t) for r,t in cover) for p in permutations(range(4)))
    check(minimum==1,"permutation cover")
    return score.real


def main():
    u = fixture()
    c.check_unitary(u)
    J = symplectic_form(3)
    check(e.multiply(e.multiply(u,J),transpose(u))==J,"symplectic fixture")
    h,lam = birkhoff(u)
    w = [lam[i]*lam[j]*(1 if i==j else 2) for i,j in PAIRS]
    coefficients(u,w)
    a = [Q(2,1),Q(-1,3),Q(4,-2)]
    b = [Q(3,-1),Q(1,2),Q(-2,-1)]
    db = e.diagonal([z for z in b for _ in range(2)])
    m = e.multiply(e.multiply(u,db),e.adjoint(u))
    for i in S:
        m[i][i] += a[i//2]
    skew = e.multiply(m,J)
    check(all(skew[i][j]==-skew[j][i] for i in S for j in S),"complex skew matrix")
    p = [c.O for _ in P]
    for k,perm in enumerate(P):
        for i,j in enumerate(perm):
            p[k] *= a[i]+b[j]
    root = sum(weight*val for weight,val in zip(lam,p))
    check(root==pfaffian(skew) and root*root==e.determinant(m),"Pfaffian and determinant")
    radius = h*h/1920
    v = e.multiply(c.rotation(6,0,2,radius/4),u)
    c.check_unitary(v)
    check(max(c.abs2(x-y) for row,base in zip(v,u) for x,y in zip(row,base))<=radius**2,
          "within entrywise ball")
    # A block with unequal singular values cannot be a left/right unitary
    # change of basis of a scalar times a quaternionic unit.
    split = False
    for r in range(3):
        for s in range(3):
            block_v = [[v[2*r+i][2*s+j] for j in range(2)] for i in range(2)]
            gram = e.multiply(block_v,e.adjoint(block_v))
            split |= gram[0][0]!=gram[1][1] or gram[0][1]!=Q(0)
    check(split,"perturbation must escape block-gauged symplectic class")
    repaired = [mass+new-old for mass,new,old in zip(w,signed(v),signed(u))]
    check(min(repaired)>=h*h/8 and sum(repaired)==1,"perturbed positivity")
    coefficients(v,repaired)
    print(f"Six-site Pfaffian law passes. Block minimum {h}; certified radius {radius}.")
    print(f"All 141 coefficients match on both fixtures. Four-pair cover {four_pair_limit()} < 1.")


if __name__ == "__main__":
    main()
