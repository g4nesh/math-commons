"""Classify and solve spectrum-independent weights for two (2,2,2) spectra."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd
from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location

INDICES = tuple(range(6))
PAIRS = tuple(combinations_with_replacement(INDICES,2))
PARITY = (1,1,1,-1,-1,-1)


def parity(p):
    return (-1)**sum(x>y for i,x in enumerate(p) for y in p[i+1:])


PERMS = tuple(sorted(permutations(range(3)),key=lambda p: -parity(p)))


def check(ok, reason):
    if not ok:
        raise ValueError(reason)


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def rank(rows):
    pivots = {}
    for row in rows:
        row = list(row)
        for j,pivot in pivots.items():
            if row[j]:
                a,b = pivot[j],row[j]
                row = [a*x-b*y for x,y in zip(row,pivot)]
                g = gcd(*row)
                if g>1:
                    row = [x//g for x in row]
        j = next((i for i,x in enumerate(row) if x),None)
        if j is not None:
            pivots[j] = row
    return len(pivots)


def monomial_rows():
    rows = {}
    for col,(i,j) in enumerate(PAIRS):
        assignment = PERMS[i]+PERMS[j]
        for take_b in product((False,True),repeat=6):
            a = tuple(sum(not take_b[v] and v%3==r for v in INDICES) for r in range(3))
            b = tuple(sum(take_b[v] and assignment[v]==r for v in INDICES) for r in range(3))
            rows.setdefault(a+b,[0]*21)[col] += 1
    return rows


def kernel():
    return [[int(i==j==v) if i==j else PARITY[i]*PARITY[j]*int(v in (i,j))
             for i,j in PAIRS] for v in INDICES]


def construct_law(q):
    """Exact difference-constraint solver; None means no common law for these data."""
    q = tuple(map(F,q))
    check(len(q)==21 and sum(q)==1, "expected 21 normalized rational weights")
    arcs = []
    for k,(i,j) in enumerate(PAIRS):
        if i==j:
            arcs.append((2*i,2*i+1,2*q[k]))
        else:
            # Even nodes encode +u, odd nodes -u; arc source->target enforces
            # potential(target) - potential(source) <= length.
            sign = int(PARITY[i]==PARITY[j])
            arcs.extend([(2*j+1-sign,2*i+sign,q[k]), (2*i+1-sign,2*j+sign,q[k])])
    distance = [F(0)]*12
    for _ in range(12):
        changed = False
        for source,target,length in arcs:
            if distance[source]+length < distance[target]:
                distance[target] = distance[source]+length
                changed = True
        if not changed:
            break
    else:
        return None
    u = [(distance[2*i]-distance[2*i+1])/2 for i in INDICES]
    weights = [q[k]+(u[i] if i==j else PARITY[i]*PARITY[j]*(u[i]+u[j]))
               for k,(i,j) in enumerate(PAIRS)]
    check(min(weights)>=0 and sum(weights)==1, "decoded law")
    return weights


def simple_cycles():
    rays = set()
    def visit(path):
        i,color = path[-1]
        for j in INDICES:
            if (PARITY[i]==PARITY[j]) != (color==0):
                continue
            step = (j,1-color)
            if step == path[0]:
                sequence = path+[step]
                count = Counter(tuple(sorted((a[0],b[0]))) for a,b in zip(sequence,sequence[1:]))
                values = [count[i,j]*(2 if i==j else 1) for i,j in PAIRS]
                divisor = gcd(*values)
                rays.add(tuple(x//divisor for x in values))
            elif step > path[0] and step not in path:
                visit(path+[step])
    for i in INDICES:
        for color in (0,1):
            visit([(i,color)])
    return rays


def symmetry_group():
    def inverse(p):
        return tuple(p.index(j) for j in range(3))
    index_perms = set()
    for rho,gamma,transpose in product(PERMS,PERMS,(False,True)):
        action = []
        for p in PERMS:
            base = inverse(p) if transpose else p
            q = tuple(rho[base[inverse(gamma)[i]]] for i in range(3))
            action.append(PERMS.index(q))
        index_perms.add(tuple(action))
    target = {p for p in permutations(INDICES) if all(
        PARITY[i]*PARITY[j]==PARITY[p[i]]*PARITY[p[j]] for i,j in PAIRS)}
    check(index_perms==target and len(target)==72,"physical symmetry group")
    return [tuple(PAIRS.index(tuple(sorted((p[i],p[j])))) for i,j in PAIRS) for p in index_perms]


def primitive(row):
    divisor = gcd(*row)
    return tuple(v//divisor for v in row)


def classify(rows):
    walks = simple_cycles()
    sub = {v:frozenset(i for i,x in enumerate(v) if x) for v in walks}
    supports = set(sub.values())
    minimal = {s for s in supports if not any(t<s for t in supports)}
    rays = {v for v,s in sub.items() if s in minimal}
    check(len(walks)==573 and len(rays)==189 and len(minimal)==189,"extreme supports")
    check(all(dot(v,d)==0 for v in walks for d in kernel()),"projected flow balance")
    group = symmetry_group()
    orbits = Counter(min(tuple(v[i] for i in p) for p in group) for v in rays)
    check(sorted(orbits.values())==[9,18,18,36,36,36,36],"seven orbits")
    published = [
        [(2,2,1),(2,5,1),(5,5,1)],
        [(2,2,2),(2,4,1),(2,5,1),(4,5,1)],
        [(1,2,1),(1,5,1),(2,4,1),(4,5,1)],
        [(1,2,1),(1,5,1),(2,4,1),(4,4,1),(5,5,1)],
        [(1,2,2),(1,5,2),(2,3,1),(2,4,1),(3,4,1),(5,5,2)],
        [(0,2,1),(0,5,1),(1,1,1),(1,4,1),(2,3,1),(3,4,1),(5,5,1)],
        [(0,2,1),(0,5,1),(1,2,1),(1,4,1),(2,3,2),(3,3,2),(4,5,1)]]
    reps = []
    for record in published:
        counts = {(i,j):v for i,j,v in record}
        reps.append(tuple(counts.get(e,0) for e in PAIRS))
    check(set(reps)==set(orbits),"written representatives")
    direct = {primitive(r) for r in rows.values()}
    check(len(direct)==71 and len(direct & rays)==45,"monomial inequalities")
    for i,v in enumerate(reps):
        orbit = {tuple(v[k] for k in p) for p in group}
        check(len(orbit)==orbits[v] and (orbit<=direct if i<3 else orbit.isdisjoint(direct)),
              "whole-orbit classification")
        if i>=3:
            q = [F(not x) for x in v]
            q[next(k for k,x in enumerate(v) if x)] = -F(1,4)
            total = sum(q)
            q = [x/total for x in q]
            check(dot(q,v)<0 and all(dot(q,z)>0 for z in rays if z!=v),"irredundancy witness")
            check(all(dot(q,r)>=0 for r in rows.values()),"abstract monomial positivity")
            check(construct_law(q) is None,"negative cycle detection")
    return rays


def physical_checks(rows,rays):
    source = Path(__file__).resolve().parents[1]/"paired-matching-boundary"/"check.py"
    spec = spec_from_file_location("paired_fixture",source)
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    labels = (0,0,1,1,2,2)
    def table(c):
        return tuple(m.histogram(c,(2*i,2*i+1)) for i in range(3))
    index = {table(tuple(y for i in range(3) for y in (PERMS[a][i],PERMS[b][i]))):k
             for k,(a,b) in enumerate(PAIRS)}
    check(len(index)==21,"contingency tables")
    matrix = m.fourier_entries()
    for i in (4,5):
        for k in (1,3,5):
            matrix[i][k] = m.times(matrix[i][k],(0,1))
    for i,j,t in [(0,2,F(2,5)),(1,5,F(1,2)),(2,4,F(6,11))]:
        matrix = m.mixed_rows(matrix,i,j,t)
        squared,signed = m.eisenstein_data(matrix)
        q = [F(0)]*21
        for coloring,value in signed.items():
            q[index[table(coloring)]] += value
        check(min(q)<0 and min(dot(q,y) for y in rays)>0,"strict complex fixture")
        law = construct_law(q)
        check(law is not None,"expected a nonnegative law")
        physical = Counter()
        for s in m.SUBSETS:
            exp_a = tuple(2-k for k in m.histogram(labels,s))
            for tcols in combinations(INDICES,len(s)):
                key = exp_a+m.histogram(labels,tcols)
                physical[key] += squared(s,tcols)
        check(set(physical)==set(rows) and all(dot(q,r)==dot(law,r)==physical[key]
              for key,r in rows.items()),"direct grouped Cauchy–Binet")
    print("Three complex fixtures: 924 squared minors each; negative starting masses repaired exactly.")


def main():
    rows = monomial_rows()
    check(len(rows)==141 and rank(rows.values())==15 and rank(kernel())==6,"complete kernel")
    check(all(dot(d,r)==0 for d in kernel() for r in rows.values()),"kernel identity")
    rays = classify(rows)
    physical_checks(rows,rays)
    print("21 products, 6 independent repairs; 573 cycles give 189 rays in 7 orbits.")
    print("45 ray inequalities are grouped squared minors; the other 4 universal orbit signs have separate Grassmann certificates (see double-paired-five-term/GRASSMANN.md).")


if __name__ == "__main__":
    main()
