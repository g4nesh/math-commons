"""Exact finite checks for the five alternating triple correlations."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location

SITES = tuple(range(6))
EDGES = tuple(combinations(SITES,2))
TRIPLES = tuple(combinations(SITES,3))
SUBSETS = tuple(s for k in range(7) for s in combinations(SITES,k))
FIVE = ((0,1,2), (0,1,3), (0,1,4), (0,2,3), (0,2,4))
COLORS = (0,1,2)


def check(ok, description):
    if not ok:
        raise ValueError(description)


def rank(rows):
    """Fraction-free elimination with primitive integer rows."""
    pivots = {}
    for row in rows:
        row = list(row)
        for j,pivot in pivots.items():
            if row[j]:
                a,b = pivot[j],row[j]
                row = [a*x-b*y for x,y in zip(row,pivot)]
                g = gcd(*row)
                if g > 1:
                    row = [x//g for x in row]
        j = next((i for i,x in enumerate(row) if x), None)
        if j is not None:
            pivots[j] = row
    return len(pivots)


def hist(c,s):
    return tuple(sum(c[i] == j for i in s) for j in COLORS)


def alt(h):
    x,y,z = h
    return (x-y)*(y-z)*(x-z)//2


def parity(p):
    return (-1)**sum(x > y for i,x in enumerate(p) for y in p[i+1:])


def coloring(m,p):
    c = [None]*6
    for (i,j),label in zip(m,p):
        c[i] = c[j] = label
    return tuple(c)


def full_histograms(law):
    values = defaultdict(F)
    for c,mass in law.items():
        for s in SUBSETS:
            values[s,hist(c,s)] += mass
    return values


def main():
    source = Path(__file__).resolve().parents[1]/"paired-matching-boundary"/"check.py"
    spec = spec_from_file_location("paired_matching", source)
    paired = module_from_spec(spec)
    spec.loader.exec_module(paired)

    matchings = tuple(paired.matchings(SITES))
    colorings = [coloring(m,p) for m in matchings for p in permutations(COLORS)]
    check(len(set(colorings)) == 90, "enumeration of labeled matchings")
    incidence = [[int(set(e) <= set(s)) for s in TRIPLES] for e in EDGES]
    coordinates = [[int(s == t) for s in TRIPLES] for t in FIVE]
    check(rank(incidence) == 15 and rank(incidence + coordinates) == 20, "five coordinates")
    colored_edges = [[int(c[i] == c[j] == label) for c in colorings]
                     for label in COLORS for i,j in EDGES]
    full = [[int(hist(c,s) == h) for c in colorings]
            for s in SUBSETS for h in sorted({hist(c,s) for c in colorings})]
    triples = [[alt(hist(c,s)) for c in colorings] for s in TRIPLES]
    reduced = colored_edges + [triples[TRIPLES.index(s)] for s in FIVE]
    check(rank(colored_edges) == 38 and rank(reduced) == 43, "reduced rank")
    check(len(full) == 358 and rank(full) == rank(full+reduced) == 43, "full equivalence")
    frame = [[alt(hist(coloring(m,COLORS),s)) for m in matchings] for s in TRIPLES]
    for s,row in zip(TRIPLES,frame):
        for t,other in zip(TRIPLES,frame):
            gram = sum(x*y for x,y in zip(row,other))
            check(gram == (-9,3,-3,9)[len(set(s)&set(t))], "frame Gram coefficient")
    for i,s in enumerate(TRIPLES):
        for j,m in enumerate(matchings):
            for p in permutations(COLORS):
                h = hist(coloring(m,p),s)
                check(alt(h) == parity(p)*frame[i][j], "orientation sign")
                check(alt(h) == (0,1,-1)[(h[1]+2*h[2]) % 3], "cube-root sign")
    check(rank(frame) == 5, "frame rank")
    # Vary orientation at fixed positive matching marginals, then reconstruct
    # the entire target from the biases and a fresh matching decomposition.
    base = [F(k,120) for k in range(1,16)]
    biases = [w * F((k%5)-2,3) for k,w in enumerate(base)]
    target_law = {coloring(m,p): (w+parity(p)*b)/6
                  for m,w,b in zip(matchings,base,biases) for p in permutations(COLORS)}
    target = full_histograms(target_law)
    delta = [sum(target[s,h]*alt(h) for h in {hist(c,s) for c in colorings}) for s in TRIPLES]
    check(any(delta) and all(x == sum(b*v for b,v in zip(row,biases))
                              for x,row in zip(delta,frame)), "nonzero target biases")
    # Algebraic frame correction, allowing signed intermediate weights.
    corrections = [sum(v*w for v,w in zip(col,delta))/216 for col in zip(*frame)]
    check(all(6*sum(v*w for v,w in zip(row,corrections)) == d
              for row,d in zip(frame,delta)), "frame right inverse")
    x = paired.edge_vector(target)
    t = 1-sum(map(abs,biases))
    residual = {e:(x[e]-sum(abs(v) for m,v in zip(matchings,biases) if e in m))/t for e in EDGES}
    recovered = paired.decompose_matching(residual)
    law = {coloring(m,p): (t*recovered.get(m,0)+abs(b)+parity(p)*b)/6
           for m,b in zip(matchings,biases) for p in permutations(COLORS)}
    paired.check_law(law,target)
    # Check the principal-minor expression against the 924 original quantum minors.
    matrix = paired.fourier_entries()
    for i in (4,5):
        for k in (1,3,5):
            matrix[i][k] = paired.times(matrix[i][k],(0,1))
    matrix = paired.mixed_rows(matrix, 0,2,F(1,7))
    matrix = paired.mixed_rows(matrix, 1,5,F(1,8))
    roots = [(1,0),(0,1),(-1,-1)]
    def conjugate(z):
        a,b = z
        return a-b,-b
    def entry(i,j):
        summands = [paired.times(paired.times(matrix[i][k],roots[k//2]),conjugate(matrix[j][k]))
                    for k in SITES]
        return tuple(sum(coeff)/3 for coeff in zip(*summands))
    unitary = [[entry(i,j) for j in SITES] for i in SITES]
    squared,_ = paired.eisenstein_data(matrix)
    physical = paired.grouped(squared,(0,0,1,1,2,2))
    moments = []
    for s in TRIPLES:
        coefficient = F(0)
        for p in permutations(s):
            term = (F(1),F(0))
            for i,j in zip(s,p):
                term = paired.times(term,unitary[i][j])
            coefficient += parity(p)*term[1]
        target_a = sum(w*alt(h) for (rows,h),w in physical.items() if rows == s)
        check(coefficient == target_a, "principal minor of the order-three unitary")
        moments.append(target_a)
    check(any(moments), "complex fixture exercises a nonzero triple")
    print("45 colored edges have rank 38; five triple signs recover the full rank 43.")
    print("Frame Gram and orientation identities hold on 15 matchings and all 90 colorings.")
    print("Cyclic reconstruction matches all 358 coefficients. All 20 physical cube-root minors agree.")


if __name__ == "__main__":
    main()
