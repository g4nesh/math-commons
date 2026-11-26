"""Replay six-site paired-spectrum laws with exact rational arithmetic."""
from collections import defaultdict, deque
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations

SITES = tuple(range(6))
SUBSETS = tuple(s for k in range(7) for s in combinations(SITES, k))
EDGES = tuple(combinations(SITES, 2))
TRIPLES = tuple(combinations(SITES, 3))
COLORS = tuple(range(3))
LABELS = tuple(permutations(COLORS))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def histogram(labels, rows):
    return tuple(sum(labels[i] == j for i in rows) for j in COLORS)


def matchings(sites):
    if not sites:
        yield ()
        return
    i, *others = sites
    for j in others:
        for rest in matchings(tuple(k for k in others if k != j)):
            yield ((i, j),) + rest


def decompose_matching(x):
    """Peel matching vertices while retaining every edge and odd-cut inequality."""
    remainder = dict(x)
    total = F(1)
    law = {}
    while total:
        slack = {s: total - sum(remainder[e] for e in combinations(s, 2)) for s in TRIPLES}
        best, weight = None, F(0)
        for matching in matchings(SITES):
            limits = [remainder[e] for e in matching]
            for s in TRIPLES:
                if all(i not in s or j not in s for i, j in matching):
                    limits.append(slack[s])
            amount = min(limits)
            if amount > weight:
                best, weight = matching, amount
        require(weight > 0, "no feasible perfect matching in residual face")
        law[best] = law.get(best, 0) + weight
        for e in best:
            remainder[e] -= weight
        total -= weight
    require(all(v == 0 for v in remainder.values()), "matching mass remains")
    return law


def uniform_labels(matching_law):
    law = defaultdict(F)
    for matching, weight in matching_law.items():
        for perm in LABELS:
            coloring = [None] * 6
            for (i, j), color in zip(matching, perm):
                coloring[i] = coloring[j] = color
            law[tuple(coloring)] += weight / 6
    return law


def matrix_product(left, right):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*right)] for row in left]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def rotations(sites):
    matrix = identity(6)
    for k, (i, j) in enumerate(combinations(sites, 2)):
        t = F(k + 1, 2*k + 5)
        c, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
        matrix[i], matrix[j] = ([c*x+s*y for x,y in zip(matrix[i], matrix[j])],
                                [-s*x+c*y for x,y in zip(matrix[i], matrix[j])])
    return matrix


def real_minor_squared(matrix):
    @cache
    def det(rows, cols):
        if not rows:
            return F(1)
        return sum((-1)**j * matrix[rows[0]][c] * det(rows[1:], cols[:j]+cols[j+1:])
                   for j, c in enumerate(cols))
    return lambda r,c: det(r,c)**2


def birkhoff(matrix, sites):
    remaining = [[matrix[i][j]**2 for j in sites] for i in sites]
    law = {}
    for indices in permutations(range(len(sites))):
        weight = min(remaining[i][j] for i, j in enumerate(indices))
        if weight:
            p = list(SITES)
            for i, j in enumerate(indices):
                p[sites[i]] = sites[j]
                remaining[i][j] -= weight
            law[tuple(p)] = weight
    require(sum(law.values()) == 1, "three-site Birkhoff law")
    return law


def orient_edges(edges, demands):
    residual = defaultdict(dict)

    def add(u, v, cap):
        residual[u][v] = cap
        residual[v][u] = 0

    for edge, weight in edges.items():
        add("source", edge, weight)
        for i in edge:
            add(edge, i, F(1))
    for i, weight in demands.items():
        add(i, "sink", weight)
    while True:
        parent = {"source": None}
        queue = deque(parent)
        while queue and "sink" not in parent:
            u = queue.popleft()
            for v, mass in residual[u].items():
                if mass > 0 and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if "sink" not in parent:
            break
        path = []
        v = "sink"
        while v != "source":
            u = parent[v]
            path.append((u, v))
            v = u
        amount = min(residual[u][v] for u, v in path)
        for u, v in path:
            residual[u][v] -= amount
            residual[v][u] += amount
    require(all(residual[i]["sink"] == 0 for i in demands), "endpoint flow failed")
    return {(e,i): residual[i][e] for e in edges for i in e}


def four_site_law(matrix, sites, labels):
    active = {c: tuple(i for i in sites if labels[i] == c) for c in set(labels[i] for i in sites)}
    squared = real_minor_squared(matrix)
    law = defaultdict(F)
    if len(active) == 1:
        return {labels: F(1)}
    if len(active) == 2:
        c, other = active
        for selected in combinations(sites, len(active[c])):
            result = list(labels)
            for i in sites:
                result[i] = c if i in selected else other
            law[tuple(result)] += squared(selected, active[c])
        return law
    u, v = sorted(c for c in active if len(active[c]) == 1)
    repeated = next(c for c in active if c not in (u,v))
    cols = tuple(sorted(active[u] + active[v]))
    weights = {e: squared(e, cols) for e in combinations(sites, 2)}
    flow = orient_edges(weights, {i: matrix[i][active[u][0]]**2 for i in sites})
    for ((i,j), first), mass in flow.items():
        result = list(labels)
        for site in sites:
            result[site] = repeated
        second = j if first == i else i
        result[first], result[second] = u, v
        law[tuple(result)] += mass
    return law


def boundary_fixture(three_first):
    small, large = (0,1,2), (2,3,4,5)
    a, b = rotations(small), rotations(large)
    matrix = matrix_product(a, b) if three_first else matrix_product(b, a)
    # With the large factor on the left, two rows miss one column of each color.
    labels = (0,0,1,1,2,2) if three_first else (0,1,2,0,1,2)
    law = defaultdict(F)
    for p, w in birkhoff(a, small).items():
        conditional_labels = labels if three_first else tuple(labels[p[i]] for i in SITES)
        for coloring, mass in four_site_law(b, large, conditional_labels).items():
            result = tuple(coloring[p[i]] for i in SITES) if three_first else coloring
            law[result] += w * mass
    if not three_first:
        for c in COLORS:
            pair = tuple(i for i in SITES if labels[i] == c)
            matrix = matrix_product(matrix, rotations(pair))
        require(all(x != 0 for row in matrix for x in row), "hidden zero rectangle fixture")
    require(matrix_product(matrix, list(zip(*matrix))) == identity(6), "fixture orthogonality")
    return real_minor_squared(matrix), labels, law


def times(z, w):
    """Product in Q[zeta_3], in the basis 1,zeta_3."""
    a,b = z
    c,d = w
    return (a*c-b*d, a*d+b*c-b*d)


def fourier_entries():
    rotations2 = [(F(3,5), F(4,5)), (F(5,13), F(12,13)), (F(7,25), F(24,25))]
    factors = [((c,s),(-s,c)) for c,s in rotations2]
    roots = [(1,0), (0,1), (-1,-1)]
    return [[times((factors[i//2][i%2][j%2],0), roots[(i//2 * (j//2)) % 3])
             for j in SITES] for i in SITES]


def eisenstein_data(matrix):
    """All entries have a common omitted 1/sqrt(3) normalization."""
    @cache
    def det(rows, cols):
        if not rows:
            return (F(1),F(0))
        a = b = F(0)
        for j,c in enumerate(cols):
            x,y = times(matrix[rows[0]][c], det(rows[1:], cols[:j]+cols[j+1:]))
            a += (-1)**j * x
            b += (-1)**j * y
        return a,b

    def squared(rows, cols):
        a,b = det(rows,cols)
        return (a*a-a*b+b*b) / 3**len(rows)

    a,b = det(SITES,SITES)
    conjugate_det = (a-b,-b)
    signed = {}
    for matching in matchings(SITES):
        for order in LABELS:
            pairs = [matching[i] for i in order]
            seq = sum(pairs, ())
            sign = (-1)**sum(seq[i] > seq[j] for i,j in EDGES)
            z = conjugate_det
            coloring = [None]*6
            for j,pair in enumerate(pairs):
                z = times(z, det(pair, (2*j,2*j+1)))
                for i in pair:
                    coloring[i] = j
            signed[tuple(coloring)] = sign*(z[0]-z[1]/2)/3**6
    return squared, signed


def mixed_rows(matrix, i, j, t):
    c,s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    rows = matrix.copy()
    rows[i] = [tuple(c*a+s*b for a,b in zip(x,y)) for x,y in zip(matrix[i],matrix[j])]
    rows[j] = [tuple(-s*a+c*b for a,b in zip(x,y)) for x,y in zip(matrix[i],matrix[j])]
    return rows


def grouped(squared, labels):
    out = defaultdict(F)
    count = 0
    for s in SUBSETS:
        for t in combinations(SITES, len(s)):
            out[s, histogram(labels,t)] += squared(s,t)
            count += 1
    require(count == 924 and len(out) == 358, "coverage of grouped coefficients")
    require(all(sum(w for (rows,_),w in out.items() if rows == s) == 1 for s in SUBSETS),
            "minor distributions are not normalized")
    return out


def check_law(law, target):
    require(min(law.values()) >= 0 and sum(law.values()) == 1, "invalid coloring probabilities")
    require(all(histogram(c,SITES) == (2,2,2) for c in law), "invalid coloring")
    reconstructed = defaultdict(F)
    for c, weight in law.items():
        for s in SUBSETS:
            reconstructed[s,histogram(c,s)] += weight
    require(all(reconstructed[key] == v for key,v in target.items()), "grouped coloring mismatch")


def edge_vector(target):
    return {e: sum(target[e,tuple(2 if c == j else 0 for c in COLORS)] for j in COLORS) for e in EDGES}


def check_matching(target):
    x = edge_vector(target)
    matching = decompose_matching(x)
    averaged = {(s,h): sum(target[s,tuple(h[j] for j in p)] for p in LABELS)/6
                for s,h in target}
    check_law(uniform_labels(matching), averaged)
    triple_slack = [1-sum(x[e] for e in combinations(s,2)) for s in TRIPLES]
    return min(x.values()), min(triple_slack), averaged == target


def repair(target, signed):
    orbits = defaultdict(list)
    for c in signed:
        pairs = tuple(sorted(tuple(i for i in SITES if c[i] == j) for j in COLORS))
        orbits[pairs].append(c)
    floors = {m: sum(signed[c] for c in orbit)-6*min(signed[c] for c in orbit)
              for m,orbit in orbits.items()}
    residual = edge_vector(target)
    for m,floor in floors.items():
        for e in m:
            residual[e] -= floor
    total = 1-sum(floors.values())
    require(total >= 0 and min(residual.values()) >= 0, "orbit-floor edge constraints")
    require(all(sum(residual[e] for e in combinations(s,2)) <= total for s in TRIPLES),
            "orbit-floor cut constraints")
    law = {m:F(0) for m in orbits}
    if total:
        law.update({m:w*total for m,w in decompose_matching(
            {e:x/total for e,x in residual.items()}).items()})
    full = {}
    for m,orbit in orbits.items():
        average = sum(signed[c] for c in orbit)/6
        for c in orbit:
            full[c] = (law[m]+floors[m])/6 + signed[c]-average
    check_law(full, target)
    return total


def main():
    for order, name in [(False, "zero edge"), (True, "tight triple")]:
        squared, labels, law = boundary_fixture(order)
        target = grouped(squared, labels)
        check_law(law, target)
        edge, cut, symmetric = check_matching(target)
        require((cut == 0 if order else edge == 0) and not symmetric, "wrong boundary fixture")
        print(f"{name}: full and averaged laws pass; minimum edge {edge}, minimum cut slack {cut}")
    labels = (0,0,1,1,2,2)
    matrix = fourier_entries()
    squared, signed = eisenstein_data(matrix)
    target = grouped(squared, labels)
    edge, cut, symmetric = check_matching(target)
    require(edge > 0 and cut > 0 and symmetric, "expected symmetric interior")
    require(min(F(1,2), 5*edge, 5*cut/2)/5400 == F(2,528125), "certified ball radius")
    print(f"Fourier family: full law passes; minimum edge {edge}, minimum cut slack {cut}")
    require(repair(target, signed) == 1, "symmetric residual mass")
    for i,j,t in [(0,2,F(1,200)), (1,5,F(1,300))]:
        matrix = mixed_rows(matrix, i,j,t)
    squared,signed = eisenstein_data(matrix)
    target = grouped(squared,labels)
    edge,cut,symmetric = check_matching(target)
    require(edge > 0 and cut > 0 and not symmetric, "expected nonsymmetric interior")
    require(sum(x < 0 for x in signed.values()) == 12, "signed-law obstruction fixture")
    residual_mass = repair(target,signed)
    require(0 < residual_mass < 1, "expected nontrivial orbit correction")
    print("Nonsymmetric complex perturbation: full law passes the orbit-floor test.")
    print("Each fixture: 924 squared minors, 358 grouped coefficients, all 64 row subsets.")


if __name__ == "__main__":
    main()
