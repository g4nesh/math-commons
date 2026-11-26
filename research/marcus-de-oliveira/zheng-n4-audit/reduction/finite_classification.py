#!/usr/bin/env python3
"""Independent exact re-derivation of the finite classification in Section 4 of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (Proposition 4.1, Lemma 4.2,
the 214 orbits, Appendices D and E, Table 5).

Everything here is the auditor's own code, integer/rational arithmetic only
(python-flint for exact integer linear algebra). Nothing from the supplement is
executed or imported. The supplement is used only for task 2 (and the optional
cross-check of task 7), as JSON data read with json.load.

usage: finite_classification.py [SUPP_DIR] [--json OUT.json] [--lane-npy facets4.npy]
       SUPP_DIR defaults to the environment variable MOC4_SUPP. Without it, task 2 is
       reported as SKIPPED (and the run does not count as complete).
Exit status: 0 if every claim line is PASS, 1 otherwise.
"""
import itertools, json, os, sys, time
from math import gcd, comb
from fractions import Fraction
import flint

T0 = time.time()
FAILS = []
RESULTS = {}


def claim(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    if not ok:
        FAILS.append(name)
    print(f"[{tag}] {name}" + (f" :: {detail}" if detail else ""), flush=True)
    return ok


def note(msg):
    print(f"[info] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Permutations, indicator embedding
# ---------------------------------------------------------------------------
PERMS = list(itertools.permutations(range(4)))          # lexicographic, 0-based images
IDX = {p: i for i, p in enumerate(PERMS)}
N = 24
assert PERMS[0] == (0, 1, 2, 3)


def oneline(p):
    return "".join(str(x + 1) for x in p)


SING = [frozenset([i]) for i in range(4)]
PAIR = [frozenset(s) for s in itertools.combinations(range(4), 2)]
COORDS = [(I, J) for I in SING for J in SING] + [(I, J) for I in PAIR for J in PAIR]  # 16 + 36


def image(p, I):
    return frozenset(p[i] for i in I)


E = [[1 if image(p, I) == J else 0 for (I, J) in COORDS] for p in PERMS]   # 24 x 52


def rank(rows):
    return flint.fmpz_mat(rows).rank() if rows else 0


def primitive(v):
    g = 0
    for x in v:
        g = gcd(g, int(x))
    return tuple(int(x) // g for x in v) if g else tuple(int(x) for x in v)


def bits(v):
    m = 0
    for i, x in enumerate(v):
        if x == 0:
            m |= 1 << i
    return m


def popc(m):
    return bin(m).count("1")


def setof(m):
    return [i for i in range(N) if m >> i & 1]


# ---------------------------------------------------------------------------
# Task 0: dimension of the function space and the affine hull
# ---------------------------------------------------------------------------
E1 = [[1] + r for r in E]
claim("span{chi_IJ} + constants on S_4 has dimension 14", rank(E1) == 14, f"rank={rank(E1)}")
claim("constants already lie in span{chi_IJ} (rank of the 52 functions alone is 14)", rank(E) == 14)
claim("affine hull of the 24 indicator points has dimension 13",
      rank([[a - b for a, b in zip(E[s], E[0])] for s in range(1, N)]) == 13)

# a basis V (24 x 14) of the function space W: greedy choice of columns of [1 | E]
cols = []
for c in range(53):
    trial = cols + [c]
    if rank([[E1[s][k] for k in trial] for s in range(N)]) == len(trial):
        cols = trial
    if len(cols) == 14:
        break
V = [[E1[s][k] for k in cols] for s in range(N)]
VM = flint.fmpz_mat(V)
claim("chosen basis V of W has rank 14", VM.rank() == 14, f"columns {cols}")


def in_W(y):
    return rank([V[s] + [y[s]] for s in range(N)]) == 14


# ---------------------------------------------------------------------------
# The matrix P of the paper re-derived from the matching basis (A.1), and Q ~ conv{e_sigma}
# ---------------------------------------------------------------------------
MATCH = [((1, 2), (3, 4), (5, 6), (7, 8)), ((1, 2), (3, 4), (5, 8), (6, 7)), ((1, 2), (3, 6), (4, 5), (7, 8)),
         ((1, 2), (3, 8), (4, 5), (6, 7)), ((1, 2), (3, 8), (4, 7), (5, 6)), ((1, 4), (2, 3), (5, 6), (7, 8)),
         ((1, 4), (2, 3), (5, 8), (6, 7)), ((1, 6), (2, 3), (4, 5), (7, 8)), ((1, 6), (2, 5), (3, 4), (7, 8)),
         ((1, 8), (2, 3), (4, 5), (6, 7)), ((1, 8), (2, 3), (4, 7), (5, 6)), ((1, 8), (2, 5), (3, 4), (6, 7)),
         ((1, 8), (2, 7), (3, 4), (5, 6)), ((1, 8), (2, 7), (3, 6), (4, 5))]   # transcribed from (A.1), page 35


def polymul(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = tuple(x + y for x, y in zip(ma, mb))
            out[m] = out.get(m, 0) + ca * cb
    return {m: c for m, c in out.items() if c}


def bracket(i, j):          # [ij] = det(w_i, w_j), w = (z, 1): z_i - z_j  (labels 1-based)
    ei = tuple(1 if k == i - 1 else 0 for k in range(8))
    ej = tuple(1 if k == j - 1 else 0 for k in range(8))
    return {ei: 1, ej: -1}


def prodpoly(pairs):
    p = {(0,) * 8: 1}
    for i, j in pairs:
        p = polymul(p, bracket(i, j))
    return p


Bpol = [prodpoly(m) for m in MATCH]
Ppol = [prodpoly([(i + 1, 5 + p[i]) for i in range(4)]) for p in PERMS]
mons = sorted(set().union(*[set(b) for b in Bpol], *[set(q) for q in Ppol]))
Amat = flint.fmpq_mat([[Bpol[a].get(m, 0) for a in range(14)] for m in mons])
AtA = Amat.transpose() * Amat
P = []
okP = True
for q in Ppol:
    rhs = flint.fmpq_mat([[q.get(m, 0)] for m in mons])
    x = AtA.solve(Amat.transpose() * rhs)
    okP &= (Amat * x == rhs)
    row = [x[a, 0] for a in range(14)]
    okP &= all(r.q == 1 for r in row)
    P.append([int(r.p) for r in row])
claim("(4.5): each p_sigma is an integer combination of B_0..B_13 (P re-derived)", okP)
claim("P has rank 14 and last column all ones (P kappa = 1)",
      rank(P) == 14 and all(r[13] == 1 for r in P))
SUB = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18]
claim("Appendix A: 14x14 submatrix of P on rows {0,1,2,3,4,6,7,8,9,10,12,13,16,18} has det 1",
      flint.fmpz_mat([P[i] for i in SUB]).det() == 1)
claim("column space of P = column space of [1|E]: Q is affinely isomorphic to conv{e_sigma}",
      rank([P[s] + E1[s] for s in range(N)]) == 14)

# ---------------------------------------------------------------------------
# Candidate generators: n_{I,J} and h_{pi,eps}
# ---------------------------------------------------------------------------


def n_vals(I, J):
    return tuple(1 if image(p, I) == J else 0 for p in PERMS)


def h_vals(pi, eps):
    out = []
    for s in PERMS:
        t = sum(1 for i in range(4) if s[i] == pi[i])
        S12, S13 = image(s, {0, 1}), image(s, {0, 2})
        d = (S12 == frozenset({pi[0], pi[2]})) + (S12 == frozenset({pi[1], pi[3]})) \
            - (S13 == frozenset({pi[0], pi[1]})) - (S13 == frozenset({pi[2], pi[3]}))
        assert (t + eps * d) % 2 == 0
        out.append((t + eps * d) // 2)
    return tuple(out)


CAND = {}     # value vector -> (type, spec)
for I in SING:
    for J in SING:
        CAND[n_vals(I, J)] = ("n1", f"{min(I)+1},{min(J)+1}")
for I in PAIR:
    for J in PAIR:
        if 0 in I:                       # representative of {I, I^c}
            v = n_vals(I, J)
            assert n_vals(frozenset(range(4)) - I, frozenset(range(4)) - J) == v
            CAND[v] = ("n2", "".join(str(x + 1) for x in sorted(I)) + "," + "".join(str(x + 1) for x in sorted(J)))
for pi in PERMS:
    for eps in (1, -1):
        CAND[h_vals(pi, eps)] = ("h", oneline(pi) + ("," + ("+" if eps == 1 else "-")))
cnt = {}
for v, (t, _) in CAND.items():
    cnt[t] = cnt.get(t, 0) + 1
claim("candidate list: 16 + 18 + 48 = 82 pairwise distinct value vectors",
      len(CAND) == 82 and cnt == {"n1": 16, "n2": 18, "h": 48}, str(cnt))
claim("every candidate value vector lies in W (is a v(c) for some c in V_R)", all(in_W(v) for v in CAND))
claim("every candidate is >= 0 at all 24 vertices", all(min(v) >= 0 for v in CAND))
pat = {}
for v, (t, _) in CAND.items():
    pat.setdefault(t, set()).add((v.count(0), v.count(1), v.count(2)))
claim("zero/one/two counts: n1 (18,6,0), n2 (20,4,0), h (13,10,1)",
      pat == {"n1": {(18, 6, 0)}, "n2": {(20, 4, 0)}, "h": {(13, 10, 1)}}, str(pat))
claim("each candidate's zero set has rank 13 (a hyperplane of aff Q)",
      all(rank([V[s] for s in range(N) if v[s] == 0]) == 13 for v in CAND))
claim("the 82 zero sets are pairwise distinct", len({bits(v) for v in CAND}) == 82)

# ---------------------------------------------------------------------------
# Task 1a: facets by exact double description (cone {c : V c >= 0} in R^14)
# ---------------------------------------------------------------------------


def double_description(A, prune_dim=None):
    """Extreme rays of the pointed cone {u : A u >= 0}, returned as primitive integer
    value vectors y = A u (len(A) entries). A: list of integer rows of rank k = ncols.
    Combinatorial adjacency test; if prune_dim is given (cone known full dimensional
    of that dimension at every stage), also uses |Z| >= prune_dim - 2."""
    m, k = len(A), len(A[0])
    R = []
    for i in range(m):
        if any(A[i]) and rank([A[r] for r in R] + [A[i]]) == len(R) + 1:
            R.append(i)
        if len(R) == k:
            break
    assert len(R) == k
    inv = flint.fmpq_mat([A[r] for r in R]).inv()
    rays = []
    for j in range(k):
        u = [inv[r, j] for r in range(k)]
        den = 1
        for x in u:
            den = den * int(x.q) // gcd(den, int(x.q))
        ui = [int(x * den) for x in u]
        rays.append(primitive([sum(a * b for a, b in zip(A[i], ui)) for i in range(m)]))
    K = 0
    for r in R:
        K |= 1 << r
    for i in range(m):
        if K >> i & 1 or not any(A[i]):
            K |= 1 << i
            continue
        pos = [y for y in rays if y[i] > 0]
        neg = [y for y in rays if y[i] < 0]
        zer = [y for y in rays if y[i] == 0]
        zs = [bits(y) & K for y in rays]
        new = []
        for p in pos:
            zp = bits(p) & K
            for q in neg:
                Z = zp & bits(q)
                if prune_dim is not None and popc(Z) < prune_dim - 2:
                    continue
                adj = True
                for y, zy in zip(rays, zs):
                    if y is p or y is q:
                        continue
                    if Z & zy == Z:
                        adj = False
                        break
                if adj:
                    new.append(primitive([p[i] * b - q[i] * a for a, b in zip(p, q)]))
        rays = pos + zer + new
        K |= 1 << i
    return sorted(set(rays))


t = time.time()
DD = double_description(V, prune_dim=14)
note(f"double description: {len(DD)} rays in {time.time()-t:.1f}s")
claim("DD: the cone {c : v_sigma(c) >= 0} has exactly 82 extreme rays (= facets of Q)", len(DD) == 82)
claim("DD rays coincide with the 82 candidates (16 n1 + 18 n2 + 48 h)", set(DD) == set(CAND))

# ---------------------------------------------------------------------------
# Task 1b: completeness certificate A -- exhaustive search over all 13-subsets
# ---------------------------------------------------------------------------
t = time.time()
found = set()
n_indep = 0
Vint = V
for S in itertools.combinations(range(N), 13):
    X, nul = flint.fmpz_mat([Vint[s] for s in S]).nullspace()
    if nul != 1:
        continue
    n_indep += 1
    c = [int(X[r, 0]) for r in range(14)]
    sgn = 0
    ok = True
    for s in range(N):
        val = sum(a * b for a, b in zip(Vint[s], c))
        if val:
            if sgn == 0:
                sgn = 1 if val > 0 else -1
            elif (val > 0) != (sgn > 0):
                ok = False
                break
    if ok:
        y = primitive([sgn * sum(a * b for a, b in zip(Vint[s], c)) for s in range(N)])
        found.add(y)
note(f"exhaustive 13-subsets: {comb(24,13)} subsets, {n_indep} affinely independent, "
     f"{len(found)} distinct supporting hyperplanes, {time.time()-t:.1f}s")
claim("exhaustive check over all C(24,13) vertex 13-subsets finds exactly the 82 candidates",
      found == set(CAND))

# ---------------------------------------------------------------------------
# Symmetry group G (order 1152) acting on the 24 vertices
# ---------------------------------------------------------------------------


def compose(a, b):          # (a o b)(x) = a(b(x))
    return tuple(a[b[x]] for x in range(4))


def inverse(a):
    r = [0] * 4
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)


gens = []
for rho in [(1, 0, 2, 3), (1, 2, 3, 0)]:
    gens.append(tuple(IDX[compose(rho, s)] for s in PERMS))        # sigma -> rho sigma
    gens.append(tuple(IDX[compose(s, rho)] for s in PERMS))        # sigma -> sigma rho'
gens.append(tuple(IDX[inverse(s)] for s in PERMS))                 # sigma -> sigma^{-1}
ident = tuple(range(N))
G = {ident}
frontier = [ident]
while frontier:
    nf = []
    for g in frontier:
        for h in gens:
            gh = tuple(h[g[x]] for x in range(N))
            if gh not in G:
                G.add(gh)
                nf.append(gh)
    frontier = nf
G = sorted(G)
explicit = set()
for rho in PERMS:
    for rho2 in PERMS:
        explicit.add(tuple(IDX[compose(compose(rho, s), rho2)] for s in PERMS))
        explicit.add(tuple(IDX[compose(compose(rho, inverse(s)), rho2)] for s in PERMS))
claim("G generated as a permutation group of the 24 vertices has order 1152 and equals "
      "{s -> r s r', s -> r s^-1 r'}", len(G) == 1152 and set(G) == explicit)


def act_vals(g, y):         # (g.y)_{g(s)} = y_s
    out = [0] * N
    for s in range(N):
        out[g[s]] = y[s]
    return tuple(out)


def act_set(g, m):
    r = 0
    for s in range(N):
        if m >> s & 1:
            r |= 1 << g[s]
    return r


FSET = set(CAND)
okG = all(act_vals(g, y) in FSET for g in G for y in CAND)
okType = all(CAND[act_vals(g, y)][0][0] == CAND[y][0][0] for g in G for y in CAND)
claim("G permutes the 82 facets (G acts by coordinate permutations of the indicator embedding)", okG)
claim("G maps minor vectors to minor vectors and h-vectors to h-vectors", okType)
forb = {}
for y in CAND:
    forb.setdefault(frozenset(act_vals(g, y) for g in G), []).append(y)
note("facet orbits under G: sizes " + str(sorted(len(o) for o in forb)))

# ---------------------------------------------------------------------------
# Task 1c: completeness certificate B -- ridge closure + connectivity of the facet graph
# ---------------------------------------------------------------------------
t = time.time()
FB = {y: bits(y) for y in CAND}


def ridges_of(y):
    Zl = [s for s in range(N) if y[s] == 0]
    if len(Zl) == 13:
        return [sum(1 << s for s in Zl if s != r) for r in Zl]
    out = set()
    for S in itertools.combinations(Zl, 12):
        X, nul = flint.fmpz_mat([V[s] for s in S]).nullspace()
        if nul != 2:
            continue
        # a functional vanishing on S that is not a multiple of y on Z(y)
        for col in range(2):
            c = [int(X[r, col]) for r in range(14)]
            g = [sum(a * b for a, b in zip(V[s], c)) for s in Zl]
            if any(g):
                break
        if min(g) >= 0 or max(g) <= 0:
            out.add(sum(1 << s for s, gv in zip(Zl, g) if gv == 0))
    return sorted(out)


RID = {}
for orb in forb:
    rep = min(orb)
    rr = ridges_of(rep)
    for y in orb:
        g = next(g for g in G if act_vals(g, rep) == y)
        RID[y] = [act_set(g, m) for m in rr]
ok_r, edges, nr = True, {}, 0
for y in CAND:
    for m in RID[y]:
        nr += 1
        cont = [z for z in CAND if FB[z] & m == m]
        if len(cont) != 2 or y not in cont:
            ok_r = False
        else:
            other = cont[0] if cont[1] == y else cont[1]
            edges.setdefault(y, set()).add(other)
seen, stack = {DD[0]}, [DD[0]]
while stack:
    a = stack.pop()
    for b in edges.get(a, ()):
        if b not in seen:
            seen.add(b)
            stack.append(b)
ridge_counts = {}
for y in CAND:
    ridge_counts.setdefault(CAND[y][0], set()).add(len(RID[y]))
note(f"ridges per facet by type {ridge_counts}; {nr} facet-ridge incidences, {time.time()-t:.1f}s")
claim("every ridge of every listed facet lies in exactly two listed facets", ok_r)
claim("the facet-ridge graph on the 82 listed facets is connected", len(seen) == 82)

# ---------------------------------------------------------------------------
# Task 3: Appendix D numbering
# ---------------------------------------------------------------------------
APPD = ["1,4", "12,34", "13,34", "14,34", "1,3", "12,24", "4,1", "3,1", "12,23", "14,23",
        "13,24", "13,23", "14,24", "2,1", "1,2", "12,14", "2,4", "3421,+", "3412,-", "2413,-",
        "2431,+", "4,2", "3412,+", "14,12", "4312,+", "3142,-", "4132,-", "3,2", "3421,-", "13,12",
        "4321,-", "4123,+", "3124,+", "12,13", "4321,+", "4312,-", "2,3", "2341,+", "2314,-", "14,13",
        "4213,+", "4,3", "4123,-", "2413,+", "2143,-", "13,14", "3241,-", "3,4", "3142,+", "2143,+",
        "2341,-", "1423,-", "1342,+", "13,13", "4231,-", "3,3", "4132,+", "2431,-", "2134,+", "14,14",
        "3214,+", "4,4", "3124,-", "2134,-", "2314,+", "1432,+", "1324,-", "2,2", "4213,-", "3241,+",
        "4231,+", "3214,-", "12,12", "1243,+", "1432,-", "1423,+", "1234,-", "1342,-", "1234,+", "1243,-", "1324,+",
        "1,1"]          # transcribed from Appendix D, page 59 (checked against the page image)
SPEC2V = {sp: v for v, (t, sp) in CAND.items()}
claim("Appendix D lists 82 specs, each one of the candidates, bijectively",
      len(APPD) == 82 and len(set(APPD)) == 82 and all(s in SPEC2V for s in APPD))
F = [SPEC2V[s] for s in APPD]
FTYPE = [CAND[v][0] for v in F]
MINOR = [t != "h" for t in FTYPE]
FZ = [bits(v) for v in F]

# ---------------------------------------------------------------------------
# Task 2: replay of the supplement's completeness certificate (data only)
# ---------------------------------------------------------------------------
args = [a for a in sys.argv[1:]]
out_json = None
lane_npy = None
if "--json" in args:
    i = args.index("--json"); out_json = args[i + 1]; del args[i:i + 2]
if "--lane-npy" in args:
    i = args.index("--lane-npy"); lane_npy = args[i + 1]; del args[i:i + 2]
supp = args[0] if args else os.environ.get("MOC4_SUPP")
task2 = {}
if not supp:
    claim("task 2: supplement location given (argv or MOC4_SUPP)", False, "SKIPPED")
else:
    D = json.load(open(os.path.join(supp, "data", "root_permutation_polytope_complete_20261022.json")))
    FD = json.load(open(os.path.join(supp, "data", "facet-definitions.json")))
    claim("supp data: permutations_1_based is the lexicographic order",
          [tuple(x - 1 for x in p) for p in D["permutations_1_based"]] == PERMS)
    claim("supp data: matching_basis_1_based equals (A.1) as transcribed",
          [tuple(tuple(e) for e in m) for m in D["matching_basis_1_based"]] == MATCH)
    claim("supp data: 'vertices' equals P re-derived above from (4.5)", D["vertices"] == P)
    Fd = D["facets"]
    okf = all([sum(a * b for a, b in zip(P[s], f["normal"])) for s in range(N)] == f["values"] and
              f["zero_vertices"] == [s for s in range(N) if f["values"][s] == 0] for f in Fd)
    claim("supp data: each facet's 'values' = P*normal and 'zero_vertices' consistent", okf)
    claim("supp data: the 82 facet value vectors equal the auditor's 82 (as a set, primitive)",
          {primitive(f["values"]) for f in Fd} == set(CAND))
    same_order = [primitive(f["values"]) for f in Fd] == F
    note(f"supp data: facet order in the JSON equals Appendix D order: {same_order}")
    defs = sorted(FD["definitions"], key=lambda d: d["index"])
    claim("facet-definitions.json agrees with the Appendix D transcription",
          [d["index"] for d in defs] == list(range(82)) and [d["spec"] for d in defs] == APPD and
          all((d["kind"] == "h") == (FTYPE[d["index"]] == "h") for d in defs))
    S = D["boundary_simplices"]
    Sset = {tuple(s) for s in S}
    claim("6846 boundary simplices, distinct, sorted 13-element vertex sets",
          len(S) == 6846 and len(Sset) == 6846 and all(len(s) == 13 and list(s) == sorted(set(s)) for s in S))
    one = []
    for s in S:
        m = sum(1 << x for x in s)
        one.append([j for j in range(82) if FZ[j] & m == m])
    claim("each simplex lies in exactly one of the 82 hyperplanes", all(len(o) == 1 for o in one))
    per = [0] * 82
    for o in one:
        if o:
            per[o[0]] += 1
    claim("every listed facet contains at least one simplex", min(per) >= 1,
          f"simplices per facet: min {min(per)}, max {max(per)}")
    sumv = [sum(P[s][k] for s in range(N)) for k in range(13)]
    W24 = [[24 * P[s][k] - sumv[k] for k in range(13)] for s in range(N)]
    dets = [int(flint.fmpz_mat([W24[x] for x in s]).det()) for s in S]
    claim("every det[24(v_i - b)] (13x13, sorted vertex order) is nonzero", all(d != 0 for d in dets))
    npos = sum(d > 0 for d in dets)
    note(f"determinant signs in sorted order: {npos} positive, {len(dets)-npos} negative; "
         f"|det| values: {sorted(set(abs(d) for d in dets))[:8]}...")
    bd = {}
    for s, d in zip(S, dets):
        sg = 1 if d > 0 else -1
        for k in range(13):
            f = s[:k] + s[k + 1:]
            bd[tuple(f)] = bd.get(tuple(f), 0) + sg * (-1) ** k
    nz = sum(1 for v in bd.values() if v)
    claim("boundary of the oriented chain z = sum sgn(det) [sorted simplex] is zero", nz == 0,
          f"{len(bd)} distinct 12-faces, {nz} with nonzero coefficient")
    claim("number of distinct boundary 12-faces is 44,499", len(bd) == 44499, str(len(bd)))
    inc = {}
    for s in S:
        for k in range(13):
            f = tuple(s[:k] + s[k + 1:])
            inc[f] = inc.get(f, 0) + 1
    note("faces by number of incident simplices: " + str(sorted({v: sum(1 for w in inc.values() if w == v)
                                                                for v in inc.values()}.items())))
    # unoriented control: the chain with all coefficients +1 is NOT a cycle (orientation is essential)
    bd0 = {}
    for s in S:
        for k in range(13):
            f = tuple(s[:k] + s[k + 1:])
            bd0[f] = bd0.get(f, 0) + (-1) ** k
    note(f"control: with all signs +1 in sorted order, {sum(1 for v in bd0.values() if v)} faces fail to cancel")
    task2 = {"simplices": len(S), "faces": len(bd), "nonzero_boundary": nz, "positive_dets": npos,
             "json_facet_order_equals_appendix_D": same_order, "simplices_per_facet": per}

# ---------------------------------------------------------------------------
# Task 4: orbits of G on six-element subsets
# ---------------------------------------------------------------------------
t = time.time()
orbit_of = {}
orbits = []
for S6 in itertools.combinations(range(N), 6):
    m = sum(1 << x for x in S6)
    if m in orbit_of:
        continue
    orb = {act_set(g, m) for g in G}
    k = len(orbits)
    for x in orb:
        orbit_of[x] = k
    rep = min(tuple(setof(x)) for x in orb)
    orbits.append({"rep": list(rep), "size": len(orb)})
note(f"orbit enumeration {time.time()-t:.1f}s")
claim("C(24,6) = 134596 six-subsets", len(orbit_of) == 134596)
claim("G has exactly 214 orbits on six-subsets", len(orbits) == 214, str(len(orbits)))
claim("orbit sizes sum to 134596 and divide 1152",
      sum(o["size"] for o in orbits) == 134596 and all(1152 % o["size"] == 0 for o in orbits))


def FT_of(m):
    return [j for j in range(82) if FZ[j] & m == m]


for o in orbits:
    m = sum(1 << x for x in o["rep"])
    o["F_T"] = FT_of(m)
    o["n_minor"] = sum(1 for j in o["F_T"] if MINOR[j])
    o["n_nonminor"] = sum(1 for j in o["F_T"] if not MINOR[j])
    o["minor_only"] = o["n_nonminor"] == 0
    o["nonminor_generators"] = [j for j in o["F_T"] if not MINOR[j]]
# invariance check of the minor-only flag inside each orbit (sample: all members)
inv_ok = True
for k, o in enumerate(orbits):
    m = sum(1 << x for x in o["rep"])
    for g in G[:: 37]:
        mm = act_set(g, m)
        f2 = FT_of(mm)
        if len(f2) != len(o["F_T"]) or sum(1 for j in f2 if not MINOR[j]) != o["n_nonminor"]:
            inv_ok = False
claim("|F_T| and number of non-minor generators are constant on orbits (checked on 32 group elements per orbit)", inv_ok)
nmo = sum(o["minor_only"] for o in orbits)
claim("124 orbits are minor-only and 90 contain a non-minor generator", nmo == 124 and 214 - nmo == 90,
      f"{nmo} minor-only, {214-nmo} non-minor")
APPE = [[0, 1, 2, 3, 6, 7], [0, 1, 2, 3, 6, 8], [0, 1, 2, 3, 7, 8], [0, 1, 4, 5, 7, 8], [1, 2, 3, 4, 6, 8],
        [0, 1, 2, 4, 6, 8], [0, 1, 2, 4, 7, 8], [0, 2, 5, 6, 7, 8], [0, 1, 2, 3, 11, 12], [0, 1, 3, 5, 7, 8],
        [0, 1, 3, 5, 6, 8], [0, 2, 3, 5, 9, 12], [1, 2, 3, 4, 8, 12], [0, 2, 3, 5, 11, 12], [1, 2, 3, 4, 6, 7],
        [0, 2, 5, 6, 9, 12], [0, 1, 4, 5, 9, 12], [0, 1, 3, 6, 7, 8], [0, 1, 4, 6, 7, 8], [0, 1, 5, 6, 7, 8],
        [0, 3, 5, 9, 11, 12], [2, 3, 4, 8, 10, 12], [2, 3, 5, 9, 11, 12], [2, 3, 4, 6, 7, 8], [3, 4, 7, 8, 10, 12],
        [2, 3, 5, 7, 8, 12], [2, 3, 4, 7, 8, 12], [2, 3, 5, 6, 9, 12], [0, 4, 5, 6, 7, 8], [2, 5, 9, 12, 13, 14],
        [0, 1, 5, 6, 9, 12], [0, 2, 3, 7, 11, 12], [0, 3, 6, 7, 11, 12], [2, 3, 5, 9, 10, 12], [2, 3, 4, 6, 9, 12],
        [2, 3, 4, 6, 11, 12], [0, 5, 6, 7, 9, 12], [1, 2, 3, 8, 11, 12], [2, 3, 5, 7, 9, 12], [0, 1, 5, 9, 14, 18],
        [0, 1, 3, 6, 11, 12], [1, 2, 3, 6, 7, 8], [1, 2, 4, 6, 7, 8], [1, 4, 5, 6, 7, 8], [1, 2, 3, 7, 11, 12],
        [3, 4, 5, 9, 10, 12], [1, 3, 5, 6, 7, 8], [0, 2, 5, 9, 12, 13], [0, 3, 5, 6, 9, 12], [1, 3, 5, 6, 10, 11],
        [1, 3, 5, 10, 11, 12], [1, 3, 5, 10, 11, 14], [0, 1, 4, 6, 9, 12], [2, 3, 5, 6, 11, 12], [1, 2, 3, 6, 11, 12],
        [3, 4, 5, 8, 9, 12], [0, 1, 5, 7, 9, 12], [0, 3, 5, 7, 9, 12], [3, 4, 5, 6, 9, 12], [2, 3, 5, 7, 11, 12],
        [1, 2, 3, 9, 11, 12], [3, 5, 7, 8, 9, 12], [1, 2, 3, 10, 11, 12], [0, 3, 5, 9, 14, 18], [0, 4, 5, 6, 9, 12],
        [2, 4, 5, 6, 9, 12], [2, 4, 5, 7, 9, 12], [1, 4, 5, 7, 9, 12], [0, 2, 5, 7, 9, 12], [1, 3, 4, 8, 9, 12],
        [1, 3, 4, 10, 11, 12], [1, 3, 4, 7, 8, 12], [0, 4, 5, 7, 9, 12], [3, 4, 5, 7, 8, 12], [3, 5, 7, 8, 10, 12],
        [0, 1, 3, 10, 11, 12], [0, 1, 3, 10, 11, 16], [1, 4, 5, 6, 9, 12], [1, 3, 4, 8, 11, 12], [1, 3, 4, 6, 8, 12],
        [1, 2, 10, 11, 12, 13], [2, 5, 9, 10, 12, 13], [3, 5, 6, 9, 14, 18], [2, 5, 9, 11, 12, 13], [0, 5, 7, 9, 14, 18],
        [4, 5, 8, 9, 12, 13], [4, 5, 7, 9, 12, 14], [1, 5, 6, 10, 13, 14], [2, 5, 9, 10, 13, 14], [1, 5, 6, 9, 14, 18]]
# transcribed from Appendix E, page 60 (T_1..T_90)
claim("Appendix E: 90 six-element sets transcribed", len(APPE) == 90 and all(len(set(T)) == 6 for T in APPE))
Torb = [orbit_of[sum(1 << x for x in T)] for T in APPE]
claim("T_1..T_90 lie in 90 distinct orbits", len(set(Torb)) == 90)
nonminor_orbits = {k for k, o in enumerate(orbits) if not o["minor_only"]}
claim("the orbits of T_1..T_90 are exactly the 90 non-minor orbits", set(Torb) == nonminor_orbits)
nu_of_orbit = {k: nu + 1 for nu, k in enumerate(Torb)}
cases = []
for nu, T in enumerate(APPE, 1):
    m = sum(1 << x for x in T)
    ft = FT_of(m)
    cases.append({"nu": nu, "T": T, "orbit": Torb[nu - 1], "orbit_rep": orbits[Torb[nu - 1]]["rep"],
                  "orbit_size": orbits[Torb[nu - 1]]["size"], "F_T": ft,
                  "generators": [{"j": j, "spec": APPD[j], "type": FTYPE[j]} for j in ft],
                  "n_minor": sum(MINOR[j] for j in ft), "n_nonminor": sum(not MINOR[j] for j in ft)})
hist = {}
for c in cases:
    hist[c["n_nonminor"]] = hist.get(c["n_nonminor"], 0) + 1
note("non-minor generators per case (count: number of cases): " + str(sorted(hist.items())))
note("orbit sizes of the 90 cases: " + str(sorted({c['orbit_size'] for c in cases})))

# ---------------------------------------------------------------------------
# Task 5: seven-element supersets and Table 5
# ---------------------------------------------------------------------------
for c in cases:
    T = c["T"]
    fine, further = [], []
    for j in range(N):
        if j in T:
            continue
        m = sum(1 << x for x in T) | (1 << j)
        ft = FT_of(m)
        if all(MINOR[k] for k in ft):
            fine.append(j)
        else:
            other = []
            for x in T:
                k6 = orbit_of[m & ~(1 << x)]
                other.append(nu_of_orbit.get(k6))
            further.append({"j": j, "nonminor": [k for k in ft if not MINOR[k]], "F": ft,
                            "T_nu_of_six_subsets_T_minus_x_plus_j": other})
    c["sup7_minor_only"] = fine
    c["sup7_needs_argument"] = further
claim("logical point: every 6-subset of a contact set lies in one of the 214 orbits, and N_S is contained in N_T for T in S",
      len(orbit_of) == comb(24, 6))
# Table 3 (page 31): the seven new contacts j of T_85 and the symmetry reductions used there
T85 = APPE[84]
TAB3 = [(1, 40, "1234", "1234", 0, [0, 1, 5, 9, 14, 18]), (4, 55, "1243", "1234", 0, [0, 4, 5, 7, 9, 18]),
        (6, 40, "1234", "2143", 0, [5, 6, 7, 9, 14, 18]), (11, 55, "2143", "1243", 1, [0, 5, 7, 9, 11, 14]),
        (12, 55, "1234", "2134", 1, [0, 5, 7, 12, 14, 18]), (15, 55, "2134", "2143", 0, [0, 7, 9, 14, 15, 18]),
        (17, 81, "1324", "2341", 1, [0, 5, 7, 9, 17, 18])]
ok3 = True
for j, nu, r1, r2, ep, img in TAB3:
    rho = tuple(int(ch) - 1 for ch in r1); rho2 = tuple(int(ch) - 1 for ch in r2)
    mapped = sorted(IDX[compose(compose(rho, inverse(PERMS[s]) if ep else PERMS[s]), rho2)] for s in APPE[nu - 1])
    ok3 &= mapped == img and j in img and set(img) <= set(T85) | {j}
claim("Table 3 (Section 5): each A_{rho,rho',eps}(T_nu) equals the stated image, a six-subset of T_85+{j} containing j", ok3)
pos_gaps = {2, 3, 8, 10, 13, 16, 19, 20, 21, 22, 23}
claim("Section 5.2: the 11 gaps declared positive plus the 7 new contacts of Table 3 are the 18 non-contacts of T_85",
      sorted(pos_gaps | {r[0] for r in TAB3}) == sorted(set(range(N)) - set(T85)) and not pos_gaps & {r[0] for r in TAB3})
c85 = cases[84]
note("T_85 + j with a non-minor generator: j in " + str([x["j"] for x in c85["sup7_needs_argument"]]) +
     "; minor-only: " + str(c85["sup7_minor_only"]))

TAB5 = [[0, 1, 2, 3, 6, 15, 17], [0, 1, 2, 3, 6, 15, 22], [0, 1, 2, 3, 6, 17, 19], [0, 1, 2, 3, 6, 17, 22],
        [0, 1, 2, 6, 9, 11, 20], [0, 1, 2, 6, 9, 11, 23], [0, 1, 2, 7, 9, 10, 14], [0, 1, 2, 7, 9, 10, 16],
        [0, 1, 2, 7, 9, 10, 17], [0, 1, 2, 7, 9, 10, 23], [0, 1, 2, 7, 9, 16, 17], [0, 1, 2, 7, 9, 16, 23],
        [0, 1, 2, 7, 10, 16, 17], [0, 1, 2, 7, 10, 16, 23], [0, 1, 2, 7, 14, 16, 23], [0, 1, 2, 9, 10, 14, 17],
        [0, 1, 2, 9, 10, 16, 17], [0, 1, 3, 10, 11, 12, 16], [0, 1, 4, 6, 9, 12, 17], [0, 1, 5, 7, 9, 12, 17],
        [0, 1, 6, 9, 16, 17, 22], [0, 2, 3, 7, 11, 12, 23], [0, 3, 5, 6, 9, 12, 17], [0, 3, 7, 11, 12, 16, 20],
        [1, 3, 4, 8, 11, 12, 15], [1, 3, 5, 6, 7, 8, 14], [1, 4, 5, 7, 9, 12, 14]]   # Table 5, pages 62-63
aux = []
for j, S7 in enumerate(TAB5, 1):
    m7 = sum(1 << x for x in S7)
    subs = []
    for x in S7:
        k = orbit_of[m7 & ~(1 << x)]
        subs.append({"drop": x, "orbit": k, "minor_only": orbits[k]["minor_only"], "T_nu": nu_of_orbit.get(k)})
    ft = FT_of(m7)
    aux.append({"j": j, "S": S7, "subsets": subs, "F_S": ft, "nonminor": [k for k in ft if not MINOR[k]],
                "T_nu_contained": sorted({s["T_nu"] for s in subs if s["T_nu"]})})
claim("Table 5: 27 seven-element sets, each contains a six-subset from a non-minor orbit",
      len(aux) == 27 and all(a["T_nu_contained"] for a in aux))
claim("Table 5: each auxiliary seven-set itself has a non-minor generator in N_S",
      all(a["nonminor"] for a in aux), str([a["j"] for a in aux if not a["nonminor"]]))

# ---------------------------------------------------------------------------
# Task 6: Lemma 4.2 by independent double description on each N_T
# ---------------------------------------------------------------------------
t = time.time()
ok6, dims = True, {}
for o in orbits:
    T = o["rep"]
    X, nul = flint.fmpz_mat([V[s] for s in T]).nullspace()
    K = [[int(X[r, c]) for c in range(nul)] for r in range(14)]
    A = [[sum(V[s][r] * K[r][c] for r in range(14)) for c in range(nul)] for s in range(N)]
    rays = double_description(A) if nul else []
    rays = {r for r in rays if any(r)}
    want = {F[j] for j in o["F_T"]}
    o["dim_W_T"] = nul
    o["dim_N_T"] = rank([list(F[j]) for j in o["F_T"]]) if o["F_T"] else 0
    o["dd_rays"] = len(rays)
    if rays != want:
        ok6 = False
    dims[o["dim_N_T"]] = dims.get(o["dim_N_T"], 0) + 1
note(f"Lemma 4.2 DD on 214 representatives: {time.time()-t:.1f}s; dim N_T histogram {sorted(dims.items())}")
claim("Lemma 4.2: for all 214 representatives the extreme rays of N_T (own DD) are exactly {f_j : j in F_T}", ok6)

# ---------------------------------------------------------------------------
# Task 7: cross-check with the facet data of the minors/moc lane (numeric data, read only)
# ---------------------------------------------------------------------------
task7 = {}
if lane_npy:
    import numpy as np
    Fl = np.load(lane_npy, allow_pickle=False)
    PAIRS3 = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    VR = []
    for p in PERMS:
        M1 = [[1 if p[i] == j else 0 for j in range(4)] for i in range(4)]
        Nn = [[(image(p, S) == frozenset(Tt)) + (image(p, S) == frozenset(Tc)) for (Tt, Tc) in PAIRS3]
              for (S, Sc) in PAIRS3]
        VR.append([M1[i][j] for i in range(3) for j in range(3)] + [Nn[i][j] for i in range(2) for j in range(2)])
    claim("lane reduced coordinates: [VERTS|1] has rank 14 and the same column space as [1|E] (exact)",
          rank([r + [1] for r in VR]) == 14 and rank([VR[s] + [1] + E1[s] for s in range(N)]) == 14)
    vals = np.array(VR, dtype=float) @ Fl[:, :-1].T + Fl[:, -1]
    zs = set()
    maxdev = 0.0
    for col in range(vals.shape[1]):
        v = -vals[:, col]
        z = np.abs(v) < 1e-9
        m = sum(1 << s for s in range(N) if z[s])
        zs.add(m)
        match = [y for y in CAND if FB[y] == m]
        if match:
            y = np.array(match[0], dtype=float)
            scale = v[~z].min() / y[~z][y[~z] > 0].min() if (y[~z] > 0).any() else 1
            maxdev = max(maxdev, float(np.abs(v / scale - y).max()))
    claim("lane facets4.npy: 82 rows whose zero sets are exactly the 82 facets (float data, tol 1e-9)",
          Fl.shape[0] == 82 and zs == set(FB.values()), f"rows={Fl.shape[0]}, max rel. deviation of values {maxdev:.2e}")
    task7 = {"rows": int(Fl.shape[0]), "same_zero_sets": zs == set(FB.values()), "max_value_deviation": maxdev}
else:
    note("task 7 skipped (no --lane-npy)")

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
if out_json:
    res = {
        "facets": [{"j": j, "spec": APPD[j], "type": FTYPE[j], "values": list(F[j]),
                    "zero_set": setof(FZ[j])} for j in range(82)],
        "permutations_one_line": [oneline(p) for p in PERMS],
        "orbits": [{"k": k, "rep": o["rep"], "size": o["size"], "F_T": o["F_T"], "minor_only": o["minor_only"],
                    "n_nonminor": o["n_nonminor"], "dim_N_T": o["dim_N_T"], "T_nu": nu_of_orbit.get(k)}
                   for k, o in enumerate(orbits)],
        "cases": cases,
        "table5": aux,
        "task2": task2,
        "task7": task7,
        "fails": FAILS,
    }
    with open(out_json, "w") as fh:
        json.dump(res, fh, indent=1)
print(f"total time {time.time()-T0:.1f}s; {len(FAILS)} failing claims", flush=True)
sys.exit(1 if FAILS else 0)
