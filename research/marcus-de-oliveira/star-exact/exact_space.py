"""Exact description of torus-invariant, conjugation-symmetric functions of bidegree (3,3) on U(4).

W   : real vector space with basis  e_{p,q} = Re( conj(m_p) * m_q ),  p <= q degree-3 monomials in the entries of U
      with the same torus weight (row multiset, column multiset).  dim W = 1488.
K   : subspace spanned by  Re( rho * conj(m_p') * m_q' )  where rho runs over homogeneous unitarity relations
      (row/column orthogonality, equality of row and column norms) and m_p', m_q' over degree-2 monomials.
      Every element of K vanishes identically on U(4).
ev_U: evaluation at a unitary U, a linear functional on W that vanishes on K.

If  rank K = r_K  and the evaluations at rational unitaries have rank r_E  with  r_K + r_E = dim W,  then the restriction
of W to U(4) has dimension exactly r_E, K is the full kernel of the restriction, and an element of W vanishes on U(4)
iff it vanishes at r_E rational unitaries whose evaluations are independent.  All arithmetic here is exact."""
import itertools, sys, time, random
from fractions import Fraction as Fr
import flint

DEG = 3
def weight(m): 
    r = [0] * 4; c = [0] * 4
    for v in m: r[v // 4] += 1; c[v % 4] += 1
    return (tuple(r), tuple(c))
MONS = list(itertools.combinations_with_replacement(range(16), DEG))
MIDX = {m: i for i, m in enumerate(MONS)}
BLOCKS = {}
for i, m in enumerate(MONS): BLOCKS.setdefault(weight(m), []).append(i)
BLOCK_LIST = list(BLOCKS.values())
PAIRS = [(b[i], b[j]) for b in BLOCK_LIST for i in range(len(b)) for j in range(i, len(b))]
PIDX = {pq: k for k, pq in enumerate(PAIRS)}
NW = len(PAIRS)
MONS2 = list(itertools.combinations_with_replacement(range(16), DEG - 1))
def pair_index(P, Q):
    p, q = MIDX[P], MIDX[Q]
    return PIDX[(p, q) if p <= q else (q, p)]

def relations():
    """list of homogeneous bidegree-(1,1) forms vanishing on U(4): each a list of (coef, holomorphic var, antiholomorphic var)"""
    rel = []
    for i in range(4):
        for j in range(4):
            if i != j:
                rel.append([(1, 4 * i + k, 4 * j + k) for k in range(4)])          # rows i, j orthogonal
                rel.append([(1, 4 * k + i, 4 * k + j) for k in range(4)])          # columns i, j orthogonal
    for i in range(3):
        rel.append([(1, 4 * i + k, 4 * i + k) for k in range(4)] + [(-1, 4 * (i + 1) + k, 4 * (i + 1) + k) for k in range(4)])
        rel.append([(1, 4 * k + i, 4 * k + i) for k in range(4)] + [(-1, 4 * k + i + 1, 4 * k + i + 1) for k in range(4)])
    rel.append([(1, k, k) for k in range(4)] + [(-1, 4 * k, 4 * k) for k in range(4)])  # row-1 norm = column-1 norm
    return rel

def kernel_generators():
    """sparse integer rows (dict col -> coef) spanning K"""
    w2 = [weight(m) for m in MONS2]; rows = []
    for rel in relations():
        c0, h0, a0 = rel[0]
        for qi, mq in enumerate(MONS2):
            for pi, mp in enumerate(MONS2):
                # weight of (u_h * m_q) must equal weight of (u_a * m_p)
                if weight(tuple(sorted(mq + (h0,)))) != weight(tuple(sorted(mp + (a0,)))): continue
                row = {}
                for c, h, a in rel:
                    Q = tuple(sorted(mq + (h,))); P = tuple(sorted(mp + (a,)))
                    assert weight(Q) == weight(P)
                    k = pair_index(P, Q); row[k] = row.get(k, 0) + c
                row = {k: v for k, v in row.items() if v}
                if row: rows.append(row)
    return rows

# ---------- exact rational unitaries ----------
TRIPLES = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (4, 3, 5), (12, 5, 13), (15, 8, 17)]
class CQ:
    """Gaussian rational"""
    __slots__ = ('re', 'im')
    def __init__(self, re, im=0): self.re = Fr(re); self.im = Fr(im)
    def __mul__(self, o): return CQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    def __add__(self, o): return CQ(self.re + o.re, self.im + o.im)
    def __sub__(self, o): return CQ(self.re - o.re, self.im - o.im)
    def conj(self): return CQ(self.re, -self.im)
    def abs2(self): return self.re * self.re + self.im * self.im
def rational_unitary(rng, ngivens=9):
    U = [[CQ(1 if i == j else 0) for j in range(4)] for i in range(4)]
    for _ in range(ngivens):
        i, j = rng.sample(range(4), 2)
        a, b, h = rng.choice(TRIPLES); c, s = Fr(a, h), Fr(b, h)
        x, y, z = rng.choice(TRIPLES); ph = CQ(Fr(x, z), Fr(y, z))            # unit modulus
        # right-multiply by the unitary acting on columns i, j:  [[c, -s*conj(ph)], [s*ph, c]]
        for r in range(4):
            ui, uj = U[r][i], U[r][j]
            U[r][i] = ui * CQ(c) + uj * (ph * CQ(s))
            U[r][j] = uj * CQ(c) - ui * (ph.conj() * CQ(s))
    return U
def is_unitary(U):
    for i in range(4):
        for j in range(4):
            s = CQ(0)
            for k in range(4): s = s + U[i][k] * U[j][k].conj()
            if s.re != (1 if i == j else 0) or s.im != 0: return False
    return True
def monomial_values(U):
    flat = [U[v // 4][v % 4] for v in range(16)]; out = []
    for m in MONS:
        z = CQ(1)
        for v in m: z = z * flat[v]
        out.append(z)
    return out
def evaluation_row(U):
    mv = monomial_values(U)
    return [mv[p].re * mv[q].re + mv[p].im * mv[q].im for (p, q) in PAIRS]      # Re(conj(m_p) m_q)
def det2(U, S, T):
    return U[S[0]][T[0]] * U[S[1]][T[1]] - U[S[0]][T[1]] * U[S[1]][T[0]]

if __name__ == '__main__':
    t0 = time.time()
    print("degree-3 monomials: %d ; weight blocks: %d ; dim W = %d" % (len(MONS), len(BLOCK_LIST), NW))
    gens = kernel_generators()
    print("generators of K: %d  [%.0fs]" % (len(gens), time.time() - t0), flush=True)
    M = flint.fmpz_mat(len(gens), NW)
    for r, row in enumerate(gens):
        for k, v in row.items(): M[r, k] = v
    rK = M.rank()
    print("exact rank of K: %d  (dim W - rank K = %d)  [%.0fs]" % (rK, NW - rK, time.time() - t0), flush=True)
    rng = random.Random(20261125); nU = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    Us = [rational_unitary(rng) for _ in range(nU)]
    assert all(is_unitary(U) for U in Us)
    rows = [evaluation_row(U) for U in Us]
    print("exact evaluations at %d rational unitaries done  [%.0fs]" % (nU, time.time() - t0), flush=True)
    # every generator of K vanishes at every unitary (exact check on a subset of unitaries, all generators)
    bad = 0
    for row in rows[:12]:
        for g in gens:
            if sum(v * row[k] for k, v in g.items()) != 0: bad += 1
    print("generators of K not vanishing at the first 12 rational unitaries: %d" % bad)
    E = flint.fmpq_mat(nU, NW, [flint.fmpq(x.numerator, x.denominator) for row in rows for x in row])
    rE = E.rank()
    print("exact rank of the evaluation functionals: %d  [%.0fs]" % (rE, time.time() - t0))
    print("rank K + rank E = %d ; dim W = %d ; equal: %s" % (rK + rE, NW, rK + rE == NW))
