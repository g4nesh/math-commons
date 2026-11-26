"""Amplitude spaces as a family linear in s = t/(1+t).

Input (committed, exact): the complement polynomials C_b (those vanishing on the zero set) and the amplitude rows found
separately at several t (degree 3: ../star-exact/amplitudes_t*.json ; degree 4: inside ../star-degree4/certificate4_t*.json).
In each block the amplitude space at t is the set of w C_b with w orthogonal to at most one null direction n_b(t).

Claim checked here in exact arithmetic: scaled so that its largest entry is 1, n_b = n0 + s*n1 with the same n0, n1 at
every t available (fitted on two values of s, confirmed on the others).  Then
        A_b(s) = A0 + s*A1 ,   rows  c_e - n_e(s) c_J   (e != J, J the position of the entry 1)
spans the amplitude space at every available s (exact rank check).  Output: family<deg>.json with integer A0, A1."""
import sys, json
from fractions import Fraction as Fr
import flint
def fm(rows): return flint.fmpq_mat(len(rows), len(rows[0]), [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for r in rows for x in r])
def rank(rows): return fm(rows).rank() if rows else 0
def null_direction(C, A):
    """the vector n (unique up to scale) with W n = 0 where A = W C ; None if A and C have the same row space"""
    rc, ra = len(C), len(A)
    assert rank(C) == rc and rank(A) == ra and rank(C + A) == rc, "amplitude rows must lie in the complement space"
    if ra == rc: return None
    assert ra == rc - 1, "more than one null direction in a block: not implemented"
    # W from A = W C : reduce [C^T | A^T] exactly (C has full row rank, so the first rc columns are the pivots)
    n = len(C[0]); R, rk = fm([[C[i][j] for i in range(rc)] + [A[k][j] for k in range(ra)] for j in range(n)]).rref()
    assert rk == rc and all(R[i, i] == 1 for i in range(rc))
    W = [[Fr(int(R[i, rc + k].p), int(R[i, rc + k].q)) for i in range(rc)] for k in range(ra)]
    for row, w in zip(A, W): assert all(sum(w[i] * C[i][j] for i in range(rc)) == row[j] for j in range(len(row)))
    R, rk = fm(W).rref(); assert rk == ra
    piv = [next(c for c in range(rc) if R[r, c] != 0) for r in range(ra)]; f = next(c for c in range(rc) if c not in piv)
    n = [Fr(0)] * rc; n[f] = Fr(1)
    for r, c in enumerate(piv): n[c] = -Fr(int(R[r, f].p), int(R[r, f].q))
    big = max(n, key=abs); return [x / big for x in n]
def linear_family(comp, amps):
    """comp: list of blocks (integer rows); amps: {s: list of blocks}.  Returns per block (A0, A1) with integer entries."""
    ss = sorted(amps); out = []; kinds = {}
    for bi, C in enumerate(comp):
        if not C:
            assert all(not amps[s][bi] for s in ss); out.append(([], [])); continue
        ns = {s: null_direction(C, amps[s][bi]) for s in ss}
        if all(v is None for v in ns.values()): out.append(([list(map(int, r)) for r in C], [[0] * len(C[0]) for _ in C])); kinds['none'] = kinds.get('none', 0) + 1; continue
        assert all(v is not None for v in ns.values()), "a null direction present at some t only"
        s1, s2 = ss[0], ss[1]; n1 = [(x2 - x1) / (s2 - s1) for x1, x2 in zip(ns[s1], ns[s2])]; n0 = [x1 - s1 * y for x1, y in zip(ns[s1], n1)]
        for s in ss: assert [a + s * b for a, b in zip(n0, n1)] == ns[s], "null direction is not linear in s"
        assert all(x.denominator == 1 for x in n0 + n1)
        J = next(i for i in range(len(C)) if n0[i] == 1 and n1[i] == 0)
        A0 = [[int(C[e][j] - n0[e] * C[J][j]) for j in range(len(C[0]))] for e in range(len(C)) if e != J]
        A1 = [[int(-n1[e] * C[J][j]) for j in range(len(C[0]))] for e in range(len(C)) if e != J]
        for s in ss:
            As = [[a + s * b for a, b in zip(r0, r1)] for r0, r1 in zip(A0, A1)]
            assert rank(As) == len(As) == len(amps[s][bi]) and rank(As + amps[s][bi]) == len(As), "family does not span the amplitude space"
        key = 'coordinate dropped' if not any(n1) else 'pair tied with ratio s'; kinds[key] = kinds.get(key, 0) + 1
        out.append((A0, A1))
    return out, kinds
if __name__ == '__main__':
    deg = int(sys.argv[1])
    if deg == 3:
        comp = json.load(open('../star-exact/complement_t1.json'))['complement']; amps = {}
        for name in ('1.5', '2', '3', '5'):
            d = json.load(open('../star-exact/amplitudes_t%s.json' % name)); t = Fr(d['t']); amps[t / (1 + t)] = [[[Fr(q) for q in row] for row in a] for a in d['amplitudes']]
    else:
        comp = json.load(open('../star-degree4/complement4.json'))['complement']; amps = {}
        for name in ('1', '1-2', '1-3', '2-3', '1-5'):
            d = json.load(open('../star-degree4/certificate4_t%s.json' % name)); t = Fr(d['t']); amps[t / (1 + t)] = [[[Fr(q) for q in row] for row in a] for a in d['amplitudes']]
    comp = [[[Fr(q) for q in row] for row in c] for c in comp]
    fam, kinds = linear_family(comp, amps)
    print("degree %d: values of s used: %s" % (deg, ", ".join(str(s) for s in sorted(amps))))
    print("blocks by kind:", kinds, "; amplitudes:", sum(len(a0) for a0, a1 in fam), "; of which depend on s:", sum(1 for a0, a1 in fam for r in a1 if any(r)))
    json.dump({'degree': deg, 'A0': [a0 for a0, a1 in fam], 'A1': [a1 for a0, a1 in fam]}, open('family%d.json' % deg, 'w'), separators=(',', ':'))
    print("null directions are linear in s and the family spans every amplitude space checked: True ; written family%d.json" % deg)
