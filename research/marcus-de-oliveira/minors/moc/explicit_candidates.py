"""Test explicit 'sequential determinantal' constructions for 3-value pattern certificates.
Target: for all S and count vectors q, P(counts of kappa on S = q) = G(S,q) = sum_{T: counts(T)=q} |det U[S,T]|^2."""
import sys, itertools, numpy as np
sys.path.insert(0, 'moc')
from moclib import haar_unitary

def target_stats(U, blocks):
    n = U.shape[0]; lab = np.empty(n, int)
    for j, B in enumerate(blocks): lab[list(B)] = j
    st = {}
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            for T in itertools.combinations(range(n), k):
                q = tuple(np.bincount(lab[list(T)], minlength=3))
                st[(S, q)] = st.get((S, q), 0) + abs(np.linalg.det(U[np.ix_(S, T)]))**2
    return st

def stats_of(lam, n):
    st = {}
    for col, p in lam.items():
        for k in range(1, n + 1):
            for S in itertools.combinations(range(n), k):
                q = tuple(np.bincount([col[i] for i in S], minlength=3))
                st[(S, q)] = st.get((S, q), 0) + p
    return st

def seq_candidate(U, blocks, order, mode):
    """order = (first colour, second, third).  Sample X_first ~ DPP(U[:,J_first]); then on Y = rest,
    2-colour with the DPP of the subspace spanned by U[Y, J_second] (mode 'span') ."""
    n = U.shape[0]; a, b, c = order; Ja, Jb, Jc = blocks[a], blocks[b], blocks[c]
    lam = {}
    for Xa in itertools.combinations(range(n), len(Ja)):
        pa = abs(np.linalg.det(U[np.ix_(Xa, Ja)]))**2
        if pa < 1e-15: continue
        Y = [i for i in range(n) if i not in Xa]
        B = U[np.ix_(Y, Jb)]
        gram = np.linalg.det(B.conj().T @ B).real
        for Xb_idx in itertools.combinations(range(len(Y)), len(Jb)):
            Xb = [Y[t] for t in Xb_idx]
            pb = abs(np.linalg.det(U[np.ix_(Xb, Jb)]))**2 / gram
            col = [None] * n
            for i in Xa: col[i] = a
            for i in Xb: col[i] = b
            for i in Y:
                if i not in Xb: col[i] = c
            lam[tuple(col)] = lam.get(tuple(col), 0) + pa * pb
    return lam

def err(lam, tgt, n):
    st = stats_of(lam, n); keys = set(tgt) | set(st)
    return max(abs(st.get(k, 0) - tgt.get(k, 0)) for k in keys), abs(sum(lam.values()) - 1)

rng = np.random.default_rng(1)
for pattern in [(1, 1, 1), (2, 1, 1), (2, 2, 1), (2, 2, 2)]:
    n = sum(pattern); blocks = []; s = 0
    for p in pattern: blocks.append(list(range(s, s + p))); s += p
    worst = {}
    for trial in range(3):
        U = haar_unitary(n, rng); tgt = target_stats(U, blocks)
        avg = {}
        for order in itertools.permutations(range(3)):
            lam = seq_candidate(U, blocks, order, 'span')
            e, tot = err(lam, tgt, n)
            worst[order] = max(worst.get(order, 0), e)
            for k, v in lam.items(): avg[k] = avg.get(k, 0) + v / 6
        e, tot = err(avg, tgt, n); worst['avg'] = max(worst.get('avg', 0), e)
    print("pattern", pattern, " max stat error per colour order:", {str(k): float('%.2e' % v) for k, v in worst.items()})
