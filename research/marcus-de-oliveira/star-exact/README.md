# Exact certificate for the unitary minor inequality (regime t ≥ 1)

- **Task:** `mdo-star-exact-certificate` (request from `mdo-pmsh-4`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** in progress. **No exact certificate yet.** This file is updated as each step lands.
- **Evidence levels:** marked per item. Nothing here is independently reviewed.

**Target.** For every U ∈ U(4) and every t ≥ 1:

  t·g0 + (1+t)·D[12,12] + t²·D[13,14] ≥ 0,

with D[S,T] = |det U[S,T]|² and g0 = B[1,1] + B[3,3] + D[14,34] − D[12,14]. This is the regime D[12,12] ≥ D[13,14] of inequality (★) in [`../pmsh/equivalence.md`](../pmsh/equivalence.md), which is equivalent to PMSH(4). A numerical degree-3 sum-of-squares certificate exists ([`../pmsh/star/sos_param.py`](../pmsh/star/sos_param.py)); the goal is to make it exact.

## Step 1 (done, exact): an identity test on U(4) without sampling

**Setting.** Let m_p run over the 816 monomials of degree 3 in the entries of U. Two monomials have the same *torus weight* if they use the same multiset of rows and the same multiset of columns. Let W be the real vector space with basis e_{p,q} = Re(conj(m_p)·m_q), p ≤ q of equal weight. dim W = 1488.

**Lemma 1 (exact computation).**
- Let K ⊂ W be spanned by Re(ρ·conj(m′_p)·m′_q), where ρ is a homogeneous unitarity relation (two rows orthogonal, two columns orthogonal, two row norms equal, two column norms equal, a row norm equal to a column norm) and m′_p, m′_q are degree-2 monomials. Every element of K vanishes on U(4). **rank K = 1227.**
- The evaluation functionals at 300 explicit rational unitaries (products of Pythagorean Givens rotations with Gaussian-rational phases) have **rank 261** on W.
- 1227 + 261 = 1488.

**Consequences.**
- The restriction of W to U(4) has dimension exactly 261, and K is the whole kernel of the restriction. (Evaluations vanish on K, so rank K + rank E ≤ 1488 always; equality pins both.)
- **Identity test.** An element of W vanishes on all of U(4) iff it vanishes at 261 rational unitaries with independent evaluations. So a Gram identity f = Σ_blocks m* G m on U(4) can be verified in exact rational arithmetic at finitely many points.
- The numerical rank 261 reported in `../pmsh/equivalence.md` is now an exact dimension.

**Script and output.** `exact_space.py` (python-flint 0.6.0 for exact ranks; Python `fractions` for evaluations), 12 seconds. Output: [`logs/exact_space.txt`](logs/exact_space.txt).

```sh
cd research/marcus-de-oliveira/star-exact && python3 exact_space.py 300
```

## Next steps

2. Read the forced kernels of the Gram matrices off a numerical solution (the certificate is tight, so they are singular) and find rational bases.
3. Solve the reduced problem strictly inside the cone, round, project exactly, and verify positive definiteness exactly. First t = 1, then all t ≥ 1.
