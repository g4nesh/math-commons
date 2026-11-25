# Exact certificate for the unitary minor inequality (regime t ≥ 1)

- **Task:** `mdo-star-exact-certificate` (request from `mdo-pmsh-4`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** in progress. **Negative result so far: at t = 1 no exact degree-3 certificate exists** (numerical evidence with controls, Step 2). This file is updated as each step lands.
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

## Step 2 (done, numerical with controls): at t = 1 the degree-3 certificate cannot be exact

**Claim.** f := g0 + 2·D[12,12] + D[13,14] is **not** a degree-3 sum of squares on U(4), although f ≥ 0 there and f + 10⁻⁸ is one numerically. The earlier "certificate" in `../pmsh/equivalence.md` (slack about −10⁻⁷) was a near-miss, not a certificate with slack 0.

**Argument.**
1. **Zeros force kernels.** If f = Σ_b m_b* G_b m_b on U(4) with G_b ⪰ 0, then at every zero U₀ of f the vectors Re m_b(U₀), Im m_b(U₀) lie in ker G_b.
2. **The zero set is large.** f vanishes exactly where g0 = D[12,12] = D[13,14] = 0. Sampling by local minimization (697 runs, then a walk along the zero set, 1,894 points, polished by Newton on the exact gradient) shows many strata; 22 vanishing patterns appear in the first sample alone (`zeros.py`, `zero_walk.py`, `zero_components.py`, `polish.py`).
3. **Forced kernels.** Over the polished zeros the monomial vectors span 494 of the 816 dimensions, with singular values either ≥ 10⁻² or at noise level. The complementary 322-dimensional spaces C_b (cubics vanishing on the zero set) have **integer** row-reduced bases (rounding moved entries by 1.5·10⁻¹⁴), saved in `complement_t1.json` (`reduce.py`).
4. **Reduced problem.** Any exact certificate has the form f = Σ_b (C_b m_b)* S_b (C_b m_b) with S_b ⪰ 0: 482 unknowns, 186 independent conditions. f lies in the linear range of this map (residual 3·10⁻¹⁴), but the best achievable minimum eigenvalue is **λ = −1.69·10⁻⁵ < 0**, with the identity imposed on a rational sample and, independently, on a well-conditioned Haar sample (`reduced_sdp.py`, `reduced_haar.py`).

**Controls (`control.py`).**
- Positive control: F := Σ_b |C_b m_b|² gives λ = +1.0000, as it must.
- Mixture: f + ε·F has λ = −1.6·10⁻⁵, −6.9·10⁻⁶, +8.3·10⁻⁵, +9.8·10⁻⁴ for ε = 10⁻⁶, 10⁻⁵, 10⁻⁴, 10⁻³. So f enters the reduced cone only after adding about 1.7·10⁻⁵·F.
- Fresh zeros: at 60 newly sampled and polished zeros every complement polynomial vanishes to 2.8·10⁻¹⁵.
- Tight tolerances: the unreduced bound stays at about −7·10⁻⁹ (f scaled to maximum 1) when the solver tolerance is tightened to 10⁻¹² (`tight.py`).

**Evidence level and limits.** Numerical, not a proof of non-existence. A proof would need exact zeros spanning the 494 kernel directions and an exact dual certificate for the reduced problem (margin about 10⁻⁵, so this looks feasible).

**What it means.**
- The request's premise fails at t = 1: there is nothing exact to round at degree 3.
- Plausible reason: the inequality (★) is the envelope of the t-family, and it is tight along a whole zero set, so f_t is flat to higher order in some directions near that set. Fixed-degree sums of squares then approximate but do not attain.

## Next steps

5. Same test for t > 1 at degree 3 (the zero set is the same).
6. Same test at degree 4 for t = 1.
