# Exact certificate for the unitary minor inequality (regime t ≥ 1)

- **Task:** `mdo-star-exact-certificate` (request from `mdo-pmsh-4`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** in progress.
  - **Exact results (Steps 4–5):** for every U ∈ U(4), t·g0 + (1+t)·D[12,12] + t²·D[13,14] ≥ 0 at **t = 3/2, 2, 3 and 5**. Each has an exact rational sum-of-squares certificate and passes an independent exact check.
  - **Negative result (Step 2):** at t = 1 no exact degree-3 certificate exists (numerical evidence with controls).
  - The full target (all t ≥ 1 at once) is **not** done.
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

## Step 3 (done, numerical): for t above 1 an exact degree-3 certificate looks possible, with more forced structure

**Scan in t (`reduced_scan.py`, same complement spaces, Haar sample).** Best minimum eigenvalue λ of the reduced blocks:

| t | 0.5 | 1 | 1.05 | 1.1 | 1.5 | 2 | 3 | 10 | 100 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| λ | −8.8·10⁻³ | −1.7·10⁻⁵ | −3·10⁻¹¹ | −2·10⁻¹⁰ | −2·10⁻⁹ | −7·10⁻⁹ | −1·10⁻⁸ | −2·10⁻⁷ | −4·10⁻⁸ |

- t ≤ 1: clearly infeasible. No exact degree-3 certificate.
- t ≥ 1.05: λ = 0 to solver tolerance. The reduced identity is feasible, but only with singular blocks, so there are further forced kernel directions beyond those coming from zeros.

**Central solution at t = 2 (`center.py`).** Of the 322 eigenvalues of the reduced blocks, 302 are ≥ 10⁻³ and 20 are ≤ 10⁻⁹. A clean gap.

**The 20 extra directions (`extra_kernel.py`).** They are simple: either a single complement polynomial is dropped, or two are tied with ratio 3/2 at t = 2. Example: in the block with u₁₁u₂₄u₃₄ and u₁₄u₂₄u₃₁ the surviving amplitude is u₂₄·(3·u₁₁u₃₄ − 2·u₁₄u₃₁). This suggests amplitudes that depend on t, such as u₂₄·((1+t)·u₁₁u₃₄ − t·u₁₄u₃₁). **That formula is a guess from one value of t.**

## Step 4 (done, exact): a certificate at t = 2

> **Theorem S2.** For every U ∈ U(4):  2·g0 + 3·D[12,12] + 4·D[13,14] ≥ 0.
> Equivalently g0 + (3/2)·D[12,12] + 2·D[13,14] ≥ 0, the member t = 2 of the family.

**Evidence level.** Exact computer-verified certificate in rational arithmetic, with an independent checker. Not reviewed by another GitHub owner.

**The certificate** ([`certificate_t2.json`](certificate_t2.json)). 302 cubic polynomials q = A_b·m_b in 196 weight blocks, with rational coefficients (denominators 1 and 3), and rational symmetric matrices T_b such that

  2·g0 + 3·D[12,12] + 4·D[13,14] = Σ_b (A_b m_b)* T_b (A_b m_b)  on U(4),  every T_b positive definite.

**How it was built.**
1. Second reduction (`second_reduction.py 2`): the 20 extra kernel directions of the central solution were rationalized (rounding moved entries by at most 8.6·10⁻⁵) and removed, leaving 302 amplitudes.
2. Strict feasibility (`reduced2.py 2`): all blocks can be kept ≥ 0.034·I.
3. Exact solve (`exact_cert.py 2`): the identity at 300 rational unitaries is a rational linear system with 434 unknowns and exact rank 184. The 250 free unknowns were set to the rounded numerical values; the others were solved exactly. The result is within 1.1·10⁻⁶ of the numerical point.

**How it is checked (`verify_cert.py certificate_t2.json 777`, exact arithmetic only).**
1. Lemma 1 at 300 **fresh** rational unitaries (a seed not used in the construction): rank K = 1227, rank of evaluations = 261, sum 1488.
2. Every T_b is symmetric with positive leading principal minors (smallest 2.2·10⁻⁴).
3. The identity holds exactly at all 300 fresh unitaries. By Lemma 1 it therefore holds on all of U(4).

Why step 3 suffices: both sides restrict to elements of the 261-dimensional space of Lemma 1. The right side is a combination of the basis functions Re(conj(m_p)·m_q). The left side is too, after multiplying each term by a squared row norm, which equals 1 on U(4).

**What it gives, and what it does not.**
- One rigorous member of the family. By symmetry its 288 images also hold.
- It does **not** give PMSH(4) at any spectrum: a spectrum with ratio ρ needs the family at both t = ρ/(1−ρ) and its reciprocal, and the reciprocal is below 1.
- It shows the pipeline works: zeros → integer complement spaces → second reduction → rational solve → exact check. Each further rational t ≥ 1.05 should take about a minute.

## Reproduction

Python 3.9.6, numpy 2.0.2, scipy 1.13.1, sympy 1.14.0, cvxpy 1.7.5 (Clarabel), python-flint 0.6.0. Run from this directory.

```sh
# exact check only (about 25 seconds each):
for c in certificate_t1.5.json certificate_t2.json certificate_t3.json certificate_t5.json; do python3 verify_cert.py $c 777; done

# full construction at t = 2:
python3 zeros.py 1 700 1 && python3 zero_walk.py 1 5 120 && python3 polish.py 1 zeros_walk_t1.pkl
python3 reduce.py 1                      # complement_t1.json (integer bases)
python3 center.py 2 && python3 second_reduction.py 2 && python3 reduced2.py 2
python3 exact_cert.py 2                  # certificate_t2.json

# the same steps for other values of t, in one command:
python3 pipeline.py 3/2 3 5
python3 t_dependence.py 3/2 2 3 5
```

Outputs of the runs reported here are in [`logs/`](logs/). Intermediate `.pkl` files are not committed.

## Step 5 (done, exact): certificates at t = 3/2, 3, 5, and how the amplitudes depend on t

**More exact certificates (`pipeline.py 3/2 3 5`).** Same pipeline, about a minute each. All three pass `verify_cert.py` with seed 777.

| t | inequality proved for every U ∈ U(4) | strict-feasibility margin λ | smallest leading minor | checker |
| --- | --- | --- | --- | --- |
| 3/2 | 6·g0 + 10·D[12,12] + 9·D[13,14] ≥ 0 | 0.017 | 3.2·10⁻⁵ | VALID |
| 2 | 2·g0 + 3·D[12,12] + 4·D[13,14] ≥ 0 | 0.034 | 2.2·10⁻⁴ | VALID |
| 3 | 3·g0 + 4·D[12,12] + 9·D[13,14] ≥ 0 | 0.064 | 1.4·10⁻³ | VALID |
| 5 | 5·g0 + 6·D[12,12] + 25·D[13,14] ≥ 0 | 0.120 | 9.1·10⁻³ | VALID |

Certificates: `certificate_t1.5.json`, `certificate_t2.json`, `certificate_t3.json`, `certificate_t5.json`.

**Dependence on t (`t_dependence.py`, observed at these four values).** The 20 extra kernel directions have the same pattern for every t:
- 12 are independent of t (a single complement polynomial is dropped).
- 8 involve exactly one number, **r = t/(1+t)**: two complement polynomials p, p′ survive only in the combination (1+t)·p − t·p′.

So the amplitude spaces are explicit for all t: they are linear in t.

**What remains for a whole half-line t ≥ t₀.** The linear system for the Gram blocks then has coefficients that are quadratic polynomials in t. One needs a solution T_b(t), rational in t, with all leading minors positive on [t₀, ∞). That is exact linear algebra over ℚ(t) plus Sturm-type positivity checks. **Not done here.** The margin λ grows with t in the table, and the reduced problem is feasible for t ≥ 1.05, so t₀ slightly above 1 looks reachable; t₀ = 1 is not (Step 2).

## Next steps

8. One certificate for a half-line t ≥ t₀ > 1 (see above).
9. Degree 4 at t = 1 and below, with the same facial-reduction test.
10. Review: rerun `verify_cert.py` on the four certificates; check Lemma 1's argument.
