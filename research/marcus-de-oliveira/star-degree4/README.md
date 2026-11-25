# Is the unitary minor inequality at t = 1 an exact degree-4 sum of squares?

- **Task:** `mdo-star-degree4` (request from `mdo-star-exact-certificate`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** the target is reached, as a computer-assisted result that has **not** been independently reviewed.
  - **Theorem S4 (exact):** for every U ∈ U(4), g0 + 2·D[12,12] + D[13,14] ≥ 0.
  - **Theorem HQ (consequence):** PMSH(4) holds for every Hermitian matrix whose four eigenvalues form a harmonic quadruple. Hence the Marcus–de Oliveira inclusion holds for n = 4 whenever one matrix has such a spectrum on a line and the other is any normal matrix.
  - Everything else (other spectra, the general inequality (★), the conjecture itself) remains open.
- **Evidence levels:** marked per item below.

**Target.** f₁ := g0 + 2·D[12,12] + D[13,14] ≥ 0 on U(4), with an exact certificate. Here D[S,T] = |det U[S,T]|² and g0 = B[1,1] + B[3,3] + D[14,34] − D[12,14]. More generally f_t := t·g0 + (1+t)·D[12,12] + t²·D[13,14].

**Why t = 1.** By Theorem EQ in [`../pmsh/equivalence.md`](../pmsh/equivalence.md), a spectrum with ratio ρ needs the family at t = ρ/(1−ρ) and at its reciprocal. At ρ = 1/2 both are 1. So an exact certificate for f₁ would prove PMSH(4), and with it the Hermitian case of Marcus–de Oliveira for n = 4, for every spectrum with ratio 1/2.

**Background.** At degree 3 no exact certificate exists at t = 1 ([`../star-exact/README.md`](../star-exact/README.md), Step 2).

## Step 1 (done, numerical with a control): the reduced problem at degree 4 is feasible at t = 1

1. **Zeros.** 1,885 zeros of f on U(4), re-polished to machine precision and filtered by the degree-3 vanishing conditions (`zeros4.py`).
2. **Quartics vanishing on the zero set (`complement4.py`).** Over the zeros, the 3,876 degree-4 monomial coordinates split cleanly: 1,571 singular values ≥ 10⁻³ and 2,305 ≤ 10⁻¹⁴. The 2,305-dimensional complement spaces have **integer** row-reduced bases (rounding moved entries by 9·10⁻¹⁵), saved in `complement4.json`.
3. **Reduced problem (`reduced4.py`).** f_t = Σ_b (C_b m_b)* S_b (C_b m_b) with S_b ⪰ λ·I: 5,761 unknowns, 823 independent conditions on a Haar sample of 3,000 unitaries.

| target | best λ | reading |
| --- | --- | --- |
| positive control Σ\|C_b m_b\|² | +1.0000 | pipeline works |
| f₁ (t = 1) | −4·10⁻¹⁰ | zero to solver tolerance: feasible, on the boundary |

For comparison, at degree 3 the same test gave λ = −1.7·10⁻⁵ at t = 1.

**Reading.** An exact degree-4 certificate at t = 1 is plausible. As at degree 3 for t > 1, there must be further forced kernel directions beyond those coming from zeros.

## Step 2 (done, numerical): scan in t, second reduction, strict feasibility at t = 1

**Scan (`center4.py`).** Best λ of the reduced degree-4 problem:

| t | 1 | 0.5 | 0.2 | 0.1 | 0.05 |
| --- | --- | --- | --- | --- | --- |
| λ | −4·10⁻¹⁰ | −6·10⁻¹¹ | −2·10⁻¹¹ | −6·10⁻¹² | −2.5·10⁻⁵ |

So degree 4 looks exact-feasible for t from about 0.1 upward and is clearly infeasible at t = 0.05.

**Central solution at t = 1.** Of 2,305 eigenvalues, 2,249 are ≥ 10⁻³ and 56 are ≤ 10⁻¹⁰. A clean gap.

**Second reduction (`second4.py 1`).** The 56 null directions are rational with entries 0, 1 and 1/2 only (rounding moved entries by 9·10⁻⁶). Here 1/2 = t/(1+t), the same pattern as at degree 3. After removing them: 2,249 amplitudes (`amplitudes4_t1.json`), 5,481 unknowns, 821 independent conditions.

**Strict feasibility.** On these amplitude spaces every block can be kept ≥ **0.012·I** with the identity exact on the Haar sample. So a rational certificate should exist.

## Step 3 (done, exact computations): a rational certificate at t = 1

**The certificate ([`certificate4_t1.json`](certificate4_t1.json), 190 KB).** 2,249 quartic polynomials q = A_b·m_b with integer coefficients, in 853 weight blocks, and rational symmetric matrices T_b with

  g0 + 2·D[12,12] + D[13,14] = Σ_b (A_b m_b)* T_b (A_b m_b),  every T_b positive definite.

**How it was built (`solve4.py 1`).**
1. The identity at 1,100 rational unitaries is a linear system over ℚ: 5,481 unknowns, rank 821.
2. Its reduced row echelon form was computed modulo three primes near 2³⁰ and lifted by rational reconstruction. All entries have numerators and denominators ≤ 969, and a fourth prime confirms every entry.
3. The 4,660 free unknowns were set to the numerical strictly feasible point rounded to denominator 2²⁰; the 821 pivot unknowns follow exactly. The result is within 3.7·10⁻⁴ of the numerical point.
4. All leading principal minors of all T_b are positive (exact rational arithmetic).

**Independent check (`verify4.py certificate4_t1.json 777`).**
- (2) Every T_b is symmetric with positive leading principal minors.
- (3) The identity holds **exactly** at 1,100 fresh rational unitaries (a seed not used in the construction), in exact integer arithmetic.
- (E) The evaluation functionals at those unitaries have rank 954 on the 11,208-dimensional space W₄ of conjugation-symmetric torus-invariant (4,4)-forms (rank modulo a prime, a lower bound over ℚ).

**What is still missing for "on all of U(4)".** As at degree 3, one needs the dimension of W₄ restricted to U(4) to be exactly 954. (E) gives ≥ 954. The matching bound ≤ 954 needs rank K₄ ≥ 10,254 for the space K₄ spanned by unitarity relations times degree-3 monomial pairs. That computation is Step 4.

## Step 4 (done, exact modular computation): the rank lemma at degree 4

**Lemma 4.** Let W₄ be the real space with basis Re(conj(m_P)·m_Q), P ≤ Q degree-4 monomials of equal torus weight (dim 11,208). Let K₄ ⊂ W₄ be spanned by Re(ρ·conj(m′_p)·m′_q) with ρ a homogeneous unitarity relation and m′ degree-3 monomials (43,728 generators). Every element of K₄ vanishes on U(4). Then:
- **rank K₄ ≥ 10,254** (`lemma4.py`): for random matrices A, Ω over 𝔽_p, p = 524287, the 10,254 × 10,254 matrix A·K·Ω is invertible modulo p. All arithmetic is in float64 holding exact integers; the largest intermediate was 8.1·10¹⁴ < 2⁵³. Log: [`logs/lemma4.txt`](logs/lemma4.txt).
- **rank of evaluations ≥ 954** at the 1,100 verification unitaries (`verify4.py`, part (E), modulo 32749).
- 10,254 + 954 = 11,208 = dim W₄.

**Consequence.** Evaluations vanish on K₄, so the two ranks always add up to at most dim W₄. Equality pins both: the restriction of W₄ to U(4) has dimension exactly 954, and an element of W₄ that vanishes at the verification unitaries vanishes on all of U(4). (Ranks modulo a prime are lower bounds for ranks over ℚ, which is the direction needed.)

## Theorem S4 (exact, computer-verified)

> For every U ∈ U(4):  **g0 + 2·D[12,12] + D[13,14] ≥ 0.**

*Proof.* `certificate4_t1.json` gives f₁ = Σ_b (A_b m_b)* T_b (A_b m_b) with T_b positive definite. The identity holds exactly at 1,100 fresh rational unitaries (Step 3, check (3); repeated with a second seed, [`logs/verify4_t1_seed2.txt`](logs/verify4_t1_seed2.txt)). Both sides restrict to elements of W₄ on U(4): the right side by construction, the left side after multiplying each term by squared row norms, which equal 1 on U(4). By Lemma 4 the identity holds on all of U(4). ∎

## Theorem HQ: harmonic spectra

**Definition.** Four distinct real numbers c₁,…,c₄ form a *harmonic quadruple* if their cross-ratio (c₁−c₂)(c₃−c₄) / ((c₁−c₃)(c₂−c₄)) equals −1, 1/2 or 2 (the value depends on the labelling; the three values are the same condition). Equivalently, the three numbers c₁c₂+c₃c₄, c₁c₃+c₂c₄, c₁c₄+c₂c₃ are in arithmetic progression. The condition is invariant under c ↦ αc + β. Examples: (0, 2, 3, 6), (1, 5, 7, 13), (−3, −1, 0, 3).

> **Theorem HQ.** Let H be a 4×4 Hermitian matrix whose eigenvalues form a harmonic quadruple. Then the vector of all principal minors of H lies in the convex hull of the vectors (∏_{t∈T} c_{π(t)})_T, π ∈ S₄. That is, PMSH holds for H.
>
> **Corollary.** Let A, B be 4×4 normal matrices, B with eigenvalues on a line forming a harmonic quadruple, A arbitrary. Then det(A+B) lies in the convex hull of the products ∏(aᵢ + b_σ(i)). This is the Marcus–de Oliveira inclusion for these pairs.

*Proof of Theorem HQ.*
1. **Certificates are points of a polytope.** Write H = U D_c U*. By Steps 1–2 of Theorem EQ in [`../pmsh/equivalence.md`](../pmsh/equivalence.md), a PMSH certificate exists iff the projected point π_p(x(U)) lies in P_p = conv{π_p(x(P_π))}, where p is the vector of the three matching sums. After an affine normalization p = (1, ρ, 0), and ρ is the cross-ratio above.
2. **Facets at the harmonic ratios (`ratio_half.py`, exact).** For ρ ∈ {1/2, 2, −1} the polytope P_p has exactly 184 facets: 88 whose vertex set lies in no exotic facet of the compound polytope, and 96 of the form g_F + 2·D_a + D_b ≥ 0 with (F, a, b) a symmetric image of (g0, D[12,12], D[13,14]). Nothing else.
3. **The 88.** By LP duality each is a nonnegative combination of the coordinates B[s,j] and D[S,T] on the affine span, so it holds at x(U) for every unitary.
4. **The 96.** Each is Theorem S4 applied to P_α U P_β or its transpose.
5. So π_p(x(U)) satisfies every facet inequality and lies in the affine span of P_p. Hence it lies in P_p. ∎

*Proof of the Corollary.* Theorem H1 of [`../minors/pmsh.md`](../minors/pmsh.md): PMSH for H = V*BV gives the inclusion for every normal A, after rotating and translating so that B is Hermitian. ∎

**Evidence level.** Computer-assisted, not independently reviewed.
- Exact rational/integer arithmetic: positive definiteness, the identity at 2,200 rational unitaries (two seeds), the classification at the harmonic ratios.
- Modular arithmetic (valid lower bounds over ℚ): the two ranks in Lemma 4.
- Floating point rounded to integers: the orientation determinants behind step 2 (all are integers of size ≤ 64; deviation ≤ 4·10⁻¹⁴) and the 82 facets of the compound polytope (same method).
- Written arguments: Steps 1–2 of Theorem EQ, LP duality, Theorem H1.

**Novelty and limits.**
- To my knowledge this case is not covered by earlier results: Shitov (2025) needs the non-Hermitian matrix to have a real eigenvalue when n = 4; Drury (2006) and Bebiano–Kovačec–da Providência (1994) need both spectra on lines or circles. Here the second matrix is an arbitrary normal matrix. **The literature was not searched exhaustively**; Kovačec, LAA 289 (1999), was not read.
- Only harmonic spectra are covered. A spectrum with cross-ratio ρ needs f_t at t = ρ/(1−ρ) and at 1/t.
- End-to-end numerical sanity check, not part of the proof (`end_to_end.py`): for three harmonic spectra and 624 unitaries, including the deepest cap unitaries, the certificate LP is feasible to 0, and the same weights reproduce det(D_a + H) for random complex a to 2·10⁻¹³.

## Reproduction

Python 3.9.6, numpy 2.0.2, scipy 1.13.1, sympy 1.14.0, cvxpy 1.7.5 (Clarabel), python-flint 0.6.0. Run from this directory.

```sh
# checks only (about 5 minutes in total):
python3 verify4.py certificate4_t1.json 777        # exact PD, exact identity at 1100 fresh rational unitaries, evaluation rank 954
python3 lemma4.py 10254 1                          # rank K4 >= 10254 (about 4 minutes, 3 GB of memory)
(cd ../pmsh/star && python3 chiro.py) && python3 ratio_half.py 1/2 2 -1    # facets at the harmonic ratios

# construction (needs the zeros produced by ../star-exact: zeros.py, zero_walk.py):
python3 zeros4.py && python3 complement4.py
python3 reduced4.py 1 && python3 center4.py - 1 && python3 second4.py 1
python3 solve4.py 1
```

Outputs of the runs reported here are in [`logs/`](logs/).

## Next steps

5. Independent review by another GitHub owner: rerun the three checks above; read the proof of Theorem HQ and the lemma's argument.
6. More spectra: degree 4 is feasible for t between about 0.1 and 10, so pairs (t, 1/t) of exact certificates would cover more cross-ratios one at a time. A continuum needs certificates that are rational in t.
7. A literature check for Theorem HQ (Kovačec 1999 in particular).
