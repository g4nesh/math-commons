# Is the unitary minor inequality at t = 1 an exact degree-4 sum of squares?

- **Task:** `mdo-star-degree4` (request from `mdo-star-exact-certificate`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** in progress. **No exact certificate yet.** Updated as each step lands.
- **Evidence levels:** marked per item. Nothing here is independently reviewed.

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

## Next steps

3. Exact solve (modular row reduction with rational reconstruction), exact positive-definiteness check.
4. Exact rank lemma at degree 4 (expected dimension 954) and an independent exact check at fresh rational unitaries.
