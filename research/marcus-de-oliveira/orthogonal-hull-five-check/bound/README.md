# Independent bound: Φ ≤ 1.1371 on O(5)

- **Context:** independent check of task `mdo-orthogonal-hull-five`. The functional Φ, the unitary V and the claim r₅(U(5)) ⊄ conv r₅(O(5)) are due to the agent `ganesh-claude-minors` (GitHub owner `g4nesh`). This folder re-derives the key inequality with its own certificate, its own generator and its own exact checker. It does not use that task's certificate or code.
- **Check by:** agent `ishaan-air` (GitHub owner `IshaanAyaan`), November 25, 2026.
- **Status:** computer-assisted proof with an exact rational certificate. **Not independently reviewed.**
- **Levels:** EXACT = exact arithmetic or a complete proof; NUMERIC = floating point only.

## Results

> **Theorem B′ (EXACT, computer-assisted).** For every real orthogonal 5×5 matrix O, Φ(O) ≤ 1.13702917… < 1.1371.

> **Proposition A (EXACT, re-checked here).** V is unitary and Φ(V) = 37/25 = 1.48.

> **Corollary (EXACT).** r₅(V) ∉ conv r₅(O(5)), with margin 37/25 − 1.1371 > 0.34. The bound 13/10 claimed in `mdo-orthogonal-hull-five` is confirmed with room to spare.

Φ is the functional of `mdo-orthogonal-hull-five`: for a 5×5 matrix M, with indices in ℤ₅,

  Φ(M) = − Σ_v |m_vv|² + Σ_v ( |det M[{v+1,v−1},{v+2,v−2}]|² + |det M[{v+2,v−2},{v+1,v−1}]|² ) − Σ_v Σ_{d=1,2} ( |det M[{v,v+d},{v,v−d}]|² + |det M[{v,v−d},{v,v+d}]|² ).

## Proof of Theorem B′

Write x_ij (0 ≤ i, j ≤ 4) for the entries of a real 5×5 matrix. Φ is a polynomial of degree 4 in the x_ij. Let g_k be the 30 orthogonality relations Σ_j x_ij x_kj − δ_ik (i ≤ k) and Σ_i x_ij x_il − δ_jl (j ≤ l).

[`certificate.json`](certificate.json) contains a rational number c = 1137021/1000000, 346 rational symmetric matrices G_b with monomial vectors v_b, and 30 rational polynomials h_k of degree at most 4. Define the polynomial

  R(x) = c − Φ(x) − Σ_b v_b(x)ᵀ G_b v_b(x) − Σ_k h_k(x) g_k(x).

[`check_certificate.py`](check_certificate.py) establishes, in exact rational arithmetic:

1. every G_b is positive definite (exact LDLᵀ elimination, all pivots positive);
2. R has 4,804 nonzero coefficients, and the sum of their absolute values is ‖R‖₁ = 8.17·10⁻⁶ (largest coefficient 2.4·10⁻⁷);
3. B := c + ‖R‖₁ = 1.13702917… satisfies B < 1.1371 < 13/10 < 37/25.

Now let O be orthogonal. Then g_k(O) = 0 for all k. Each row of O is a unit vector, so |O_ij| ≤ 1, every monomial has absolute value at most 1 at O, and |R(O)| ≤ ‖R‖₁. Each v_b(O)ᵀ G_b v_b(O) ≥ 0. Hence

  Φ(O) = c − Σ_b v_b(O)ᵀ G_b v_b(O) − R(O) ≤ c + ‖R‖₁ = B. ∎

The relations g_k hold on both components of O(5), so no separate treatment of det O = −1 is needed. (It would also follow from Φ(−O) = Φ(O), since −O has the opposite determinant when n = 5.)

**What the checker builds itself.** It constructs Φ from the definition by 1×1 and 2×2 determinants (and confirms the 35 terms are distinct minors), constructs the 30 relations from their definitions, reads only c, the blocks (v_b, G_b) and the multiplier coefficients from the file, and expands everything. The soundness argument above does not depend on how the blocks were chosen, so the checker does not need to trust the block structure. As a sanity test it prints Φ on the 120 permutation matrices: 1 on 60 and −5 on 60 (EXACT), as stated in `mdo-orthogonal-hull-five`. Tampering tests: flipping the sign of one Gram diagonal entry makes the checker stop (block not positive definite); perturbing one multiplier coefficient by 0.1 raises B to 1.637; replacing c by 1 leaves B at 1.13703, because the difference moves into R.

**How the certificate was built** ([`generate_certificate.py`](generate_certificate.py), not needed for verification). Φ and the relations are invariant under sign changes of single rows and columns, so the monomials of degree ≤ 3 split into 346 blocks by the parities of their row and column degrees (sizes 42, 26, 13, 6, 5, 2). A semidefinite program (cvxpy with Clarabel) minimises c subject to the coefficient identity, with all G_b ⪰ 10⁻⁵·I so that rounding cannot destroy positivity. Its optimum is 1.1370208. Gram entries and multiplier coefficients are rounded to integers over 10⁹ and c is rounded up to 1137021/1000000; nothing else is corrected, and the rounding error is exactly R.

## Connection with the alternating sum E_S

**Lemma (EXACT, complete proof).** For O ∈ O(n) and index sets I, J of equal size,

  det O[I,J]² = det(O) · Σ_{σ : σ(I)=J} sgn(σ) Π_i O_{i,σ(i)}.

*Proof.* The generalised Laplace expansion along the rows I groups the terms of det O by σ(I) = J; the group for J equals (−1)^{ΣI+ΣJ} det O[I,J] det O[Iᶜ,Jᶜ]. Jacobi's complementary-minor identity with O⁻¹ = Oᵀ gives det O[Iᶜ,Jᶜ] = (−1)^{ΣI+ΣJ} det(O) det O[I,J]. Multiply and use det(O)² = 1. ∎

Put v_σ = Φ(P_σ) for the permutation matrix P_σ (entries P_{i,σ(i)} = 1), let S = {σ : Φ(P_σ) = −5} and E_S(O) = Σ_{σ∈S} sgn(σ) Π_i O_{i,σ(i)}. Since Φ is linear in the squared minors and r(P_σ)_{I,J} = 1 exactly when σ(I) = J, the lemma gives on all of O(5)

  Φ(O) = det(O) · Σ_σ sgn(σ) v_σ Π_i O_{i,σ(i)} = det(O)(det O − 6 E_S(O)) = 1 − 6 det(O) E_S(O).

So Theorem B′ is equivalent to E_S(O) ≥ −(B − 1)/6 ≈ −0.02284 on SO(5). The set S has 60 elements, 30 even and 30 odd; it contains the identity and the transposition (0 1) (EXACT: Φ(I) = −5 and Φ(P_(01)) = −5 by hand and by code). The identity Φ = 1 − 6 det(O) E_S was also tested on 200 Haar orthogonal matrices with largest error 7·10⁻¹⁵ (NUMERIC, [`logs/sanity_numeric.txt`](logs/sanity_numeric.txt)).

## Numerical observations (not part of the proof)

| Observation | Level |
| --- | --- |
| Φ(V) = 1.48 and V unitary to 2·10⁻¹⁶ (the exact check is `exact_phi_at_V.py`) | NUMERIC |
| 20 local maximisations of Φ over O(5) all end at 1.0000000000; 20,000 Haar samples give at most 0.9918 | NUMERIC |
| The degree-6 programme above has optimum 1.13702 with margin 10⁻⁵ (consistent with the 1.1365 reported without margin in `mdo-orthogonal-hull-five`) | NUMERIC |
| Restricting to SO(5) and adding the relations O = cof(O) and det O = 1 at the same degree lowers the optimum only to 1.13654 (exploratory run, not committed) | NUMERIC |

The sharp inequality Φ ≤ 1 on O(5) remains open as a theorem. Degree-6 certificates cannot reach it.

## Files

| File | Content |
| --- | --- |
| `certificate.json` | c, 346 Gram blocks, 8,310 multiplier coefficients; integers over the common denominator 10⁹ (0.7 MB) |
| `check_certificate.py` | exact checker, Python standard library only, about 3 seconds |
| `generate_certificate.py` | builds the certificate (numpy, scipy, cvxpy, Clarabel), about 35 seconds |
| `exact_phi_at_V.py` | Proposition A in exact arithmetic (sympy) |
| `sanity_numeric.py` | floating-point checks: Φ(V), the E_S identity, local maxima |
| `logs/` | outputs of the four scripts |

## Reproduction

```sh
python3 check_certificate.py        # Theorem B', exact; Python 3.9+, no packages
python3 exact_phi_at_V.py           # Proposition A, exact; needs sympy
python3 generate_certificate.py     # rebuilds certificate.json
python3 sanity_numeric.py           # numerics
```

Versions used: checker on Python 3.11.7; generator on Python 3.14 with numpy 2.5.5, scipy 1.18.1, cvxpy 1.9.3, Clarabel 0.11.1; sympy 1.14.0.
