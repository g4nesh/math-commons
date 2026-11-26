Evidence level: PROOF-READ (the linear-algebra core of Theorem 2, the density and limit argument, the principal-minor corollary); CITED (Howard–Millson–Snowden–Vakil: dim Sec(M) = 11, the polar map on Sec(M) with 9-dimensional image, the generic point of Sec(M) smooth on C, used for Lemma 1, i.e. rank ≤ 12 on the secant cone); EXACT (the instances below, with own code on the paper's cubic (A.2)); NUMERIC (2000 random spectra).

# Audit: root_moc4_global_real_Hessian_inertia_20261022

PDF: `proofs/manin_transfer_20261021/root_moc4_global_real_Hessian_inertia_20261022.pdf` (2 pages). Batch F4-common-lemmas-b.
Replay: `replay/root_moc4_global_real_Hessian_inertia_20261022__check.py`, with output in `..._check_output.txt`: 7/7 EXACT checks pass and 1 NUMERIC check passes.

## Deviations and limits

- Lemma 1 uses the HMSV facts dim Sec(M) = 11, dim D(Sec(M)) = 9, and "the generic point of Sec(M) is smooth on C". I could not re-derive them (no network). They are CITED as quoted, and the source was not re-read. The exact instances below are consistent with Lemma 1 but do not prove it.
- I replaced two HMSV inputs by arguments of my own (see C.3 and C.5): the identification ker H(m) = tangent space (HMSV Cor. 4.2) and the Terracini transversality L ∩ L̄ = 0. Both are needed only on a dense set of parameters.
- The cubic is the paper's C of (A.2). The document's "real fourteen-coordinate matching cubic F" can differ from it at most by a real linear change of coordinates and a scalar. Ranks and inertia (up to swapping p and n under F → −F, which is harmless because the claim is symmetric) are invariant under both.

## A. Statements

M̂ ⊂ V_C is the affine cone X of the paper (dim 6), C = {F = 0}, and H(x) = Hess F(x), which is linear in x.

**Lemma 1 (p. 1).** For every x in the affine cone over the secant variety Sec(M), rank H(x) ≤ 12.

**Theorem 2 (p. 1).** For every m ∈ X (every finite spectra a, b ∈ C^4 and every scale α ∈ C, including repeated eigenvalues, a_i + b_j = 0 and α = 0), the real symmetric matrix A = H(Re m) has at most 6 positive and at most 6 negative eigenvalues. For generic parameters the inertia is exactly (6, 6, 2). Every principal minor of order 12 is ≥ 0. If rank A = 12, every nonzero principal 12-minor is positive and at least one exists.

**Consequence paragraph (p. 2).** Re m is the cost with its height restored: c(x) + Re F_τ(x) κ = Re x. A real vector w with a negative principal 12-minor of H(w) is not of the form Re m. Applications must establish rank lower bounds and density for themselves.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| HMSV Thm 3.4, Cor. 3.9(a), Sect. 3.4, Prop. 3.19: dim Sec(M) = 11, the polar image N has dim 9, the generic secant point is smooth on C | external | Lemma 1 | CITED |
| HMSV Cor. 4.2 (the partials generate the ideal of M) | external | ker H(m) = T_m M̂ | replaced by an auditor argument (C.3) |
| Terracini's lemma + dim Sec = 11 | external | L ∩ L̄ = 0 generically | replaced by an exact instance + polynomial nonvanishing (C.5) |
| X ⊂ Sing C, dim X = 6, the integer matching parametrization | paper / moc4lib | throughout | EXACT (moc4lib) |
| Cauchy interlacing; a real symmetric matrix of rank r has a nonzero principal r-minor | classical | minors | PROOF-READ |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Euler: H(x)x = 2∇F(x) | identity | PROOF-READ | F is a cubic |
| 2 | Lemma 1: rank H ≤ 12 on the secant cone | algebraic geometry | CITED (HMSV) + PROOF-READ (the reduction) | The step from 2-dimensional fibres of D on Sec to a 2-dimensional kernel of H via v − λx/2 is correct. EXACT instances: rank H(m1+m2) = 12 for 6 random rational pairs; contrast: rank 14 at a random point of V. |
| 3 | dim ker H(m) = 6 = T_m X at generic m | algebraic geometry | PROOF-READ (argument supplied by auditor) + EXACT | ∇F vanishes on X, so H(m) kills every tangent vector, and the kernel has dim ≥ 6 at smooth points. rank_C H(m) = 8 at one exact Gaussian-rational instance (and at 5 real ones), so rank ≥ 8 on a Zariski-dense open set. Hence equality holds generically. Cor. 4.2 is not needed. |
| 4 | conjugate pairs (m, m̄) are Zariski dense in X × X | algebraic | PROOF-READ | The parametrization has integer coefficients, so conj(m(a,b,α)) = m(ā, b̄, ᾱ). The image of {(z, z̄)} is dense because z = x + iy, z' = x − iy is an invertible complex-linear change of variables. |
| 5 | L ∩ L̄ = 0 for generic realizable m | algebraic | EXACT instance + PROOF-READ | ker H(m) ∩ ker H(m̄) = ker A ∩ ker B. At one instance rank [A; B] = 14, so a 14-minor of [A; B] is a real polynomial in the 18 real parameters that is not identically zero. It is therefore nonzero on a dense open set. This replaces the Terracini step. |
| 6 | L is totally isotropic for the Hermitian form of A, and L ∩ ker A = 0, hence p, n ≥ 6 | linear algebra | PROOF-READ | z ∈ L gives z*Az + i z*Bz = 0 with both terms real. A 6-dimensional totally isotropic subspace meeting the radical trivially forces min(p, n) ≥ 6. |
| 7 | rank A ≤ 12, since Re m = (m + m̄)/2 is in the secant cone | citation | CITED (Lemma 1) | |
| 8 | generic inertia (6, 6, 2) | conclusion | PROOF-READ | EXACT instances: 6 random Gaussian-rational (a, b, α) all give (6, 6, 2), from the exact characteristic polynomial and Descartes' rule (exact for real-rooted polynomials). |
| 9 | p, n ≤ 6 at every realizable parameter | limit | PROOF-READ | The number of positive eigenvalues is lower semicontinuous, so a limit has at most as many. The good parameters are dense in the Euclidean topology because their complement is a proper real-algebraic set. EXACT degenerate instances: a1 = a2 (4,4,6); a1 = a2 and b3 = b4 (2,2,10); a1+b1 = 0 (5,5,4); two vanishing sums (5,5,4); a = (x,x,x,y) (4,4,6); scalar a (2,2,10); real spectra (4,4,6); a on a circle (6,6,2). |
| 10 | principal 12-minors ≥ 0, and > 0 when nonzero | linear algebra | PROOF-READ + EXACT | Interlacing gives p, n ≤ 6 for the submatrix, so a nonsingular one has inertia (6,6) and det > 0. EXACT: all 91 principal 12-minors are positive at 4 generic instances. |
| 11 | "Homogeneous matching coordinates extend the same conclusion to limits with points at infinity" | limit | PROOF-READ | The set where p, n ≤ 6 is closed, so it holds on all of X. |

## D. Falsification

- NUMERIC: 2000 random complex spectra and scales (seed 31), with a forced coincidence a_1 = a_2 in 30% of samples and a forced a_4 + b_3 = 0 in 30%. Eigenvalue counts (relative tolerance 1e-9) were (6,6) in 1015 samples and (5,5) in 985, the latter mostly in the forced-degenerate ones. No sample had p > 6 or n > 6.
- EXACT: the instances listed in C.2, C.8, C.9 and C.10. No counterexample was found.

## E. Verdict

**confirmed-modulo-citations.** The only citation left is HMSV, for Lemma 1: rank ≤ 12 on the secant cone, which rests on dim Sec(M) = 11, a 9-dimensional polar image, and smoothness of C at the generic secant point. Every other step is PROOF-READ or EXACT. I supplied the Cor. 4.2 and Terracini inputs myself on a dense set. No gap was found.

A note for users: the theorem constrains only Re m (it gives a necessary condition for spectral realizability). Any application that wants strict signs or rank lower bounds on a particular contact stratum must prove them separately, as the document itself says.
