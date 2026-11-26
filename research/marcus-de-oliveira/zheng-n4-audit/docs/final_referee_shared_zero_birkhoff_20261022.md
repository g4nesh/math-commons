# Audit: final_referee_shared_zero_birkhoff_20261022

Evidence level: EXACT for all finite identities, and for (1) as a polynomial identity on 8 exact unitaries (including s = 0 and s = 1) and on all 16 relabeled shared pairs; PROOF-READ for the general argument; CITED only for the classical Birkhoff and Caratheodory theorems.

Document: `proofs/manin_transfer_20261021/final_referee_shared_zero_birkhoff_20261022.pdf`, 3 pages, "A shared spectral value in order four: a Birkhoff construction".
Replay: `replay/final_referee_shared_zero_birkhoff_20261022__check.py`, output in `replay/final_referee_shared_zero_birkhoff_20261022__check_output.txt` (all EXACT and NUMERIC checks pass).

## Deviations and limits

- The general theorem is a statement about every U in U(4). The generic algebra is replayed symbolically: the marginal identities (2) imply (6), and the permutation-matrix identity holds for all 24 sigma. The unitary-specific facts (Cauchy-Binet marginals, the complementary-minor identity) are PROOF-READ in general and replayed EXACTLY on instances: 5 Gaussian-rational unitaries, 1 rational SO(4), one with s = 0 and one with s = 1.
- Birkhoff's theorem and Caratheodory's theorem are classical and used as stated. The face dimension 8 was computed exactly.
- This document matches the paper's Theorem B.1 and its proof (pp. 37-40) step for step. The paper's B.1 is therefore covered by the same checks.
- Theorem 2 (failure of a uniform law for n >= 5) is not used for MOC(4). It was still checked.

## A. Statements

- **Theorem 1 (p.1)** = paper Theorem B.1. For every U in U(4) there is a probability law (p_sigma) on S_4 that depends only on U, such that for all d, a1, a2, a3, c1, c2, c3 in C, with a0 = c0 = d,
  det(diag(a) - U diag(c) U^*) = sum_sigma p_sigma prod_i (a_i - c_sigma(i)).
  The law has at most 10 nonzero weights and satisfies sum_{sigma(0)=0} p_sigma = |u00|^2. Corollary: MOC(4) holds for any normal A, B such that some eigenvalue of A is the negative of some eigenvalue of B. In the paper's notation, with c = -b, this means a_i + b_j = 0 for some i, j. All other eigenvalues are arbitrary complex numbers, and repeated eigenvalues are included. By the relabeling (checked exactly for all 16 pairs (i0, j0)), the shared pair can be any pair.
- **Theorem 2 (p.3).** "One law for all remaining labels" holds for n <= 4 and fails for every n >= 5: take U = O ⊕ I_{n-4} with the displayed O, and the shared pair a0 = c4 = 0. This is not a MOC statement; the document notes that pointwise MOC is not disproved.

No contact sets appear.

## B. Dependencies

- Birkhoff-von Neumann theorem (classical; external).
- Caratheodory's theorem (classical; external). It is used only for the bound of at most 10 weights.
- Cauchy-Binet and the unitary complementary-minor identity |det U[I,J]| = |det U[I^c,J^c]| (classical; PROOF-READ).
- No other document.

## C. Steps (Theorem 1)

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | translate by d to get a0 = c0 = 0 | identity | PROOF-READ | both sides depend only on differences |
| 2 | row marginal sum_j m_ij = s - \|u_i0\|^2 | identity (Cauchy-Binet) | PROOF-READ + EXACT on 8 U | Gram determinant of rows 0, i restricted to N |
| 3 | column marginal sum_i m_ij = s - \|u_0j\|^2 | identity | PROOF-READ + EXACT on 8 U | \|u0k col_l - u0l col_k\|^2 = \|u0k\|^2 + \|u0l\|^2 |
| 4 | total = 2s | identity | EXACT on instances, PROOF-READ | |
| 5 | W in (3) is doubly stochastic | identity | PROOF-READ + EXACT | |
| 6 | Birkhoff decomposition; W00 = 0 forces pi(0) != 0 | citation + logic | CITED (Birkhoff) / PROOF-READ | exact greedy decomposition computed |
| 7 | phi is an involution on {pi(0) != 0}, and m(phi(pi)) is the lower-right 3x3 block of P_pi | finite | EXACT | all 18 pi |
| 8 | (4): sum p_sigma m_ij(sigma) = m_ij; sum_{sigma(0)=0} p = 1 - s = \|u00\|^2 | finite | EXACT on instances, PROOF-READ | |
| 9 | (5): principal-minor expansion plus Cauchy-Binet with a0 = c0 = 0 | identity | PROOF-READ; EXACT as a consequence on instances | which (I,J) survive: 0 in I and 0 not in J; the 3x3 term uses complementary minors |
| 10 | (2) + (5) => (6) with the displayed H_ij | identity | EXACT | symbolic in 9 free m_ij and 6 labels |
| 11 | (6) at P_sigma: prod (a_i - c_sigma(i)) = sum H_ij m_ij(sigma) | finite | EXACT | all 24 sigma; m(sigma) = squared minors of P_sigma also checked |
| 12 | averaging gives (1) for all six labels at once | logic | PROOF-READ + EXACT | (1) replayed as a polynomial identity on 8 unitaries |
| 13 | s = 0: the determinant vanishes identically; any sigma fixing 0 works | logic | PROOF-READ + EXACT | instance diag(1, V) |
| 14 | face {W00 = 0} has dimension 8, so at most 9 Birkhoff terms, +1 gives <= 10 | Caratheodory | CITED / dimension EXACT | rank of 17 difference vectors = 8 |
| 15 | MOC consequence: A = V D_a V^*, B = W D_b W^*, U = V^* W, c = -b, shared index ordered first | logic | PROOF-READ | |

## C'. Steps (Theorem 2)

| # | claim | status | how checked |
|---|---|---|---|
| 1 | n <= 3: Birkhoff on (\|u_ij\|^2); order-2 minors in order 3 are complementary entries | PROOF-READ | |
| 2 | O orthogonal; minor values 4/49, 144/2401 (x2), 400/2401 (x2), d_BB = 1; G(U) = -158/2401 | EXACT | det O = 1 |
| 3 | G(sigma) >= 0 for every sigma in S_n | PROOF-READ (general n) + EXACT (n = 4, 5, 6, 7, exhaustive) | the combinatorial argument was read and is correct |
| 4 | each minor used is (-1)^\|J\| times the coefficient of a surviving monomial under a0 = c4 = 0 | EXACT (n = 5) | coefficient extraction from the symbolic determinant |
| 5 | contradiction by averaging | PROOF-READ | |

## D. Falsification

- (1) with a floating-point law on 300 Haar-random complex unitaries, 5 random label sets each (seed 7): max relative error 2.5e-14. NUMERIC.
- The exact instances above are exact, not numeric.

## E. Verdict

**confirmed-modulo-citations.** The only citations are the classical Birkhoff and Caratheodory theorems, whose hypotheses I checked. Every other step was replayed exactly or read and found correct. No gaps.

Note for the orchestrator: the displayed O is not entrywise equal to `O_STAR` in the shared library, which stores the paper's matrix from Example 1.3. Rows 2-4 differ. The two produce the same value -158/2401, but I did not check whether they are related by a symmetry. Nothing in this audit depends on that.
