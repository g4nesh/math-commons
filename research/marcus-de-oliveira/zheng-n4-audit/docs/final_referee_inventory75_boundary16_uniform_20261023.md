Evidence level: EXACT for the source, the zero-coefficient reductions, the atom orbit memberships, C0, Q0, (3), the Hessian-determinant factorization (4) and both positive chart arrays (own sympy + python-flint code, 0 failures); PROOF-READ for the logic; CITED for the earlier contact cases, the SO(4) atoms (foundations) and radial necessity; NUMERIC only for supplementary tests.

# final_referee_inventory75_boundary16_uniform_20261023 (Table 5 set j = 18)

PDF: `proofs/manin_transfer_20261021/final_referee_inventory75_boundary16_uniform_20261023.pdf`, 3 pages. This is a branch document, the seven-contact input S16 = {0, 1, 3, 10, 11, 12, 16} for paper T_76 (and listed for T_78).
Replay: `replay/final_referee_inventory75_boundary16_uniform_20261023__check.py` (with `--numeric`; about 15 s for the exact part).

## Deviations and limits

- No deviation. I computed G myself: the exact 14 × 14 determinant of Hess C(w0 + zκ), with C the paper's cubic (A.2), by fraction-free Bareiss elimination in python-flint at d = 1. This is legitimate because both sides of (4) are homogeneous of the same degree. I then extracted the square root and substituted the charts. The supplement's G, and its scripts and JSON, were not used.
- deg_z det Hess C(w0 + zκ) = 3. So the identity holds in the form Q0³·det = −4·(scale)·Q0·G², which is the same as the document's Q0⁴ det = −4d²Q0²G² after multiplying by Q0. The constant is −65536 = −4·2¹⁴, because my matrix is built from 2m to clear the halves of f30.
- The orbit memberships of the atoms were computed with the group action on generator indices. The atoms themselves are proved in the foundations (so4-atoms.md, 3.2 old-H and 3.3 chain). They are cited, not re-proved.

## A. Statements

**Theorem 1 (p. 1).** For every x ∈ X with P c(x) ≥ 0 and (P c(x))_σ = 0 for all σ ∈ S16 = {0,1,3,10,11,12,16}, Φ_{c(x)}(O) ≥ 0 for all O ∈ O(4), and hence for U(4) via Theorem 1.2. Contacts are non-strict. Spectra are arbitrary: the proof never uses an eight-distinct assumption, only x ∈ X through radial necessity. S16 is literally Table 5, j = 18 (EXACT).

Intermediate claims:
- (1) (chain): b ≥ d²(1/q + 1/a) is sufficient.
- (2) (conic sum of two old-H atoms): b ≥ d + d²/(a+j) is sufficient.
- (3): S·C0 − aj·Q0 > 0 for positive coefficients.
- (4): Q0⁴ det Hess C(w0 − (C0/Q0)κ) = −4d²Q0²G².
- G > 0 on both chart domains (6) and (7).

## B. Dependencies

| kind | ref | used for | match checked |
|---|---|---|---|
| paper-case | old T39 = T_40 {0,1,5,9,14,18} (stage 1) | q = 0 (extra tie 20) | yes, T_40 lies in the zero set of the remaining rays |
| paper-case | old T28 = T_29 {0,4,5,6,7,8} (stage 1) | j = 0 (tie 14), c = 0 (tie 13) | yes |
| paper-case | old T43 = T_44 {1,4,5,6,7,8} (stage 1) | e = 0 (tie 6) | yes |
| paper-case | old T30 = T_31 {0,1,5,6,9,12} (stage 1) | b = 0 (tie 7) | yes |
| paper-case | old T10 = T_11 {0,1,3,5,6,8} (stage 3) | a = 0 (tie 5) | yes; stage 3 < stage 4 of T(7)_18 |
| common-lemma | chain atom (foundations so4-atoms.md 3.3), tuple (centre f10; leaves f0, f15; bad f30) | (1) | yes, in the orbit of (2; 0,12; 70) |
| common-lemma | old-H atom (so4-atoms.md 3.2), tuples (aux f15, main f10; f30) and (aux f8, main f10; f30) | (2) | yes, in the orbit of (12, 2; 70) with the printed roles; the reversed roles are not in the orbit |
| common-lemma | radial necessity (C(u) = 0 and Hess C(u) v = 0 for x = u + iv ∈ X) | "spectral realizability forces det Hess = 0" | yes: the document's own derivation (p.2) is the same argument |
| core-doc | paper (3.12): grad C(κ) = 0, C(κ) = 0 | C(w0 + zκ) is affine in z | EXACT in lib; deg_z = 1 replayed |
| core-doc | Lemma 4.2 (complete source) | decomposition | yes |
| external | Theorem 1.2 | unitary version | yes |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | source = f0,f8,f2,f12,f30,f10,f15 for q,j,c,e,d,b,a; exclusive gaps 20,14,13,6,2,7,5; only f30 non-minor (p.1) | finite | EXACT | lib |
| 2 | d = 0 gives a minor sum; each other zero coefficient gives an eight-contact face containing an earlier case (p.1) | citation | EXACT + CITED | zero set of the remaining 6 rays contains the named case |
| 3 | (1) is the chain atom (p.1) | citation | EXACT + CITED | orbit membership with roles |
| 4 | (2) as a conic sum: split d = d_a + d_j, centre portions d_a + d_a²/a and d_j + d_j²/j (p.1) | identity | EXACT | sympy |
| 5 | C(w0 + zκ) = C0 + zQ0 with the displayed C0, Q0 (p.2) | identity | EXACT | ratio 1 to the paper cubic |
| 6 | (3) and S(S+d) − aj = a² + aj + j² + d(a+j) > 0 (p.2) | identity | EXACT | |
| 7 | realizable means C(u) = 0 and det Hess C(u) = 0 for u = w0 + zκ, z = Re F_{π*} (p.2) | argument | PROOF-READ (= radial necessity) | |
| 8 | Q0 ≠ 0 by (3), so z = −C0/Q0 (p.2) | argument | PROOF-READ | |
| 9 | (4): factorization with a perfect square G of degree 8 with 465 terms (p.2) | identity | EXACT | own determinant (6462 terms, checked against a Fraction determinant at 3 random points); square root by flint; 465 terms |
| 10 | failure of (1) and (2) implies (5); charts (6) (q ≥ 1) and (7) (q < 1) are onto supersets of the failure region (p.2) | argument | PROOF-READ | B/(1+B) maps (0,∞) onto (0,1) |
| 11 | chart numerators: 1432 and 2103 terms, every coefficient strictly positive; pure-A terms 8A² (+8A³ + 2A⁴) and A² (+2A³ + A⁴); denominators ((1+Q)A(1+B))² and (1+Q)³((A+J)(1+B))² are positive (p.2) | identity | EXACT | so G > 0 on the whole failure domain, Q = 0 included |
| 12 | contradiction with G = 0 (p.3) | argument | PROOF-READ | |

## D. Falsification

- 4000 random points of the failure domain (seed 20261125, d = 1). min(−det Hess C(w0 − (C0/Q0)κ)) = 5.65 > 0, consistent with G ≠ 0.
- Six numerical realizable points with contacts ⊇ S16 (lib stratum search, margin 0). All have d = g2 ≈ 0 (|d| < 4e-11) and much larger contact sets. C(Re x) ≈ 0, as expected. min Φ ≥ −9e-13. This agrees with the theorem's mechanism: no realizable point with d > 0 in the failure domain was found. It is not evidence of more than that.

## E. Verdict

**confirmed-modulo-citations.** Citations: T_40, T_29, T_44, T_31, T_11 (all earlier in Table 6), the chain and old-H atoms (foundations, established), radial necessity (established), Lemma 4.2, Theorem 1.2. No gap. The decisive inequality is exact: an exact square, with positive-coefficient chart numerators.
