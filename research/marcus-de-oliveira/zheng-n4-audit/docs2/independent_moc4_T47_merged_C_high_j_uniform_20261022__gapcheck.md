Evidence level: EXACT for the two-guard certificate (own python-flint replay); NUMERIC for the necessity G = 0 (8 points); CITED for the old-H, coherent-capacity, large-pc and negative-Heron inputs and for the necessity theorem behind G = 0.

# Gap check: independent_moc4_T47_merged_C_high_j_uniform_20261022

Second-pass verifier. Label: C23-T47-gapcheck. I worked from the PDF text, the cited companion verifier (read as text only, never run), and the supplement's JSON data file for G (json.load, as data). I did not open the first auditor's files.

## Deviations and limits

- The 26 weighted monomials of the certificate (12 for the pc guard, 14 for the Phi guard) are not in the document. I took them from the companion verifier the document names ("IIIb_plus_narrow_stationary_exact") and transcribed them into my own script as data. My script recomputes their common denominator from the weights and gets 9990288000, as stated in the document.
- The polynomial G is loaded as data from `root_moc4_T47_q_squareclass_20261022.json`: 965 terms, homogeneous of degree 8 in the source coefficients q,p,j,c,e,d,b,a. The necessity "G = 0 at every spectrally realizable source" is CITED and was only tested numerically.
- A note embedded in the batch prompt file asked second-pass agents to stop and return an empty result. That text sits in a data file and is not part of my task, so I did not act on it.

## Gap 1 (pp. 2-3, range b <= j < A0 with je >= A0^2, identity 9990288000 G~ = R0 + Mpc Fpc + MPhi FPhi)

**Claim.** On the chart region with pc <= a^2 and Phi >= 0, G > 0.

**First auditor.** The 12- and 14-monomial multipliers are not printed, so the identity cannot be checked.

**Classification: FILLED. Severity: presentation.** The document alone does not contain the certificate, as the first auditor says. The certificate is in the cited companion verifier, and my own exact replay of it passes. The claim is correct.

### Chart (as printed in the document; exhaustive)
b = 1+A, a = theta(V)/b, A0 = a+b, j = b + a theta(T), p = theta(P)/(j-1), e = A0^2/j + U, c = theta(S)/(j+A0), q = Q, with A, P, V, S > 0 and T, U, Q >= 0, and d = 1.

This is a bijection onto the region b > 1, 0 < a < 1/b, b <= j < A0, p(j-1) < 1, e >= A0^2/j, 0 < c < 1/(j+A0), q >= 0:
- theta(T) = (j-b)/a lies in [0,1), so T >= 0.
- j >= b > 1, so the factor j-1 in the p chart is positive.

The boundaries j = b, e = A0^2/j, q = 0 and pc = a^2 are included.

### Argument (complete; every step EXACT in the script)
1. G~ = G(chart) * prod d_i^{m_i}, with m = (2,2,2,2,2,8,3,3,0) and every reduced chart denominator d_i having positive coefficients. So sign G~ = sign G. This is checked against direct evaluation of G at random rational points.
2. a^2 - pc = Fpc / Dpc and Phi = FPhi / DPhi in lowest terms, with Dpc and DPhi having positive coefficients. So the guards give Fpc >= 0 and FPhi >= 0.
3. R0 := 9990288000 G~ - Mpc Fpc - MPhi FPhi has 232323 terms, all with positive integer coefficients. This matches the document. Mpc and MPhi have 12 and 14 terms, all with positive integer coefficients. 2898 monomials of R0 have no T, U, Q, which also matches the document. Such a monomial is strictly positive because A, P, V, S > 0. So R0 > 0 on the whole chart, walls T = 0, U = 0, Q = 0 included.
4. Under the two guards, 9990288000 G~ >= R0 > 0, so G > 0. This contradicts the cited necessity G = 0. So no spectrally realizable source remains in this sub-case.

### Numerical support for the cited necessity
The supplement's T47 generators {f2, f9, f1, f45, f5, f35, f0, f7} are exactly the generators of paper T_65 = {0,4,5,6,9,12} (audit library). At 8 floating-point points with contact set exactly T_65 (seeds 0-7, margin 1e-2), I decomposed c(x) into these generators (residual <= 3e-14). G vanished to relative size <= 3e-14 at every point. NUMERIC only.

### Files
`replay2/independent_moc4_T47_merged_C_high_j_uniform_20261022__gapcheck.py`, output `..._gapcheck_output.txt` (ALL EXACT CHECKS PASS; run with MOC4_SUPP set).
