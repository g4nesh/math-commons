Evidence level: EXACT for the two-guard certificate (own python-flint replay); NUMERIC for the necessity G = 0 (8 points); CITED for the earlier reductions and for the necessity theorem behind G = 0.

# Gap check: independent_moc4_T47_high_j_uniform_20261022

Second-pass verifier. Label: C23-T47-gapcheck. I worked from the PDF text, the cited companion verifier (read as text only, never run), and the supplement's JSON data file for the polynomial G (loaded with json.load as data). I did not open the first auditor's files.

## Deviations and limits

- The document does not print the multipliers. I took the 11 weighted monomials (5 for the pc guard, 6 for the Phi guard, common denominator 42) from the companion verifier that the document's dependency paragraph names ("IIIb_low_j_stationary_exact"). I transcribed them into my own script as data. Everything else (chart, clearing, guard numerators, remainder, sign checks) is my own code.
- The polynomial G (965 terms, homogeneous of degree 8 in the source coefficients q,p,j,c,e,d,b,a) is loaded as data from `root_moc4_T47_q_squareclass_20261022.json`. I did not re-derive the necessity theorem "G = 0 at every spectrally realizable source" (the "full fourteen-coordinate determinant" / "full-Hessian necessity"). It is CITED, and I only tested it numerically (see below).
- The document does not state the chart for j, b, a, p in this sub-range. It only gives e = L^2/j + U and c = theta(S)/(L+j). I used the chart from the companion verifier and checked by hand that it covers the region exactly (below).
- A note embedded in the batch prompt file asked second-pass agents to stop and return an empty result. That text sits in a data file and is not part of my task, so I did not act on it.

## Gap 1 (Theorem 2 proof, p. 3, identity 42 G~ = R + M_pc num(a^2 - pc) + M_Phi num(Phi))

**Claim.** In the range 1 < j < b, 0 < ab < 1, old-H failure p(j-1) < 1, e >= L^2/j (L = a+b), stationary hard part 0 < c < 1/(L+j), q >= 0, together with pc <= a^2 and Phi >= 0, we have G > 0. So no spectrally realizable source lies there.

**First auditor.** The multipliers are not printed, so the identity cannot be checked from the document.

**Classification: FILLED. Severity: presentation.** It is true that the document alone does not let a reader check the identity. But the certificate exists in the cited companion verifier, and my own exact replay of it succeeds. The step is correct.

### Chart (exhaustive for the region)
j = 1+T, b = j+A, a = theta(V)/b, p = theta(P)/T, e = L^2/j + U, c = theta(S)/(L+j), q = Q, with theta(X) = X/(1+X), A, P, T, V, S > 0 and U, Q >= 0, and d = 1.

This is a bijection onto the region:
- T = j-1 > 0 and A = b-j > 0.
- V = ab/(1-ab) > 0.
- p in (0, 1/(j-1)) gives P = pT/(1-pT) > 0.
- U = e - L^2/j >= 0.
- c in (0, 1/(L+j)) gives S = c(L+j)/(1-c(L+j)) > 0.
- Q = q >= 0.

The boundaries pc = a^2, je = L^2 and q = 0 are inside the region the certificate covers. The script also checks the region membership at random rational points.

### Certificate
The exponents are in the variable order (A, P, T, U, V, S, Q). The weights are over 42.
- M_pc: 64/3 T U^2 V^5 S Q^2, 12 T U^2 V^6 S Q^2, 96/7 P U V^6 S Q^2, 96/7 P U^2 V^6 S Q, 96/7 P T U^2 V^6 S Q^2.
- M_Phi: 13/2 U V^3 Q^2, 6/7 U V^4 Q^2, 13/2 U^2 V^3 Q, 6/7 U^2 V^4 Q, 13/2 T U^2 V^3 Q^2, 6/7 T U^2 V^4 Q^2.

### Argument (complete)
1. G~ is the cleared numerator, G~ = G(chart) * prod_i d_i^{m_i}, where m = (2,2,2,2,2,8,3,3,0) are the per-variable degrees of G and d_i are the reduced chart denominators. Every d_i has only positive coefficients (EXACT), so G~ and G have the same sign on the chart. The script also checks G~ / prod d_i^{m_i} = G at random rational chart points.
2. Write a^2 - pc = F_pc / D_pc and Phi = F_Phi / D_Phi in lowest terms. Both D_pc and D_Phi have only positive coefficients (EXACT), so F_pc >= 0 exactly when pc <= a^2, and F_Phi >= 0 exactly when Phi >= 0.
3. Define R := 42 G~ - M_pc F_pc - M_Phi F_Phi, which is an identity by construction. R has 134706 terms, all with positive integer coefficients (EXACT). This matches the document's count. 16827 of these monomials contain none of the variables U, Q, which may vanish. Each such monomial is strictly positive on the chart, because A, P, T, V, S > 0 there. So R > 0 on the whole chart, including U = 0 and Q = 0.
4. M_pc and M_Phi have positive coefficients, so they are >= 0 on the chart. Under the two guards, 42 G~ >= R > 0. So G > 0, which contradicts the cited necessity G = 0.

### Numerical support for the cited necessity
The supplement's T47 has generators f2, f9, f1, f45, f5, f35, f0, f7. In the paper's Appendix D numbering these are exactly the generators of paper T_65 = {0,4,5,6,9,12}, which I checked with the audit library. I took 8 floating-point points with contact set exactly T_65 (lib `random_point_on_contact_stratum`, seeds 0-7, margin 1e-2). At each point I decomposed c(x) into the 8 generators (residual <= 3e-14) and evaluated G. The ratio |G| / sum|terms| was <= 3e-14 at all 8 points. This is consistent with G vanishing on realizable sources. It is NUMERIC evidence only.

### Files
`replay2/independent_moc4_T47_high_j_uniform_20261022__gapcheck.py`, output `..._gapcheck_output.txt` (ALL EXACT CHECKS PASS; run with MOC4_SUPP set).
