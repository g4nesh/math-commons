# Certificate replay for root_moc4_T37_uniform_20261022 (paper case T_38)

Evidence level: EXACT for the rebuilt radial factor G9, all 17 chart positivity certificates (15 direct charts plus the two lower-b charts with their factor identities), denominator positivity, and the square absorptions; PROOF-READ (argument written out below, auditor-supplied in part) for the covering of the hard region by the 17 charts, with an EXACT grid sanity test; CITED for the capacities (2)-(5), the four zero-face cases T_4, T_8, T_11, T_46 and Theorem 1.2. These citations are unchanged from the first audit.

Auditor label cert-T37. This continues the first audit (`docs/root_moc4_T37_uniform_20261022.md`), which returned could-not-check because the chart certificates had not been replayed. Replay: `replay2/root_moc4_T37_uniform_20261022__certificates.py`. Output: `replay2/root_moc4_T37_uniform_20261022__certificates_output.txt` (33 PASS, 0 FAIL, about 40 s).

## Deviations and limits

- None in scope. Every certificate the document names was rebuilt and checked. G9 was rebuilt from the paper's cubic C (A.2) via lib `cubic_C` and the generator vectors lib `F_VECS`. It was not taken from the supplement's tables. The supplement's stored degree-9 factor was loaded only as DATA, for an equality comparison after the rebuild.
- For the square absorptions I used my own coefficient allocation, not the supplement's stored one. Each negative coefficient m sits at monomial x^e. I took all of the neighbour mass R at x^(e+dir) and L = m^2/(4R) at x^(e-dir). This allocation differs from the document's (for example, the document uses 1/2 J^2(J+1)U^2C^2Q^2(B-1)^2 in the first chart). Both are valid. Validity only needs m^2 <= 4LR with L, R taken from available nonnegative coefficients, and the remainder nonnegative.
- I did not re-audit the capacities (2)-(5), the zero-face cases, or the radial necessity / Q != 0 steps. The first auditor did those (EXACT / PROOF-READ / CITED). I re-read the radial necessity polarization argument and agree (see section E).
- The supplement's sha256 digests of the cleared polynomials are not available as data (the checker's JSON report is not in the copy). So each cleared polynomial is matched to the document by term count and by sign pattern, not by hash.

## A. What each certificate asserts

Notation: the source is g = q f1 + p f2 + j f15 + c f9 + e f10 + d f56 + b f0 + a f53 (document eq. (1)), with d = 1 and q, j, e, a > 0, p, c, b >= 0. The hard region H is the set where the capacities (2)-(5) all fail:
- pj < 1;
- c = 0 or p < 1 + 1/c;
- j < 1 + 1/(a+e);
- b = 0 or j < 1/e + 1/b.

G9 is defined by Q^4 det Hess C(g - (C/Q) kappa) = -4 Q^2 G9^2, normalized so that its d^9 coefficient is +1. Here C(g + h kappa) = C + hQ.

Each chart certificate has the following form. A rational map Phi from a parameter domain in (T,J,U,A,C,B,Q) into (q,p,j,c,e,b,a) is given. Write each coordinate as n_i/den_i. The cleared polynomial is sum_m g_m prod_i n_i^(m_i) den_i^(deg_i - m_i), where deg_i is the degree of G9|_{d=1} in variable i: (q,p,j,c,e,b,a) degrees (2,2,2,2,2,3,2). This equals G9(Phi) * prod den_i^deg_i. The assertion is that after the stated binomial-square absorptions every coefficient is >= 0, and that some strictly positive monomial uses only variables that are strictly positive on the domain. Together with den_i > 0 on the domain, this gives G9 > 0 on Phi(domain).

## B. How it was rebuilt

1. X ⊂ Sing C sanity check: grad C = 0 at 5 random rational matching vectors (EXACT).
2. w = g + h kappa in the paper's 14 coordinates. The generator coordinates are integers, so the common denominator is 1. C(w) is affine in h: C has 19 terms and Q has 14.
3. The 14x14 Hessian of C at w was formed (linear entries in q,p,j,c,e,d,b,a,h). Its determinant was computed by fraction-free Bareiss elimination in fmpz_mpoly (own code; 19509 terms; degree 3 in h).
4. E = sum_k det_k (-C)^k Q^(4-k) = Q^4 det Hess(g - (C/Q) kappa). Q^2 divides E exactly, and E/Q^2 = -4 G^2 with G an integer polynomial (fmpz_mpoly sqrt checked by squaring back).
5. Normalized with d^9 coefficient +1, G9 has 1154 terms and is homogeneous of degree 9. It equals the negative of the supplement's stored degree-9 factor, coefficient for coefficient. This means the supplement's "skew cubic F" coordinates coincide with the paper's (A.2) coordinates here.
6. Charts: the substitutions were rewritten by hand in lowest terms (listed in the script). Every denominator has positive coefficients and a monomial in domain-strict variables, so it is positive on the domain.

## C. Results per chart (all EXACT)

| # | chart | domain | cleared terms (doc) | negatives (doc) | squares, direction | strict witness |
|---|---|---|---|---|---|---|
| 1 | both small, u0<1: p=T/(1+T), j=J/(1+J), e=U/((1+U)j) | T,C,B>=0; J,U,A,Q>0 | 3186 (3186) | 2 (2) | 2, B | U^2 (coeff 1) |
| 2 | both small, u0=1: e=1/j | same | 1309 (1309) | 2 (2) | 2, B | constant 1 |
| 3 | both small, u0>1, b>=j: e=(1+U)/j, b=j+(e/U-j)B/(1+B) | same | 12522 (12522) | 0 | - | J^3U^3 |
| 4 | j>1, u0<1: j=1+J, p=T/((1+T)j), e=u0/j, a=(1/J-e)A/(1+A) | same | 12416 (12416) | 0 | - | 3A^2 |
| 5 | j>1, u0=1 | same | 2564 (2564) | 0 | - | 6A^2 |
| 6 | j>1, u0>1 simplex: D=1+U+A, e=(1+U/(JD))/j, a=A/(JDj), b=(JD+U)B/(jU(1+B)) | same | 31706 (31706) | 20 (20) | 20, TB | 3JU^4A^2 (as stated) |
| 7 | p>1, u0<1: p=1+J, j=T/((1+T)p), c=C/(J(1+C)) | C,B>=0; T,J,U,A,Q>0 | 6354 (6354) | 0 | - | 2J^2U^2 |
| 8 | p>1, u0=1 | same | 2943 (2943) | 0 | - | 2J^2 |
| 9 | p>1, u0>1: b=eB/((u0-1)(1+B)) | same | 23198 (23198) | 0 | - | 2T^3J^2U^3 |
| 10 | p=1, u0<1: j=J/(1+J) | C,B>=0; J,U,A,Q>0 | 1077 (1077) | 0 | - | 2U^2 |
| 11 | p=1, u0=1 | same | 450 (450) | 0 | - | constant 2 |
| 12 | p=1, u0>1 (b capped) | same | 2843 (2843) | 0 | - | 2J^3U^3 |
| 13 | j=1, u0<1: p=T/(1+T), e=U/(1+U) | T,C,B>=0; U,A,Q>0 | 833 (833) | 0 | - | constant 2 |
| 14 | j=1, u0=1 | same | 296 (296) | 0 | - | constant 10 |
| 15 | j=1, u0>1, b>=1: e=1+U, b=1+B/(U(1+B)) | same | 1359 (1359) | 0 | - | 36U^3 |
| 16 | both small, u0>1, b<j: b=jB/(1+B) | T,U(*),A,C,B,Q>=0; J>0 | Z: 7134 (7134) | 57 (signed) | factor identity | Z >= 10J^6(J+1) |
| 17 | j=1, u0>1, b<1: e=1+U, b=B/(1+B) | as 13 | Z1: 906 | - | factor identity | Z1 >= 10 |

(*) U > 0 on the chart; the bound does not need it.

Chart 16 (EXACT). Z0 is the A^0 part of Z (2330 terms), and ZA = Z - Z0 has nonnegative coefficients. flint factors Z0 as (J+1) * P5 * F11, with 2, 40 and 331 terms. R11 = F11 - UCQ(J+1)^2(TJB+TB+T-J)^2 has 309 terms, all positive. The identity Z = ZA + (J+1) P5 (R11 + UCQ(J+1)^2(TJB+TB+T-J)^2) holds exactly. P5 has positive coefficients including 2J^2, and R11 includes 5J^4. So Z >= (J+1)·2J^2·5J^4 > 0 for J > 0 and all other variables >= 0.

Chart 17 (EXACT). Z1 was computed by direct substitution. Z has J-degree 7, and Z1 = [J^7]Z. P5 and R11 have J-degree exactly 2 and 4, so taking the J^7 coefficient of the factorization gives Z1 = Alead + Plead (Rlead + UCQ(TB-1)^2). All three factors are coefficientwise nonnegative, with constant terms 2 (Plead) and 5 (Rlead). So Z1 >= 10.

Totals. 17 certificates were replayed, with 24 negative coefficients absorbed by squares (2 + 2 + 20), and every term count agrees with the document.

## D. Covering argument (written out; the grid test is EXACT)

Assume a point of H. Since pj < 1, exactly one of the following holds:
- p < 1 and j < 1;
- j > 1 (then p < 1/j < 1);
- p > 1 (then j < 1);
- p = 1 (then j < 1);
- j = 1 (then p < 1).

p = j = 1 is excluded by (2). Put u0 = je and split into u0 < 1, = 1, > 1. In each chamber the conditions of H reduce as follows.
- (3) fails automatically when p <= 1. When p > 1 it is c < 1/(p-1) = 1/J.
- (4) fails automatically when j <= 1. When j > 1 it is a + e < 1/J.
- (5) fails automatically when u0 <= 1, since then 1/e >= j. When u0 > 1 it is b < e/(u0-1).

Inverse maps, so each chart is onto its chamber. For each case the inverse sends a point of the chamber into the chart domain:
- p in [0,1) gives T = p/(1-p) in [0,∞); j in (0,1) gives J = j/(1-j) > 0; u0 in (0,1) gives U = u0/(1-u0) > 0.
- j > 1: T = pj/(1-pj).
- p > 1: T = jp/(1-jp) > 0 and C = cJ/(1-cJ) >= 0, using cJ < 1.
- Capped variables x in [0,cap) are written x = cap·B/(1+B), with B = s/(1-s) and s = x/cap in [0,1).
- Chart 3: cap e/U = (1+U)/(jU) > j because j < 1. The two pieces b in [j, e/U) (chart 3) and [0, j) (chart 16) cover [0, cap). Similarly, at j = 1, [1, 1+1/U) (chart 15) and [0, 1) (chart 17) cover [0, e/(e-1)).
- j > 1 with u0 <= 1: e = u0/j <= 1/j < 1/J, so a in (0, 1/J - e) is covered by chart 4 or 5, with A = s/(1-s) > 0.
- j > 1 with u0 > 1 (chart 6): x = je - 1 > 0, y = ja > 0, s = J(x+y) in (0,1) (equivalent to a + e < 1/J). Then D = 1/(1-s), U = JDx, A = JDy satisfies D = 1+U+A, and the b-cap becomes (JD+U)/(jU) = e/(u0-1).

Boundary values are included: p = 0, c = 0, b = 0, b = j (chart 3), b = 1 at j = 1 (chart 15), and the layers u0 = 1, p = 1, j = 1. Each chart's strict witness avoids T where p = 0 is allowed, and avoids C and B. The other zero resources (q, j, e, a = 0) are the four closed zero faces, which are CITED.

Grid sanity test (EXACT). I took 3000 rational points of H on a grid that includes the layers p = 1, j = 1, u0 = 1, b = j and the zeros p = 0, c = 0, b = 0. For each point I applied the inverse map of the chart chosen by the argument above. In every case the chart parameters lie in that chart's domain, and the chart map returns exactly the original point. All 17 charts were hit, the smallest with 18 points. G9 > 0 holds by exact evaluation at all 3000 points.

## E. Re-examination of the first auditor's PROOF-READ steps that depend on the certificates

- Step 9 (Hessian identity): it was UNCHECKED as a polynomial identity and is now EXACT. Q^2 divides E = Q^4 det Hess exactly, and E/Q^2 = -4 G9^2 with the constant -4. This agrees with the stored factorization (factors of degree 2 and 9, each with multiplicity 2, constant -4).
- Step 10 (case split) and step 11 (chart coverage): confirmed. The written covering argument is above, and the inverse maps were tested exactly.
- Steps 12 and 13 (chart certificates, lower-b factorization): now EXACT.
- Step 14 (conclusion): PROOF-READ, and I agree. On a realizable source, w = Re x lies in g + R kappa. C(w) = 0 and det Hess C(w) = 0 hold by polarization: grad C(x) = grad C(x̄) = 0, and Hess C(x+x̄) annihilates x - x̄, or x itself by Euler's identity when x is real. Q != 0 follows from failure of (4) (first auditor, EXACT). So w = g - (C/Q) kappa and G9(g) = 0. But G9 > 0 on all of H, so H contains no realizable source.

## Verdict

confirmed-modulo-citations. No certificate fails and none is missing. The remaining citations are the capacities (2)-(5) (SOC, old-H one-pair and chain inequalities via canonical role tuples), the zero-face cases paper T_46, T_4, T_8, T_11 (supplement T45, T4, T7, T10), Theorem 1.2, and X ⊂ Sing C (paper App. A, replayed in lib and spot-checked here).
