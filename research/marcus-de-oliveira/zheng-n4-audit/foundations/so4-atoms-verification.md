Evidence level: EXACT (polynomial identities, covariance and orbit checks by my own code in `so4_atoms_verify.py`, 0 failures) plus PROOF-READ (inequality steps written out below) for SOC, old-H and four, each proved from scratch. Chain: EXACT re-check of the supplement's Gram identity plus PROOF-READ of the Cauchy-Schwarz step; my own from-scratch chain proof is incomplete (GAP; only a NUMERIC multiplier certificate). NUMERIC attack: 0 violations. Two sharpness remarks in `so4-atoms.md` are wrong.

## Deviations and limits

- **Chain is not proved independently.** I reduced it to a cleaner pointwise statement, but I could not close that statement by hand, and no SDP solver is installed. I therefore re-checked the supplement's (PB) Gram identity EXACTLY with my own code. The cubic features, the congruence T and the four sign maps are read from the supplement JSON with `json.load`; nothing is executed. To learn the data layout (which index is which, and the mod-N convention) I read `so4_atoms.py` section A5. Everything else in step 1 was derived without that file.
- **Step 1 independence.** For SOC, old-H and four I wrote my own proofs before reading Sections 3.1 to 3.4 of `so4-atoms.md`; the comparison in section 2 came afterwards.
- **Wrong variants.** For chain, factor 1.8 is not a failing variant: the chain stays nonnegative at factor 1.8, and the numerical threshold lies between 1.550 and 1.575 (section 5). For SOC the pair is symmetric, so "roles swapped" means the pair (1, 9) put with a wrong non-minor facet k.
- **Spot-check.** I checked 16 cited tuples against the PDF text extracts, not all 58.
- **Compute.** Run through RUNJOB in six sections, about 9 minutes in total.

## 0. Setting and reduction method

- O = L_p R_{conj q} Pi, with p, q unit quaternions and Pi a fixed signed permutation (det +1). This covers SO(4).
- P_i = p_i^2, Q_i = q_i^2, D_i = P_i - Q_i, so sum D = 0. 2x2 minors have bidegree (2,2) and Phi has bidegree (4,4).
- **Exact reduction.** Each identity is bihomogenised with |p|^2 and |q|^2 to a common bidegree, for example D_i -> P_i|q|^2 - Q_i|p|^2. It is then compared as a polynomial over Q in python-flint. On S^3 x S^3 the homogenised and original forms agree, so a polynomial identity of the homogenised forms proves the identity on SO(4).
- For the chain certificate, the convention is the supplement's: reduce modulo N = |p|^2 - |q|^2 by the normal form p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2. This is exact because N is monic in p0.

## 1. Own proofs (representatives; transfer to orbits in section 3)

### 1.1 SOC: Phi_22 + 2|x1 x9| >= 0 (Pi = I)

EXACT (E1):
- x1 = D0 + D1 and x9 = D0 + D3;
- Phi_22 = 2K - 2D0D2 - 2D1D3, with K = (p0p2 - p1p3)^2 + (q0q2 - q1q3)^2;
- x1 x9 = D1D3 - D0D2, using sum D = 0.

So Phi_22 + 2|x1x9| = 2K - 4 min(D0D2, D1D3). It remains to show min(D0D2, D1D3) <= K/2. PROOF-READ:
- If either product is <= 0, the claim holds because K >= 0.
- Otherwise both products are > 0. Since sum D = 0, one pair is positive and the other negative. Swapping p and q negates D and fixes K, so take D0, D2 > 0 > D1, D3.
- Put a = |p0p2| - |q0q2| and b = |q1q3| - |p1p3|. Both are > 0, because P0 > Q0, P2 > Q2, Q1 > P1 and Q3 > P3.
- AM-GM gives P0Q2 + Q0P2 >= 2 sqrt(P0P2Q0Q2), so D0D2 <= a^2. In the same way D1D3 <= b^2.
- Put X = |p0p2| - |p1p3| and X' = |q0q2| - |q1q3|. Then K >= X^2 + X'^2 by the reverse triangle inequality, and X - X' = a + b.
- So K/2 >= (X^2 + X'^2)/2 >= (X - X')^2/4 = (a+b)^2/4 >= ab >= min(a^2, b^2) >= min(D0D2, D1D3).

No division occurs, so zero coordinates are covered.

### 1.2 old-H: Phi_2 + Phi_70 + 2|x12 x2| >= 0 (aux f12, main f2; Pi = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])

EXACT (E2):
- x2 = D0 + D2 and x12 = -(D0 + D3);
- Phi_2 = x2^2;
- with K' = (p0q3 - p2q1)^2 + (p1q2 - p3q0)^2, Phi_2 + Phi_70 + 2 sgn·x2x12 equals 2K' - 4D0D1 for one sign and 2K' - 4D2D3 for the other.

So the left side is 2K' - 4 min(D0D1, D2D3). PROOF-READ:
- If the minimum is <= 0 we are done.
- Otherwise sum D = 0 forces opposite signs on {0,1} and {2,3}. K' is invariant under p <-> q (EXACT), so take D0, D1 > 0 > D2, D3.
- Put r = |p|, t = |q| componentwise, X = r0t3 - r2t1 and Y = r1t2 - r3t0. Then K' >= X^2 + Y^2, and X + Y = (r0t3 - t0r3) + (r1t2 - t1r2).
- Aczél, an EXACT identity: (r0t3 - t0r3)^2 - (r0^2 - t0^2)(t3^2 - r3^2) = (r0r3 - t0t3)^2.
- With r0 > t0 and t3 > r3 this gives r0t3 - t0r3 >= sqrt(D0|D3|). Similarly r1t2 - t1r2 >= sqrt(D1|D2|).
- So K'/2 >= (X+Y)^2/4 >= sqrt(D0 D1 D2 D3) >= min(D0D1, D2D3).

### 1.3 four: Phi_54 + Phi_10 + 2(|x0x6| + |x0x10| + |x6x10|) >= 0 (Pi = I)

EXACT (E3), with m03 = p0q3 - p3q0 and m12 = p1q2 - p2q1:
- x0 = m03 - m12, x6 = -(m03 + m12) and x10 = -(D0 + D2);
- Phi_10 = x10^2;
- Phi_54 = x0^2 + x6^2 - XY, with X = (P0 - P3) - (Q1 - Q2) and Y = (P1 - P2) - (Q0 - Q3).

Put, with E = P + Q and L = D:
- lam0 = L0 - L3, lam1 = L1 - L2;
- eps0 = E0 - E3, eps1 = E1 - E2;
- sigma = (eps0 - eps1)/2.

Then (EXACT) XY = (lam0 + lam1)^2/4 - sigma^2 and x10 = (lam0 - lam1)/2 on sum D = 0. Also |x0| + |x6| = 2M with M = max(|m03|, |m12|). So the left side equals

  4M^2 + 2M|lam0 - lam1| + sigma^2 - lam0 lam1.

PROOF-READ:
- **Pair fact.** With z = p_i + i q_i and w = p_j + i q_j, EXACT: |z^2 - w^2|^2 = (|z|^2 - |w|^2)^2 + 4 Im(conj(z) w)^2. Since lam = Re(z^2 - w^2), this gives lam_k^2 <= eps_k^2 + 4m_k^2.
- **Sign lemma.** If lam > 2|m|, then sign(eps) = sign(Lambda) != 0, where Lambda = L_i + L_j.
  - Put s = z + w = S e^{i alpha} and d = z - w = R e^{i beta}. Then lam = SR cos(alpha + beta), eps = SR cos(alpha - beta), 2m = SR sin(alpha - beta) and Lambda = (S^2 cos 2alpha + R^2 cos 2beta)/2.
  - Let A = alpha + beta and B = alpha - beta. From cos A > |sin B| we get cos B != 0.
  - If cos B > 0, then |A| < pi/2 - |B|, so cos(A ± B) > 0 and Lambda > 0.
  - If cos B < 0, replace B by B - pi; then cos(A ± B) < 0 and Lambda < 0.
- **Main argument.** If lam0 lam1 <= 0 there is nothing to prove. The expression is invariant under lam -> -lam (p <-> q) and under swapping the pairs, so take lam0 >= lam1 > 0, and put t = 2M >= 2|m_k|.
  - If lam1 <= t: lam0 lam1 = lam1^2 + lam1(lam0 - lam1) <= t^2 + t(lam0 - lam1).
  - If lam1 > t: by the sign lemma, sign eps_k = sign Lambda_k. Since Lambda_0 + Lambda_1 = sum L = 0, we get eps0 eps1 < 0, so sigma^2 = (|eps0| + |eps1|)^2/4 >= |eps0||eps1| >= sqrt((lam0^2 - t^2)(lam1^2 - t^2)). Because (lam0 - t)(lam1 + t) >= (lam0 + t)(lam1 - t), this is >= (lam0 + t)(lam1 - t) = lam0 lam1 - t^2 - t(lam0 - lam1).

  Both cases give the claim. This proof does not use the T5 certificate.

### 1.4 chain: Phi_70 + 2|x2|(|x0| + |x12|) >= 0

- **Own reduction (EXACT, E2).** In the old-H coordinates, Phi_70 = 2K' - (D1 - D3)(D0 - D2), x12 = -(D0 + D3) and x2 = D0 + D2. Equivalently, chain = R - x2^2 + 2|x2||x0|, where R >= 0 is the old-H left side.
- **What I could not finish (GAP).** I found a stronger pointwise statement that implies the chain, but I have not proved it:

    F = Phi_70 - 2 O33 x0 x2 - (A11 + B11) x2 x12 >= 0, with O33 = -p0q1 - p1q0 + p2q3 + p3q2, A11 = P0 + P1 - P2 - P3 and B11 = Q0 + Q1 - Q2 - Q3.

  Since |O33| <= 1 and |A11 + B11| <= 2, F >= 0 gives the chain at once. F is genuinely of degree 6 on S^3 x S^3, so it does not reduce to a quartic. It is NUMERIC only: min -4e-16 over 1000 multistarts. Searches for a quadratic multiplier in p alone or in q alone failed (LP over samples).
- **Supplement certificate (EXACT re-check, E4).** In the supplement's model x -> p x q, the following hold modulo N:
  - 2Phi_k = |p|^2 |u_k|^2 for k in {0, 2, 12};
  - Phi_2 = x2^2, Phi_12 = x12^2 and Phi_0 = |p|^2|q|^2 x0^2;
  - (PB) 2Phi_70 = |p|^2(|u_70|^2 - 2<u0,u2> - 2<u2,u12>).

  Here u_k in R^8 is built from T^{-1} times the cubic features under the four sign maps.
- **Cauchy-Schwarz step (PROOF-READ).** On unit p, q we have |u_k| = sqrt2 |x_k|. So 2Phi_70 + 4|x2|(|x0| + |x12|) >= |u_70|^2 >= 0.

## 2. Comparison with Sections 3.1 to 3.4 of `so4-atoms.md`

- **3.1 SOC.** Correct. Same identity as mine, written as Phi_22 = 2(a0A2 + ...) - 4p0p1p2p3 - 4q0q1q2q3. I expanded mine to it, and they agree. Their increment lemma plays the role of my AM-GM step plus the triangle inequality, and step 5 is ordered correctly. No missing absolute value, and no factor error. Zero boundaries are handled because there is no division.
- **3.2 old-H.** Correct. They use Hadamard coordinates and I use Pi; the structure is the same: a two-case split with min(v0v2, v1v3) against my min(D0D1, D2D3). Their regrouping in step 2 checks out term by term. The coefficient step c = d + c' with t c' >= d^2 is right, including c = d, which forces d = 0.
- **3.3 chain.** The logic is correct, and I re-verified the identity exactly with independent code. The proof rests entirely on certificate data that no PDF derives (their Deviations note says the same). The coefficient step (AM-GM per leaf) is right.
- **3.4 four.** A replay of the T5 certificate. "Z PSD iff ab(c-d) >= d^2(a+b+c+d)" is used only in the direction "condition implies PSD", and that direction is right. The pointwise form then follows by their Cauchy-Schwarz minimisation. My proof in 1.3 is an independent route to the same pointwise form.
- **Errors found (Sections 4 and 5 only, not in the proofs).**
  - (a) "the chain stays nonnegative with factor 0.5" is false. Phi_70 + |x2|(|x0| + |x12|) takes the value -3.2703e-4 at an exact rational point of SO(4), an EXACT counterexample (section 5).
  - (b) "Chain reaches 0 only at degenerate points" is false. The chain is 0 at points with |x2| about 0.3 and x0 = x12 = Phi_70 = 0 (NUMERIC, section 5).

  Neither affects any atom.
- **Cited chain tuples.** 5 of the 12 are printed with the two leaves in the opposite order to every orbit image: (2;39,67;70), (5;23,6;54), (2;5,6;46), (59;53,55;54) and (59;53,67;54). The chain form is symmetric in its leaves, so this is harmless, but "printed roles" holds only for the unordered leaf pair.

## 3. Orbits (EXACT, section `orbits`)

- **Covariance.** For a generic symbolic 4x4 O and every g in G (1152 elements): w_{g.sigma}(g.O) = eps(g) w_sigma(O), where eps(g) = det(g.O)/det(O) = ±1 (576 elements of each sign). Since F_VALUES[g.j] = g.F_VALUES[j], this gives Phi_{g.j}(g.O) = eps(g) Phi_j(O), so Phi^{O4} = Phi/det is invariant. G maps minors to minors.
- **Minors.** minor_{g.j}(g.O) equals ± minor_j(O) or ± its complementary minor; transposition and row/column relabelling can land on the complement. By Jacobi on SO(4) (EXACT in the quaternion model), |complementary minor| = |minor|, so |x| is preserved on O(4).
- **det -1.** For J = diag(1,1,1,-1): Phi_c(OJ) = -Phi_c(O) for all 82 c, and |minors| are unchanged. So for det O = -1, Phi^{O4}_c(O) = Phi_c(OJ) with OJ in SO(4). Each atom therefore holds on O(4) in the Phi^{O4} form. Since g permutes O(4), the atom for the tuple g.tau holds on O(4), and on SO(4) in particular.
- **Orbit sizes.**

| family | representative | ordered orbit | role symmetry | distinct atoms |
|---|---|---|---|---|
| SOC | (1, 9; 22) | 288 | pair swap in orbit | 144 |
| old-H | (12, 2; 70) | 288 | reversed order never in orbit | 288 |
| chain | (2; 0, 12; 70) | 576 | leaf swap never in orbit | 576 |
| four | (0, 6; 10; 54) | 576 | a <-> b swap in orbit | 288 |

  (3,15;26) lies in the old-H orbit. All 58 cited tuples lie in their orbits with the printed roles (chain up to leaf order, see section 2). For no cited old-H tuple is the reversed order in the orbit.
- **PDF spot-check (16 tuples, against `txt/proofs`).**
  - SOC (1,11;17), T3 eq. (1).
  - T42 p.1: old-H (3,5;54) and (23,5;54); chains with centre f5 and endpoints (f0,f3) and (f23,f6); SOC (5,2;54).
  - bad70 p.1 coefficient roles: (PB) gives chain (2;0,12;70), (CQ) gives (2;39,67;70), and (H) gives old-H with main f2 and aux f12 or f39.
  - T5: (0,6;10;54), with f0 = (A-B)^2, f6 = (A+B)^2, f10 = D^2 and f54 = H.
  - T34 p.2: old-H (9,5;68), chain (5;7,9;68), and four (7,0;1;68) and (7,67;53;68).
  - T38: four (0,7,12,40).

  All agree with the inventory's roles.

## 4. Coefficient forms (PROOF-READ)

For each family, "safe" means that the combination is >= 0 on SO(4). Write X = |x_.|, and use Phi_j = x_j^2 for minors.

- **SOC: ab >= d^2.** a x_i^2 + b x_j^2 + d Phi_k >= aX^2 + bY^2 - 2dXY >= 0, since the form is PSD when a, b >= 0 and ab >= d^2. If d = 0 it is trivial. If a = 0 or b = 0 then d = 0.
- **old-H: t(c-d) >= d^2.** t X_t^2 + c X_c^2 + d Phi_k >= t X_t^2 + (c-d) X_c^2 - 2d X_t X_c >= 0. For d > 0 the condition forces t > 0 and c > d.
- **chain: a >= d^2/p + d^2/b.** p X_l^2 + b X_m^2 + a X_c^2 + d Phi_k >= (a - d^2/p - d^2/b) X_c^2 >= 0, by 2d X_c X_l <= p X_l^2 + (d^2/p) X_c^2 and the same for the other leaf. Zero p or b with d > 0 reads as +inf, which is excluded. d = 0 is trivial.
- **four: ab(c-d) >= d^2(a+b+c+d).**
  - Put a' = a + d and b' = b + d. Then a X^2 + b Y^2 + c Z^2 + d Phi_d >= a' X^2 + b' Y^2 + c Z^2 - d(X+Y+Z)^2.
  - By Cauchy-Schwarz, (X+Y+Z)^2 <= (1/a' + 1/b' + 1/c)(a'X^2 + b'Y^2 + cZ^2). So the form is >= 0 when d(1/a' + 1/b' + 1/c) <= 1.
  - Multiplying by a'b'c, that condition equals a'b'c - d(b'c + a'c + a'b') = ab(c-d) - d^2(a+b+c+d) >= 0 (expanded by hand and checked).
  - For d > 0 the condition forces c > d and a, b > 0, so there is no division by 0. d = 0 is trivial.

## 5. Attack (NUMERIC; wrong variants EXACT at rational points)

- **Per family.** 2000 Nelder-Mead multistarts on the representative, plus 30 random other orbit members with 70 starts each. The minimiser was re-evaluated with mpmath at 40 digits.

| family | representative min (mpmath) | worst of 30 members (mpmath) |
|---|---|---|
| SOC | 2.9e-28 | 1.7e-27 |
| old-H | 6.3e-29 | 2.5e-28 |
| chain | 9.7e-26 | 1.3e-28 |
| four | 7.1e-26 | 1.6e-29 |

- **Wrong variants.** The minimiser is rounded to an exact rational point of SO(4) (rational unit quaternions) and the value is computed exactly with Fractions.

| family | variant | float min | EXACT value at rational point |
|---|---|---|---|
| SOC | pair (1,9) with wrong facet k = 28 | -4.478e-2 | -4.478e-2 |
| SOC | factor 1.8 | -2.153e-5 | -2.153e-5 |
| SOC | 2 x1 x9 without absolute value | -0.5 | -0.5 |
| old-H | aux <-> main, (2,12;70) | -4.481e-2 | -4.481e-2 |
| old-H | factor 1.8 | -8.861e-6 | -8.861e-6 |
| old-H | 2 x12 x2 without absolute value | -0.277 | -0.277 |
| chain | centre <-> leaf, (0;2,12;70) | -5.156e-2 | -5.156e-2 |
| chain | factor 1.8 | 1.4e-30 (no failure) | - |
| chain | x2 x0 without absolute values | -0.5 | -0.5 |
| chain | half compensation (factor 1) | -3.270e-4 | -3.270e-4 |
| four | b <-> c, (0,10;6;54) | -3.646e-2 | -3.646e-2 |
| four | factor 1.8 | -2.221e-5 | -2.221e-5 |
| four | x0 x6 without absolute value | -0.333 | -0.333 |

Chain sharpness (NUMERIC):
- By bisection, Phi_70 + lam|x2|(|x0| + |x12|) has negative minima for lam <= 1.550. None were found for lam >= 1.575. So the critical compensation factor is about 0.78, not "between 0.3 and 0.5" as in `so4-atoms.md`.
- Equality of the true chain (lam = 2) holds at a point with x2 = 0.3175, where x0, x12 and Phi_70 are all below 2e-21 (mpmath agrees). This is a nondegenerate equality point.
- The multiplier form F >= 0 of 1.4: min -4.1e-16 over 1000 multistarts.

## 6. Verdicts

- **SOC: CERTIFIED.** Own proof, exact identities. Agrees with 3.1.
- **old-H: CERTIFIED.** Own proof, exact identities. Agrees with 3.2.
- **four: CERTIFIED.** Own proof, independent of the T5 certificate, with exact identities. Agrees with 3.4.
- **chain: CERTIFIED**, on the supplement's Gram identity, which I re-checked EXACTLY with my own code. It is not derived independently: my own reduction to F >= 0 is NUMERIC only (GAP).
- **Orbit transfer: CERTIFIED** (EXACT covariance, including det -1 and complementary minors). Orbit sizes 144 / 288 / 576 / 288 are confirmed.
- **Errors.** Two sharpness remarks in `so4-atoms.md` Sections 4 and 5 are false (chain at factor 0.5; degenerate-only equality). There is also a leaf-order nuance in 5 cited chain tuples. No atom fails.
