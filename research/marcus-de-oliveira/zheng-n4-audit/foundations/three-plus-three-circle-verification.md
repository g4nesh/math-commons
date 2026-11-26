Evidence level: EXACT for the quaternion identities behind the master atom, the twelve atom identities, the dual vertex lists (own double description), the replay of all 158 chamber certificates against my own F_sigma, and membership tests on 3,900 random rational configurations (148,241 profile extreme rays); PROOF-READ for the reduction steps of the document, the closure argument, and my own proof of the inequality (E1); CITED for "every nonzero x in X is alpha times a matching vector of eight points" (paper Sect. 3) and the classical surjection S^3 x S^3 -> SO(4); NUMERIC for the attack. Verdict: CERTIFIED.

# Verification of the three-plus-three generalized-circle theorem (sixcircle-verify)

Object checked: `foundations/three-plus-three-circle.md` (statement, proof, coverage table) and its claims. My code is `foundations/six_circle_verify.py` (written from the definitions; it imports only `lib/moc4lib.py` and reads three supplement JSON files with `json.load`). Output: `foundations/six_circle_verify_output.txt`. I read `six_circle_replay.py` only after my code was written, to compare.

## Deviations and limits

- (E1) and the master atom: I did not re-read the sibling gap-check note's argument before deriving my own. My proof (Section 3) uses different coordinates and a different final lemma. It proves (E1) and, with it, Phi_m >= 0.
- Chamber uniformity rests on the supplement's 158 stored certificates (polynomials with nonnegative coefficients). I did not regenerate them. I replayed each one exactly, as data, against my own F_sigma, profile columns and own vertex lists. My independent exact sampling (Section 4) is evidence on finitely many points, not a second uniform proof.
- The atom list itself (12 vectors) is taken from the supplement JSON as data. Its validity is proved here, not assumed.
- The numeric attack minimizes over SO(4) by multistart BFGS. A positive result is evidence only.
- Not checked: the supplement's hidden proof note (not public), the 27,210-simplex cycle (not needed), and the separate one-spectrum theorems (Appendix B/C, shared-zero) that the citing documents also use.
- Compute: all stages ran through the job limiter, about 17 minutes in total (atoms 5 s, dd 1 s, chambers 22 s, random 241 s, numeric 382 + 338 s). The numeric stage ran as two seeded batches of 10,000 because of the 600 s per-job limit.

## 1. My own proof plan (written before reading the document's proof)

1. Write z_i = a_i and z_{4+j} = -b_j, so [i,4+j] = z_i - z_{4+j}. Every matching uses each of the eight vectors once, so a Moebius map T (acting linearly on lifts) multiplies the matching vector by det(T)^4 times the rescaling factors. The element x of V_C is unchanged once alpha absorbs this scalar. So c(x), P c(x) and Phi_{c(x)} do not change. Choose T to send the circle to R ∪ {∞}.
2. With six real points, x is bilinear in the lifts of the two free points. So Re x lies in the real span of four real vectors X_ab, and c(x) lies in L = span(X_11, X_12, X_21, X_22, kappa), of dimension at most 5.
3. Prove the stronger claim that every c in L with Pc >= 0 has Phi_c >= 0 (or lies in a cone of valid vectors). This is a linear statement for each six-point configuration, depends only on the configuration modulo real Moebius maps, and splits by the cyclic colour pattern: AAACCC, AACACC (= AACCAC up to reflection), ACACAC.
4. Handle degenerations by continuity, since the conclusion is a closed condition.

The document follows exactly this route (Steps 1, 3, 5, 6). For step 4 it uses the closedness of a cone condition rather than a perturbation of c. Both work (Section 2 (vii)).

## 2. Step-by-step check of the document

(i) Moebius step (Step 1). PROOF-READ, correct. [Tu,Tv] = det T [u,v], and every matching monomial has degree one in each of the eight lifts. So m(T v) = det(T)^4 m(v), and rescaling lift k by lambda_k multiplies m by prod lambda_k. x = alpha m(points) is then alpha' m(T points) with alpha' = alpha / (det(T)^4 prod lambda_k). Since x itself is unchanged, the matching vector, the gap vector for any pi*, the contact set and Phi are all unchanged. A point at infinity is the lift (1,0), so the brackets are still defined (for example [∞, c] = 1 in the lib convention), with no exceptional case. Minor slip: Step 3 writes "Re(alpha x)" where x already contains alpha. Immaterial.

(ii) Profile space L (Step 3). PROOF-READ and EXACT, correct. A_4 and C_4 each occur in exactly one bracket of every matching, so x is bilinear in their lifts. For real Γ the four coefficient vectors X_ab are real. Hence c(x) = Re x - Re F_{pi*} kappa lies in L(Γ), and dim L <= 5. Own exact checks:
- The document's formula F_sigma = -C_sigma(2)(1 - C_sigma(3))(zeta - C_sigma(4)) for A = (∞, 0, 1, zeta) agrees with my bracket computation.
- The four coefficient columns are P-images of real vectors identically in s, t, u.
- In 3,900 random configurations (Gaussian-rational circle frames, from the span of Re and Im of the four complex coefficient vectors plus kappa) the exact rank of L was 5 every time. Rank at most 5 is the concyclicity consequence. Rank equal to 5 shows no hidden further reduction.
- The claim that L is invariant under real Moebius maps is correct (the lifts of the free points transform linearly and range over all of C^2).

(iii) Validity of the atoms. Re-derived from scratch: EXACT identities plus my own PROOF-READ inequality, Section 3. The master atom is m = f70 + 2 f2 + f12, with f2 = n_{13,34}, f12 = n_{14,24}, f70 = h_{4231,+} (checked against the lib's Appendix D table). Phi_m >= 0 and (E1) hold on all of SO(4).

(iv) Symmetry images. EXACT for all twelve stored atoms, plus a PROOF-READ general argument.
- For each stored atom c_k I found a row permutation r, a column permutation c and a transpose flag. With S1, S2 signed permutation matrices, det S1 det S2 = +1 (one column sign flip when sgn(r) sgn(c) = -1), the identity Phi_{c_k}(O) = Phi_m(S1 O^(T) S2) holds as a polynomial identity in the 16 symbolic entries of O. Hence Phi_{c_k} >= 0 on SO(4).
- The general covariance argument of Step 0 is correct. For sigma -> rho sigma rho', the Leibniz sum picks up sgn(rho) sgn(rho'), and flipping one column restores both the sign and det = 1, because every Leibniz product uses exactly one entry of that column. Inversion corresponds to O -> O^T.
- Relabelling the A's (or C's) is the G-action sigma -> sigma rho^{-1} (or rho sigma), so g.K consists of minors (G permutes them; lib, EXACT) and valid atom images.

(v) Cone membership by dual vertices (Step 4). EXACT, correct.
- The sum of the 34 minor vectors is 7 kappa and the minors have rank 14. So kappa is interior to K_0 ⊂ K, every nonzero w in K* has w(kappa) > 0, and R, R' are polytopes whose vertices are the normalized extreme rays of K*.
- By the bipolar theorem (K is finitely generated, hence closed), c ∈ K iff w(c) >= 0 at every vertex. The test is an equivalence.
- My own exact double description (combinatorial adjacency test, integer arithmetic) gives exactly 72 vertices for R and exactly 86 for R'. These equal the stored lists as sets, and the relaxation file's 72 facet normals are the same set. Every vertex of R' is feasible with an active set of rank 13.
- The atoms cut off exactly two vertices of R (indices 45 and 65 of the stored 72-list).
- NUMERIC cross-check of my DD code: 500 random LPs over each polytope reach the minimum over my vertex list.

(vi) Chamber decomposition and certificates (Step 5). PROOF-READ and EXACT, correct.
- Exhaustion: for six distinct points with the A's sent to (∞, 0, 1), the C's lie in the three arcs, and the C counts per arc form (3,0,0), (2,1,0) or (1,1,1). The Moebius maps permuting {0, 1, ∞} together with A-relabelling act on the arcs as the full S_3 (each arc is named by the opposite A). This group is simply transitive on ordered pairs of distinct arcs, so (2,1,0) can always be put as two C's in (0,1) and one in (1,∞). Sorting the C's inside an arc is a column relabelling. Orientation-reversing Moebius maps are harmless, since arcs map to arcs.
- The three parametrizations by x, y, z > 0 are bijections onto the open chambers (checked by hand).
- Certificates: I replayed all 48 + 48 + 62 = 158 records against my own F_sigma and my own value-space form of each vertex functional. Every coefficient of every W_sigma and of R is a nonnegative rational, and R is a nonzero polynomial. The kappa identity sum_sigma W_sigma = R and the four column identities sum W_sigma v_sigma(X_k) = R w(X_k) hold exactly as rational functions of x, y, z after substituting the chamber charts.
- The record targets are exactly the 48 / 48 / 62 non-ordinary vertices of my own DD lists. A negative control (a record tested against a wrong vertex) fails in each chamber.
- Signs on the whole open chamber: R and all W_sigma have nonnegative coefficients and R ≠ 0, so R > 0 and W_sigma >= 0 on the open orthant. The chart denominators 1+x+y+z, 1+x+y and 1+y are positive there. So w(c) = sum (W_sigma/R) (Pc)_sigma >= 0 for every c in L with Pc >= 0, uniformly on each chamber.
- In chambers I and II the minors alone suffice (vertex list R). In chamber III the atoms are needed (see Section 4: 1,065 of my exact extreme rays lie outside the minor cone).

(vii) Closure argument (Step 6). PROOF-READ, correct.
- Work with lifts in (S^1)^6 (unit vectors of R^2) after the circle has been mapped to RP^1. Every degenerate Γ (coincidences, collisions A_i = C_j, points at ∞) is the limit of generic configurations in the same space, and is therefore still six points on one circle. So every Γ lies in the closure of one of the finitely many classes Ω_g.
- For fixed g and a vertex w of the dual of g.K, consider the condition that there is lambda in the simplex Δ_23 with (sum lambda_sigma v_sigma - w)(X_ab(Γ)) = 0 for all a, b. This condition is closed in Γ: the X_ab are continuous, Δ_23 is compact, and the condition is invariant under rescaling lifts. It is equivalent to w|_L ∈ cone(v_sigma|_L), because kappa ∈ L forces the weights to sum to w(kappa) = 1. It implies w(c) >= 0 for every c ∈ L(Γ) with Pc >= 0, which is the only direction needed.
- Since c(x) ∈ L(Γ) also holds for degenerate Γ, every configuration in the hypothesis is covered.
- My exact tests confirm this on degenerate configurations, including the "alt+AC" family (C_1 = A_1 adjacent on the circle). There the random tie-break sends the configuration to chamber II or chamber III, so the same gap vector must lie in K_0 and in g.K. All passed.
- An alternative (my plan, step 4) also works: approximate c by c_n + eps_n kappa in L(Γ_n).

(viii) From cone membership to Phi >= 0. Correct. Phi_c is linear in c, minors give det O[I,J]^2 >= 0 ((4.2), lib), and the atom images are valid by (iii) and (iv). So c(x) ∈ g.K implies Phi_{c(x)} >= 0 on SO(4).

No error found. The one notational slip (Step 3, "Re(alpha x)") needs no repair.

## 3. My proof of (E1) and of Phi_m >= 0

Coordinates: O(p,q) is the matrix of x ↦ p x conj(q) on H = R^4 (basis 1, i, j, k). Its entries are bilinear in p and q. EXACT: O^T O = |p|^2 |q|^2 I and det O = |p|^4 |q|^4. By the classical fact (CITED) that every element of SO(4) has this form with |p| = |q| = 1, it suffices to prove the inequalities on S^3 x S^3.

Change coordinates by the orthogonal Hadamard matrix H = (1/2)[[1,1,1,1],[1,-1,-1,1],[1,1,-1,-1],[1,-1,1,-1]], with P = Hp and Q = Hq. I found H as the common eigenbasis of the two commuting quadratic forms 2(p0p3 + p1p2) and 2(p0p1 + p2p3), which are the two rotation-matrix entries appearing in x2 and x3. Put v_i = P_i^2|Q|^2 - Q_i^2|P|^2, delta = P0Q2 + P1Q3 and Delta = P2Q0 + P3Q1. Own sympy code checks these as polynomial identities (bihomogeneous of degree (4,4)):

- Phi_{f2} = x2^2 and Phi_{f12} = x3^2, where x2 = det O[13,34] and x3 = det O[14,24];
- x2 = v0 + v1, x3 = v1 + v3, and v0 + v1 + v2 + v3 = 0 identically;
- **Phi(f2+f70) = 2(delta^2 + Delta^2)|P|^2|Q|^2 - 2(v0 v3 + v1 v2)**;
- x2 x3 = v0 v3 - v1 v2;
- Phi_m = Phi(f2+f70) + x2^2 + x3^2 = 2(delta^2 + Delta^2) + (v0 - v3)^2 - 4 v1 v2 (on unit p, q).

This is the document's identity up to renaming v2 ↔ v3 (a reordering of the Hadamard rows).

Proof of (E1), Phi(f2+f70) + 2|x2 x3| >= 0 (PROOF-READ, mine). Take |P| = |Q| = 1. Using x2 x3 = v0v3 - v1v2, a two-case computation gives

  -(v0v3 + v1v2) + |x2 x3| = -2 min(v1v2, v0v3),

so (E1) is equivalent to delta^2 + Delta^2 >= 2 min(v1v2, v0v3). This is trivial unless v1v2 > 0 and v0v3 > 0. Since sum v_i = 0, the pairs then have opposite signs. The swap P ↔ Q maps v to -v and exchanges delta and Delta, so assume v1, v2 > 0 > v0, v3. Put a_i = |P_i| and b_i = |Q_i|, and A = v1, B = v2, C = -v0, D = -v3 (all positive), with A + B = C + D.

Then a1^2 = b1^2 + A and b3^2 = a3^2 + D. By Cauchy-Schwarz, sqrt((r+u)(s+w)) >= sqrt(rs) + sqrt(uw), so a1 b3 - a3 b1 >= sqrt(AD). Likewise a2 b0 - a0 b2 >= sqrt(BC). With the triangle inequality,

  |delta| + |Delta| >= (a1 b3 - a0 b2) + (a2 b0 - a3 b1) >= sqrt(AD) + sqrt(BC),

so delta^2 + Delta^2 >= (|delta| + |Delta|)^2 / 2 >= (sqrt(AD) + sqrt(BC))^2 / 2.

Normalize A + B = C + D = 1 and write A, B = (1 ± sin φ)/2 and C, D = (1 ± sin ψ)/2 with φ, ψ ∈ [-π/2, π/2]. Then (sqrt(AD) + sqrt(BC))^2 = (1 + cos(φ+ψ))/2 = cos^2((φ+ψ)/2), while 4AB = cos^2 φ and 4CD = cos^2 ψ. Since |φ+ψ|/2 <= max(|φ|, |ψ|) <= π/2 and cos^2 is decreasing in |t| on [0, π/2],

  cos^2((φ+ψ)/2) >= min(cos^2 φ, cos^2 ψ) = 4 min(AB, CD),

which is exactly delta^2 + Delta^2 >= 2 min(v1v2, v0v3). ∎

Consequences: Phi_m = Phi(f2+f70) + x2^2 + x3^2 >= (|x2| - |x3|)^2 >= 0, and every G-image of m is valid by (iv). Equality occurs (for example P = (0, 1, 1, 0)/√2 and Q = (1, 0, 0, 1)/√2 give Phi_m = 0), so the atom is tight. NUMERIC sanity: on 200,000 random (p,q), min(Phi(f2+f70) + 2|x2x3|) = 6.7e-4 > 0; the trigonometric lemma holds on 10^6 random points.

## 4. Own computations (`six_circle_verify.py`)

| item | level | result |
|---|---|---|
| quaternion identities above; 12 atoms in im P, each Phi_m(S1 O^(T) S2) identically | EXACT | pass |
| own DD: R has 72 vertices, R' has 86; equal to the stored lists; rank-13 active sets | EXACT | pass |
| 500 LPs per polytope versus own vertex list | NUMERIC | pass |
| own F_sigma = document formula; columns in im P identically | EXACT | pass |
| 158 certificates replayed (nonneg. coefficients, R ≠ 0, kappa + 4 column identities), targets = all non-ordinary vertices; negative controls rejected | EXACT | pass |
| 75 rational chamber points (25 per chamber): all 3,727 profile extreme rays in K_0 (I, II) or K (III) | EXACT | pass; 11 chamber-III rays need atoms |
| 3,900 random rational configurations (13 kinds × 300), 3 Gaussian-rational directions each, reference = argmin, every prescribed relabelling (all 6 in chamber III): gap vector in g.K | EXACT | 0 failures |
| strong form on all 3,900: exact DD of the profile cone, 148,241 extreme rays, each in g.K | EXACT | 0 failures; 1,065 chamber-III rays need atoms |
| attack: 20,000 random six-on-a-circle configurations, multistart BFGS on SO(4), worst cases in mpmath | NUMERIC | no negative value; worst +3.8e-20, matching in mpmath |

About the random families:
- The 13 kinds are: generic; two A's equal; two C's equal; an A equal to a C; a point at ∞; a free point equal to one of the six; both free points real; collision plus ∞; two coincidences; free points on the circle (eight concyclic points); forced alternating pattern; alternating with an adjacent A = C; alternating with ∞ and a free point equal to one of the six.
- About two thirds use a random Gaussian-integer Moebius frame, so the circle is a genuine circle and not R.
- The chamber is read from the cyclic order of the preimages on RP^1, with ties broken by a random perturbation, as the closure argument prescribes. Chambers by cyclic order: I 1,006, II 1,885, III 1,009.
- As the document observed, realizable gap vectors almost never need atoms: 0 of the chamber-III gap vectors did. The extreme-ray test is the real stress test.

NUMERIC attack result. Setup:
- 20,000 configurations in two seeded batches of 10,000.
- Six points are Cauchy-random reals mapped onto a random circle by a random complex Moebius map. Every 5th configuration has a near-coincidence at 1e-3, every 7th has a free point equal to one of the six, and every 11th has both free points real.
- The free slots (I, J) are random, the direction is a random unit complex number, and the reference is the argmin, so the minimum gap is 0. Gaps are normalized to max 1.
- Phi was screened on 3,000 fixed random SO(4) points. BFGS (quaternion chart) then ran from the 2 best of them, plus 2 random starts whenever the screened minimum was below 0.02. In total there were 47,652 local minimizations, taking 720 s.

Results:
- The worst minimum found over SO(4) is +3.8e-20 (batch 1) and +5.1e-20 (batch 0). Since Phi = 0 at the signed permutation matrix of the reference, the true minimum is 0, and no negative value was found.
- The 10 worst points of each batch, re-evaluated in mpmath at 50 digits, give the same values (worst +3.83e-20).

## 5. Coverage of the citing documents (PROOF-READ)

I read the invoking passage in four documents.

| citing document (supplement) | chart | six real points when v ∈ R | theorem covers the use? |
|---|---|---|---|
| root_moc4_T4 (paper T_4), p. 2 | A = (z,0,1,x), C = (w,v,∞,tv); x, t shown real first | A1..A3 = 0,1,x and C1..C3 = v,∞,tv | yes |
| final_referee_moc4_T12 (paper T_13), p. 1 | A = (z,1,0,x), C = (w,v,∞,tv); t then x shown real | 1,0,x and v,∞,tv | yes |
| root_moc4_T17 (paper T_18), pp. 1-2 | A = (x,1,1+k(v-1),A3), C = (∞,0,v,w); x, k real | x,1,1+k(v-1) and ∞,0,v | yes |
| root_moc4_T7 (paper T_8), pp. 1-2 | A = (∞,v,xv,z), C = (r,1,0,w); x, r real | ∞,v,xv and r,1,0 | yes |

In each use:
- the eight labels are distinct;
- gaps are >= 0 with the stated contacts (extra contacts allowed);
- the reference permutation is pi_0 or pi_1 and the direction is a complex kappa;
- the conclusion needed is Phi >= 0 on SO(4), extended to U(4) by the paper's reduction (Theorem 1.2 / Corollary 2.2, outside this theorem).

The theorem as stated (three A's and three C's on one generalized circle, everything else arbitrary) covers exactly these uses. None of them needs strict contact sets or the stronger cone form.

Separate dependence, not covered by this theorem: in root_moc4_T4 (p. 2), the sentence about repeated values in one spectrum ("at most three distinct values ... covered by the generalized-circle theorem") concerns one spectrum on a circle with the other spectrum arbitrary. That is a Theorem C.1 / Appendix B type statement, not the 3+3 theorem. T_12, T_17 and T_7 likewise send repeated spectra and shared values to separate theorems. This agrees with the document's own caveat.

## 6. Verdict

**CERTIFIED.** Every step of `three-plus-three-circle.md` was checked and my own computations agree. No repair is needed (one harmless notational slip in Step 3). The theorem holds in the stated strong form for all configurations: coincidences, points at infinity, every (I, J), every direction and every reference. It rests on:
1. paper Sect. 3: every nonzero x ∈ X is alpha times a matching vector of eight points of P^1 (CITED);
2. (4.2) for the minors (lib, PROOF-READ there);
3. the classical surjection S^3 x S^3 -> SO(4) (CITED), my EXACT identities, and my proof of (E1);
4. the supplement's 158 stored chamber certificates, used as data and replayed EXACT by my own code, with vertex completeness from my own DD.
