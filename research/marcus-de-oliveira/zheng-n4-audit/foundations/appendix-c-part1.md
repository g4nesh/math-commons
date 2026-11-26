Evidence level: mixed. EXACT for the coordinate identities (C.3), (C.4), (C.5), (C.6), the Moebius rule, X in Sing(C), and for rank, inertia and (C.8) at explicit exact points. PROOF-READ for Lemmas C.2, C.3 and C.7, for Corollary C.5 (with my own proof that avoids HMSV) and for the linear algebra of Proposition C.6. CITED (HMSV, an external source I did not re-read) for the closed embedding of the Gale quotient and for part of Proposition C.4. NUMERIC for the generic statements, tested on random samples. Each step in the tables below carries its own level.

# Appendix C audit, part 1: skeleton of Theorem C.1; Lemmas C.2, C.3; Sections C.1-C.3

- **Auditor:** `appC1`.
- **Paper:** D. Zheng, "The Marcus-de Oliveira conjecture for n = 4", 24 Oct 2026, Appendix C.
  - My part is pp. 42-50: Theorem C.1 statement, Lemmas C.2 and C.3, Sections C.1-C.3.
  - Pp. 50-58 (Sections C.4-C.7 and the final proof) were read for the architecture. Part 2 (`appendix-c-part2.md`) audits them.
- **Supplement notes consulted (text only):**
  - the "one circle" note (12 pp.), which is the same argument with other names: r(U) = rho(U), e = kappa, F = d, Lemmas 4-8 = Lemmas C.7-C.11;
  - the "principal minors" notes (two versions), which give a *different* proof of the same theorem by a principal-minor route (their Corollary 11);
  - the 2-page "global real Hessian inertia" note, which corresponds to Proposition C.6.
- **Own code:** `appC1_checks.py`. It reads nothing from the supplement and rebuilds everything from (A.1), (A.2), [ij] = z_i - z_j and the Cauchy-Binet formula. Its captured output is `appC1_checks_output.txt`.

## Deviations and limits

- **HMSV not re-read.** HMSV is Howard-Millson-Snowden-Vakil, arXiv:1103.5255, the paper's reference [15]; I had no network. Where possible I replaced each HMSV input by my own argument or an exact computation. What remains cited is listed in Section 4.
- **Theorem 1.2 (unitary to two orthogonals)** is CITED; it is audited elsewhere. I checked only the linear step that carries it to the 14 coordinates, (C.7).
- **Proposition C.4 is about general points.** My evidence is:
  - an exact check at one rational configuration;
  - numeric checks at 40 random complex configurations;
  - a numeric dimension count.

  This supports the statement; it does not replace a proof. The written proof is a sketch: see gap G2.
- **Section F of the script crosses into part 2.** It is an end-to-end numeric test of the critical-point mechanism, which needs Sections C.3-C.6 together. I report it because it tests the interface between the two parts.
- **Numeric tests** use the fixed numpy seeds printed in the output, with the tolerances printed on each line.

## 1. Skeleton of the proof of Theorem C.1 (whole appendix)

Notation:
- z = (a1..a4, -b1..-b4) are the eight parameters, m = m(z) is in X, and theta is in C.
- y + i w = theta m with y, w real.
- f(U) = <rho(U), y> = Re(theta det(D_a + U D_b U*)).

Dependency list:

1. **(C.3)** For every V in GL4(C), det(D_a + V D_b V^{-1}) = <rho(V), m(a, -b)>. Here rho(V) is the coefficient vector, in the basis B_j, of the polynomial with coefficients det V[I,J] det V^{-1}[J,I]. Also <rho, kappa> = 1.
   - Feeds: steps 2, 5 and 12; the Moebius step 0; the identification rho(U) <-> squared minors.
2. **(C.7)** From Theorem 1.2 [main text, CITED], rho(U) = lambda rho(O-) + (1 - lambda) rho(O+). So max over U(4) of f equals max over O(4), and changing the sign of one row puts the maximiser O in SO(4).
   - The maximum is taken over U(4); SO(4) suffices because rho depends only on squared minors.
   - Theorem 1.2 is used only here.
3. **(C.5)** C(u + h kappa) = C(u) + h q(u), with q = u0 u6 - u1 u5.
4. **(C.6) is where the hypothesis enters.** If a is real, then q(Re(theta m)) = 0. This is the only use of the generalized-circle hypothesis. (The rest of the argument also lets a be real; see 0 below.)
   - The block (m0, m1; m5, m6) equals ([12][34], [14][23])^t ([56][78], [58][67]), so q(y) = 0 for all theta exactly when [12][34]/[14][23] is real.
   - That ratio is the cross-ratio of a1..a4, which is real exactly when a lies on a generalized circle. Checked EXACT for real a and NUMERIC for concyclic a.
5. **Sigma2-hat is contained in {C = 0} and contains Re x for x in X.**
   - Re x lies on the secant line through x and conj(x), and X is defined over R.
   - Sigma2-hat lies in {C = 0} because C(s x + t x') = 0 for x, x' in Sing C. This one-line argument needs only X in Sing(C), which is EXACT; the paper cites HMSV for it.
6. **Proposition C.4:** (M8')* = Sigma2. Used in three places:
   - the inclusion (M8')* in Sigma2 is used in (C.8) and in step 10 of part 2;
   - full equality is used for biduality, Con(M8')^t = Con(Sigma2), in step 11 of part 2.
7. **Lemma C.7 + Section C.3 + (C.8):** C|K_O = alpha(O) det S_O, where:
   - K_O is the 7-dimensional space of "centered critical costs" at a stable orthogonal octad [I4, O];
   - S_O is the 3x3 phase Hessian.
8. **Lemma C.9 (part 2):** alpha(O) != 0 when the octad is stable and d(O) != 0. It uses (C.9), (C.10), Plücker-divisor purity and an explicit rational witness.
9. **Lemma C.8 / C.10 (part 2):** when d(O) = 0, common permutation weights exist for all squared minors, via Lam's double-dimer identity.
10. **Lemma C.11 (part 2):** on a dense open set of (a real, b, theta):
    - y is a smooth point of Sigma2-hat;
    - Sigma2-hat meets y + C kappa only at y.

    This uses (C.2): det H restricted to {C = 0} equals gamma0 s^2 (HMSV, CITED).
11. **Lemma C.3 (Drury-Cload; the paper proves n = 4)** finishes: det = 0 implies 0 is in the hull. It uses Lemma C.2 for n = 3.
12. **Assembly** (proof of Thm C.1, pp. 56-58). Take O maximising f on U(4), then reduce it to SO(4) by step 2. There are three cases:
    - **Decomposable O:** Birkhoff block weights (Lemma C.2) settle it.
    - **d(O) = 0:** Lemma C.10 settles it.
    - **Otherwise:** set c = y - h kappa, which lies in K_O. Then:
      1. C(c) = 0, by steps 3, 4 and 5.
      2. det S_O(c) = 0, by (C.8) and Lemma C.9.
      3. ([rho(O)], [c]) lies in Con(M8'), by Lemma C.7.
      4. [c] lies in Sigma2, by Proposition C.4.
      5. h = 0, by the uniqueness in Lemma C.11.
      6. rho(O) annihilates T_y Sigma2-hat, which contains w, by biduality; so theta det(D_a + O D_b O^t) = 0.
      7. Lemma C.3 then gives max_sigma Re(theta z_sigma) >= 0 = h.

    **Degenerate cases** (repeated eigenvalues, a_i + b_j = 0, parameters outside the open set) are closed by continuity. Both support functions are continuous because U(4) is compact and the sigma-points are finitely many, and the inequality is a closed condition.
0. **From a generalized circle to the real line.** Take a Moebius f with its pole on the circle and away from the eight parameters, so that the circle maps to R.
   - Then m(f(z)) = (ad - bc)^4 / prod_j (c z_j + d) * m(z), EXACT.
   - By (C.3), for each fixed U the determinant and all 24 sigma-points are multiplied by one common nonzero complex factor. Hull membership is invariant under a common complex scaling.
   - So MOC is not Moebius-invariant in general; what is invariant is the triple (determinant at U, sigma-points), for fixed U, up to a common scalar. The new spectra are a'_i = f(a_i) and b'_j = -f(-b_j).
   - If b is the concyclic spectrum, swap A and B. For a line, use an affine f.

   Checked: correct (EXACT identity, PROOF-READ argument).

## 2. Referee table for my part

| # | Place | Statement | Justification in paper | My check | Status |
|---|---|---|---|---|---|
| 1 | Lemma C.2, p. 43 | MOC for n <= 3 | (2.1) + complementary-minor identity for unitary 3x3 + Birkhoff | identity rederived; NUMERIC 500 samples (err 5e-15; det always inside hull) | PROOF-READ, correct (known result [2]) |
| 2 | Lemma C.3 statement, p. 43 | normal 4x4 A, B, det(A+B) = 0 => 0 in conv{z_sigma} | Drury-Cload [11], proof given for n = 4 | see rows 3-6 | proved as written |
| 3 | p. 43, kernel vector | v unit with Av = Cv, C = -B; p_i, q_j squared coordinates; sum p a = sum q e, sum p\|a\|^2 = sum q\|e\|^2 | <v,Av> = <v,Cv>, \|\|Av\|\|=\|\|Cv\|\| | rederived; NUMERIC (residual < 1e-8, 1800 cases) | PROOF-READ |
| 4 | p. 43-44, extreme point | at most 5 positive weights (2 normalizations + 3 real moment equations) | basic feasible solution | rederived; NUMERIC: max support 5 | PROOF-READ |
| 5 | p. 44, single support | p_i0 = 1 => sum q\|e - a_i0\|^2 = 0 => some a_i0 + b_j = 0 | expansion | rederived | PROOF-READ |
| 6 | p. 44, 3x3 reduction | equal Gram matrices of (x, D_alpha x), (y, D_gamma y) give unitary W; D_alpha - W D_gamma W* singular; Lemma C.2; multiply by the 4th factor | standard | rederived; NUMERIC replay of the construction (W residual 8e-9; 3x3 det 1e-14; 3x3 hull margin <= -2e-4) | PROOF-READ, complete |
| 7 | p. 44, M8' embedding | Gale-invariant ring generated in degree 1 => M8/<g> embeds in P(V*) | HMSV Props 8.1, 8.2 | not re-read | CITED |
| 8 | p. 45, (C.2) | Sing C = M8; Sigma2 in C; G(Sigma2 \ M8) dense in M8'; det H\|C = gamma0 s^2 | HMSV Thm 3.4, Prop 3.13 proof, Prop 3.19 | Sigma2 in {C=0}: own one-line proof given X in Sing C (EXACT). Sing C = M8 needs HMSV Cor 4.2 (ideal generation). G(Sigma2) has dimension 9 (cone 10): NUMERIC (rank 10). det H\|C = gamma0 s^2: CITED (part 2 has consistency checks) | partly EXACT, rest CITED |
| 9 | p. 45, (C.3) | det(D_a + U D_b U*) = <rho(U), m(a,-b)>, <rho,kappa> = 1, <rho,n_IJ> = \|det U[I,J]\|^2 | minor expansion | EXACT: det(V) rt(V) satisfies all 56 linear relations of span(B) as polynomials in 16 entries (so it holds for all V in GL4(C)); NUMERIC on 200 Haar U (err 1e-14) | EXACT |
| 10 | p. 45, (C.4) | G(m(z)+m(w)) = det A0 det A1 rho(A0^{-1}A1) | "seventy pairings agree" | EXACT: all 70 pairings agree as polynomials in 16 variables, with exactly the paper's normalization (no stray constant or sign) | EXACT |
| 11 | p. 46, (C.5) | C(u + h kappa) = C(u) + h q(u) | inspection of (A.2) | EXACT | EXACT |
| 12 | p. 46, (C.6) | q(Re theta m(a,-b)) = 0, a real | rank-one block | EXACT (symbolic in a, Re b, Im b, theta); also holds for concyclic a (NUMERIC 1e-17), fails for generic a (1e-1) | EXACT |
| 13 | p. 46, (C.7) | Theorem 1.2 in rho coordinates; max over U(4) attained on O(4) | rho linear in squared minors | EXACT: the 70 functionals n_IJ have rank 14, so rho is a linear image of r(U); Theorem 1.2 CITED | EXACT modulo Thm 1.2 |
| 14 | p. 46-47, Prop C.4 | (M8')* = Sigma2 | elliptic-curve sketch + HMSV | see Section 3 | modulo citations; sketch (G2) |
| 15 | p. 47-48, Cor C.5 | rank H(c) <= 12 on Sigma2-hat; rank H(Re x) <= 12 | dominance of G on Sigma2 | own proof below (uses only (C.4) and dimension counting); NUMERIC: rank 12 on 300 secant points and on 300 Re x | proved (own argument); paper's version modulo HMSV |
| 16 | p. 48-49, Prop C.6 | n+, n- <= 6 for H(Re x), x in X; generic inertia (6,6,2); 12x12 principal minors >= 0 | isotropic-subspace argument + density | linear algebra rederived (correct); EXACT: rank H(x) = 8, rank [H(Re x); H(Im x)] = 14 and inertia (6,6,2) at 3 Gaussian-rational points; NUMERIC 600 generic + 1500 degenerate samples: n+, n- <= 6 always | proved, modulo "every nonzero x in X is lambda m(z)" (GIT, standard) |
| 17 | p. 49, Lemma C.7 | injectivity of B; conormal limit | cotangent generators; Grassmannian limit | PROOF-READ, correct. NUMERIC on M8: rank[rho0, A, B] = 13 at 30 points (20 real, 10 complex O) | proved as written (abstract lemma) |
| 18 | p. 50, stability | indecomposable O => [I4,O] stable; trivial stabilizer; diagonal self-duality unique | Hilbert-Mumford for points + orthogonality | PROOF-READ, correct (also over C with the bilinear form) | PROOF-READ |
| 19 | p. 50, tangent decomposition | M8 smooth, dimension 9 at stable fixed points; eigenspaces 6 (+1) and 3 (-1); gauge {D1 + O D2 O^t} has dimension 7 | differential of V -> V^{-T} is X -> -X^t on X = dV O^{-1} | rederived; NUMERIC: rho even in symmetric directions (exactly 0), gauge in Hessian kernel; EXACT at 2 rational O | PROOF-READ + EXACT at points |
| 20 | p. 50, dim K_O = 7; c -> S_O(c) onto Sym^2 | from Lemma C.7 | NUMERIC: rank 7 and rank 6 at 30 points; EXACT at the 2 rational O | proved |
| 21 | p. 50, (C.8) | C\|K_O = alpha(O) det S_O | det S_O = 0 => [c] in Sigma2 (Lemma C.7 + Prop C.4) => C(c) = 0; det of symmetric 3x3 irreducible; equal degrees | logic PROOF-READ (correct given Prop C.4). NUMERIC: relative residual 3e-13 at 30 points. EXACT as cubic forms in 7 variables at one rational O with all entries nonzero and at the Lemma C.9 matrix | proved modulo Prop C.4 |
| 22 | p. 51, alpha as section of a line bundle | changing lift/basis multiplies by units; pure codimension one zero set | standard | PROOF-READ (see G4) | PROOF-READ |

### Own proof of Corollary C.5 (replaces the HMSV dominance citation)

1. Every point of a dense subset of Sigma2-hat has the form m(z) + m(w), because s m(z) = m(s^{1/4} z).
2. By (C.4), G(m(z) + m(w)) = det A0 det A1 rho(V). Here rho(V) is invariant under V -> D1 V D2, so it is a function on GL4 modulo the 7-dimensional torus, a 9-dimensional space.
3. So the image of G on Sigma2-hat lies in a cone of dimension at most 10.
4. Hence at every smooth point c of Sigma2-hat, whose dimension 12 was checked numerically, dG_c = H(c) restricted to T_c has rank at most 10, so ker H(c) has dimension at least 2.
5. By closedness, rank H(c) <= 12 on all of Sigma2-hat.

An explicit kernel vector also exists: H(s x + t x')(s x - t x') = 0 for x, x' in X, since G(s x + t x') = 2 s t Q(x, x'); NUMERIC 4e-17. This argument uses only (C.4), which is EXACT, and dimension counting. It does not need G(Sigma2) to be dense in M8'.

### Proposition C.6, linear-algebra part (rederived)

- **Isotropy.** H(x) = A + iB with A = H(Re x) and B = H(Im x), because H is linear in its argument. For v, w in Lambda = ker H(x), Av = -iBv gives w*Av = w*Bv = 0, so Lambda is totally A-isotropic.
- **Trivial intersection.** Lambda ∩ ker A is contained in ker A ∩ ker B = ker H(x) ∩ ker H(conj x).
- **Inertia bound.** A totally isotropic L with L ∩ ker A = 0 satisfies dim L <= min(n+, n-), since L meets (negative space ⊕ ker A) only in 0.
- **What this needs, and how it is checked:**
  - dim Lambda >= 6. This holds at every x in X, because T_x X ⊂ ker H(x) as G vanishes on X. It does not need HMSV.
  - Transversality at one point (EXACT at 3 points: rank [A; B] = 14).
  - Real-algebraic density. This is correct as written: the conjugate pairs (t, conj t) are Zariski dense, and the complement of a proper real algebraic set is dense.
- **Upper bound.** Eigenvalue continuity gives n+, n- <= 6 everywhere.
- **Principal minors.** Cauchy interlacing gives nonnegative 12x12 principal minors.

## 3. Proposition C.4 in detail

**What the proof uses:**
- the elliptic quartic E through a general octad n in P^3;
- the curve f(Pic^2 E) in M8, of degree 4, spanning a P^3 called W_n;
- the "paired secant" lines joining f(A) and f(L ⊗ A^{-1}), which (C.4) maps under G to the point n.

**Order of the argument:**
1. c^t H(c) u = 2 G(c)^t u = 0 for u in T Sigma2. Together with surjectivity of dG, this makes general secant points tangent hyperplanes to M8'. This is the inclusion Sigma2 ⊂ (M8')*, and it is clean.
2. The conormal space at n is a P^3 containing the paired secants, hence equal to W_n.
3. Secants of an elliptic normal quartic fill P^3, so W_n ⊂ Sigma2. This gives (M8')* ⊂ Sigma2.

**Checks:**
- **NUMERIC, 40 random complex V:**
  - the affine tangent of the cone over M8' at rho(V) has dimension 10, so the conormal space L(V) is 4-dimensional (a P^3);
  - on L(V): C = 0 (relative 5e-17) and rank H = 12;
  - the 14 quadrics grad C restricted to L(V) span exactly a 2-dimensional pencil, so X ∩ P(L) is a complete intersection of two quadrics, a quartic curve;
  - for 3 quadrics of the pencil through the 8 points [I, V] and both rulings, the octads f(A) built directly from the rulings have m(f(A)) in L(V) (distance 3e-14);
  - G(m(ruling 1) + m(ruling 2)) is parallel to rho(V) (1 - |cos| = 4e-16).

  So the paper's description of the conormal space is reproduced exactly.
- **NUMERIC, dimension of the dual:** the union of the conormal spaces has affine dimension 12 (Jacobian rank 12 with a clear gap: 0.17 against 2e-9). So dim (M8')* = 11 = dim Sigma2. Since both are irreducible, the inclusion (M8')* ⊂ Sigma2 then forces equality.
- **EXACT at the rational V = [[3,-2,-1,3],[-1,-3,1,2],[2,2,1,-3],[2,-1,2,2]]:** tangent rank 10; dim L = 4; C restricted to L vanishes identically; the restricted quadrics have rank 2; rank H(c) = 12 at a point of L.

Status: the statement is well supported. The written proof is a sketch (G2), and it relies on HMSV for the closed embedding of M8' and for (C.2).

## 4. Gaps and remarks (page, equation)

None of these is an error in a statement. G1-G4 are places where the text is incomplete or leans on citations.

- **G1 (p. 45, (C.2), and p. 44 embedding): CITED inputs that remain.**
  - M8/<g> ⊂ P(V*) is a closed embedding (HMSV Props 8.1, 8.2). Lemma C.7 needs it at the fixed points.
  - Sing(C) = M8, i.e. the Jacobian ideal is the ideal of X (HMSV Cor 4.2). It is used only through dim ker H(x) = 6 at general x, which I verified exactly at 3 points; and the Prop C.6 argument needs only dim >= 6.
  - det H|C = gamma0 s^2, used by part 2's Lemma C.11.
  - Everything else in (C.2) that my part uses is proved here independently: Sigma2 ⊂ {C = 0}; rank <= 12 on Sigma2; dim Sigma2 = 11 (numeric).
- **G2 (pp. 46-47, Prop C.4): the proof is a sketch.** Three claims are asserted without detail:
  - "general points on the resulting paired secant lines form the general fiber of G";
  - "the polar map has surjective differential at a general point";
  - "the same argument ... gives Sigma2 ⊂ (M8')*".

  The proof also uses HMSV's G(Sigma2 \ M8) = M8'. My numerics reproduce every intermediate claim at random points and the decisive dimension count. Proving (M8')* ⊂ Sigma2 rigorously needs only: the conormal space at general n equals the span W_n of an elliptic quartic. That is checked exactly at one point and numerically at 40, but not proved for general n beyond the sketch. Status: modulo citations, with an expository gap; no evidence against it.
- **G3 (p. 44): "the projective image is isomorphic to the finite quotient".** Lemma C.7 assumes a closed embedding of the local quotient Y/iota near the fixed point, not only a bijective image. HMSV Prop 8.2 (as quoted) gives generation of the invariant ring in degree 1, which implies that Proj of the invariant ring is embedded. Fine, but CITED.
- **G4 (pp. 50-51): alpha as a section of a line bundle on the stable fixed locus.**
  - The gluing across charts is asserted briefly.
  - The irreducibility of the stable fixed locus, which purity needs, is not stated. Part 2 notes the same.
  - Expository only.
- **Remark (p. 43 overview vs. the proof): only first-order criticality is used.** The chain therefore predicts that at every critical point of f on SO(4) with [I, O] stable and d(O) != 0, f(O) = 0 and det(D_a + O D_b O^t) = 0.
  - Script section F tested this: 40 random instances, 480 converged critical points. Of these, 465 were decomposable, 15 were indecomposable but on the Cauchy divisor, and 0 were in the stable d != 0 branch. So nothing contradicts the chain, but the branch is not exercised.
  - Part 2 found the same: 0 such points in 1200.
  - This suggests that for generic data such critical points do not exist at all. That is consistent with the proof, which allows them only at value 0.
- **Remark: Lemma C.3's proof is complete and correct.** It is a self-contained n = 4 version of the Drury-Cload argument.
- **Remark: an alternative route.** The supplement's "principal minors" note proves a stronger statement: one probability vector for all b, Hermitian a. Its Corollary 11 implies Theorem C.1 independently of Sections C.2-C.6. It is not used by the paper; I did not audit it. Part 2 tested its claim numerically.
- **Discrepancies with the supplement (cosmetic):**
  - The principal-minors note uses [ij] = z_j - z_i. This is harmless, since all products have degree 4.
  - The one-circle note cites HMSV "Theorems 3.4, 3.15, 3.18; Prop 3.19" where the paper cites "Thm 3.4, proof of Prop 3.13, Prop 3.19". I cannot tell which is right without the source.
  - The one-circle note says "G(Z) is dense in X"; the paper writes "closure of G(Sigma2 \ M8) = M8'". These are equivalent.
  - The cubic printed in the one-circle note is identical to (A.2) (EXACT).
- **Exact (C.8) points.** The Cayley matrix from 3K = [[0,4,3,2],[-4,0,3,-3],[-3,-3,0,-3],[-2,3,3,0]] has all entries nonzero and d = -400/131769; there alpha = -10673289/928287872. The Lemma C.9 matrix gives alpha = -729/80000.
- **How to read the principal-minor line in the output.** For degenerate spectra of rank 10 or less (for example a1 = a2), the printed "min normalized 12x12 principal minor -1.0" is rounding noise divided by a maximum that is itself zero. All 12x12 minors vanish there, because the rank is below 12.
- **Exact alpha values depend on normalization.** At the Lemma C.9 example matrix (O_00 = 0, stable, d = -16/729) I get alpha = -729/80000 with my phase basis; part 2 gets 729/16 with theirs. Alpha is defined only up to a unit (the square of a basis-change determinant) and the normalization of the lift, so only "alpha != 0" is meaningful. Both computations agree on that.

## 5. Interface to part 2 (Sections C.4-C.7)

**Part 1 delivers to part 2:**
- **(C.3)** for all V (EXACT), and **(C.7)** (EXACT, modulo Theorem 1.2).
- **(C.5), (C.6)** (EXACT). The Moebius rule (EXACT) for the last paragraph of the proof.
- **Sigma2-hat ⊂ {C = 0}** and Re x ∈ Sigma2-hat (proved).
- **Proposition C.4** (M8')* = Sigma2. The inclusion ⊂ is used in (C.8) and step 10; equality is used for biduality, step 11. Status: modulo HMSV, sketch G2, strong numeric support.
- **Lemma C.7 + Section C.3** (stability, dim K_O = 7, surjectivity of S_O) and **(C.8)**, with alpha a section whose zero set is chart-independent and of pure codimension 1. Status: proved modulo Prop C.4 and G3/G4; (C.8) EXACT at 2 matrices, NUMERIC at 30.
- **Lemma C.3** (complete) and **Lemma C.2**.
- **Corollary C.5 and Proposition C.6** (rank <= 12; n+, n- <= 6; 12x12 principal minors >= 0). These are not used by part 2's assembly; they are used by the ninety case proofs.

**Part 1 assumes from part 2:** nothing. Section C.3 defines alpha; Lemma C.9 (part 2) is what makes it nonzero.

**From outside:** HMSV ((C.2), closed embedding of M8'), Theorem 1.2 (main text), projective biduality (standard).

## 6. Verdicts

| Item | Verdict |
|---|---|
| Lemma C.2 | proved as written (classical) |
| Lemma C.3 | proved as written; the n = 4 argument is complete; no counterexample in 1800 random cases or under adversarial search |
| Proposition C.4 | correct modulo HMSV citations; the written proof is a sketch (G2); every intermediate claim reproduced numerically, and the key structure exactly at one point |
| Corollary C.5 | proved (own argument from (C.4), which is EXACT); paper's version modulo HMSV dominance |
| Proposition C.6 | proved, modulo GIT surjectivity (standard); exact inertia (6,6,2) at 3 points; n+, n- <= 6 in every numeric sample including degenerate spectra |
| Lemma C.7 | proved as written; injectivity checked numerically on M8 at 30 points |
| (C.8) | proved modulo Prop C.4; EXACT at 2 matrices |

No error found in pp. 42-51.
