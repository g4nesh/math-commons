Evidence level: mixed. EXACT for the 30 finite computations in `appC_fill_checks.py` (all pass; output in `appC_fill_checks_output.txt`; check IDs A0-D3 are cited below). PROOF-READ for my own arguments written out below, which turn those finite checks into proofs. CITED for exactly three inputs: Theorem 1.2 of the main text (audited elsewhere); hypothesis (H3) of Section 2, taken from HMSV's degree-one generation theorem as the paper quotes it; and standard textbook facts (projective biduality, Krull's principal ideal theorem, Chevalley's theorem, GIT and Luna slices). The usage map in Section 7 is PROOF-READ of citation sentences only. No case mathematics is audited there.

# Appendix C: the remaining soft spots, closed or classified

**Auditor:** `appCfill`.

**Inputs:**
- Paper Appendix C, pp. 42-58.
- `appendix-c-part1.md` and `appendix-c-part2.md`, with their scripts and outputs.
- `appendix-b.md`.
- `../deps/dependency-map.*`.
- The supplement's proof texts, used only for the usage map.

## Deviations and limits

- **Definition of M8'.** M8' is the Zariski closure in P(V_C^*) of the image of V -> [rho~(V)], with rho~(V) = det(V) rho(V) built from the complementary products of [I4, V] (Section 1.1). This is the variety the paper calls the projective image of M8-bar/<g>, and it is the variety whose conormal variety Lemma C.7 produces. I do not use HMSV's assertion that this image is isomorphic to the quotient, except through (H3).
- **Sources not re-read.** HMSV, Lam, Postnikov and Drury-Cload were not re-read; I had no network. Where an item depended on one of them, I replaced it by a proof or stated precisely what is still cited.
- **Genericity arguments.** Every "general point" statement is proved the same way. First I show the condition is Zariski-open on an irreducible parameter space. Then one exact rational point shows the open set is nonempty. Both halves are written out; neither is a numeric claim.
- **(H3) is not proved for all stable fixed points.** Section 2 says what would be needed.
- **The usage map was built with help.** A helper agent read the supplement texts with grep and sed only. I spot-checked 8 of its locations: T_7 p.2, T_24 p.3, T_4 p.3, T_10 p.1-2, T_1 Prop. 10, T_45 conic note p.1, T_42 bad54 p.3, T_65 negative_Phi p.1. All 8 were confirmed.
  - The map records what each document says it uses. It does not judge whether that use is correct.
  - Two notes may not be on the final proof path, and I did not settle which: T_45 `a_Heron`, and the two T_47 Heron notes.
- **The script.** `appC_fill_checks.py` reads nothing from the supplement. It imports `../lib/moc4lib.py` for (A.1), P, kappa, C and n_{I,J}. It is deterministic and runs in about 15 s.

## 0. Summary

| Item | Class | Result |
|---|---|---|
| G2: Proposition C.4, (M8')^* = Sigma2 | (a) | Complete proof (Section 1). It also proves dim Sigma2 = 11 exactly. Checks A0-A9. |
| G3: closed embedding used in Lemma C.7 | (b) | Used only through (H3): at a stable fixed point, rank[rho0, A, B] = 13. HMSV Props. 8.1-8.2, as quoted, imply it. Not proved here for every point. |
| G4: alpha as a section; irreducibility of the stable fixed locus | (a), given (H3) | The zero locus of alpha has an intrinsic definition, so no gluing is needed. F is irreducible, via OG+(4,8). Purity follows by Krull (Section 3). |
| Part 2, p. 53: Pluecker divisors irreducible, one S8-orbit | (a) | Full proof (Section 3.3). |
| G1: (C.2) from HMSV | no longer needed for Theorem C.1 | Lemma C.11 is reproved without (C.2) (Section 4; checks B1-B6). The rest of (C.2) is unused or elementary. Four supplement case documents still use (C.2)'s set equality; for them it stays (b). |
| Part 2, p. 55: density step in Lemma C.10 | (a) | Baire argument, plus (C.13) and (C.14) proved as identities (D0-D3). |
| Part 2, p. 54: Lam's double-dimer theorem (one-colour boundary hypothesis) | (a), replaced | Explicit nonnegative weights for every totally positive Vandermonde point (C1-C4). No external theorem is needed. |

**Theorem C.1 now rests on** Theorem 1.2, (H3) and standard facts (Section 6).

**Main usage-map findings (Section 7):**
- Theorem C.1 with four distinct points of one spectrum on a line is indispensable in exactly three cases: T_7, T_11 and T_24. Only the real-spectrum case is needed.
- Twelve other cases invoke a different "three-plus-three" or "six-circle" statement. It is not Theorem C.1, and I could not locate a proof of it anywhere in the paper or the supplement.
- Four cases use the full HMSV set equality {C = det Hess C = 0} = Sigma2-hat. That is more than Corollary C.5 gives.

## 1. G2: Proposition C.4, complete proof (class (a))

### 1.1 Objects

For V in GL4(C), define rho~(V) in V_C^* by

  <rho~(V), n_{I,J}> = (-1)^{sum I + sum J} det V[I,J] det V[I^c,J^c]   for all 70 pairs (I,J), |I| = |J|.

These conditions are consistent, as a polynomial identity in the 16 entries (A0). Moreover

  <rho~(V), m(a,-b)> = det(D_a V + V D_b)

holds identically in a, b and V, 24 variables in all (A0). This is (C.3) multiplied by det V, so rho~ = det(V) rho. It also holds for every V, invertible or not, so rho~ is a polynomial map C^16 -> V_C^*, homogeneous of degree 4.

Notation:
- M-hat is the closure of rho~(GL4). It is the affine cone over M8'.
- Sigma2-hat is the closure of {x + x' : x, x' in X}.
- N(V) = {xi in V_C : <d rho~_V(delta), xi> = 0 for all delta} is the annihilator of im d rho~_V.
- By Euler's identity rho~(V) lies in im d rho~_V. Hence every xi in N(V) also annihilates rho~(V).

**Lemma 1.1 (dimension).** rank d rho~_V <= 10 for all V, with equality on a dense open set U0. Hence dim M-hat = 10 and dim M8' = 9.

*Proof.*
1. A1 gives the identity rho~(D1 V D2) = det(D1) det(D2) rho~(V). So d rho~_V maps the space {D1 V + V D2} into C rho~(V).
2. Suppose all entries of V are nonzero. Then (D1, D2) -> D1 V + V D2 has kernel {(tI, -tI)}: the condition (D1)_ii + (D2)_jj = 0 for all i, j forces this. So that space has dimension 7, and ker d rho~_V has dimension at least 6.
3. Hence rank d rho~_V <= 10 on a dense set. The condition rank <= 10 is closed, so it holds everywhere.
4. Equality holds at V0 (A4). The image is a cone, since rho~(tV) = t^4 rho~(V). Its dimension is the generic rank, which is 10.

**Lemma 1.2 (conormal criterion).** Let z be in C^8. If the matrix ((z_i - z_{4+j}) V_ij), with i, j = 1..4, has rank at most 2, then m(z) lies in N(V).

*Proof.*
1. With a = (z_1..z_4) and b = -(z_5..z_8), the matrix is L(V) = D_a V + V D_b.
2. By A0, V' -> <rho~(V'), m(z)> equals det(L(V')), with L linear.
3. The function det vanishes to second order at every matrix of rank at most n - 2: all its first partials are cofactors.
4. So the value and the differential of V' -> <rho~(V'), m(z)> vanish at V', which is the claim.

**Lemma 1.3 (Segre configurations).** Let z, w be in C^8. Let A be the 4x8 matrix with columns (1, z_k, w_k, z_k w_k)^t, split as A = [A0 | A1] with A0 invertible, and put V = A0^{-1} A1. Then rank((z_i - z_{4+j}) V_ij) <= 2, and the same holds with w in place of z.

*Proof.*
1. Scale the columns by z: A0 D_{z1..4} has rows (z_k, z_k^2, z_k w_k, z_k^2 w_k).
2. Rows 1 and 3 of A0 D_z are rows 2 and 4 of A0. So A0 D_z = S A0 + E0, where S is the fixed matrix with S_12 = S_34 = 1 and E0 vanishes outside rows 2 and 4.
3. In the same way, A1 D_{z5..8} = S A1 + E1, with the same S.
4. Then A0 (D_z V - V D_z') = (S A0 + E0) V - S A1 - E1 = E0 V - E1, which vanishes outside rows 2 and 4. So it has rank at most 2.
5. Its (i,j) entry is (z_i - z_{4+j}) V_ij, up to the invertible factor A0.
6. For w, permute rows 2 and 3 of A first. This changes neither V nor the argument.

**Proposition 1.4.** There is a nonempty Zariski-open set Omega in C^16 such that for (z,w) in Omega, with V = V(z,w), we have rank d rho~_V = 10 and N(V) is contained in Sigma2-hat.

*Proof.* Let Omega be defined by five conditions:
- (i) det A0 != 0.
- (ii) rank d rho~_V = 10.
- (iii) M_z = ((z_i - z_{4+j}) V_ij) has rank exactly 2, and so does M_w. They have rank at most 2 everywhere, by Lemma 1.3.
- (iv) Transversality. Let l and r run over bases of the left and right kernels of M_z. The 4 linear forms dz -> l^t ((dz_i - dz_{4+j}) V_ij) r are independent. The same holds for M_w.
- (v) The 14 x 10 matrix [m(z), m(w), dm_z(T_z), dm_w(T_w)] has rank 4. Here T_z and T_w are the kernels of the forms in (iv).

*Openness.*
1. Restrict to the open set where one fixed 2x2 minor of M_z, and one of M_w, is nonzero. There the kernels have bases polynomial in (z,w), by cofactors.
2. Restrict further to where a fixed 4x4 minor of the 4x8 system in (iv) is nonzero. There T_z and T_w have bases rational in (z,w).
3. On that open set every condition says that a matrix rational in (z,w) has full rank, which is an open condition.

*Nonemptiness.* All five conditions hold at the rational point z = (0,1,3,-2,5,-4,7,2), w = (2,-1,4,3,-3,6,1,-5) (A2-A8).

*Main step.* Fix (z,w) in Omega.
1. By (iii) and (iv), Inc_V = {z' : rank((z'_i - z'_{4+j}) V_ij) <= 2} is a 4-dimensional complex submanifold near z. Indeed, near a rank-2 matrix the rank <= 2 locus is a smooth submanifold of codimension 4. Its normal space is detected by the pairings l^t (.) r, and (iv) says the affine map z' -> M_{z'} is transverse to it. The same holds near w.
2. Define Phi(z', w', s, t) = s m(z') + t m(w') on Inc_V x Inc_V x C^2.
3. By Lemma 1.2, Phi takes values in N(V). Since X is a cone, the values also lie in Sigma2-hat.
4. By (v), dPhi at (z, w, 1, 1) has rank 4. By (ii), dim N(V) = 4. So Phi is a submersion onto N(V) there.
5. Its image therefore contains a nonempty open subset of N(V). Sigma2-hat ∩ N(V) is Zariski closed in the irreducible N(V) and contains an open set, so it equals N(V). ∎

**Theorem 1.5.** (M8')^* = Sigma2, and dim Sigma2 = 11.

*Proof of the inclusion (M8')^* ⊂ Sigma2.*
1. At the test point, (z,w,lambda) -> lambda rho~(V(z,w)) has differential rank 10 (A5). So its restriction to Omega x C^* has dense constructible image in the 10-dimensional M-hat.
2. By Chevalley's theorem, that image contains a dense open subset U1 of M-hat. Let D0 = U1 ∩ (M-hat)_smooth, minus 0.
3. Take p in D0, so p = lambda rho~(V) with V = V(z,w) and (z,w) in Omega. Then im d rho~_V has dimension 10 = dim T_p M-hat and lies inside T_p M-hat, so the two are equal.
4. Hence the conormal fibre of M8' at [p] is P(N(V)), which is contained in Sigma2 by Proposition 1.4.
5. Con(M8') is the closure of the P^3-bundle over the smooth part of M8'. Its restriction over P(D0) is dense in it.
6. So (M8')^* = pr2(Con(M8')), the closure of pr2 of that restriction, is contained in Sigma2.

*Proof of the dimension.*
1. Let E = {(V, xi) : V in U0, xi in N(V)}. It is a rank-4 vector bundle over U0, so it is smooth and irreducible of dimension 20.
2. Over the dense open set of V with rho~(V) smooth on M-hat, N(V) is the conormal fibre. So pr2(E) is contained in the affine cone over (M8')^*.
3. At (V0, xi0) the tangent space of E has dimension 20 = dim E, and its projection has rank 12 (A9). By lower semicontinuity of rank, and the constant rank theorem at a nearby point of maximal rank, the affine cone over (M8')^* has dimension at least 12.
4. So 11 <= dim (M8')^* <= dim Sigma2 <= 11. The last bound holds because Sigma2-hat is the closure of X + X and dim X = 6.
5. Both varieties are irreducible, as images of irreducible varieties, so (M8')^* = Sigma2. ∎

**Consequences for the paper.**
- **The statement holds as written.** The written sketch's three asserted steps are not needed: the general fibre of G, surjectivity of the polar differential, and the reverse inclusion via G. Neither is HMSV's density of G(Sigma2 \ M8).
- **All three uses are covered.** The inclusion is what (C.8) and step 9 of the assembly need. The equality is what biduality needs in step 11.
- **Corollary C.5 is complete.** Part 1's own proof needed dim Sigma2-hat = 12, which it had only numerically; Theorem 1.5 now proves it.
- **The paper's elliptic-quartic picture is consistent.** Inc_V maps to a curve in X ∩ P(N(V)), but the proof above does not use it.
- **Remark on Lemma C.7.**
  - The closed-embedding hypothesis enters Lemma C.7 only through the injectivity of B.
  - Once B is injective, the limit matrix [rho0, A, 2B(v . -)] has rank 10. So for small epsilon the image near rho(0, epsilon v) is an immersed 10-dimensional germ of M-hat.
  - At smooth points of M-hat on that germ, which are dense on it, the conormal directions lie in Con(M8'), and Con(M8') is closed.
  - So Lemma C.7's conclusion holds for the image variety M8' under (H3) alone.

## 2. G3: the embedding needed in Lemma C.7 (class (b))

**What is used.** At a stable Gale-fixed point p = [I4, O], take the local expansion rho(x, z) = rho0 + A x + B(z . z) + ... of Lemma C.7. The hypothesis is

  (H3)  the 13 vectors rho0, A e_1..A e_6, B(e_a . e_b) (a <= b) are linearly independent.

Equivalently: dim K_O = 7 and c -> S_O(c) maps K_O onto Sym^2(C^3)^*.

**Where it is used:**
1. In Lemma C.7 at the maximiser O, which is stable with d(O) != 0. This gives the rank-10 limit and hence ([rho(O)], [c]) in Con(M8').
2. In (C.8) at the same O. The identity C|K_O = alpha det S_O presupposes S_O onto.
3. In the purity argument of Lemma C.9 (Section 3). There it is needed on an open subset F° of the stable fixed locus that contains the C.9 example point. The local equation alpha needs S_p onto.

**Why the citation gives it.**
1. HMSV Props. 8.1-8.2, as quoted on p. 44, say the Gale-invariant ring R^g is generated in degree one.
2. Then M8-bar/<g> = Proj(R^g) maps to P^13 by a closed embedding. This is standard: a graded ring generated in degree one embeds its Proj.
3. A closed embedding makes the 13 dehomogenised coordinates generate the maximal ideal of the local ring of the quotient at p.
4. In equivariant coordinates that ring is C[[x, z_a z_b]]. Its cotangent space is minimally spanned by x_1..x_6 and z_a z_b (a <= b), 12 classes in all.
5. Their linear parts span it exactly when (H3) holds. The paper's argument on p. 49 is correct.

**A weaker statement would suffice, but I could not prove it.** (H3) is only needed at the maximisers that actually occur, plus at one point on each Pluecker divisor. One could try a genericity route instead: show that for (a, b, theta) in a dense open set, no critical point lies in the bad locus {(H3) fails}, a proper closed subset of F. That needs a dimension count for the critical-point condition, and the count can fail where the rank of A drops. I did not complete it.

**Evidence that (H3) holds:**
- Exact rank 13: part 1, at 2 rational O.
- Exact (C.8) with alpha != 0: part 2, at 4 matrices, including the boundary example of Lemma C.9.
- Numeric: 30 points (part 1).

**Status.** CITED. In my judgment it is a correct citation, not a gap. If (H3) failed at some stable O with d(O) != 0, Lemma C.7 and (C.8) would fail there. The assembly would then not cover parameters whose maximiser is that O.

## 3. G4 and the Pluecker divisors (class (a), given (H3))

### 3.1 The stable fixed locus F is irreducible

1. **Every stable self-associated configuration comes from an isotropic subspace.** Let Lambda in Gr(4,8) be a row space representing a stable self-associated configuration. Its Gale dual is Lambda^perp. So Lambda^perp = Lambda D for some diagonal D, and Lambda is isotropic for diag(d_k).
2. **All d_k are nonzero.** If d_k = 0, then e_k is in (Lambda^perp)^perp = Lambda. Then the other seven points lie in a plane. A plane may contain at most six points of a stable configuration, so this is impossible.
3. **Normalise the form.** Scaling the columns by sqrt(d_k) makes Lambda isotropic for q0 = sum x_k^2. This does not change the configuration.
4. **So F is an image.** F is the image of the stable part of OG(4,8), taken for q0.
5. **One component suffices.** A sign change of one column swaps the two components of OG(4,8) and fixes the configuration. So F is the image of the stable part of OG+(4,8).
6. **Conclusion.** OG+(4,8) is irreducible, so F is irreducible.

This also shows that F ∩ {Pluecker_1234 != 0} is the paper's chart {[I4, O] : O in SO(4,C), stable}, since O O^t = I is isotropy for q0 after scaling the O-columns by i.

### 3.2 The zero locus of alpha, intrinsically

Let F° be the open set of F where (H3) holds; it is all of F if HMSV holds. For p in F°:
- K_p is the 7-dimensional space of c that annihilate the affine tangent space, at rho(p), of the image of F, together with rho(p).
- Put Z = {p in F° : C|K_p ≡ 0}.

This definition involves no chart, no lift and no basis. The paper's gluing remarks (p. 51) are correct but no longer needed.

**Z is locally the zero set of one regular function.**
1. Near p, choose a regular frame of K, and equivariant coordinates.
2. det S is then a cubic form on K whose coefficients are regular functions. It is not identically zero, because S_p is onto.
3. Fix one coefficient of det S that does not vanish at p.
4. alpha is the matching coefficient of C|K divided by that coefficient. It is regular, and Z = {alpha = 0} near p, by (C.8).

**Z is a divisor.** F° is irreducible, and alpha is not identically zero (case 1 of Lemma C.9: the exact identities (C.5), (C.10), (C.9)). So by Krull's principal ideal theorem every component of Z has codimension 1. Smoothness of F is not needed for this.

**Z is S8-invariant.**
1. A relabelling sigma acts linearly on V_C, with m(z o sigma) = T_sigma m(z). It maps the image of F to itself, and K_{sigma p} to the corresponding image of K_p.
2. C(T_sigma c) = sgn(sigma) C(c). This is EXACT for adjacent transpositions (`moc4lib` test).
3. Hence C|K_{sigma p} ≡ 0 if and only if C|K_p ≡ 0.

### 3.3 Pluecker divisors (the part-2 item on p. 53)

For a 4-set S of the 8 labels, let D_S = {p in F : the points labelled S are dependent}.

**Claim.** Each nonempty D_S is irreducible of dimension 5, and D_S = D_{S^c}. S8 permutes the D_S transitively.

*Proof.*
1. **Dependence as an intersection.** For Lambda in OG+(4,8): Delta_S(Lambda) = 0 if and only if Lambda contains a nonzero vector supported on S^c, that is, Lambda ∩ W0 != 0 with W0 = span(e_j : j in S^c). W0 is nondegenerate for q0.
2. **Complements.** Lambda ∩ W0 != 0 if and only if Lambda + W0^perp != C^8. Take orthogonal complements and use Lambda^perp = Lambda: this says Lambda ∩ W0^perp != 0, and W0^perp = span(e_S). So D_S = D_{S^c}.
3. **The incidence variety I.** Let I = {(l, Lambda) : l an isotropic line in W0, l ⊂ Lambda, Lambda in OG+(4,8)}.
   - The isotropic lines of W0 form a smooth quadric surface, which is irreducible.
   - The maximal isotropic Lambda of the given family containing l correspond to the maximal isotropic subspaces of one family in l^perp/l. That space is 6-dimensional and nondegenerate, and those subspaces form OG+(3,6) ≅ P^3, which is irreducible.
   - I is closed in a product of projective varieties, and its projection to the quadric is surjective with irreducible equidimensional fibres. So I is irreducible.
4. **The divisor itself.** pr2(I) = {Delta_S = 0} ∩ OG+(4,8). It is irreducible. It is cut out by one equation on the irreducible 6-dimensional OG+(4,8), so it has pure dimension 5.
5. **Pass to F.** The stable part of pr2(I) maps onto D_S, so D_S is irreducible.
6. **Transitivity.** sigma in S8 permutes coordinates and preserves q0, so it maps {Delta_S = 0} to {Delta_{sigma S} = 0}. If sigma swaps the two components of OG(4,8), compose with one column sign change, which leaves the configuration unchanged. S8 is transitive on 4-sets. ∎

**Completion of Lemma C.9.**
1. Case 1 shows alpha != 0 at chart points where all Pluecker coordinates are nonzero and d != 0.
2. Points outside the chart lie on D_{1234}.
3. So every component of Z is either some D_S or lies in the closure of {d = 0} ∩ chart.
4. Suppose Z contained some D_S. By S8-invariance and transitivity, Z would contain the divisor through the C.9 example point. But alpha != 0 there: it is stable, d = -16/729, and alpha = 729/16 exactly (part 2). That is a contradiction.
5. Hence Z ∩ chart is contained in {d = 0}. So alpha(O) != 0 for every stable [I4, O] with d(O) != 0. ∎

This needs the example point to lie in F°. Its exact (C.8) value with alpha != 0 is consistent with that, and HMSV gives it anyway.

## 4. G1: (C.2) is no longer an input of Theorem C.1

**Where each part of (C.2) was used, and what replaces it:**

| Part of (C.2) | Used in Appendix C | Replacement |
|---|---|---|
| Sing C = M8, i.e. X = {grad C = 0} | Prop. C.6, "ker H(x) = T_x X" | Only the inclusion T_x X ⊂ ker H(x) is needed, from X ⊂ Sing C (EXACT by substitution, `moc4lib`). Part 1 already shows this. |
| Sigma2 ⊂ C | assembly, step 7 | X ⊂ Sing C gives C(s x + t x') = 0. Elementary. |
| closure of G(Sigma2 \ M8) = M8' | Prop. C.4; paper's proof of C.5 | Not needed: Section 1; part 1's proof of C.5 uses only G(Sigma2-hat) ⊂ M-hat (from (C.4) EXACT and Lemma 1.1). |
| dim Sigma2 = 11 | C.5, C.11 | Theorem 1.5. |
| det H restricted to C = gamma0 s^2 | Lemma C.11 only | Section 4.1. |

### 4.1 Lemma C.11 without s

Set W = {0,1,5,6}, K = the other ten coordinates, K' = K \ {13}, and q = u0 u6 - u1 u5.

**Lemma 4.1.** For every u in V_C,

  det(H(u) + h H(kappa)) = p0(u) + p1(u) h + p2(u) h^2 - 2 q(u) det H(u)[K',K'] h^3.

*Proof.*
1. **Structure of C in u13.** dC/du13 = q and d^2C/du13^2 = 0 (B1). So H(kappa) = Hess q, which is supported on W x W. Also H(u) e13 = Hess(q) u is supported on W. In particular H(u)_{k,13} = 0 for every k in K, 13 included.
2. **No h^4 term.** Expand the determinant by multilinearity in the rows W. The h^4 term is det(Hess q)_{WW} times det H(u)[K,K]. The column 13 of H(u)[K,K] is zero, so this term vanishes.
3. **The h^3 term.** It is a sum over 3-sets S = W \ {s}, T = W \ {t} of ± det(Hess q)[S,T] det H(u)[K ∪ {s}, K ∪ {t}]. The first factor is nonzero only when t is the partner of s, the pairs being 0↔6 and 1↔5.
4. **Expand the 11x11 minor.** Expand along column 13: its only nonzero entry is d_s q, in row s. Then expand along row 13: its only nonzero entry is d_t q, in column t. The minor becomes ± d_s q · d_t q · det H(u)[K',K'].
5. **So p3 = det H(u)[K',K'] (alpha u0 u6 + beta u1 u5)** for absolute constants alpha and beta.
6. **The constants.** Solving at two random rational u and checking at four more gives alpha = -2 and beta = 2 (B3). All 6 points also have degree <= 3 (B2). ∎

**Proposition 4.2 (uniqueness on the kappa-line).** Let a be real and y = Re(theta m(a, -b)). Suppose p2(y) != 0. Then y is the only point of Sigma2-hat on y + C kappa.

*Proof.*
1. y lies in Sigma2-hat.
2. q(y) = 0, by (C.6), which is EXACT.
3. Corollary C.5 gives rank H <= 12 on Sigma2-hat. So det H(y) = 0, and adj H(y) = 0, which gives p1 = tr(adj H(y) H(kappa)) = 0.
4. By Lemma 4.1 and step 2, det H(y + h kappa) = p2(y) h^2.
5. For h != 0 this is nonzero, so rank H(y + h kappa) = 14 and y + h kappa is not in Sigma2-hat. ∎

p2(y) = det H(y + kappa) is a real polynomial in (a, Re b, Im b, Re theta, Im theta). It is nonzero at the paper's example, where it equals -114105877589169227887017984; there rank H(y) = 12 and p0 = p1 = p3 = 0 (B4, B5). So it is nonzero on a dense open set.

**Proposition 4.3 (smoothness of y).** Let Z = V(C, all 13x13 minors of H). Sigma2-hat is contained in Z, by Corollary C.5. Suppose dC_y, together with the differentials at y of three fixed 13x13 minors of H, spans a space of rank at least 2. Then y is a smooth point of Sigma2-hat.

*Proof.*
1. The Zariski tangent space T_y Z lies in the common kernel of these functionals, so dim T_y Z <= 12.
2. Z contains Sigma2-hat, which is irreducible of dimension 12 (Theorem 1.5) and passes through y. So dim O_{Z,y} >= 12.
3. Hence O_{Z,y} is regular of dimension 12, and in particular a domain.
4. The prime ideal of Sigma2-hat in O_{Z,y} has a quotient of dimension 12, so that prime is zero.
5. Hence Sigma2-hat = Z near y, and y is smooth on Sigma2-hat. ∎

**The condition holds on a dense open set.**
- The three minors are fixed index sets: the first three (dropped row, dropped column) pairs in lexicographic order whose block has rank 12 at the example.
- The condition says a 2x2 minor of a 4x14 matrix that is polynomial in y does not vanish. So it is open in the parameters.
- At the example the rank is exactly 2 (B6). So the condition holds on a dense open set.
- G(y) != 0 is not needed.

**Conclusion.** Lemma C.11's two conclusions (y is the unique intersection with y + C kappa, and y is a smooth point) now hold on a nonempty open dense parameter set, and its proof uses no part of (C.2). The paper's statement "deg s <= 1 on lines in {C = 0, q = 0}" is not needed.

### 4.2 What is still cited from (C.2), and by whom

- **The set equality {C = 0, det Hess C = 0} = Sigma2-hat (HMSV)** is invoked by four supplement documents to apply rank <= 12 at points not known to lie on Sigma2-hat:
  - T_42: `final_referee_T41_bad54`, p.3.
  - T_45: `blind_moc4_T44_conic_branch_exclusion`, p.1; `independent_moc4_T44_a_Heron_branch`, p.1.
  - T_47: `final_referee_T46_negative_j_Heron`, p.2; `final_referee_T46_two_Heron_necessity`, p.2.
  - T_65: `root_moc4_T47_negative_Phi`, p.1; `final_referee_T47_negative_p_Heron`, p.1.
- **Neither my results nor Corollary C.5 prove that equality.** I only show Z = Sigma2-hat near the specific points y. For these cases it stays CITED (class (b)) unless the case auditors find a local replacement. Some sibling documents (T_40, T_42 `negative_Heron_region` and `bad17`) prove the rank bound by their own bordered-minor identities.
- **X = {grad C = 0} (HMSV ideal generation)** is used, according to the helper's reading, by the T_85, T_87, T_89 and T_90 documents. This is not Corollary C.5 or Proposition C.6. It is still CITED.

## 5. Lemma C.10: Lam's theorem and the density step (class (a))

### 5.1 Nonnegative weights without Lam

Positions 1..8 on R carry nodes t_1 < ... < t_8. Write Vand(t_S) = prod over i < j in S of (t_j - t_i) for a set S of positions. A noncrossing perfect matching pi is compatible with S if each of its pairs joins S and S^c.

**Proposition 5.1.** There are 14 polynomials w_pi in the gaps s_k = t_{k+1} - t_k, one for each noncrossing perfect matching pi of 8 points on a line, with nonnegative coefficients and not identically zero, such that for all 70 four-sets S

  Vand(t_S) Vand(t_{S^c}) = sum over pi compatible with S of w_pi.

*Proof.* The incidence matrix has rank 14 (C2). Solve on 14 independent rows. Then check all 70 equations as polynomial identities in s_1..s_7 (C3). Every coefficient is >= 0; the weights have 64 to 3144 terms (C4). ∎

So w_pi > 0 whenever all s_k > 0.

**Use in Lemma C.10 (entries nonzero, real O in SO(4) with d(O) = 0).**

1. **Twisted-cubic coordinates.** Lemma C.8 and the alternation step (PROOF-READ in part 2) give the following. In the t-order, the columns of [I4 | O] are A v(t_k) c_k, with A, c_k and t_k real. The coordinate columns S0 occupy the odd positions, or the even ones.
2. **Moduli of minors.** As points of Gr(4,8), with columns reordered, the two Pluecker vectors agree up to a scalar and position-dependent signs. Delta_S(A[v c]) = det A · prod over S of c_k · Vand(t_S). For S = (S0 \ I) ∪ J, |Delta_S([I | O]) Delta_{S^c}([I | O])| = det O[I,J]^2, because complementary minors of a special orthogonal matrix agree up to sign. Hence

   det O[I,J]^2 = Vand(t_S) Vand(t_{S^c}) / (Vand(t_{S0}) Vand(t_{S0^c})).

   No sign bookkeeping and no total positivity is needed.
3. **Weights.** Every noncrossing matching of a line joins odd to even positions. So every pi is compatible with S0 and defines a bijection sigma_pi from coordinate labels to O-labels.
4. **Compatibility.** A pair (i, j') is cut by S = (S0 \ I) ∪ J exactly when (i in I) is equivalent to (j in J). So pi is compatible with S exactly when sigma_pi(I) = J.
5. **Conclusion.** p_sigma = w_pi / (Vand(t_{S0}) Vand(t_{S0^c})) are probabilities. They satisfy sum over sigma(I) = J of p_sigma = det O[I,J]^2 for all I, J.

This replaces (C.12) and the citation of Lam and Postnikov entirely. The one-colour boundary issue noted by part 2 disappears.

### 5.2 The density step (p. 55)

**Identities used:**
- **(C.13)**: d(O) = O01 O02 O03 O10 O20 O30 on {O00 = 0} ∩ SO(4).
  - It is proved as a rational identity on the hypersurface T00 = 0 of the quaternion parametrisation T = L(1,a,b,c) R(1,d,e,f), solved for f (D0, D1).
  - The parametrisation omits only q0 p0 q3 = 0, which has dimension 4 inside the real 5-manifold {O00 = 0}.
  - So the identity extends to all of {O00 = 0} ∩ SO(4) by continuity.
- **(C.14)**: dd(E_i O) = -a b x_j x_k y_i z_i for i = 1, 2, 3. The common constant is -1 (D2, D3).
  - It is proved with a rational circle parametrisation of (a,b) and Euler-Rodrigues coordinates for [x y z] in SO(3).
  - Both are dense parametrisations, so the identity extends by continuity.
  - Part 1 and part 2 had these only at 12 and 15 samples.

**Corrected density argument.** Let O be indecomposable with zero entries and d(O) = 0. Then O is a smooth point of the real hypersurface H = {d = 0}, by (C.14) and the case analysis in the text.

Suppose that near O every point of H has some zero entry.
1. By Baire, some zero-entry divisor E_ij contains an open piece of H near a point O' close to O.
2. O' is still indecomposable: entries that are nonzero at O stay nonzero nearby. It also lies on H and on E_ij.
3. Near O', E_ij is smooth, since its differential e_i ∧ v_j is nonzero. So H and E_ij coincide there.
4. By (C.13), transported to (i,j), the product of the other six entries of row i and column j vanishes on an open piece of E_ij.
5. The ring of germs of real-analytic functions at a point of a connected manifold has no zero divisors. So one factor O_kl vanishes on an open piece of E_ij.
6. Then E_ij and E_kl coincide near an indecomposable point. That contradicts the independence of the differentials of two zero entries (p. 55, PROOF-READ by part 2).

Hence points with all entries nonzero are dense in H near O. The limit argument in the compact simplex finishes the proof. ∎

## 6. What Theorem C.1 rests on after this work

| Step | Status |
|---|---|
| (C.3), (C.4), (C.5), (C.6), Möbius rule | EXACT (parts 1/2; A0 here) |
| Lemma C.2, Lemma C.3 | proved (part 1) |
| Proposition C.4 | proved (Section 1) |
| Corollary C.5 | proved (part 1's argument plus Theorem 1.5) |
| Proposition C.6 | proved (part 1) |
| Lemma C.7 | proved, given (H3) |
| (C.8) | proved, given (H3), Lemma C.7 and Proposition C.4 |
| Lemma C.8 | proved ((C.9) EXACT) |
| Lemma C.9 | proved, given (H3) (Section 3) |
| Lemma C.10 | proved (Section 5, Lemma C.8, Birkhoff) |
| Lemma C.11 | proved (Section 4) |
| Assembly | PROOF-READ (part 2) |

**Remaining external inputs:**
- Theorem 1.2 (main text).
- (H3), from HMSV degree-one generation.
- Projective biduality in characteristic 0.
- Standard GIT: stable points with trivial stabiliser are smooth points of M8-bar, and equivariant local coordinates exist at Gale-fixed points (Luna slice).

None of these is in doubt in my judgment. Only (H3) is a specialised citation.

The supplement's `principal_minors` note claims a second proof of Theorem C.1 (Corollaries 10 and 12) that avoids HMSV. I did not audit it. Part 2 tested its central claim numerically.

## 7. Where Theorem C.1, Corollary C.5 and Proposition C.6 are used

**In the paper itself:**
- Proposition 1.4(i), p. 5, is Theorem C.1.
- p. 27 (Group IV discussion): repeated eigenvalues and a_i + b_j = 0 are "handled by Proposition 1.4". This is class (i).
- p. 27 derives C(u) = 0 and Hess C(u) v = 0 for x = u + iv in X by elementary means.
- Appendix F, p. 61, points to Corollary C.5 and Proposition C.6 for the rank and inertia restrictions.
- Section 5 (T_85) makes no use of Appendix C.
- Proposition C.4 is used only inside the proof of Theorem C.1. C.5 and C.6 do not depend on it.

**Classes in the supplement** (paper numbering via `../deps/dependency-map.md`):

- **(A) A full four-point spectrum on a line; Theorem C.1 indispensable.**
  - Which cases:
    - T_7: `root_moc4_T6`, p.2. "If v is real, then z is real, so the entire spectrum A lies on the real projective line."
    - T_11, both branches: `root_moc4_T10_B_branch`, p.2; `final_referee_T10_A_branch`, p.1.
    - T_24: `independent_moc4_T23_uniform`, p.3; `T23_harmonic_reduction`, p.2. Also the auxiliary T^(7)_6 (`T23_g13`, p.2), which feeds T_24.
  - In all of them the spectrum is real, so only the real case of Theorem C.1 (the core of Appendix C) is used. Lemma B.2 cannot replace it.
  - Downstream through the boundary graph, 51 cases have T_7, T_11 or T_24 among their transitive prerequisites: T_8, T_12, T_13, T_14, T_16, T_30, T_32, T_33, T_38, T_41, T_45-T_56, T_59-T_71, T_73-T_87 and T_89.
- **(B) A different statement, the "three-plus-three" / "six-circle" theorem.** It says that three eigenvalues of A and three of C lie on one generalized circle.
  - Which cases: T_4 (p.3, "previously proved three-plus-three generalized-circle theorem"); T_8 with T^(7)_1; T_9; T_10 (p.1 lists "the 3A+3C common-circle theorem" among proved inputs); T_12; T_13; T_15; T_16; T_17; T_18; T_19; T_20.
  - This is not Theorem C.1: in every case the off-circle points A_0 and C_0 are free.
  - No proof of it was found in the paper or in any of the 208 texts. T_1's Proposition 10 (`six_contact_reduction`, pp. 6-7) is an argument specific to T_1's contact set: two contacts give Re L = 0 and add contacts.
  - The helper notes a plausible per-case repair: two contacts sigma, tau with sigma(0) = tau(0) = 0 give the same entry-face argument. It was not checked.
  - **Status: GAP in the citations of these twelve case documents.** It lies outside Appendix C and is flagged for the case auditors.
- **(C) Corollary C.5 / Proposition C.6.**
  - Only det Hess C(Re x) = 0 is used, which is replaceable by the elementary fact (Hess C(u) v = 0 with v = Im x; if v = 0 then T_x X ⊂ ker): T_47 (`T46_whole_assembly`, p.2), T_53 (`inventory52_whole`, p.2), T_89 (`inventory88`, p.1), and cores B, C, D (T^(7)_8, T^(7)_11, T^(7)_12).
  - Inertia, or 12x12 principal minors >= 0, at Re x (Prop. C.6, which is proved; the elementary facts give only rank <= 13 and do not suffice): T_49 (`T48_negative_Phi`), T^(7)_19 (→ T_53), T^(7)_12 (→ T_51), T^(7)_22 (→ T_32, T_33), T^(7)_24 (→ T_89, T_90), T^(7)_25 (→ T_79; uses the positive-index bound). These documents prove their own rank bound by bordered minors.
  - The HMSV set equality, which goes beyond C.5 (Section 4.2): T_42, T_45, T_47, T_65. These stay CITED.
- **(D) Repeated eigenvalues or degenerate configurations only (class (i)); replaceable by Lemma B.2 and Theorem B.1.**
  - Which cases: T_1, T_2, the general reduction, and the collision sentences of T_32, T_33, T_45, T_49, T_50, T_53, T_54, T_58, T_59, T_61, T_64, T_66-T_71, T_73, T_75, T_76, T_78, T_80 and T_81.
  - The collision lines in T_3-T_6, T_9-T_20 and T_24 are also of this class.
  - A spectrum with a repeated eigenvalue is covered directly by Lemma B.2. A cross-spectrum common point is covered by Theorem B.1.

**Answer to "which of the ninety cases would be affected if Theorem C.1 or Proposition C.4 were unproved":**
- Directly: T_7, T_11 and T_24, through the real-spectrum case.
- Conditionally through the boundary graph: the 51 cases listed under (A).
- The class (i) uses are not affected, since Lemma B.2 and Theorem B.1 cover them.

Independently of Appendix C, case (B) needs attention in its own right. Its 12 root cases have 49 downstream cases; the union of (A), (B) and the four HMSV-equality cases with all their downstream cases covers 65 of the 90.

## 8. Files

- `appC_fill_checks.py` (sections A-D; deterministic; about 15 s).
- `appC_fill_checks_output.txt`: 30 PASS, 0 FAIL.
