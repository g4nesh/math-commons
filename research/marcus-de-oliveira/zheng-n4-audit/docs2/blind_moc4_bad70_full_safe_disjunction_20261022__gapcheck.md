Evidence level: EXACT for the key polynomial identities (own sympy code, replay2/blind_moc4_bad70_full_safe_disjunction_20261022__gapcheck.py); PROOF-READ (argument supplied by this auditor) for the inequality that fills the gap; NUMERIC cross-checks only as corroboration.

# Second-pass gap check: blind_moc4_bad70_full_safe_disjunction_20261022 (paper case T_34)

Label: C13-T32-T34-gapcheck. Document: "The contact inequality for T34", 3 pages, contact set T = {2, 3, 5, 9, 10, 12}, generators (f0, f1, f2, f12, f39, f67, f70, f72) (replayed: equals `generators(T)`).

## Deviations and limits

- Only the one reported gap (p.1, Section 1, the sufficiency of (H)) was examined. Sections 2-4 (radial equation, R9, Bernstein expansions) and the (PB)/(CQ) Gram claims were not audited here.
- The first auditor's files were not opened. Nothing from the supplement was executed or imported; I read the T_34 text extract and, for orientation only, the supplement note "A scalar compensation inequality for rotations" (final_referee_T2_scalar_compensation). The argument below does not use that note; it is self-contained.
- One standard external fact is used: every O in SO(4) is the map x -> p x conj(q) on the quaternions for some unit quaternions p, q (the classical double cover S^3 x S^3 -> SO(4)).

## The gap

Claim (p.1, Section 1): if d > 0, b + c > 0 and a >= d + d^2/(b+c), then Phi_O(g) >= 0 for every O in SO(4), for g as in (1).

The document's argument is this. Split d in the ratio b : c, use Cauchy's inequality for the bookkeeping, and then appeal to "the two endpoint inequalities". These are

 (EP_b) Phi_O(a' f2 + b f12 + d' f70) >= 0 whenever b > 0 and a' >= d' + d'^2/b,
 (EP_c) the same with (f39, c) in place of (f12, b).

They are never stated as lemmas and never proved. The phrase "scalar compensation and quaternion Gram identities give..." does not supply them. The (PB) and (CQ) Gram forms do not imply them either, because those need p > 0 or q > 0 and here p = q = 0 is allowed.

The bookkeeping itself is correct. With d_b = db/M, d_c = dc/M and M = b + c, the thresholds satisfy (d_b + d_b^2/b) + (d_c + d_c^2/c) = d + d^2/M.

Reduction (mine). Write x2 = det O[13,34], x3 = det O[14,24] and x4 = det O[14,13]. By (4.2), Phi(f2) = x2^2, Phi(f12) = x3^2 and Phi(f39) = x4^2 (replayed EXACT). Minimizing over b > 0 gives inf_b [(1 + 1/b) x2^2 + b x3^2] = x2^2 + 2|x2 x3|. Hence the family (EP_b) for d' = 1 is equivalent to the single pointwise inequality

 (E1) Phi_O(f2 + f70) + 2 |x2 x3| >= 0 for all O in SO(4),

and (EP_c) is equivalent to (E2), the same inequality with x4 in place of x3. Conversely, (H) at (b, c -> 0) implies (E1). So the claim under dispute is true if and only if (E1) and (E2) hold.

(E1) is not a triviality. Phi(f2 + f70) is itself negative somewhere (NUMERIC min about -0.0255), and Phi(f70) has min about -0.0674. Both absolute-value terms are needed: Phi(f2+f70) - 2|x2 x3| reaches about -0.277.

## Classification: FILLED (the step is true; the complete argument follows)

### Step 1. Coordinates (EXACT identities)

Let O = O(p, q) be x -> p x conj(q). Write P = W p and Q = W q, where W = (1/2)[[1,1,1,1],[1,-1,-1,1],[1,-1,1,-1],[1,1,-1,-1]] is orthogonal. Since W is orthogonal, (P, Q) ranges over S^3 x S^3, and every O in SO(4) arises. Put a_i = P_i^2 and A_i = Q_i^2, so on the unit spheres sum a_i = sum A_i = 1. Put v_i = a_i - A_i, so sum v_i = 0. Finally put

 delta = P0 Q3 + P1 Q2 and Delta = P2 Q1 + P3 Q0.

The script checks the following as polynomial identities in the 8 variables, homogenized with |P|^2 and |Q|^2. It also checks O^T O = |P|^2 |Q|^2 I and det O = (|P|^2 |Q|^2)^2.

- x2 = v0 + v1, x3 = v1 + v2, and hence x2 x3 = v0 v2 - v1 v3.
- KEY: Phi_O(f2 + f70) = 2 [ delta^2 + Delta^2 - (v0 v2 + v1 v3) ].

### Step 2. Reduction of (E1)

By Step 1, Phi(f2+f70) + 2|x2 x3| = 2 [ delta^2 + Delta^2 - (v0 v2 + v1 v3) + |v0 v2 - v1 v3| ]. This equals 2 [ delta^2 + Delta^2 - 2 min(v0 v2, v1 v3) ]. If the minimum is <= 0, we are done.

Otherwise v0 v2 > 0 and v1 v3 > 0. Since sum v_i = 0, the pairs (v0, v2) and (v1, v3) have opposite strict signs.

### Step 3. Lemma (square-root increment)

For r, s, u, v >= 0 we have sqrt((r+u)(s+v)) - sqrt(rs) >= sqrt(uv). Proof: (r+u)(s+v) - (sqrt(rs) + sqrt(uv))^2 = rv + us - 2 sqrt(rsuv) >= 0 by AM-GM, and both sides are nonnegative.

### Step 4. Case A: v0, v2 < 0 < v1, v3

Put alpha = -v0, gamma = -v2, M = v1 and N = v3, all > 0. Then A0 = a0 + alpha, A2 = a2 + gamma, a1 = A1 + M and a3 = A3 + N.

By the triangle inequality, |delta| >= |X| with X = sqrt(a1 A2) - sqrt(a0 A3), and |Delta| >= |Y| with Y = sqrt(a3 A0) - sqrt(a2 A1). Regroup the sum:

 X + Y = [sqrt((A1+M)(a2+gamma)) - sqrt(A1 a2)] + [sqrt((A3+N)(a0+alpha)) - sqrt(A3 a0)].

By the Lemma, X + Y >= sqrt(M gamma) + sqrt(N alpha) > 0. Hence

 delta^2 + Delta^2 >= X^2 + Y^2 >= (X+Y)^2 / 2 >= (sqrt(M gamma) + sqrt(N alpha))^2 / 2 >= 2 sqrt(M N alpha gamma) = 2 sqrt((v0 v2)(v1 v3)) >= 2 min(v0 v2, v1 v3).

### Step 5. Case B: v0, v2 > 0 > v1, v3

Put alpha = v0, gamma = v2, M = -v1 and N = -v3. Then a0 = A0 + alpha, a2 = A2 + gamma, A1 = a1 + M and A3 = a3 + N.

Use X' = sqrt(a0 A3) - sqrt(a1 A2) and Y' = sqrt(a2 A1) - sqrt(a3 A0). Then X' + Y' = [sqrt((A0+alpha)(a3+N)) - sqrt(A0 a3)] + [sqrt((A2+gamma)(a1+M)) - sqrt(A2 a1)] >= sqrt(alpha N) + sqrt(gamma M). The chain of Case A then gives the same bound. This proves (E1) on all of SO(4), including every zero-coordinate boundary, because no division occurs.

### Step 6. (E2) from (E1) by symmetry

Let g = (rho, rho', 0) with rho = rho' = the reversal 4321. On matrices g acts by g.O = J O J, where J is the reversal permutation matrix (an even permutation, so g.O is in SO(4)). On permutations it acts by sigma -> rho sigma rho. Then w_{g.sigma}(g.O) = w_sigma(O) directly from the definitions, so Phi_{g.c}(g.O) = Phi_c(O).

The script checks exactly that g fixes the value vectors of f2 and f70 and sends f12 to f39, and that the Leibniz-sum invariance holds at 5 rational SO(4) points. P is injective, so these are equalities of vectors. Hence Phi_O(a' f2 + c f39 + d' f70) = Phi_{JOJ}(a' f2 + c f12 + d' f70) >= 0.

### Step 7. Conclusion for (H)

(E1) gives (EP_b): a' x2^2 + b x3^2 + d' Phi(f70) >= d' [x2^2 + (d'/b) x2^2 + (b/d') x3^2 + Phi(f70)] >= d' [Phi(f2+f70) + 2|x2 x3|] >= 0. Likewise (E2) gives (EP_c).

Then:
- If b, c > 0, the document's split together with (EP_b) and (EP_c) gives Phi(a f2 + b f12 + c f39 + d f70) >= 0. The excess (a - d - d^2/M) f2 and the other minor generators contribute >= 0.
- If c = 0, (H) is (EP_b) itself. If b = 0, it is (EP_c). The document's "by continuity" is unnecessary.

So the claim under dispute holds.

### Corroboration (NUMERIC)

- min of (E1) and (E2) over SO(4): about -7e-30 and -2e-30 (40 Nelder-Mead starts each, seed 11). The bound is sharp, with equality attained.
- (H) taken at equality a = d + d^2/(b+c), d = 1, over 16 (b, c) pairs (c = 0 included) with 10 BFGS starts each: min = 2e-14.

## Severity

As written the step is a substantive missing argument. The endpoint inequalities are genuine SO(4) inequalities with equality cases, and the document's text does not supply a proof of them. It is now filled by the argument above (EXACT identities + PROOF-READ inequality).

## Side observation (not assessed)

Reading the (PB) Gram form naively, its vector entries for indices 0, 2, 3 would be +-O14, +-x2, +-x3 with fixed signs. Then Phi(f70) + 2 s1 O14 x2 + 2 s2 x2 x3 would have to be a square. Numerically it takes negative values for all four fixed sign choices (min about -1.9 or lower over 4000 random O). The same holds for the analogous (CQ) reading with O22, x2, x4. Pointwise-varying signs or different Gram vectors could still make (PB) and (CQ) correct. This was outside this check and should be looked at by whoever audits Section 1's (PB)/(CQ) claims.
