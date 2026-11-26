Evidence level: (R) EXACT + PROOF-READ by an **own, independent proof** (no Gram data used); (F1), (F2) EXACT + PROOF-READ by an own re-verification of the other auditor's Gram certificates (own code; the supplement's cubic features read as data); (F3), (F4) EXACT reductions to proved atoms. NUMERIC attacks are reported separately and are not part of any proof. 96 exact checks, 0 failures (`so4_atoms_part2a_verify_output.txt`).

# Independent verification of SO(4) atoms, part 2a

Label `atoms2a-verify`. This file verifies the claims of [so4-atoms-part2a.md](so4-atoms-part2a.md) (label `atoms2a`). The other auditor's proofs were not read before my own attempt at (R), (F1) and (F2). For (R) I found an independent proof. For (F1) and (F2) I did not, and switched to a line-by-line re-verification. Files: `so4_atoms_part2a_verify.py` (all code; deterministic seeds), `so4_atoms_part2a_verify_output.txt` (captured run).

## Deviations and limits

- **Independent proofs.** (R) was proved from scratch by an elementary argument (section 1). For (F1) (5), (6), double-SOC and (F2) (9), (10) I did not find an independent proof in reasonable time. For those I re-verified the Gram certificates with my own code (step 2 of the task). They therefore rest on the supplement's cubic features, as the other auditor's proofs do. This does not weaken the proof: the identity is checked exactly, and any exact SOS identity proves the inequality, whatever the polynomials' origin.
- **Data handling.** The feature cubics, the four signed coordinate maps, the congruence T, the source order, the budget diagonal and the relation matrix R are read with `json.load`. The polynomial strings are parsed by my own tokenizer, with no `eval`. Nothing of the supplement is executed or imported. The data's R matrix is used only for a consistency cross-check, never in a proof.
- **Convention not taken on trust.** I tried the four congruence conventions g = T F, T^{-1} F, T^T F and T^{-T} F. Only g = T^{-1} F gives an identity. Its constant is kappa = 2 on all three faces.
- **Scope.** As in the other report, "every real orthogonal matrix" means SO(4) (O(4) in the Phi/det form). Only sufficiency was checked. The documents' necessity claims for the fixed degree-six cone and their spectral arguments (radial roots, Bernstein tables, Heron sectors) are not audited here. Only the atom applications were.
- **Uses.** I read 2 to 4 applying documents per item, listed in section 6. I did not search all ~46 documents that name the face Gram families.
- **Compute.** All runs went through RUNJOB. The exact part takes about 10 s, the numeric part a few minutes.

## Verdicts

| item | verdict | basis |
|---|---|---|
| (R) rank-two, coefficient form (1) and pointwise form (2), orbit 144 | **CERTIFIED** | own proof (section 1), plus the Gram certificate re-verified (section 2) |
| (F1) T24 (5) coherent, (6) parallel, double-SOC, with both extensions | **CERTIFIED** | re-verified Gram certificates; PSD region derived by me symbolically |
| (F2) T25 Lemma 1 (9), (10), with the forced choice (12) | **CERTIFIED** | re-verified; the forcing of (12) re-derived |
| (F3) T41 bad17 (3), (4) | **CERTIFIED** | exact transport plus formula identities |
| (F4) inventory56 (1), (2), (3) | **CERTIFIED** | chain allocation identity, orbit roles, (2) = (R) |

No error was found and no repair was needed. Two remarks only: a typo in the T24 PDF and an interpretation of the T25 statement of (10) (section 6).

## 1. (R): own proof

**Statement** (final_referee_inventory56_rank_two_capacity, p.1, read from the page image):
- Coefficient form: for a, j, d > 0 and b >= 0, if b >= d/4 + d^2/(4 min(a, j)), then b Phi12 + j Phi2 + a Phi29 + d Phi35 >= 0.
- Pointwise form (2): Phi35 + Phi12/4 + sqrt(Phi12 (Phi2 + Phi29)) >= 0.
- Here f2 = n13,34, f12 = n14,24, f29 = n13,12 and f35 = h4312,-.

**Coordinates (EXACT, R1).**
- Take unit quaternions u, v. Put p = u conj(1+i) and q = v conj(j+k), and let Ohat be half the matrix of x -> p x conj(q). Then Ohat^T Ohat = |u|^2|v|^2 I and det Ohat = (|u|^2|v|^2)^2.
- As (u, v) ranges over S^3 x S^3, Ohat ranges over all of SO(4).
- With a_i = u_i^2 and A_i = v_i^2, as homogeneous polynomial identities:
  - Phi35 = 2 sum_i a_i A_i + 4(u0u1u2u3 + v0v1v2v3);
  - x12 := det Ohat[14,24] = x with x = alpha + alpha' - 1, where alpha = a0 + a3 and alpha' = A0 + A3. Phi12 = x^2.
  - Phi2 + Phi29 = 2[(u1u2 + u0u3)^2 + (v1v2 + v0v3)^2].

**Reduction to two variables per quaternion (PROOF-READ; identities EXACT, R2).**
1. Put m = |u0u3|, n = |u1u2| and similarly m', n'. Then u0u1u2u3 >= -mn and (u1u2 + u0u3)^2 >= (m - n)^2. So the left side E of (2) satisfies E >= 2 sum a_iA_i - 4mn - 4m'n' + x^2/4 + |x| sqrt(2[(m-n)^2 + (m'-n')^2]).
2. Use a0A0 + a3A3 = [(a0+a3)(A0+A3) + (a0-a3)(A0-A3)]/2 and |a0 - a3| = 2e with e = sqrt(alpha^2/4 - m^2). This gives a0A0 + a3A3 >= alpha alpha'/2 - 2ee', and the same for the (1,2) pair.
3. Write (m, e) = (alpha/2)(c1, s1) and (n, f) = (beta/2)(c2, s2), with c, s >= 0 and c^2 + s^2 = 1; the same with primes. Here beta = 1 - alpha. If alpha = 0, take c1 = 0. Then 1 - s s' >= (c^2 + c'^2)/2, since 2 - 2ss' - c^2 - c'^2 = (s - s')^2.
4. Together: E >= Qd + x^2/4 + |x| r, where
   - Qd = alpha alpha'(c1^2 + c1'^2)/2 + beta beta'(c2^2 + c2'^2)/2 - alpha beta c1c2 - alpha' beta' c1'c2',
   - y = alpha c1 - beta c2, z = alpha' c1' - beta' c2', and r = sqrt((y^2 + z^2)/2).

**The two-variable inequality (PROOF-READ; identities EXACT, R2).**
- Swapping u and v preserves everything, so assume delta := alpha' - alpha >= 0.
- Put h = |x|/2. Since r >= |z - y|/2, we get Qd + x^2/4 + |x| r >= Qd + h^2 + h(z - y).
- The right side equals (f_u + f_p)/2 (exact identity), where
  - f_u = (y - h)^2 + delta(alpha c1^2 - beta c2^2),
  - f_p = (z + h)^2 - delta(alpha' c1'^2 - beta' c2'^2).
- **f_u >= 0 on [0,1]^2.**
  - If delta = 0 this is clear. Otherwise alpha' > 0, and f_u = alpha alpha'(c1 - k/alpha')^2 + (delta/alpha')(k^2 - alpha' beta c2^2), with k = beta c2 + h >= 0.
  - So it suffices that k >= sqrt(alpha' beta) c2.
  - If beta >= alpha', this is clear.
  - If alpha' > beta, then h = (alpha' - beta)/2 and c2 <= 1. So k - sqrt(alpha' beta) c2 >= (alpha' - beta)/2 - sqrt(beta)(sqrt(alpha') - sqrt(beta)) = (sqrt(alpha') - sqrt(beta))^2/2 >= 0.
- **f_p >= 0.** f_p(alpha, alpha', c1', c2') = f_u(1 - alpha', 1 - alpha, c2', c1') (exact). The substitution preserves alpha <= alpha' and |x|.
- Hence E >= 0, which is (2).

**Coefficient form.** Let x = min(a, j) and let b0 = d/4 + d^2/(4x) be the threshold, so b - b0 >= 0. Then b Phi12 + j Phi2 + a Phi29 + d Phi35 >= d(Phi35 + Phi12/4) + (d^2/(4x)) Phi12 + x(Phi2 + Phi29). By AM-GM (exact square identity in R2) this is >= d[Phi35 + Phi12/4 + |x12| sqrt(Phi2 + Phi29)] >= 0. This holds on the whole stated region, the boundary b = b0 included.

**Orbit.** All 144 members of the orbit follow by the covariance of `so4-atoms.md` section 1. T1 checks the orbit size, 144.

**Side observation (NUMERIC only).** The pointwise form (2) is sharp in its compensating term. Halving the sqrt term gives a negative minimum (-2.7e-3), and so does replacing sqrt(Phi2 + Phi29) by |x29| alone (-1.2e-2). So both leaves are needed pointwise, not only in the Gram cone.

## 2. Gram certificates re-verified (F1, F2, and (R) a second time)

**Method (own code).**
- Model: x -> p x q, with O(p,q) the unnormalised matrix. For each face I build g = T^{-1} F per parity.
- M_ik = sum over the parities and the four sign maps S_j of g_i(S_j y) g_k(S_j y).
- Delta(G) = |p|^2 sum G_ik M_ik - 2 Phi_{cost(diag G)}, reduced to normal form mod N = |p|^2 - |q|^2. Each Phi_f is computed with lib `F_VALUES` and Leibniz weights.
- The 36 columns of Delta give a matrix with 394 normal-form rows. Its row space is my relation space.
- **Lemma (PROOF-READ).** If G is PSD with Delta(G) = 0, write G = V^T V. On unit (p, q) we have N = 0, so 2 Phi_cost = sum |V g(S_j y)|^2 >= 0. The map onto SO(4) is surjective.

**Checks (EXACT, G0-G3).**
- On all three faces the data's source order equals the printed one, and the sources are generators of the stated contact sets (lib `generators`).
- The non-minors are f35; f54 (T24 also has f70 and f1 off the Gram, both with zero coefficient); and f70.
- The relation space has dimension 10 on every face and equals the row space of the data's R.
- **inventory56:** the ten relations printed in the PDF (rank-two p.2) are independent and span exactly my relation space. This is an independent link between the data and the PDF text.
- Every construction satisfies the relations identically in its parameters, checked symbolically. Each also passes an end-to-end identity at a rational point.

**Constructions and their PSD regions (own symbolic LDL, G4).**
- (R), G4 = v v^T + d x w w^T: the cost is exactly 4x[x f2 + x f29 + (d/4 + d^2/(4x)) f12 + d f35]. G4 is a sum of rank-one terms.
- (F1) (5): the last pivot is e - [(1-x)^2/j + L x^2]. At x = 1/(1 + jL) it equals e - 1/(j + 1/L), which is >= 0 iff j >= 1/e - 1/L.
  - The other pivots are j, 1, b, (ab-1)/b, all > 0 for ab > 1.
  - j = 0 member: last pivot [e(ab-1) - (a+1)(b+1)]/(ab-1), which is >= 0 iff e >= L. This is the j = 0 case of the printed condition.
  - The printed homogeneous threshold is d times the normalised one.
- (F1) double-SOC: pivots j, e - 1/j, 1, b, a - 1/b. So it is valid for **every ab >= d^2** with ej >= d^2 (extension confirmed).
- (F1) (6), with b = r^2, a = s^2, tau = (1 - rs)/(r - s), on the order (E, D, B, A, J):
  - pivots e, 1, r^2, then 0 with a zero residual row, then j - 1/e - tau^2;
  - tau^2 = b(1-w)^2/(b-w)^2 = a(1-w)^2/(a-w)^2 with w = sqrt(ab);
  - so the bound is **symmetric and valid for every a != b** (extension confirmed), including s = 0, where tau^2 = 1/b.
  - ab < d^2 is not needed for safety.
- (F2) Gbase pivots: a, hb, q - z^2/hb, c, 1. G = K Gbase K^T is PSD iff q hb >= z^2.
  - The cost coefficients are e = e(t) of (8) and p = r^2/c + (r-z)^2/(hb+q-2z).
  - (9), z = r: p = r^2/c, and PSD iff q hb >= r^2.
  - (10), z = wr, q = w^2 r^2/hb: PSD boundary (third pivot 0), hb + q - 2z = (hb - wr)^2/hb > 0, and p = r^2/c + hb r^2 (1-w)^2/(hb - wr)^2.
  - With P = kappa V + (r/c) C and kappa free, exactly one relation fails. Solving it gives kappa = (r - z)/(hb + q - 2z), which is **the choice (12) of the PDF, forced** (G3).
- Larger e, p, j or q are minor additions (f1, f2, f5, f67 are minors). So the safe sets are upward closed.

## 3. (F3) T41 bad17 (final_referee_T41_bad17_uniform p.2; large_product p.1)

- A group element maps (59, 53, 55, 67, 54) to (1, 53, 4, 16, 17) in the printed order (T2). The roles are j_T24 = p, e_T24 = j, b_T24 = a, a_T24 = b.
- (3), A >= 1: this is exactly (5) (sympy). At A = 1 it is the double-SOC.
- (4): (1-z)^2/(S-2z) = a(1-z)^2/(a-z)^2 when b = z^2/a. This is the tau^2 of (6) with r^2 = a and s^2 = b.
  - The hypothesis a > z is guaranteed by S >= 1 + z.
  - The ab = 0 limit 1/S is the s = 0 member of (6).
- large_product (1) is the homogeneous form of (5).
- Direction: the reduced source (J <= j, c = 0) is shown safe by the atom, and the original cost adds (j - J) f53 + c f6, which are minors. This is consistent with the proved statement.

## 4. (F4) inventory56 (boundary17_uniform p.1-2; boundary17_large_product p.1)

- (1): with alpha = qj/(q+j), beta = ea/(e+a) and d1 = d alpha/(alpha+beta), d2 = d beta/(alpha+beta), we get d1 + d2 = d and d1^2(1/q + 1/j) + d2^2(1/e + 1/a) = Bch (sympy). large_product (1) is the same inequality.
  - Chains (12; 0, 2; 35) and (12; 7, 29; 35) are in the proved chain orbit with centre f12 (T3). Their leaves are disjoint.
- (3) is the four-coefficient atom (7, 0; 5; 35), in the orbit with the c slot on f5 (T3).
- (2) uses (R) with the same roles: j on f2, a on f29, b on f12, d on f35, and the same min(a, j).

## 5. Numeric attack (NUMERIC)

Method:
- Multistart BFGS over (p, q) in R^8, normalised, x -> p x conj(q).
- 25 coefficient points on each stated threshold, 20 starts each (500 per family). The pointwise form of (R) gets 500 starts.
- Boundary cases are included: a = j; e = L; ab = 1; a = 0 in (6); t = 0, w = 0 and w = 1 in (10); ab = 0 in T41 (4).
- The worst point is re-evaluated with mpmath at 50 digits.
- Weakened variants (threshold coefficient halved or dropped) test detection.

Results (seed base 20261125; family seeds 11-52):

| family | in-region min (BFGS) | mpmath 50 digits | weakened variant | weakened min |
|---|---|---|---|---|
| (R) coefficient boundary | 1.0e-17 | 1.0e-17 | b = 0 | -6.2e-2 |
| (R) pointwise (2), 500 starts | 2.2e-16 | 2.2e-16 | sqrt term x 1/2 | -2.7e-3 |
| (F1) (5) incl. e = L, ab near 1 | 5.7e-18 | 5.7e-18 | j, e halved | +5.6e-17 |
| (F1) (6), a != b, incl. a = 0 | 5.8e-18 | 5.8e-18 | j halved | +9.4e-18 |
| (F1) double-SOC, ab >= 1, ej = 1 | 9.1e-18 | 9.1e-18 | j halved | +5.7e-18 |
| (F2) (9) incl. t = 0 | 7.6e-18 | 7.6e-18 | p = 0 | -3.8e-2 |
| (F2) (10) incl. w = 0, 1 | 9.6e-18 | 9.6e-18 | p = 0 | -4.4e-2 |
| (F3) T41 (3) | 1.8e-17 | 1.8e-17 | p, j halved | +7.6e-17 |
| (F3) T41 (4) incl. ab = 0 | 5.4e-18 | 5.4e-18 | p = 0 | -3.9e-2 |
| (F4) (1) b = Bch | 2.5e-18 | 2.5e-18 | b = 0 | -3.2e-2 |
| (F4) (3) four boundary | 7.2e-18 | 7.2e-18 | c = 0 | -5.0e-2 |

- No in-region negative value was found. The minimum is 0, attained, for example, at signed permutation matrices of contacts.
- Where the weakening dropped a coefficient, the attack found clear violations. Where it only halved one (the (5), (6), double-SOC and T41 (3) rows), the minimum stayed nonnegative. So those thresholds are not pointwise sharp, which agrees with the other auditor. The PDFs claim sharpness only for the fixed degree-six Gram cone.
- Evidence only; no proof above uses it.

## 6. Uses checked (step 5)

- **(R):** inventory56_boundary17_uniform p.2, condition (2) (same roles, same min(a, j), same diagonal (0, x, 0, 0, d/4 + d^2/(4x), x, d)); inventory56_boundary17_large_product p.1 (same face and roles); inventory56_57_whole_assembly (names the rank-two certificate only as an alternative).
- **(F1):**
  - T41_bad17_uniform p.2 and large_product p.1: see section 3.
  - blind_moc4_T27 p.2: transport (59,53,55,67,54) -> (5,2,0,67,40), in the orbit with the printed order. Its (4) and (5) are (6) in both orientations (b >= 1 or a >= 1 with ab < 1 forces a != b). Its (6) is (5). ab = 1 is the double-SOC.
  - final_referee_T25 p.2: (6) and (7) are (5) and (6) under j -> f2, e -> f39, {b, a} -> {f6, f67}. Both orders are realised by group elements. The "symmetric version of (7)" is covered by the a != b extension.
- **(F2):** T25 p.5, region II (18): there q hb = bq(1 - t^2/(ab)) > 1 - t^2 >= r^2, so (9) applies (identity in G4). Region IIIb- (20) uses (9) and (21) uses (10): q of (21) equals w^2 r^2/hb with w = W/(1+W) in [0, 1) (exact), and hb > wr. Same document, so the same roles.
- **(F3), (F4):** sections 3 and 4.
- **Remarks.**
  - T24 p.2 prints the Gram diagonal as "(q, p, j, c, e, q + 1, b, a)". The data's budget diagonal reads (q, p, j, c, e, q + d, b, a), which is the d = 1 normalisation. All constructions use q = 0, so this does not matter.
  - The T25 statement of (10) ("if q = w^2 r^2/hb ...") is an equality condition. Larger q is a minor addition, or is covered by (9) once q hb >= r^2. The PDF's use in (21) is at equality.

## 7. Status

Every claimed inequality holds on its whole stated region, boundary included. (R) is proved independently of the supplement data. (F1) and (F2) are proved by exact identities that I re-derived, with the PSD regions derived by me. (F3) and (F4) reduce exactly to proved atoms with the roles as used. The extensions found by the other auditor (double-SOC for every ab >= d^2; (6) for every a != b) are confirmed.
