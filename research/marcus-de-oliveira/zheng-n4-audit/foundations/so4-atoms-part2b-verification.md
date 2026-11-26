Evidence level: EXACT for every Gram identity, closed-form threshold, allocation, orbit membership, transport and source identity listed below (own code `so4_atoms_part2b_verify.py`, written without reusing the part-2b script; captured run `so4_atoms_part2b_verify_output.txt`). PROOF-READ for the PSD and closure arguments, which are short linear-algebra steps written out below. NUMERIC (attacks over SO(4)) as supporting evidence only. CITED where a step rests on part 2a (the T24 face) or on the four core atom families (part 1 and its re-verification).

# Independent re-verification of so4-atoms-part2b (label atoms2b-verify)

Scope: the universal SO(4) inequalities and allocation lemmas listed in [so4-atoms-part2b.md](so4-atoms-part2b.md), groups G1 (T44), G3 (T46), G4 (T47), G2 (T45), G5 (allocations), checked in that order. Historical labels and paper f_j indices are used throughout (0-based, Appendix D via `moc4lib`).

## Deviations and limits

1. **Data from the supplement.** As permitted, the cubic feature basis, the four sign maps and the congruences T (T44, T46, T47) and V (T47) are read from supplement JSON as data. Nothing else is taken from it. In particular I did not use their cost columns, relation matrices R or D, or any stored hashes. I computed the Leibniz sums myself from the f_j value vectors and my own quaternion model, then derived the full relation system myself.
2. **Choice of coordinates.** I fixed the model (x -> p x q), the sign-map convention and the coordinate order of each face by testing, not by reading the other auditor's code.
   - The other three model variants leave the costs outside the span of the Gram columns.
   - The order of each face's g-coordinates was read off from which unit Gram represents which single generator.
   - This choice cannot make a false statement pass. Whatever basis is used, the identity is checked against the true Leibniz sums of the generators the document names, with the coefficients it claims.
3. **Numerics are evidence, not proof.** Each statement was attacked at 4 to 9 coefficient samples placed exactly on the threshold. Each sample had 300 batched Adam starts over S^3 x S^3, plus a BFGS polish of the best 5. The worst point was re-evaluated in mpmath at 40 digits. Seeds: numpy `default_rng(101..150)`. Pass tolerance: -1e-9 on min Phi / max|coef|.
4. **Not re-checked here:**
   - The T24-face inputs (part 2a): the coherent bounds behind the T27/T38 eighth chamber, T36 (7), T38 gen. 40, T47 IIIb_minus_swapped and T47 (3). Only their group transports were checked.
   - All radial, polynomial, Bernstein, Sturm and chart-coverage parts of the documents.
   - The cone-classification claims ("exactly", Theorem 2 of T46 core, the converse part of T46 six-ray).
5. **Faces with more generators.** The T44, T46 and T48 contact sets have 9, 10 and 10 generators (T44 adds f0, T46 adds f0 and f3, T48 adds f1 and f9). All Gram statements concern the 8-generator faces the documents fix. The reduction from the full cone to those faces is outside this file.
6. **Process.** All runs used RUNJOB. One attempted run found no free slot and did nothing.

## Method for the Gram statements (EXACT)

- **Model.** O(p,q) is the matrix of x -> p x q, which maps S^3 x S^3 onto SO(4). For each face I form the 36 Gram polynomials sum_{parity, 4 sign maps} F_a(Sy) F_b(Sy) and the 8 Leibniz sums 4 Phi_f(O(p,q)) (degree 8). Each is reduced to normal form modulo N = |p|^2 - |q|^2 (division with leading term p0^2), and the identity checked is (|p|^2+|q|^2) SOS_M = 4 Phi_c (mod N).
- **Ranks.** For T44, T46 and T47 alike, the Gram columns have rank 18, the costs rank 8, and the costs lie in the Gram span (total rank 18). The solution space {(M, c)} therefore has dimension 26, as the other auditor reported.
- **Checking a document Gram M.** Pull M back to M_F = W^T M W. Then check, symbolically in all of the document's parameters, that every one of the 18 independent relations vanishes with the document's coefficient vector c.
- **Rational-point spot check.** Each identity was also checked at one rational parameter point and one rational point of S^3 x S^3 (exact arithmetic): SOS = 2 Phi_c holds exactly. So a PSD M gives Phi_c >= 0.
- **Control.** Reading the T44 Prop. 1 d-coefficient as |D|^2 = 1 + v^2, the naive "diagonal" reading, is rejected. The coefficient read-off is not purely diagonal: the D-diagonal equals d + a on T44 and d + c on T46.

## G1. Historical T44 (face q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15)

Sources: blind_moc4_T44_coherent_chain_capacity p.1-2; T44_conic_branch_exclusion p.1; independent_moc4_T44_c_below_j p.1. My reading of each statement matches part 2b.

- **Prop. 1** (r != k, w > 0; p >= x^2+y^2, j >= r^2(1+v^2), c >= k^2(1+v^2), a, e >= v^2, b >= w^2, q >= 0, d = 1).
  - Gram of the eight stated R^3 vectors: identity verified symbolically in (r, k, v, w), with c = (0, x^2+y^2, r^2(1+v^2), k^2(1+v^2), v^2, 1, w^2, v^2). PSD: it is a Gram of real vectors. Surplus coefficients add minors.
  - **CERTIFIED.**
- **Prop. 2.**
  - The specialisation v = 0, r = sqrt j, k = -z, w = sqrt b gives (1+z^2/b)/(sqrt j+z)^2, with derivative 2(z sqrt j - b)/(b (sqrt j+z)^3). So the minimum over z in [0, sqrt c] is at z = min(sqrt c, b/sqrt j). Both branch values and their agreement at jc = b^2 are verified, as are the j = 0 case (1/b + 1/c) and the example (1/2) f1 + f53 + f12 + f7 + f18 >= 0.
  - The b = 0 direct Gram (J = (0, sqrt j, 0), P = (0, -1/sqrt j, 0), D = (1, 0, 0)) satisfies the identity with p = 1/j and d = 1.
  - **CERTIFIED.**
- **Scalar test (1).** The substitution v = t/sqrt(1-t^2), r = sqrt j sqrt(1-t^2), k = -sqrt c sqrt(1-t^2), w = sqrt b gives exactly the displayed threshold, with j and c thresholds equal to j and c. Also v^2 <= m iff t^2 <= m/(1+m). Requires J + C > 0 and b > 0, as stated. **CERTIFIED.**
- **SOC and Hsum.**
  - SOC(1,53;18), H(15,53;18) and H(8,53;18) are in the proved orbits with these roles. The reversed old-H roles (53,15;18) are not, so the role order matters.
  - Splitting d in ratio a:e puts both old-H atoms exactly at threshold, summing to j = d + d^2/(a+e).
  - **CERTIFIED.**
- **J\*.** Some g in G maps (1,53,12,7,18) to (53,1,8,7,18). The transported Prop. 2 gives exactly (1) of c_below_j with roles p <-> j, c -> e. The document's literal "(1,3,0,2) + transpose" action was not decoded; only existence of g was checked. **CERTIFIED.**

## G3. Historical T46 (face q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4)

Sources: final_referee_T46_core_coherent_capacity p.1-5; T46_six_ray_coherent_capacity p.1-4; T46_a0_all_small_product Lemma 1; root_moc4_T46_positive_a_mixed (1). Statements match part 2b.

### Core Thm 1 (q, b > 0, l = qb - beta^2 > 0, j > (1+(q+b-2beta)/l) t^2)

- **Identity.** The Gram X G0 X^T (base Q, B, J0, D; J = J0 + tD; P = lambda(Q+B) + mu J0) satisfies the identity with c = (q, p0, j, 0, 0, 1, b, 0), symbolically in (q, b, beta, t, j). Also |P|^2 = p0 = threshold (1).
- **PSD.**
  - det of the (Q, B, J0) block = l (j - (1+(q+b-2beta)/l) t^2), so the hypotheses give positive leading minors and G0 is PD.
  - The Gram of (Q+B, J0) has determinant Delta.
- **Formula (2) and the optimisation.** Verified exactly: the derivative, u - 4/s = (q-b)^2/(s l), the t\* value 1/j + (s+2k)/4 + k^2/(4(s+4)), its beta = 1-q-b form 1/j + k(k+3)/(k+4), the t = 0 threshold and derivative, and the example 399/1543.
- **Cor. 3** (beta = 1 - S, t = 0 gives 1/j + 2 - S).
  - The strict leaf condition is l > 0.
  - Closure: the only critical point of qb - (q+b-1)^2 is (2/3, 2/3), where the value is 1/3 != 0. So every point of the closed region (S-1)^2 <= qb is a limit of strict points, and the threshold is continuous. Valid coefficient vectors form a closed set (pointwise limits of Phi_c >= 0), so the closed statement follows.
  - The T47 version uses the transport [1,59,6,27,28] -> [45,5,7,0,35], verified.
- **CERTIFIED.**

### Six-ray Thm 1

- **Identity.** Base (Q, B, J0, Z); P = lambda6 S + mu6 J0; C = kc(J0 - (v/b)B); D = C - Z. The identity holds with c = (q, p\*, j, c\*, 0, 1, b, 0), symbolically in (q, b, beta, t, v, j). |P|^2 = p\* and |C|^2 = c\*.
- **PSD.** det of the (Q, B, J0) block = l(j-t^2) - (b t^2 - 2 beta t v + q v^2). So the hypothesis is exactly PD of G6. Then Delta > 0 (independence of S and J0) and j - t^2 - v^2/b > 0 (a Schur complement).
- **CERTIFIED.**

### Singular leaf (Lemma 1 of a0_all_small_product; positive_a_mixed (1))

- **Setup.** q = u^2, b = w^2, beta = -uw, v = -rt/q. The identity re-verified at this substitution (EXACT), and p\* = k + (1-eta t)^2/(j - kappa t^2), c\* = rho^2 t^2/(j - kappa t^2).
- **PSD directly, with no limit needed.**
  - The leaf block has rank one and (t, v) = (t/u)(u, -w) lies in its range, so the generalised Schur complement is j - kappa t^2 > 0.
  - Delta = (u-w)^2 (j - kappa t^2) > 0 because q != b.
  - j - t^2 - v^2/b = j - kappa t^2.
- **Rational section.** With t = eta j/(D0 + A0): j - kappa t^2 = jN/(D0+A0)^2, N - rho^2 A0/c = D0^2 + D0 A0 + A0^2 > 0, consumed p = k + D0^2/(jN). Checked in T47 roles, which is the same algebra.
- **Cor. 2 and Cor. 4 of six-ray (T47 roles):**
  - the t = eta e/kappa form;
  - the e-condition difference identity;
  - both Cor. 4 branches (t = eta e/D0 and t = 1/eta), with unused-p identities;
  - the residual-difference identity versus Cor. 3.
- **CERTIFIED.**

### Remaining T46 items

- **four(4,6;2;28):** in the proved orbit with these roles. **CERTIFIED** (reduction).
- **Theorem 2 / "exactly" classification:** UNCHECKED (not a universal inequality).

## G4. Historical T47 (face q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7)

Sources: independent_moc4_T47_PCBE_energy p.1-2; T47_IIIb_plus_large_pc p.1; T47_merged_C_high_j p.1-2. Statements match part 2b.

### Item 4 of the brief: the "spectral coordinates"

- **What the coordinates are.** In W = V^-1 T^-1, each unit Gram E_ii represents exactly the generator of its own letter for all seven minors (q, p, j, c, e, b, a). D is the non-minor slot. So these are diagonal-coefficient coordinates for the T47 face in the documents' order (Q, P, J, C, E, D, B, A).
- **Why the budget coordinates fail.** In the budget coordinates T^-1 only five unit Grams are pure, under different letters. The PCBE matrix read in T^-1 violates 18 of 18 relations.
- **Conclusion.** Every identity below is checked against the true Leibniz sums of f2, f9, f1, f45, f5, f35, f0, f7 with the documents' coefficients. The statements are therefore about the **same generators and coefficients** the T47 documents use. The basis V^-1 T^-1 is only a certificate basis (det V = 1, det T = 1/4). No issue.

### Statements

- **PCBE Thm 1.**
  - **Identity:** the free matrix satisfies the identity identically in (x, y, z, w, L, j, c, e), with c = (0, L, j, c, e, 1, L, 0).
  - **Structure:**
    - P and B coincide as vectors (equal rows, M_PP = M_BB = M_PB = L), so PSD reduces to the matrix without P.
    - The Schur complement over {B, D} equals the displayed 3x3 matrix.
    - The chosen JE entry diagonalises the Schur complement of C0.
  - **Scalar facts:**
    - identity (3);
    - simplification to K(1-x)^2 - 1 and B x^2;
    - x\* solves the crossing.
  - **Proof-read steps:**
    - (1) <=> (2) by the interval-intersection argument;
    - clamping to [0, 1];
    - the closure argument: c + eps gives strict inequalities, so admissible y, z exist, and valid coefficient vectors are closed.
  - **CERTIFIED.**
- **large_pc Lemma 1.**
  - Identity holds with c = (0, p, d^2/(a+c), c, 0, d, 0, a). The products of J are as stated.
  - det of the (P, C, A) Gram = a(pc - a^2), so the Gram is PSD iff pc >= a^2. D is orthogonal and J lies in the span.
  - **CERTIFIED.**
- **merged-C Lemma 1.**
  - All six Grams (crossing, stationary, half-energy; leaf a or b) satisfy the identity. Bases: (J, leaf, D) with E = -sqrt(e/j) J, or (J, E, leaf, D).
  - The inner products and |C|^2 values are as stated, including crossing |C|^2 = d^2(L+e)/(L(sqrt j+sqrt e)^2).
  - PSD:
    - crossing: always, since it is a vector construction;
    - stationary: needs je >= L^2;
    - half-energy: needs e >= j.
  - Proportional merge: all three bound types sum exactly to the A0 formulas, and j_i e_i - a_i^2 = delta_i^2 (je - A0^2).
  - **CERTIFIED.**
- **Atoms.** H(9,1;35), SOC(45,1;35), chain(1;7,9;35) and four(0,7;5;35) are all in their orbits with these roles. The T24 core -> T47 transport (59,53,55,67,54,2) -> (45,5,7,0,35,1) holds. **CERTIFIED** (reductions; the T24 role itself is CITED, part 2a).
- **zero_q_rational:** reduces to six-ray Cor. 3 via the verified transport [27,1,59,5,28,6] -> [0,45,5,9,35,7]. **CERTIFIED.**
- **Two of the three items part 2b left UNCHECKED, now checked:**
  - **Old stated core c >= 1/e + k:** this is core Thm 1 at t = 0, beta = -sqrt(ab), transported to T47.
    - The threshold equals 1/e + b(1-r)^2/(b-r)^2 (EXACT).
    - The Gram has a rank-one leaf block and Delta = (sqrt b - sqrt a)^2 e > 0 when a != b, so the construction is directly PSD.
    - It is used where b > 1 > a.
    - **CERTIFIED.**
  - **PCBE_low_e_q_monotone:**
    - Its use of PCBE is checked: B = 1/e + 1/(b+1) iff e <= b(b+1), and failure of (2) with c = s^2 B is exactly j < J\*. **CERTIFIED** (chart algebra).
    - The q-monotonicity of the radial G (40-thousand-term tables) is radial: UNCHECKED.
- **Still unchecked or cited:**
  - **high_j_small_b conic:** radial/conic certificates, not a universal inequality. UNCHECKED.
  - **IIIb_minus_swapped and T47 (3):** CITED (part 2a).

## G2. Historical T45 (p f0 + r f1 + a f2 + b f12 + c f59 + q f67 + d f70 + s f72)

- **(1):**
  - the contact nullspace has dimension 8;
  - the eight generators vanish on the contacts;
  - their values at (18,16,13,8,2,14,11,7) form the identity matrix;
  - f70 is the only non-minor.

  **CERTIFIED.**
- **(2)-(5):** SOC(2,59;70), chain(2;0,12;70), chain(59;0,1;70) and chain(59;67,72;70) are in the orbits with the roles the coefficient forms need (centre = the d^2/... coefficient). **CERTIFIED** (reduction). The note's own small PSD8 Grams were not rebuilt.
- **(6):** H(1,59;70) and H(72,59;70) are in the orbit (aux, main). Splitting d in ratio r:s sums exactly to c = d + d^2/(r+s). **CERTIFIED.**
- **Zero faces:** g maps T33, T6 and T14 onto (3,5,9,10,12,13), (2,3,4,5,9,10) and (3,4,5,9,10,16). Each image lies in T45 plus the added contact, and those contacts are 13 = a, 2 = c and 16 = r by (1). **CERTIFIED**; the three input families are CITED.

## G5. Allocation lemmas

- **Source identities (EXACT).** All 18 response identities of T28, T42, T38, T36 and T27 hold, and the generator lists of T28, T42, T38 and T27 are as stated. All 13 faces have rank 8, vanish on their contacts and have exactly one non-minor generator. For T36:
  - the nullspace has dimension 8;
  - delta35 + delta40 = -g17;
  - decompositions (1), (2) and (3) hold on the whole nullspace.
- **Symmetries.** The T42 22<->54 exchange, the T27 40<->68 exchange and the T36 (a,b,q,p) exchange are each realised by a stabiliser of the contact set.
- **H/chain lemma (6), and the T48 allocations.**
  - **Lemma (6), algebra:** with b_H = b eta1/beta1, b_C = b delta1/beta1 (and the same for c), the two old-H and two chain atoms at threshold use b and c exactly, have non-minor total 1, and need exactly a = RHS of (6) (symbolic).
  - **Lemma (6), roles:** all 16 role memberships pass. The tuples (b, c, p, q, a, k) are T28 (33,5,4,6,2,46), T38 gen. 70 (12,39,0,67,2,70), T42 (3,23,0,6,5,54) and T36 (10) (2,29,0,7,5,40). The T36 tuple is also the transport of the T28 face.
  - **T48 general identity (EXACT):** e-side sum = s - alpha + alpha^2/a + s^2/e, and j-side sum = (1-s) - gamma + gamma^2/Q + (1-s)^2/j. Optimising alpha gives phi_a, hence T(s) of Cor. 2.
  - **T48 allocations I, II, III, IV:** u + v = 1 and the closed forms of p\* are verified, with roles H(2,5;40), H(29,5;40), chain(5;0,2;40), chain(5;7,29;40) and chain(5;41,29;40).
  - **Verdict: CERTIFIED** (reductions).
- **Single atoms.**
  - T28 (3): H(45,1;34), chain(1;0,45;34), four(6,0;5;34).
  - T38: H(2,5;40), four(0,7;12;40), four(0,67;2;40), SOC(2,5;70).
  - T27: four(67,41;39;40), four(0,67;2;40), H(2,5;40).
  - T42: SOC(5,2;54).
  - T36 (6): SOC(45,1;35), H(5,45;35), four(7,0;5;35); plus T36 SOC(5,45;40).

  All pass. **CERTIFIED.**
- **Transports.** T24 core -> T38, T27, T36 and T47; T28 face -> T36. All exist. **CERTIFIED** (as transports).
- **Eighth chamber (T27 (10)-(11), T38 p.3).**
  - The four atom (delta a, c, p, delta) is exactly on its boundary at p(delta).
  - p(delta) is increasing on [0, 1). The derivative numerator is (a+1)(ac^2 + 2acdelta - delta^2), which is >= c + 2delta - delta^2 > 0 when ac >= 1.
  - The rescaled coherent bound equals r^2/e - r^2/E_r, and e >= E_r is the rescaled four-coefficient test.
  - **CERTIFIED as a reduction; the coherent input is CITED (part 2a).**

## Numerical attacks (NUMERIC)

All minima are >= -1e-16; see the output file. A failed control below means the weakened variant was not caught, not that a statement failed. Controls:

| control | min |
|---|---|
| T44 Prop. 1, p, j -> 0 | -0.025 |
| T44 Prop. 2, p -> 0 | -0.058 |
| T46 core, p -> 0 | -0.050 |
| T46 six-ray, p -> 0 | -0.005 |
| PCBE, p, b -> 0 | -6e-4 |
| large_pc, j -> 0 | -0.024 |
| H/chain T42, a -> 0 | -0.034 |
| SOC pj = 0.9 | -3.5e-7 |
| old-H lowered 10% | -3.5e-6 |

- **Weak control on Prop. 1.** The T44 Prop. 1 control with p alone removed stayed nonnegative on its samples (p is not tight there). It was replaced by removing p and j, which is detected at -0.025. Final run: 246 PASS, 0 FAIL.
- **Thresholds are sufficient, not sharp.** Several stated thresholds are not sharp, so removing a single coefficient is the meaningful control.

## Verdicts

| item | verdict |
|---|---|
| G1 T44 Prop. 1, Prop. 2 (incl. j = 0, c = 0, b = 0), scalar test (1), SOC/Hsum, J\* | CERTIFIED |
| G3 T46 core Thm 1, (2), t\*, Cor. 3 (closed region) | CERTIFIED |
| G3 T46 six-ray Thm 1, singular leaf / rational section, six-ray Cor. 2-4 | CERTIFIED |
| G3 four(4,6;2;28) | CERTIFIED |
| G4 T47 PCBE Thm 1, large_pc Lemma 1, merged-C Lemma 1, atoms, zero_q_rational, old-core c >= 1/e + k | CERTIFIED |
| G4 "spectral coordinates" question | same generators and coefficients; no issue |
| G2 T45 (1)-(6), zero faces | CERTIFIED (zero-face families CITED) |
| G5 identities, H/chain (6), T48 I-IV and Cor. 2 identity, single atoms, transports | CERTIFIED |
| G5 eighth chamber | CERTIFIED reduction + CITED part 2a |
| T47 IIIb_minus_swapped, T47 (3), T38/T27/T36 coherent bounds | CITED (part 2a) |

**Unchecked (complete list for this group):**
- radial, polynomial, chart and Bernstein parts of all listed documents;
- PCBE_low_e q-monotonicity of G;
- high_j_small_b conic certificates;
- the cone-classification ("exactly") statements of T46 core Thm 2 and six-ray;
- the T45 note's own PSD8 Grams for (2)-(5) (the reduction to atoms replaces them);
- the documents' literal (row, col, transpose) actions (only existence of g was checked);
- the reduction of the 9/10-generator T44, T46 and T48 cones to the 8-generator faces.

No ERROR was found. No repair was needed.
