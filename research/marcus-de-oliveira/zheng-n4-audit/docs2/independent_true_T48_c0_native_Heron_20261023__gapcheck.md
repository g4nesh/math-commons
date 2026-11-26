Evidence level: EXACT (every identity of p.2 eqs (3)-(7) and p.3 eqs (8)-(10) replayed by own sympy code from the chart definitions; one input, the representative of the c amplitude, is read as data from the cited companion JSON and then validated exactly against the rebuilt chart)

# Second-pass gap check: independent_true_T48_c0_native_Heron_20261023

Label C27-T49-gapcheck. Document: `proofs/manin_transfer_20261021/independent_true_T48_c0_native_Heron_20261023.pdf` (3 pp.).
Replay: `replay2/independent_true_T48_c0_native_Heron_20261023__gapcheck.py`, output in `replay2/independent_true_T48_c0_native_Heron_20261023__gapcheck_output.txt`.

## Deviations and limits

- The first auditor's files were not opened. Nothing from the supplement was executed or imported.
- One supplement JSON was read as data (`independent_true_T48_native_harmonic_conditional_exact_20261023.json`, keys `affine_amplitudes/c` and `curve`). This is the file the note cites for "the stated coefficient arrays". Both arrays are checked against my own rebuild before they are used (see step 2).
- The chart itself (B = (t,0,1,inf), A = (y,z,w,x), eq. (4), and the variables J, W, M, ZD) is taken from the printed text of the companion note `independent_true_T48_native_harmonic_conditional_20261023` (pp. 2-3). That note was not audited here beyond what this gap needs. The chart's coverage claims are CITED.
- The `p` amplitude formula was checked by an exact Jacobian test: the difference is a function of (r,s,u,h) on the contact surface. It is not shown to be a polynomial shift. That is what "modulo real shifts" needs, because h is real on the surface.
- The upper-Heron chain lemma and the other results used on p.3 ("The remaining branch") are outside this gap and are not examined.

## The reported gap

p.2, eqs (3)-(5). The claim is that the three real rows Re P = 0, (Im P)/y = 0 and Im(C_N conj C_D)/y = 0 are linear in (RK, R, K), that det M3 = -s N^2 m F_p T_N, that the second Cramer numerator is N_R = -s N F_p T_N R0, and that m N + v R0 = t H. The first auditor said the rows depend on an unprinted representative of the c amplitude and that the checker is not in the supplement.

## Classification: FILLED (severity: presentation)

The first auditor's facts are correct. The c amplitude is not printed in the PDF, and the note's own checker (`independent_true_T48_c0_native_Heron_20261023.py/.json`) is not in the supplement. Their conclusion does not follow, for two reasons.

1. **The step does not depend on the choice of representative.** By Appendix D, f39 is the only one of the eight canonical generators that is nonzero at permutation 17, and its value there is 1. So c = g_17 exactly, and c_amp can be taken to be H_17, which is canonical. Any other representative c' = c_amp + phi with phi real (a function of the real parameters r, s, u, h) gives Im(C'_N conj C_D) = Im(C_N conj C_D) + phi |C_D|^2 Im(1) = Im(C_N conj C_D). The third row is therefore the same polynomial for every representative with the same denominator C_D. Its zero set {Im c_amp = 0} is intrinsic. The choice of representative affects only the polynomial shape of the row, that is, whether it is linear in (RK, R, K).
2. **The representative the note uses is supplied as data and is valid.** The JSON cited by the note contains the affine c amplitude, and its denominator is J(W-1)M, which matches "denominator Z{rJ-(r-s)Z}" after W = 1 - Z/J. I checked exactly that c_json - H_17, restricted to the contact surface h = H_9(J,W), is a function of (r,s,u,h) alone. The Jacobian d(diff)/dJ * dH9/dW - d(diff)/dW * dH9/dJ vanishes identically. So c_json is a legitimate representative, and Im c_json = Im H_17 = c on the surface.

With that representative, every assertion of the step replays exactly from the definitions.

### Argument and replay (all EXACT unless stated)

Step 1 (chart). From B = (t,0,1,inf), A = (y,z,w,x) and eq. (4) of the harmonic note, I rebuilt all 24 F_sigma and H_sigma = (F_sigma - F_0)/E. The contact values H0 = 0, H3 = s-r, H5 = -r, H6 = -1 and H12 = -u come out as stated. With z = -M/ZD, the numerator of H9 - h is the curve C(J,W), and its [W^2], [J^2] and W^0 coefficients equal eq. (6) of the harmonic note. The supplied `curve` array equals my C.

Step 2 (short coordinates, Heron eq. (1)). The exclusive gaps are j = g2, b = g7, c = g17, p = g11, and c+d = g1 - g7. I checked:
- j_amp = H2 = J;
- (c+d)_amp = H1 - H7 = (J-s)/Z;
- b_amp = H7 = ((r-1)J+s)/Z - (r-s);
- (a+d)_amp and (q+d)_amp agree with the printed formulas up to the real constants u-1 and 1-u on the surface;
- p_amp agrees up to a function of (r,s,u,h), by the Jacobian test. At 4 curve points the imaginary parts agree to 12 digits (NUMERIC).

Taking imaginary parts with J = x+iy and Z = z+iyv gives each line of (1) in a few lines:
- Im((J-s)conj Z) = y(z - xv + sv) = y Dn, so d/y = Dn/N when c = 0;
- Im(((r-1)J+s) conj Z) = y((r-1)t - sv), so (b+d)/y = rt/N;
- Im(-(RZ+J+r-s)/r) = -y(Rv+1)/r;
- Im(-(u-1)(J-s)/J) = (K/s) Im(1 - s/J) = Ky/A;
- Im(R(Z+r) conj(J+r-s)) = R y F_p.

Step 3 (rows). With P(J,Z) = J C(J, 1 - Z/J), C_D = Z{rJ - (r-s)Z} and C_N = c_json C_D (both are polynomials), substitute J = x+iy, Z = z+iyv, u = 1 - K/s and h = R - r + s. All three rows are polynomials in Y = y^2. Rows 1 and 2 carry the monomials {RK, R, K, 1} with a common factor 1/s, and row 3 carries {R, K, 1}. So each row is linear in (RK, R, K). This matches the claim.

Step 4 (eq. (4)). det M3 = -s N^2 m F_p T_N, an exact polynomial identity with the normalization above. That is the printed formula exactly, which also shows the note used this normalization.

Step 5. N_R / (-s N F_p T_N) is a polynomial R0. Eq. (5), m N + v R0 = t H, holds identically, and so do both identities of (6).

Step 6. Cramer's rule in the form det(M3) * R = N_R holds without dividing, for every solution of the linear system. With s N F_p T_N != 0 (s != 0 is a chart condition, N > 0, F_p != 0 since p > 0, and T_N > 0, see below), cancelling gives N m R = R0, and (5) then gives (7), m N (Rv+1) = t H. The cleared form of (7) was also checked directly.

Step 7 (T_N > 0). r J - (r-s) Z = J M. If M = 0, then z ZD = -M = 0. Either z = 0, which is the collision A_2 = B_2, excluded by eight-distinctness, or ZD = 0. In the second case W = -s/(r-s) is real, so J = r/(W-1) is real, contradicting y = j > 0. N > 0 because Z != 0.

Step 8 (m != 0). This follows the note's p.2 argument, which I re-read. If m = 0, then R0 = 0 by step 6, so tH = 0 by (5), so H = 0 because rt > 0. If r != 1, the second identity of (6) gives Dn = -r^2 < 0, contradicting Dn = N d/y > 0. If r = 1, then v = 0, and (1) gives (q+d)/y = -1 < 0, contradicting q+d > 0. (PROOF-READ, all inputs EXACT.)

Step 9 (eq. (9)). With L as in (8) (linear in RK, R, K), I checked m (L_RK N_RK + L_R N_R + L_K N_K + L_0 det M3) = -A t Dn (m - H) det M3 identically. Since det M3 != 0 once m != 0, the Cramer values are the actual (RK, R, K), and (9) follows. Eq. (8) follows from (1) by direct expansion (L = a(B+dj) + (d-j)B with a/y = K/A - Dn/N and B/y^2 = -t(Rv+1)/N). I re-derived it by hand.

Step 10 (eq. (10) and the theorem). B/y^2 = -t(Rv+1)/N = -t^2 H/(m N^2) by (7). F6/y^6 = (L/y^3)^2 - 4(d/y)^2 B/y^2 = t^2 Dn^2 (m+H)^2/(m^2 N^4) >= 0, checked as a rational identity. Since B = (b+d)(q+d) > 0 (d > 0 and b, q >= 0), H/m < 0, so (m-H)/m > 1. Then sign L = -sign(t) = -sign(r), because rt > 0.

The step is therefore justified. The polynomials it needs come from the printed chart, and from the cited JSON for the c representative, which is itself validated against the chart. What is missing is the note's checker and a printed representative, which makes this a presentation issue only. The only mathematical subtlety, the dependence on the representative, is answered by observation 1.

## Falsification

Not needed for an exact identity. The NUMERIC agreement of the p amplitude at 4 curve points (25-digit evaluation, seed 5) was only a cross-check.
