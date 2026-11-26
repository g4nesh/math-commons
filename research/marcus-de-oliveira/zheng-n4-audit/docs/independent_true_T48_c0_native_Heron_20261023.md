Evidence level: EXACT for the short coordinate functions, the c = 0 relations (1), identities (6), the derivation of (8) and of (10) from (7) and (9); NUMERIC for (7) and (9) themselves (2363 stratum points); GAP for the three-row system (3)-(5), whose coefficients are not printed and whose checker is absent.

# Audit: independent_true_T48_c0_native_Heron_20261023 (3 pp.)

"A Heron identity and orientation condition on the T48 c = 0 face". Batch C27-T49. A branch used for paper case T_49.

## Deviations and limits

- No deviation from DOC_AUDIT.md. Own code only: `replay/independent_true_T48_c0_native_Heron_20261023__check.py` rebuilds the chart from the explicit inverse of the companion note (independently replayed in that note's audit) and the paper's definitions.
- The note's three real rows (3), with matrix M3, right side b3 and Cramer numerators, depend on a particular representative CN/CD of the affine c amplitude. That representative is "supplied by the standalone checker", which is not in the supplement, and the note does not print the rows. I therefore could not replay (4), the factorization of N_R or (5) exactly. I tested their consequences (7) and (9) numerically on points of the seven-contact stratum. I built those points independently: contact 17 holds exactly when c = 0, since g17 = c.
- The "remaining branch" paragraph (p. 3) is a remark outside Theorem 1. Its citations (the upper-Heron chain lemma, the j <= d theorem) were not audited.

## A. Statements

**Theorem 1 (p. 1).** Let x in X be a realizable source with eight distinct spectral points and contacts {0,3,5,6,9,12,17} (= Table 5 #23, EXACT), in canonical coordinates (q,p,j,c,e,d,b,a) on rays (7,5,29,39,2,40,41,0). Assume c = 0, p, j, d > 0 and q, b >= 0. With B = (b+d)(q+d), L = a(B + dj) + (d-j)B and F6 = L^2 - 4d^2j^2B, we have F6 >= 0 and sign L = -sign r, where r = -H5 is the real chart parameter of the companion note. This is a necessary condition only: it proves no support inequality, and repeated eigenvalues and a_i + b_j = 0 are excluded. The note itself says the branch r > 0, j > d is unresolved.

## B. Dependencies

| kind | ref | used for | status |
|---|---|---|---|
| other-doc | independent_true_T48_native_harmonic_conditional_20261023 (chart, explicit inverse) | coordinates | audited in this batch; the chart was replayed exactly |
| other-doc | "standalone checker" (absent) | rows of (3), (4), N_R, (5) | not available |
| other-doc | blind_moc4_T48_c0_upper_Heron_chain_20261023 | remark: r < 0, j >= d branch | CITED, not used in Theorem 1 |
| other-doc | "entire j <= d theorem" | remark | CITED, unresolved, not used in Theorem 1 |
| seven-contact | Table 5 #23 | statement | EXACT |

## C. Steps

| # | claim | kind | status | check |
|---|---|---|---|---|
| 1 | short coordinate functions (j, c+d, a+d, q+d, b, p) | identity | EXACT | each equals my dual-basis amplitude up to a real shift (0, 0, -1, -h-r+s, 1, h+r) |
| 2 | (1) on c = 0: d/y = Dn/N, (b+d)/y = rt/N, (q+d)/y = -(Rv+1)/r, (a+d)/y = K/A, p/y = R Fp/abs(J+r-s)^2 | identity | EXACT | sympy with J = x+iy, Z = z+iyv |
| 3 | Dn > 0, rt > 0, Fp != 0, N > 0 from the positivity assumptions | argument | PROOF-READ | |
| 4 | (2) TN = abs(rJ - (r-s)Z)^2 > 0, i.e. M != 0 at eight-distinct sources with j > 0 | argument | PROOF-READ | rJ - (r-s)Z = JM, and M = 0 forces J real (companion note) |
| 5 | (3) three real rows linear in (RK, R, K); (4) det M3 = -s N^2 m Fp TN; N_R = -s N Fp TN R0; (5) mN + vR0 = tH | identity | **GAP** (rows not printed, checker absent) | consequence (7) holds NUMERICALLY at 2363 points |
| 6 | (6) H + (r-1)Fp = -rm, (r-1)(Dn + r^2) - H = r^2 m | identity | EXACT | |
| 7 | (7) mN(Rv+1) = tH, valid even when m = 0 | identity | NUMERIC (derivation correct given step 5) | Cramer's rule det*x = N_x holds for every solution |
| 8 | wall m = 0 impossible | argument | PROOF-READ | uses (7), (6), Dn > 0, q+d > 0 |
| 9 | (8) numerator of L/y^3 | identity | EXACT | derived from (1) |
| 10 | (9) L/y^3 = -t Dn (m-H)/(m N^2) | identity | NUMERIC (2363 points, no failure) | |
| 11 | (10) B/y^2 = -t^2H/(mN^2), F6/y^6 = t^2 Dn^2 (m+H)^2/(m^2 N^4) | identity | EXACT given (7), (9) | |
| 12 | B > 0 gives H/m < 0, so sign L = -sign t = -sign r | argument | PROOF-READ | |

Counts: 12 steps; 5 EXACT, 4 PROOF-READ, 2 NUMERIC, 1 GAP.

## D. Falsification

2363 points of the seven-contact stratum (seed 31; W solved from Im H9 = Im H17 = 0 at random J, r, s, u). At every one, (1), (7) and (9) hold to relative 1e-6 (decomposition residual <= 4.3e-12). 33 points satisfy p, j, d > 0 and q, b >= 0. At all of them F6 >= 0 (worst normalized value 1.4e-20) and sign L = -sign r. None of the sampled supported points lies in the open branch r > 0, j > d. This is not evidence either way about that branch.

## E. Verdict

**gap-found** (one gap, substantive, believed true).

- Gap, p. 2, eqs. (3)-(5): the three real rows (matrix M3, right side b3) and the Cramer numerators are not printed. They depend on a representative CN/CD of the c amplitude that the note leaves to an absent checker. So (4), the factorization of N_R and (5) cannot be checked from the note or the supplement. Their consequences (7) and (9) hold numerically at 2363 stratum points, and the rest of the argument (wall m = 0, sign and square) is correct given them. Repair: print M3, b3 and R0, or publish the checker. The identities are then a routine exact replay.
- Scope: the theorem is a necessary condition on eight-distinct sources. It does not prove the support inequality on the c = 0 face, and the note says so.
