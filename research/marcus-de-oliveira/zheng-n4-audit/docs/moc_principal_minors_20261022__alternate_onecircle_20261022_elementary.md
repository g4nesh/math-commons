# Audit: alternate_onecircle_20261022_elementary (copy in `moc_principal_minors_20261022`)

Evidence level: EXACT for every identity in the note (replayed with own sympy code); PROOF-READ for the inequalities built on them and the case split; NUMERIC falsification passed. No CITED step.

Document: `proofs/moc_principal_minors_20261022/alternate_onecircle_20261022_elementary.pdf`, 3 pages, "An elementary compensation inequality for two rotations", 22 Oct 2026.
Replay: `replay/moc_principal_minors_20261022__alternate_onecircle_20261022_elementary__check.py` (it runs the checks of the sibling script `replay/__workspace____alternate_onecircle_20261022_elementary__check.py` unchanged), output in `replay/moc_principal_minors_20261022__alternate_onecircle_20261022_elementary__check_output.txt` (28 EXACT checks and 5 NUMERIC checks, all pass).

## Deviations and limits

- This note is a standalone lemma about SO(3) and nonnegative reals. It does not mention spectra, contact sets or the Leibniz sum. Whether and how it feeds into the paper's Appendix C (one spectrum on a generalized circle) or the "principal minors" documents that cite it (as Theorem 7 "Two rotations") was **not** audited here. The verdict covers only the statements in this note.
- The supplement has three copies of this note: `__workspace__/...`, `principal_minor_result_20261022/...` and `moc_principal_minors_20261022/...`. They differ only in one sentence of the introduction (the two copies outside `__workspace__` add a sentence saying that applications to projection or determinant problems need separate arguments). I compared them line by line on the text extracts and checked the page images of the formulas. The mathematics is identical.
- The standard fact that every rotation in SO(3) is the image of a unit quaternion under v -> q v q^* is used as classical (PROOF-READ). It is also tested numerically.

## A. Statements

- **Lemma 1 (p.1).** For real u, v, x, y >= 0: sqrt((u+x)(v+y)) - sqrt(uv) >= sqrt(xy).
- **Theorem 1 (p.1).** Let a, b, c, d, A, B, C, D be nonnegative reals with a+b+c+d = A+B+C+D. Put s = b+c-B-C, t = a+b-A-B and R = ab + cd + AB + CD - 2 sqrt(bcAD) - 2 sqrt(adBC). Then R + s^2 + |st| >= 0. Nothing is excluded: zeros and equalities are allowed.
- **Theorem 2 (p.2).** For all P, Q in SO(3), with d0 = P11 - Q11, e0 = P33 - Q33 and G = 2 d0^2 + P12^2 + Q12^2 + P21^2 + Q21^2 + 2(P32 Q32 - P23 Q23), we have G + 2|d0 e0| >= 0.

None of these is a contact-case statement, so no contact set applies.

## B. Dependencies

- Classical: the surjection from unit quaternions onto SO(3) with the displayed convention. I replayed exactly that the displayed entries equal the matrix of v -> q v q^* (in homogeneous form), which is orthogonal with det 1 on |q| = 1. Surjectivity is standard. It was also tested numerically by reconstructing a quaternion for 2000 random rotations (max error 1e-15).
- No other document, case or external theorem is used.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Lemma 1: both sides nonnegative; difference of squares is (sqrt(uy) - sqrt(vx))^2 | identity | EXACT | polynomial identity in the square-root variables |
| 2 | (a,b,c,d) -> (b,a,d,c) on both quadruples preserves R and t and sends s to -s (using equal sums) | identity | EXACT | sympy; s' + s equals the total difference |
| 3 | (a,b,c,d) -> (d,c,b,a) preserves R and s and sends t to -t | identity | EXACT | same |
| 4 | therefore one may assume s, t >= 0 | logic | PROOF-READ | R + s^2 + \|st\| is invariant, and the two maps change the signs of s and t independently |
| 5 | (3): F := R + s^2 + st = \|h-k\|^2 + 2mn with h = (sqrt(bA), sqrt(aB)), k = (sqrt(cD), sqrt(dC)), m = b-B, n = D-d | identity mod total-sum relation | EXACT | the remainder after reducing D through the sum relation is 0 |
| 6 | (4): \|h-k\|^2 = \|h-Jk\|^2 - 2(h1-h2)(k1-k2) | identity | EXACT | sympy |
| 7 | if mn >= 0 then F >= 0 by (3) | inequality | PROOF-READ | immediate |
| 8 | mn < 0 means b-B and d-D have the same strict sign, so there are two cases and they are exhaustive | case split | PROOF-READ | n = -(d-D) |
| 9 | case b<B, d<D: alpha = M+t, beta = M+s, N = M+s+t, alpha*beta = MN + st, mn = -MN | identities | EXACT | parametrized replay with d fixed by the sum relation |
| 10 | case b<B, d<D: h2-h1 >= sqrt(alpha M) and k1-k2 >= sqrt(beta N) by Lemma 1, so (h1-h2)(k1-k2) <= -sqrt(alpha beta M N) <= -MN | inequality | PROOF-READ | Lemma 1 applies with (u,x,v,y) = (A, alpha, b, M) and (C, beta, d, N), all >= 0 since s,t >= 0; alpha beta >= MN >= 0 |
| 11 | case b>B, d>D: alpha = M-t = N+s, beta = M-s = N+t, M = N+s+t, alpha*beta = MN + st, mn = -MN | identities | EXACT | parametrized replay |
| 12 | case b>B, d>D: the Lemma 1 bounds in the other direction | inequality | PROOF-READ | (u,x,v,y) = (B, M, a, alpha) and (D, N, c, beta); alpha = N+s > 0 and beta = N+t > 0 |
| 13 | (5) holds in both cases, so F >= 0 by (4) | inequality | PROOF-READ | |
| 14 | the doc's quaternion formulas equal the rotation v -> q v q^* | identity | EXACT | sympy, homogeneous form |
| 15 | every rotation, including rotations through pi, has such a representation | classical | PROOF-READ (+NUMERIC) | standard double cover SU(2) -> SO(3) |
| 16 | d0 = -2s, e0 = -2t with (a,b,c,d) = (x^2,y^2,z^2,w^2), (A,B,C,D) = (X^2,...,W^2) | identity | EXACT | holds without using the norm relation |
| 17 | G/8 = s^2 + ab + cd + AB + CD + 2(yzWX + wxYZ) | identity | EXACT | sympy |
| 18 | \|yzWX\| = sqrt(bcAD), \|wxYZ\| = sqrt(adBC), so (G + 2\|d0 e0\|)/8 >= R + s^2 + \|st\| | inequality | PROOF-READ (squares EXACT) | both quadruples sum to 1 = \|q\|^2, so Theorem 1 applies |

## D. Falsification

- Theorem 1: 200000 random points on the product of two simplices (Dirichlet with concentration 0.2/1/5, seed 20261125), minimum 3.6e-6. Nelder-Mead from 300 starts, minimum -4e-17, which is 0 up to rounding. NUMERIC.
- Theorem 2: 200000 random pairs in SO(3)^2, minimum 2.1e-3. Nelder-Mead from 300 starts, minimum -7e-16 (equality at P = Q). NUMERIC.
- Information only: G by itself reaches -0.0696, so the 2|d0 e0| term is needed. The statement is not vacuous.

## E. Verdict

**confirmed.** Every identity was replayed exactly. Every inequality step is an explicit elementary argument that I read line by line, and the case split is exhaustive. The only outside fact is the classical quaternion parametrization of SO(3). No gaps.

Not covered by this verdict: how Theorem 2 is used downstream (the "principal minors" documents and the paper's Appendix C).
