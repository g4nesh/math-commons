Evidence level: EXACT for C, Q, the discriminant factorization, the e-decomposition of G, the conic determinant and the sign/positivity deductions; EXACT pointwise (25 random rational points) for det HF = -4 d^2 G^2/Q^2; PROOF-READ for the rank/Schur/continuity logic; NUMERIC for the key consequence "P < 0 on G = 0 when Phi < 0" and for the two cited global facts; UNCHECKED: the five "quadratic norm" identities (U, W0, H9, K7, W16); CITED: the secant set equality with its rank bound, and the global inertia theorem.

# Audit: root_moc4_T47_negative_Phi_20261022 (branch lemma for paper T_65)

PDF: proofs/manin_transfer_20261021/root_moc4_T47_negative_Phi_20261022.pdf, 3 pages. Replay: `replay/root_moc4_T47_negative_Phi_20261022__check.py`.

## Deviations and limits

- "Historical T47" is paper T_65 = {0,4,5,6,9,12} (EXACT: the eight roles are exactly generators(T_65), only f35 = h4312,- non-minor).
- G (965 terms) is the supplement's polynomial, loaded from `root_moc4_T47_full_radial_20261022.json` as data. I tied it to the paper's cubic by checking det HF(w(g,-C/Q)) Q^2 = -4 d^2 G^2 exactly at 25 random rational points with my own Hessian. All structural claims about G (degree, discriminant, decomposition, conic) were then proved symbolically from that G.
- I did NOT rebuild the order-12 principal minor P (3475 terms) symbolically, nor the auxiliary polynomials U, W0, H9, K7, W16, Psi-companions; so the five norm identities on p.1 are UNCHECKED. I tested their consequence numerically (P < 0 at real roots of G = 0 with Phi < 0).
- The two global facts the proof rests on (rank HF <= 12 on the secant variety via {F = det HF = 0} = Sec; nonnegativity of all order-12 principal minors at realizable responses) are not proved in this document; they are CITED. I tested both numerically on random realizable responses.

## A. Statements

**Theorem 1 (p.1).** Let g = q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7 with all coefficients >= 0, d > 0, e > 0 and Phi := [(p-c)j + (b-c)(p+d+a) - pa]^2 - 4ac(p-b)(p+d+a) < 0. Then there is no x in X with c(x) = g (i.e. g is not the gap vector of any spectral pair, at any height). Includes q = 0, j = 0, b = 0; nothing claimed for Phi >= 0 or e = 0. Being a non-existence statement on part of the T_65 cone, it is stronger than (*) there and covers all spectra (repeated and colliding included), since it only uses that w = Re(x) lies in Sec(X).

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| {F = det HF = 0} = Sec(M) and rank HF <= 12 on Sec(M) | global fact (another supplement document) | rank <= 12 at radial points and at nearby real points of G = 0 | CITED; NUMERIC: on 200 random realizable responses the 13th singular value of HF was below 1e-9 relative in every case |
| global inertia theorem: all order-12 principal minors of HF at a realizable real matching response are >= 0 | global fact | final contradiction | CITED; NUMERIC: no negative order-12 principal minor on the same 200 samples |
| radial necessity (det HF(Re x) = 0) | established (audit item 5) | G = 0 at realizable points | matches |
| supplement exact script for the five norm identities | supplement verifier | P <= 0 on G = 0 | not executed; UNCHECKED |

## C. Steps

| # | claim (p.) | kind | status | how |
|---|---|---|---|---|
| 1 | F(w(g,h)) = C + hQ with the displayed C, Q (p.1) | identity | EXACT | own cubic (A.2), F affine in h because kappa in X |
| 2 | C < 0 under the hypotheses; Q = 0 then makes C + hQ = 0 impossible | inequality | PROOF-READ | C = -d(... + d^2 + ...) with d > 0, others >= 0 |
| 3 | det HF(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (p.1) | identity | EXACT pointwise (25 points) | own Hessian, supplement G |
| 4 | P (delete coordinates 0,5) independent of height (p.1) | identity | EXACT | HF(kappa) is supported on rows/columns 0,1,5,6 (paper (3.10)) |
| 5 | G quadratic in q, q-content 1; B^2 - 4AC0 = e^2 Phi Psi, Psi 162 terms (p.1) | identity | EXACT | sympy division |
| 6 | U = AP - p2 G, u0^2 - u1(Bp0 - p1C0) = W0^2, 2Au0 - Bu1 -+ 2AW0 = e^2 Phi H9^2, e^2 Psi K7^2 (p.1) | identity | UNCHECKED | polynomials not rebuilt |
| 7 | consequence on G = 0: P e^2 Phi H9^2 = (U - W0)^2 | algebra | PROOF-READ given step 6 | quadratic norm calculation |
| 8 | G = e^2 A2(q) + e(q+d+b) B1(q) + (q+d+b)^2 K (p.2) | identity | EXACT | sympy |
| 9 | det(2H_conic) = -2acd^2 L^2 Phi^2 (p.2) | identity | EXACT | sympy |
| 10 | Phi < 0 forces a, c > 0, p > b; L > 0 (b = j = 0 gives a square) | inequality | PROOF-READ + EXACT | Phi = X^2 - 4ac(p-b)(p+d+a); specialization replayed |
| 11 | nonsingular conic, so (G_q, G_e) != 0 at every claimed point | analytic | PROOF-READ | (q,e) -> (u,v) is a local diffeomorphism for e > 0, q+d+b > 0 |
| 12 | d^2/dt^2 det = -8 d^2 G_t^2/Q^2 < 0 at G = 0, so corank <= 2; with step-B rank bound, rank = 12 | analytic | PROOF-READ (rank bound CITED) | corank k forces vanishing to order >= k |
| 13 | P <= 0 on the dense set A H9 W16 != 0 of the real radial hypersurface near the point, P = 0 excluded since W0 != 0 | inequality | rests on step 6 (UNCHECKED); NUMERIC support | see D |
| 14 | at rank 12 with an invertible 12-block, every 12-minor is det(R_S)^2 times that block's determinant | identity | PROOF-READ | M = R^t M11^{-1}... standard Schur argument, correct |
| 15 | continuity to the point; rank 12 gives a nonzero, hence negative, 12-minor; contradicts the inertia theorem | analytic | PROOF-READ (inertia CITED) | nearby real points of G = 0 are radial, hence in Sec by the cited set equality, so rank <= 12 there too |

## D. Falsification

1. P-sign test (NUMERIC, seed 7, mpmath 40 digits): 1500 random draws of (p,j,c,e,d,b,a) in (0.01,3)^7 with j = 0 or b = 0 in 20% of draws; 324 had Phi < 0; at the 132 real roots q of G = 0 (either sign of q), the height-free minor P was negative every time (largest value -0.0202). This is the consequence of the unreplayed norm identities, and it held.
2. Cited global facts (NUMERIC, seed 3): 200 random complex spectra and phases, w = Re(alpha m). In all 200, the 13th singular value of HF(w) was negligible (rank <= 12), and no order-12 principal minor was negative beyond relative 1e-8 (minimum relative value +2.7e-13).

## E. Verdict

**confirmed-modulo-citations.** The logical structure is correct and every identity I could rebuild holds exactly. The proof depends on (i) two cited global theorems (secant set equality with rank <= 12; inertia of order-12 principal minors), which survived my numerical tests, and (ii) five integer polynomial identities that only the supplement's script checks; I could not replay them in the budget (UNCHECKED), but their consequence passed a numerical test. No gap found.
