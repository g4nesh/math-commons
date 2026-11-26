Evidence level: EXACT for every displayed identity (G9 itself, (1)-(5), the principal-minor norm identities, both gcds), all rebuilt from the paper's definitions by own code; PROOF-READ for the rank/density argument; CITED for radial necessity, Corollary C.5, Proposition C.6 and (C.2) (HMSV); NUMERIC only for a small stratum sample.

# final_referee_T46_negative_j_Heron_uniform_20261023 ("Excluding a negative Heron factor on a six-contact cone", 3 pp.)

Replay: `replay/final_referee_T46_negative_j_Heron_uniform_20261023__check.py` (output `..._check_output.txt`: 19 EXACT checks pass, 2 NUMERIC lines).

## Deviations and limits

- I did not use the document's verifier or JSON. G9 was rebuilt from scratch. I took u = q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4 + h kappa in the paper's m-coordinates and the cubic C of (A.2), computed det Hess C(u) by fraction-free elimination over Z[q,...,a,h], and eliminated h through C(u) = C0 + hQ. The paper's coordinates turn out to be the document's "original" coordinates. C(u) equals C0 + hQ exactly, with the printed signs. The document's principal minor "deleting coordinates 0, 1" is height-free in the paper's coordinates, and all of its identities hold there verbatim.
- The sign of G is a convention (only G^2 is determined). No conclusion depends on it.
- "Height-free" and "original" are the document's words. Nothing in my replay depends on the supplement's coordinate dictionaries.
- Numeric sampling of realizable points with d > 0 on this stratum is very sparse: the lib sampler mostly returns costs whose f28 coordinate is <= 0, which are sums of minors (see B). So the NUMERIC test is weak.

## A. Statements

Contact set T = {1,3,5,6,7,8}. This is paper T_47 itself (the identity element of G maps one to the other; EXACT). The ten generators are {f0,...,f6, f27, f28, f59}. The document works with the eight canonical rays (q,p,j,c,e,d,b,a) = (f27,f1,f59,f5,f2,f28,f6,f4), of which only f28 = h_{3421,-} is non-minor. The other two generators are f0 = f27+f1+f5+f2-f28 and f3 = f27+f5+f4-f28 (EXACT). So every cost in N_T is either a sum of minors or has nonnegative canonical coordinates with d > 0 (absorb min(y28, y0+y3) of f28 into f0, f3). The canonical form therefore loses nothing.

**Theorem 1.** Let x in X have c(x) = g = q f27 + ... + a f4, with all eight coefficients >= 0 and b, j, d > 0. Then H6(g) >= 0, where D = a+d, M = (c+d)(a+d+e), L = M+bD, U = c+e, and H6 = [L(p+q) - bDU]^2 - 4pqML + 4cebDM.

The statement covers every spectrum, including repeated eigenvalues: it uses only radial necessity and C.5/C.6, which hold for all x in X. It is a necessary condition, not a support inequality.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| Radial necessity: C(Re x) = 0, det Hess C(Re x) = 0 for x in X (paper p. 27; audit item (5)) | core | CITED (established) |
| Corollary C.5 (rank H(c) <= 12 on the secant cone) and Proposition C.6 (12x12 principal minors >= 0 at Re x; one > 0 if rank 12) | core | CITED (foundations/appendix-c-part1, proved there modulo GIT) |
| "Every radial point with F = det HF = 0 has rank <= 12", used at NON-realizable approximating points | core | CITED. It follows from paper (C.2) (det H restricted to {C = 0} is gamma0 s^2, and the zero locus of s on {C = 0} is Sigma_2; from HMSV) together with Corollary C.5. The document gives no reference; I resolved it to these. |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | C(u) = C0 + hQ with the printed C0, Q | identity | EXACT | flint |
| 2 | C0 < 0, so Q != 0 at realizable points; h = -C0/Q | inequality | EXACT/PROOF-READ | L >= bd > 0 |
| 3 | det Hess C(w(g,-C0/Q)) = -4 G^2/Q^2; G has 1556 terms, degree 9, quadratic in j | identity | EXACT | Bareiss; R/Q is -4 times a perfect square |
| 4 | (1) B^2 - 4AC = (q+d)^2 H6 H8, with H6 (150 terms, irreducible) and H8 (279 terms) | identity | EXACT | |
| 5 | (2) H8 = W0^2 + 4dDb^2 L Xi | identity | EXACT | |
| 6 | smooth H8 wall: (H8)_p = 2W0(q+d)L; at Xi = 0, (H8)_q = 4dDb^2L(c+e+d+b) > 0 | identity + sign | EXACT | |
| 7 | the principal 12-minor deleting coordinates 0, 1 is height-free and equals -j^2 P10, with P10 quadratic in j | identity | EXACT | |
| 8 | u0^2 - u1(B gamma - beta C) = W^2; S+ = -H8 R8^2; S- = -(q+d)^2 H6 T8^2; W, R8, T8 free of j | identity | EXACT | W found as an exact square root |
| 9 | P10 S+- == (u1 j + u0 +- W)^2 mod G | identity | EXACT | exact division by G |
| 10 | (4) gcd(A,B,C) = 1, gcd(G, G_j) = 1 | finite | EXACT | |
| 11 | (5) G_j^2 - 4AG = (q+d)^2 H6 H8 | identity | EXACT | |
| 12 | At G = 0, H6 < 0, H8 < 0, W R8 T8 != 0: P10 > 0, so the minor -j^2 P10 < 0 | inequality | PROOF-READ | from 8 and 9 |
| 13 | At such points the rank is exactly 12, so all 12x12 principal minors have the sign of that minor (<= 0) | analytic | PROOF-READ + CITED | rank >= 12 from the nonzero minor; rank <= 12 from (C.2) + C.5; for symmetric rank r, det(A_S) = det(V_S)^2 det D |
| 14 | H8 > 0 is impossible at G = 0 | case split | PROOF-READ | (5) at G = 0 gives G_j^2 = (q+d)^2 H6 H8. The document leaves this case implicit. |
| 15 | H8 = 0 case: A != 0 and grad G != 0 | analytic | PROOF-READ | checked line by line against 6 and 11 |
| 16 | density: G_j != 0, and W R8 T8 != 0, are dense in the local real hypersurface G = 0 | analytic | PROOF-READ | coprimality (4). A j-free factor dividing G would divide A, B and C. |
| 17 | limit: all 12-minors <= 0 at the original point; rank exactly 12 there (second derivative -8G_x^2/Q^2); one is negative, contradicting C.6 | analytic | PROOF-READ + CITED | |

## D. Falsification

NUMERIC. I sampled the stratum with the lib sampler, seeds 0-39 at margin 0 and seeds 100-139 at margin 0.02, giving 77 points. Only one had d > 0, and there H6 = 6.5e17 > 0. G was approximately 0 at the margin-0.02 points (relative 6e-14). At margin 0 some points are numerically degenerate (coordinates up to 1e15), with relative |G| up to 5e-2. I read this as conditioning, not a violation. I made no attempt to find H6 < 0 on the d > 0 part beyond this. Note that H6 < 0 does occur at realizable points with d < 0 (seed 114: H6 = -0.40, d = -0.04). There the cost is a minor sum, outside the theorem, so the hypothesis d > 0 matters.

## E. Verdict

**confirmed-modulo-citations.** Every algebraic identity replays exactly from the definitions, and the rank/density argument is correct. Citations: radial necessity, Cor. C.5, Prop. C.6, and (C.2) of the paper (HMSV) for rank <= 12 at non-realizable points of {C = 0, det H = 0}. Presentation remarks:

- The global input is called "separately established global secant facts" without a reference.
- The case H8 > 0 is excluded only implicitly.
