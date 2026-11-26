Evidence level: EXACT for the cubic, Q, N, R, the z-independence of H_I, the discriminant factorization Delta = F+F- and all three chart coefficient tables (own code from the paper's cubic C and generators). EXACT-pointwise (exact rational arithmetic at random points, a probabilistic identity test) for det H = -4G^2/Q^2, det H_I = -dA^2 P, the identities (6) to (8) and the bordered-minor divisibility. PROOF-READ for the real-branch argument. CITED: Prop. C.6 and Cor. C.5 (proved in foundations).

# blind_moc4_inventory52_plus17_low_products_20261023 (3 pp.)

Replay: `replay/blind_moc4_inventory52_plus17_low_products_20261023__check.py` (needs MOC4_SUPP). The run took about 3 min. Output: `..._check_output.txt`, all checks PASS.

## Deviations and limits

- G (687 terms), P9 and the two octics were read as DATA from the supplement JSON. I did not derive them independently. They were then tested against my own objects:
  - G: det H(c + z* kappa) Q^2 = -4 G^2 at 12 random rational points.
  - P9: det H_I = -d A^2 P9 at 6 points.
  - Octics: Delta = F+ F- as an exact polynomial identity, and each octic is irreducible.
- A pointwise identity test is not a symbolic proof. The polynomials have degree <= 30 and the random rationals have numerators up to 60, so Schwartz-Zippel makes a false pass very unlikely, but it is a probabilistic check.
- (7) was checked as "S± / (-d F) is a rational square" at 8 points. (8) was checked as proportionality of the B-quadratic to G at 8 points. Neither was checked symbolically.

## A. Statement

**Theorem 1 (p. 1).** Let x ∈ X have c(x) in the seven-contact face T7 = {0,1,4,6,9,12,17} (Table 5, aux set 19), with all seven canonical coefficients strictly positive. Then max(ab, bc, ca) > d^2. This is a partial exclusion, not an MOC statement. Spectra are arbitrary: the argument uses only Re x and the cubic, so degenerate spectra are covered.

## B. Dependencies

- Radial necessity: C(u) = 0, and det Hess C(u) = 0 at u = Re x (foundation (5)).
- Prop. C.6: all 12x12 principal minors of H(Re x) are >= 0 for x ∈ X. This is proved in foundations/appendix-c-part1.md, modulo standard GIT surjectivity. The statement used matches.
- The S3 stabilizer of T7, EXACT in the whole_assembly replay.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | C(c + z kappa) = -N - R B + Q z, with Q = -(AB+AC+BC+d(A+B+C)+d^2), z-free part as printed, linear in B, and N, R with positive coefficients | EXACT |
| 2 | forced height z = (N + R B)/Q | EXACT (Q < 0) |
| 3 | det H(B, z*) = -4G^2/Q^2; G has 687 terms, degree 9, degree 2 in B | EXACT-pointwise (12 points) |
| 4 | H with rows and columns 0, 5 deleted is independent of z | EXACT |
| 5 | det H_I = -d A^2 P, with P of degree 9 and degree 2 in B | EXACT-pointwise (6 points) |
| 6 | Delta = F+ F-, octics of degree 8 with 136 terms each, irreducible | EXACT |
| 6' | the norm is a square W^2 | EXACT-pointwise |
| 7 | (7) S± = -d F_eps H±^2 | EXACT-pointwise (square class, 8 points) |
| 8 | (8) P S± = (U ± W)^2 mod G | EXACT-pointwise (8 points) |
| 9 | bordered 13x13 minors on I+{0} and I+{5} vanish on G = 0 | EXACT-pointwise (4 parameter points, B symbolic) |
| 10 | chart X <= Y <= Z <= 1 covers the closed chamber a <= b <= c, bc <= 1 | PROOF-READ (D = 1, x = X, y = Y - X, z = Z - Y, W0 = 1 - Z) |
| 11 | tables gamma (9772 terms, all +, 512 x^21), R beta - N alpha (13819, all +, 256 x^21), R F* - N alpha (14549, all -, -208 x^21 A and -48 x^21 C), with F* the second octic (index 1). The other octic gives mixed signs | EXACT |
| 12 | a root B >= 0 forces alpha < 0, a unique simple root, F* < 0, hence both octics < 0 | PROOF-READ |
| 13 | where W H+ H- ≠ 0: S± > 0, so P > 0, so det H_I < 0, which contradicts Prop. C.6 directly | PROOF-READ |
| 14 | at the exceptional zeros of W H+ H-: closure of the analytic root graph; bordered minors give corank 2 on the open set, so every 12-minor is <= 0 there and at limits. The second B-derivative of det is -8 G_B^2/Q^2 ≠ 0, so corank <= 2 and some 12-minor is < 0, which contradicts Prop. C.6 | PROOF-READ. The Plücker sign claim is a standard fact for symmetric matrices of corank 2 |

## D. Falsification

See the whole_assembly report: a numeric sample on the T7 stratum, testing max(ab, bc, ca) > d^2.

## E. Verdict

**confirmed-modulo-citations.** Cited: Prop. C.6 / Cor. C.5. G, P9 and the octics were read as data and validated pointwise, not derived symbolically.
