Evidence level: EXACT. Every disputed identity was replayed by the auditor's own exact code (python-flint), starting from the paper's definitions (f_j, the cubic C of (A.2), kappa). The polynomials S, T and Omega, which exist but are not printed, were taken from the supplement's JSON data as a candidate certificate and checked.

# Gap check: blind_moc4_T44_conic_branch_exclusion_20261022 (label C18-T43-T45-gapcheck)

PDF: proofs/manin_transfer_20261021/blind_moc4_T44_conic_branch_exclusion_20261022.pdf (3 pages).

## Deviations and limits

- Psi, S, T, Omega and U are not printed in the document. Psi, G and P are determined by the document's own definitions, and I recomputed all three. For S, T and Omega I used the polynomials `H8`, `K8` and `W17` in the supplement's data file `blind_moc4_T44_principal_squareclass_20261022.json`. The file was read with `json.load` and its integer term lists were converted by my own code. Nothing of the supplement was executed or imported. I read the companion `.py` only to learn which keys hold which polynomial. I took U := A P, which is admissible because the document asks only that U be some polynomial.
- I compared my own G, P, Phi and Psi against the supplement's copies, but these comparisons are not needed for the checks.
- Scope: only identity (2) and the X = 0 octic identity are examined. The geometric use of these identities on pp. 2-3 is not judged here: the corank/rank argument, the density and continuity extension, and the appeal to the "global Hessian inertia theorem". As a by-product, the other displayed identities of Lemma 2 and the conic determinant (3) also pass exactly.

## The gap as reported

The disputed identities are on p.1, display (2), and on p.2, the identity used in the X = 0 case:

  (2a)  B^2 - 4 A C0 = e^2 Phi Psi,
  (2b)  P (e^2 Phi S^2) = (U + Omega)^2  (mod G),
  (2c)  P (Psi T^2)     = (U - Omega)^2  (mod G),
  (X0)  W^2 Psi |_{a = pc/W} = c^2 (E1 e + E0)^2.

The first auditor's objection is that Psi, S, T, Omega and U are not printed and only a companion verifier is cited.

## Classification: FILLED (severity: presentation)

All four identities hold exactly, and the polynomials they need exist as described.

### Objects rebuilt from the definitions

- **Cost and cubic restriction.** Let g = q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15 and w = g + h kappa, with kappa = e_13. Then C(w) = C + hQ exactly, with C and Q free of h (so the h^2 and h^3 terms vanish). These C and Q coincide with the supplement's C_terms and Q_terms.
- **Hessian.** H(w) is the 14 x 14 Hessian of the cubic C of (A.2) at w.
- **P.** P := det H[I,I] / (q+d)^2 with I = {2,...,13}. The minor is free of h, the division is exact, and P has degree 10 and 1743 terms, as the document says.
- **G and equation (1).** det H(g + h kappa) has h-degree 3. Substituting h = -C/Q gives the exact identity

    Q^4 det H(w(g, -C/Q)) = -4 Q^2 G^2,

  where G is a polynomial of degree 9 with 1188 terms. This verifies equation (1) with G determined up to sign. All identities below are invariant under G -> -G (with U -> -U, which swaps (2b) and (2c) up to Omega -> -Omega), so I used the supplement's sign.
- **Coefficients of G.** G = A q^2 + B q + C0, and gcd(A, B, C0) = 1.
- **Phi.** Phi was computed from the document's displayed formula, Phi = R^2 - 4dj(b-a)W{c(d+a) - ja} with W = c + d + b.

### The disputed identities (all EXACT)

- **(2a)** B^2 - 4AC0 is exactly divisible by e^2 Phi. The quotient Psi is a q-free octic with 137 terms, as stated. Psi is therefore defined by (2a) itself.
- **(2b), (2c)** Take S = H8, T = K8 (q-free octics), Omega = W17 (q-free, degree 17), all nonzero, and U = A P. Both P e^2 Phi S^2 - (U + Omega)^2 and P Psi T^2 - (U - Omega)^2 are divisible by G in Q[q,p,j,c,e,d,b,a]. A single polynomial is a Groebner basis of the principal ideal it generates, so a zero remainder is a proof of divisibility.
- **(X0)** Psi has a-degree 2, and W^2 Psi(a = pc/W) = c^2 (E1 e + E0)^2, where Z = pj - pc - dW, E1 = (p-b)Z - jd(p+W) and E0 = p(d+b)Z.

### Side claims attached to (2)

- **"q-content of G is one".** This is gcd(A, B, C0) = 1, checked above.
- **S and Omega contain no hypersurface component of G.** Both are q-free, nonzero and independent of h. Every irreducible factor of G involves q, because the q-content is 1. Over C, the zero set of such a factor projects dominantly onto the q-free coordinates. A q-free polynomial vanishing on that zero set therefore vanishes identically, which S and Omega do not. (PROOF-READ, argument written here.) The real-hypersurface version used on p.3 is part of the geometric argument, which is outside this check.

### Also checked (not disputed, EXACT)

- D R = D^2 + dj(b-a)W{c(d+a) - ja} - (d+a)(j-c) Y.
- D^2 Phi(p = ajdW/D) = (D^2 - dj(b-a)W{...})^2.
- W^2 E1^2 G(a = pc/W, e = -E0/E1) = d(d+b) Z T3 (q L1 + L0)^2. G has a-degree 2, and the substitution in e has degree 2, so the stated clearing factor is right.
- L1 + W E1 = p{pjW - (d+b)(pc+dW)}.
- G = e^2 A2 + e(q+d) B1 + (q+d)^2 K, with A2 quadratic and B1 linear in q.
- The conic determinant (3): det H = 2 j X Y Phi^2.

### Why "presentation"

The identities are correct, and the objects in them exist with exactly the stated degrees and term counts. The document's defect is that S, T, Omega and U are existential and are given only by reference to a companion verifier, and Psi, G and P are given only by their defining properties. A reader cannot check (2b) and (2c) from the PDF alone. The repair is to publish S, T and Omega with U = A P, which this replay does from the data file, and to state that Psi := (B^2 - 4AC0)/(e^2 Phi). No mathematical step is missing.

## Files

- `replay2/blind_moc4_T44_conic_branch_exclusion_20261022__gapcheck.py`: exact replay, 26 checks. It takes about 100 s, most of it on the 14 x 14 symbolic determinant.
- `replay2/blind_moc4_T44_conic_branch_exclusion_20261022__gapcheck_output.txt`: its output, ending in ALL EXACT CHECKS PASS.
