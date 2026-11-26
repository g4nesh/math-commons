Evidence level: EXACT (eq. (2), eq. (3), the three 13-minor divisibilities, all gcds and the signed square class are replayed or proved by own flint/sympy code from the paper's definitions; only the polynomial G is read as data from the supplement, and its defining identity (2) is then checked exactly). The p.2 sign argument is PROOF-READ, with a replacement certificate supplied by the auditor.

# Second-pass gap check: blind_moc4_T48_negative_Phi_uniform_20261023

Label C27-T49-gapcheck. Document: `proofs/manin_transfer_20261021/blind_moc4_T48_negative_Phi_uniform_20261023.pdf` (3 pp.).
Replay: `replay2/blind_moc4_T48_negative_Phi_uniform_20261023__gapcheck.py` (parts A, B, C), output in `replay2/blind_moc4_T48_negative_Phi_uniform_20261023__gapcheck_output.txt`.

## Deviations and limits

- The first auditor's files were not opened. Nothing from the supplement was executed or imported.
- One supplement JSON was read as data: `blind_moc4_T48_native_full_radial_20261023.json`, for G, C and Q. C and Q are compared with my own rebuild and agree exactly. G is accepted only because the exact identity (2), det H = -4G^2/Q^2, is replayed with it.
- **I did not reconstruct the note's polynomials R7, R~7 and L±.** They are not printed, and the note's verifier and its companion files (`..._negative_Phi_exact_...json`, `..._native_j_squareclass_...json`, `..._native_rank12_borders_...json`) are absent from the supplement. In their place I prove the property the argument actually uses: P·H6 (equivalently P·Φ) is a nonzero square in the function field of G. The certificate is a different one, and it is complete (see step 5). The literal congruences P S± ≡ L±^2 (mod G), with the note's specific R7 and R~7, remain unverified as printed objects.
- Run times in the output file (part B 2825 s, part C 691 s wall clock) were inflated by machine contention. Unloaded, the same computations took about 250 s and 230 s. All three parts exited 0.
- The rest of the note (smoothness, eq. (4), the second-derivative completion, the global spectrally realizable Hessian sign theorem it cites) is outside this gap and was not examined.

## The reported gap

p.2, "Height and rank charts" and "Signed quadratic classes". The claims are: three 13-minors (deleting (0,0), (5,5), (0,5)) have forced numerators divisible by G; P S± = L±^2 (mod G) with S+ = Y^2 p^2 (q+d)^2 H6 R7^2 and S- = Y^2 p^2 (q+d)^2 Φ R~7^2; and R7, R~7 and H6 are coprime to G. The first auditor said the polynomials are not printed, are absent from the supplement, and only G is supplied.

## Classification: FILLED (severity: presentation)

The first auditor's facts are right: only G (with C and Q) is in the supplement. But every fact the sign argument needs can be rebuilt from the paper's definitions and checked exactly, and I did so. Notation: Y = b+c+d, G = A j^2 + B j + C0, and P is the 12-principal of H = Hess C(w) on the indices other than {0, 5}.

### Steps (all EXACT unless stated)

1. **Inputs from the definitions.** The source is m = q f7 + p f5 + j f29 + c f39 + e f2 + d f40 + b f41 + a f0, with the paper's f_j from Appendix D (shared lib). The two profile identities hold, and the generator set of T = {0,3,5,6,9,12} is {0,1,2,5,7,9,29,39,40,41}. With w = m + h·kappa and the paper's cubic (A.2), F(w) = C + hQ is linear in h, with Q = (c+d)(d+e+j) - p(e+j). C and Q equal the supplement's arrays term by term, so the note's cubic normalization is the paper's (A.2) in the paper's coordinates.

2. **Eq. (2).** I computed det H(w) symbolically, by fraction-free elimination over Q[q,...,a,h], and substituted h = -C/Q. The cleared numerator equals -4 G^2 Q^{n-2} exactly, which is det H = -4G^2/Q^2 with the supplied G. G has 1772 terms and degree 9, and flint factorization shows it is irreducible over Q.

3. **Rank charts.**
   - P is height-free. Hess C(w) = Hess C(m) + h·H_kappa, and H_kappa = e0e6^t + e6e0^t - e1e5^t - e5e1^t lives only in rows and columns 0, 1, 5, 6, so deleting {0, 5} removes it. P has 2784 terms and has degree 2 in j.
   - The forced numerators of the three 13-minors, deleting (0,0), (5,5) and row 0 with column 5, are exactly divisible by G. The quotients have 872, 421 and 613 terms, matching the note.
   - Let S be the Schur complement of the block H[I0,I0] (I0 = the complement of {0,5}). Then P·S_00, P·S_55 and P·S_05 are, up to sign, the 13-minors deleting (5,5), (0,0) and (0,5) (or (5,0), which is the same by symmetry). So S = 0 wherever G = 0 and P Q != 0, and there rank H = 12 (PROOF-READ).
   - gcd(G,P) = gcd(G,G_j) = 1. The content of G as a polynomial in j is gcd(A, B, C0) = 1.

4. **Eq. (3).** B^2 - 4AC0 is exactly divisible by Y^2 (q+d)^2 Φ, with Φ as printed on p.1. The quotient H6 is a j-free polynomial with 159 terms, and gcd(G, H6) = gcd(G, Φ) = 1. Since H6 is j-free and G is irreducible with no j-free factor, the coprimality of H6 with G also follows from step 3.

5. **Signed square class (the disputed congruences).** Let K = Q(q,p,c,e,d,b,a)[j]/(G). This is a field because G is irreducible of degree 2 in j. Also sqrt(D) = 2Aj + B lies in K, where D = B^2 - 4AC0. Since deg_j P = 2, A·P ≡ r1 j + r0 (mod G), with r1 = A P1 - B P2 and r0 = A P0 - C0 P2. So P = (u0 + r1 sqrt D)/(2A^2) in K, with u0 = 2A r0 - B r1. For an element beta of a quadratic extension with N(beta) = n^2, the identity (beta + n)^2 = beta (Tr beta + 2n) shows that beta is a square iff Tr beta + 2n lies in k^2 ∪ D k^2 for some choice of sign of n. I computed:
   - Res_j(G,P) = (A r0^2 - B r0 r1 + C0 r1^2)/A, with 727,848 terms, is a perfect square rho^2 (flint sqrt, exact). So N(P·H6) is a square.
   - H6 (u0 + 2A rho) is a perfect square, and H6 (u0 - 2A rho) is D times a perfect square.

   Hence **P·H6 is a nonzero square in K**. Since D = Y^2 (q+d)^2 Φ H6 is a square in K, **P·Φ is a square in K** as well. These are exactly the two signed classes S+ and S- of the note, up to the square factors Y^2 p^2 (q+d)^2 R^2, which do not change the square class. As cross-checks, the same property holds exactly in Q(sqrt D) at 120 random rational specializations. At 558 real roots of G (108 of them with Φ < 0; 40-digit arithmetic, seed 1), P·H6 >= 0 always holds and P < 0 at every Φ < 0 point (NUMERIC).

6. **From the square class to the sign (PROOF-READ, argument supplied by the auditor).** From step 5 there are polynomials s0, s1 and a j-free polynomial den != 0 with den^2 · P · H6 ≡ (s0 + s1 j)^2 (mod G) over Q(q,...,a). Clearing the j-free denominators gives e·(den^2 P H6 - (s0 + s1 j)^2) = G·F, with e j-free and nonzero. At a real point of G = 0 with e·den != 0, this gives P·H6 >= 0.

   Near a smooth real point of G = 0, the real zero set is a real-analytic hypersurface. Suppose e·den vanished on an open piece of it. Then e·den would vanish on the geometric component through that point, and the point lies on only one component because it is smooth. Since e·den is defined over Q, it would then vanish on all the Galois conjugate components, that is on all of {G = 0}, so G would divide e·den. That is impossible, because e·den is j-free and G has j-degree 2 and content 1. So e·den != 0 on a dense subset, and P·H6 >= 0 holds at every real point near a smooth real point, by continuity.

   At a real point with G_j != 0, the discriminant equals G_j^2 > 0 (or B^2 > 0 when A = 0). With Φ < 0, eq. (3) gives H6 < 0, so P <= 0, and P < 0 on the dense set where P != 0 (P is coprime to G). This is the conclusion the note draws from S± < 0.

So the step is justified. The proof uses only objects rebuilt from the paper's definitions plus the supplied G, which is pinned down by (2). The note's own R7, R~7 and L± are never needed: my certificate replaces them and has no extra denominators except j-free ones, which are automatically coprime to G. The defect is that the note relies on unprinted polynomials and a verifier missing from the supplement. That is a reproducibility and presentation defect. It is not a mathematical gap.

## Not covered

- The literal identities P S± = L±^2 (mod G) with the note's R7 and R~7 (354 and 456 terms): not reproduced, because the objects are not available. The equivalent square-class statement is proved instead.
- Everything in the note outside p.2's rank-chart and signed-class paragraphs.
