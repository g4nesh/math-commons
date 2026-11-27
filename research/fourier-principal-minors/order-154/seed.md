# Order 154 · Contributor comparison evidence

**Status.** This is a kickoff proof sketch and exact reference computation, prepared on 2026-11-27. No independent-owner verdict is claimed. Cold replay workers should first freeze their own implementation and result using the [independent brief](replay.md).

**Claim under review.** Every ordinary principal minor of the complex Fourier matrix F₁₅₄ is nonzero. Let d(S) = det(ζ₁₄^(ij)) for rows and columns indexed by S ⊆ {0,…,13}. Take 154 = 14 · 11; the factors are coprime and 154 is square-free. Theorem 1.4 reduces the requested conclusion to nonvanishing of the principal minors of F₁₄ in characteristic 11. In particular the all-prime-ideals criterion suffices. [Caragea–Lee–Malikiosis–Pfander, v3, Theorem 1.4 and §3.2](https://arxiv.org/html/2505.24326v3#S1) · [Wang–Zhang, v1, Proposition 3.1](https://arxiv.org/html/2608.22119v1#S3)

**Arithmetic.** The polynomial Φ₁₄(X) = X⁶ − X⁵ + X⁴ − X³ + X² − X + 1 reduces modulo 11 to the product of X³ + 6X² + 4X + 1 and X³ + 4X² + 6X + 1. They are distinct cubics and neither has a root in the eleven-element field, hence both are irreducible. The two quotients account for all prime ideals above 11 in ℤ[ζ₁₄]. In each quotient the residue of X has exact order 14, as checked by the powers 14, 7 and 2. Thus evaluating F₁₄ at X in both fields tests every required reduction.

**Enumeration.** At each index the reference checker branches on exclusion and inclusion. Exclusion deletes the corresponding row and column. Inclusion with nonzero pivot u replaces the remaining block by its Schur complement and multiplies the selected determinant by u. This determinant identity applies to every subsequent principal subset. A zero pivot would give a vanishing minor formed by the earlier included indices and the current one, for which the program raises an exception. Both decisions on each of the 14 indices give every subset exactly once. Arithmetic uses triples of coefficients in the two cubic quotients. Inverses are checked by multiplication.

From the repository root, run the following with standard-library Python 3.8 or later:

```sh
python3 research/fourier-principal-minors/order-154/verify_seed.py
```

**Observed result.** In each field the nonzero counts, in order from size 0 to size 14, were `1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1`. This is 16,384 per field and 32,768 checks across the two fields. The seed applies Theorem 1.4 to infer P(154). The program checks the finite hypothesis; validity of the theorem application belongs in the written review. Priority has a separate [literature assignment](prior-art.md).

**Provenance.** The workspace kickoff supplied an exact finite-field elimination pilot. This comparison implementation performs a Schur-complement enumeration using the same two reductions and runs without external data, randomness or third-party dependencies. Both implementations belong to the seed side of the review. [Reference code](verify_seed.py) · [Handoff and comparison](synthesis.md)
