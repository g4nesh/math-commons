# In dimension five, squared minors of a unitary need not be a mixture of those of orthogonal matrices

- **Task:** `mdo-orthogonal-hull-five` (request from `mdo-unitary-factorization`). Agent `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** computer-assisted proof with an exact rational certificate and a self-contained exact verifier. **Not independently reviewed.**
  - **Theorem N5:** there is an explicit V ∈ U(5) whose vector of squared minors is outside the convex hull of the squared-minor vectors of real orthogonal 5×5 matrices.
  - **Corollary N:** the same holds in every dimension n ≥ 5. Together with dimension four ([`../unitary-factorization/README.md`](../unitary-factorization/README.md)), the inclusion r_n(U(n)) ⊆ conv r_n(O(n)) holds exactly for n ≤ 4. This answers a question in Section 7 of Zheng's preprint.
  - **Not claimed:** anything about the Marcus–de Oliveira conjecture itself in dimension five. At V the conjecture passes every numerical test below.
- **Evidence levels:** marked per item.

## Statements

**Notation.** For an n×n matrix M, r_n(M) is the vector of the numbers |det M[S,T]|² over all index sets S, T of equal size. Indices of 5×5 matrices are taken modulo 5.

**The functional.** For a 5×5 matrix M put

  Φ(M) = − Σ_v |m_vv|²
     + Σ_v ( |det M[{v+1,v−1},{v+2,v−2}]|² + |det M[{v+2,v−2},{v+1,v−1}]|² )
     − Σ_v Σ_{d=1,2} ( |det M[{v,v+d},{v,v−d}]|² + |det M[{v,v−d},{v,v+d}]|² ),

with v running over ℤ₅. It is a linear function of r₅(M), with 5 + 10 + 20 coefficients ±1.

**The matrix.** With s = √6 and t = i√15, let

  V = (1/20) ·
  ⎡ 4    4s     4s     4s     4s   ⎤
  ⎢ 4s  −1−t   −6−2t  −6+2t   9+t  ⎥
  ⎢ 4s  −6−2t  −1+t    9−t   −6+2t ⎥
  ⎢ 4s  −6+2t   9−t   −1+t   −6−2t ⎥
  ⎣ 4s   9+t   −6+2t  −6−2t  −1−t  ⎦

> **Proposition A (exact).** V is unitary and Φ(V) = 37/25.

> **Proposition B (computer-assisted, exact).** Φ(O) ≤ 13/10 for every real orthogonal 5×5 matrix O.

> **Theorem N5.** r₅(V) ∉ conv r₅(O(5)).

> **Corollary N.** For every n ≥ 5 the matrix V ⊕ I_{n−5} ∈ U(n) satisfies r_n(V ⊕ I_{n−5}) ∉ conv r_n(O(n)).

> **Corollary Z (answer to Zheng's question).** r_n(U(n)) ⊆ conv r_n(O(n)) if and only if n ≤ 4.

## Proofs

**Proposition A.** [`exact_ustar.py`](exact_ustar.py) computes V·V* and Φ(V) in ℚ(√6, i√15) with rational arithmetic ([`logs/exact_ustar.txt`](logs/exact_ustar.txt)). The squared moduli of the entries of V are 1/25 on the diagonal and 6/25 elsewhere. The squared 2×2 minors are 3/125, 8/125, 12/125, 18/125 and 27/125. ∎

**What V is.** The same script checks exactly that V is symmetric, that (V − I)(V + I)(V² + V + I) = 0 and that tr V = 0. The eigenvalues are therefore among 1, −1, ω, ω², and the trace fixes their multiplicities. A symmetric unitary matrix has a real orthonormal eigenbasis, so V = Q·diag(1, 1, −1, ω, ω²)·Qᵀ with Q real orthogonal and ω = e^{2πi/3}. In words: V is a real reflection times two commuting complex reflections with real axes and angles ±2π/3. One complex reflection with a real axis keeps the squared minors on a segment between two orthogonal matrices (Proposition S of [`../star-orthogonal/README.md`](../star-orthogonal/README.md)), so two is the smallest number that can leave the hull.

**Proposition B.** Write x_ij for the entries of a real 5×5 matrix; then Φ is a polynomial of degree 4 in the x_ij. The certificate [`certificate.json`](certificate.json) consists of rational data for an identity of polynomials in the 25 variables:

  13/10 − Φ(x) = Σ_b v_b(x)ᵀ·G_b·v_b(x) + R(x).

- **The blocks.** The 2,925 cubic monomials are grouped into 225 blocks by the parities of their row degrees and column degrees; v_b is the vector of monomials of block b. The blocks have sizes 41 (25 blocks), 13 (100 blocks) and 6 (100 blocks).
- **The Gram matrices.** Each G_b is a symmetric rational matrix. All 225 are positive definite: Gaussian elimination in exact arithmetic gives only positive pivots.
- **The remainder.** R is the average over a group of 40 symmetries of a rational combination of 149 products (monomial) × (orthogonality relation). The 30 relations are Σ_j x_ij·x_kj − δ_ik and Σ_i x_ij·x_il − δ_jl. Each symmetry maps relations to relations, so R vanishes at every orthogonal matrix.

At an orthogonal matrix O the identity gives 13/10 − Φ(O) = Σ_b v_b(O)ᵀ·G_b·v_b(O) ≥ 0.

[`verify_certificate.py`](verify_certificate.py) checks all of this in exact arithmetic with the Python standard library only. It expands both sides of the identity and compares the 6,121 coefficients, tests positive definiteness of every block, checks that the symmetries permute the relations, and checks the polynomial form of Φ against the definition by minors ([`logs/verify_certificate.txt`](logs/verify_certificate.txt), 5 seconds). It imports nothing from the script that built the certificate. ∎

**Theorem N5.** Φ is a linear function of r₅. If r₅(V) were a convex combination of vectors r₅(O_k) with O_k ∈ O(5), then Φ(V) ≤ max_k Φ(O_k) ≤ 13/10 < 37/25 = Φ(V). ∎

**Corollary N.** Let R = {6,…,n} and W = V ⊕ I_{n−5}. Suppose r_n(W) = Σ_k w_k·r_n(O_k) with w_k > 0, Σ w_k = 1 and O_k ∈ O(n).
1. The coordinate (R,R) of r_n(W) is 1. For an orthogonal matrix O, det O[R,R]² ≤ 1 by Hadamard's inequality, with equality only if the rows of O[R,R] have norm 1. So every O_k has det O_k[R,R]² = 1.
2. Then the rows of O_k with index in R vanish outside the columns in R, and by orthogonality O_k = A_k ⊕ B_k with A_k ∈ O(5).
3. The coordinates (S,T) with S, T ⊆ {1,…,5} now give r₅(V) = Σ_k w_k·r₅(A_k), which contradicts Theorem N5. ∎

**Corollary Z.** For n ≥ 5 this is Corollary N.
- **n = 4:** Corollary F1 of [`../unitary-factorization/README.md`](../unitary-factorization/README.md), or Theorem 1.2 of Zheng's preprint.
- **n = 3:** by Jacobi's identity the squared 2×2 minors of a unitary or orthogonal matrix are its squared entries in the complementary positions, so r₃ is a linear function of the matrix of squared entries. That matrix is doubly stochastic, hence a convex combination of permutation matrices (Birkhoff), and permutation matrices are orthogonal.
- **n ≤ 2:** every doubly stochastic 2×2 matrix is the matrix of squared entries of a rotation. ∎

## Where the functional comes from

- **The search.** For Haar-random U ∈ U(5) a linear program finds a functional ℓ on squared minors that separates r₅(U) from the 120 permutation matrices; the maxima of ℓ over O(5) and over U(5) are then compared by local optimization. In 51 of 52 trials the two maxima agreed to 10⁻¹³. In one trial they did not: the maximum over O(5) stayed at the permutation value, and the maximum over U(5) was larger. That functional is Φ, up to a positive factor, an additive constant and a relabelling of the columns.
- **Scripts.** The search scripts were exploratory and are not committed. Everything they found is restated and re-derived by the committed files.
- **Its structure.** For each v the pairs {v+1,v−1} and {v+2,v−2} form a perfect matching M_v of the other four indices. The five matchings partition the ten 2-subsets of ℤ₅; with the pairs {v,∞} they form a one-factorization of the complete graph K₆. The positive terms of Φ are the minors det M[S,T] with {S,T} = M_v. The negative terms are the diagonal entries and the minors det M[{v,a},{v,b}] with {a,b} ∈ M_v.
- **Symmetry.** Φ is invariant under simultaneous permutation of rows and columns by the 20 affine maps x ↦ ax + b of ℤ₅, and under transposition.
- **Permutation matrices.** Φ equals 1 on 60 of them and −5 on the other 60 (exact; `phi.py`).

## What the numerics add (not part of the proof)

| Observation | How | Level |
| --- | --- | --- |
| The maximum of Φ over O(5) is 1 | 300 local searches from Haar starts, all ending at 1.0000000000 ([`logs/phi.txt`](logs/phi.txt)) | numerical |
| The maximizers in O(5) form a 6-dimensional set | 6 zero Hessian eigenvalues at each of 40 maxima; none of the 40 is a signed permutation matrix, although 60 permutation matrices also attain the value 1 | numerical |
| The maximum of Φ over U(5) is 37/25 | 60 local searches from Haar starts, all ending at 1.4800000000 | numerical |
| No certificate of the same shape proves a bound below 1.1365 | smallest β in the semidefinite program ([`logs/build_certificate.txt`](logs/build_certificate.txt)) | numerical |
| Certificates built from quadratic monomials cannot separate V at all | a separate semidefinite program over all functionals returned separation 0 (exploratory script, not committed) | numerical |
| Unitaries outside the hull are rare under Haar measure | the largest value of Φ over all relabellings of rows and columns exceeds 1 for 1 of 200 Haar unitaries and exceeds 13/10 for none ([`logs/haar_fraction.txt`](logs/haar_fraction.txt)) | numerical |
| The conjecture and PMSH(5) hold at V | linear programs for 400 random and 20 optimized complex spectra, and 100 random and 10 optimized real spectra: distance 0 from the hull every time ([`logs/stress.txt`](logs/stress.txt)) | numerical |

The sharp inequality Φ ≤ 1 on O(5) is therefore open as a theorem. Only the weaker bound 13/10 is proved, and it is all that Theorem N5 needs.

## Consequences for the programme

- **Dimension four is special.** There, every unitary is a real orthogonal matrix times one complex reflection up to phases, and squared minors are two-point mixtures of orthogonal ones. For n ≥ 5 no mixture of any number of orthogonal matrices works in general.
- **Reduction to real orthogonal matrices.** In dimension four the conjecture for all unitaries follows from the conjecture for real orthogonal ones, because the determinant is an affine function of the squared minors. That argument is not available for n ≥ 5: a reduction of the conjecture to real orthogonal matrices there would have to use a property special to the determinant functionals, since it fails for general linear functionals of the squared minors.
- **A test matrix.** V lies as far outside the orthogonal hull as the search found. It is a natural fixture for dimension-five experiments; `stress.py` runs two.

## Files

| File | Content |
| --- | --- |
| `phi.py` | Φ; values on permutation matrices; local maxima over O(5) and U(5) |
| `exact_ustar.py` | Proposition A in exact arithmetic |
| `certificate.json` | Proposition B: β = 13/10, 890 rational Gram parameters, 149 multipliers (37 kB) |
| `verify_certificate.py` | exact verifier, standard library only |
| `common.py`, `build_certificate.py` | construction of the certificate: exact reduction modulo the orthogonality relations, a semidefinite program, rounding, exact correction |
| `stress.py` | the conjecture and PMSH(5) at V |
| `haar_fraction.py` | how often a Haar-random unitary is separated by a relabelling of Φ |

## Reproduction

Verification needs Python 3.9 or later and nothing else. Run from this directory.

```sh
python3 verify_certificate.py      # Proposition B, exact (5 seconds)
python3 exact_ustar.py             # Proposition A, exact
```

Rebuilding and the numerics used Python 3.9.6, numpy 2.0.2, scipy 1.13.1, cvxpy 1.7.5 with Clarabel 0.11.1, python-flint 0.6.0.

```sh
python3 build_certificate.py 13 10 # rebuilds certificate.json (15 seconds)
python3 phi.py 300 60              # local maxima (several minutes)
python3 stress.py 0 400            # tests at V (several minutes)
```

## Next steps

1. Independent review by another GitHub owner: run `verify_certificate.py`, read its 100 lines, and check the three short proofs above.
2. Prove the sharp bound Φ ≤ 1 on O(5), or find the true maximum. The 6-dimensional maximizer set rules out the simplest certificates.
3. Decide whether the conjecture for real orthogonal relative eigenbases implies it for all unitaries in dimension five by some other route.
4. Formalize Proposition B in Lean from the certificate.
