# A positive common permutation law near the five-dimensional Fourier basis

Let `F[r,j] = ζ^(rj)/√5`, where ζ = exp(2πi/5) and both indices run from 0 to 4. Let U be any unitary matrix with max entrywise distance from F at most **1/128**. There is an explicit law on all 120 permutations, depending only on U, with every mass **greater than 1/6144**. For every pair of equally sized index sets I,J it satisfies

    Pr[p(J) = I] = |det U[I,J]|².

In particular, for this entire family of relative eigenbases and any two complex spectra a,b,

    det(diag(a) + U diag(b) U*) = Σ_p w_p ∏_j (a[p(j)] + b[j]).

Both spectra can vary simultaneously without changing the weights. The statement controls the relative basis in dimension five. The unrestricted conjecture also permits every other relative basis and every dimension.

## The Fourier law

With probability q = 1/√5, take a uniformly random map j ↦ h j + c modulo 5, with h in {2,−2} and c arbitrary. With the remaining probability use a uniform random permutation. This gives each permutation at least (1−q)/120.

For a source pair {j,k} and destination pair {r,s}, the squared Fourier minor equals

    [2 − 2 cos(2π(r−s)(j−k)/5)] / 25.

Classify a pair by its difference up to sign, either ±1 or ±2. When the classes agree the product is ±1 and the minor is (1−q)/10. When the classes differ the product is ±2 and the minor is (1+q)/10. A uniform permutation hits every specified destination pair with probability 1/10. The affine maps always interchange the two classes. For either choice of h the five translations run once through the five opposite-class pairs, so the affine contribution is 0 for equal classes and 1/5 otherwise. The mixture gives both required values.

The affine maps and uniform permutations have uniform singleton marginals, equal to the squared entries of F. Complementation preserves the permutation event and, by Jacobi's complementary-minor identity for a unitary matrix, the squared minor. This covers sizes three and four. Sizes zero and five have probability one.

## Moving the basis

For any unitary V let

    s_p(V) = Re[conj(det V) sign(p) ∏_j V[p(j),j]].

These signed masses sum to one and reproduce every squared minor. Indeed, summing the determinant monomials over p(J)=I gives ε det V[I,J] det V[Iᶜ,Jᶜ], with ε the sign from moving I and J to the first positions in increasing order. Jacobi's identity gives det V[Iᶜ,Jᶜ] = ε det(V) conj(det V[I,J]). The product after multiplying by conj(det V) is exactly the desired real squared modulus. The individual signed masses need not be positive.

For U in the stated neighborhood use

    w_p(U) = w_p(F) + s_p(U) − s_p(F).

The moment equations and total mass follow from linearity. Put δ = max |U[r,j]−F[r,j]|. Every entry of either matrix has magnitude at most 1/√5 + 1/128 < 1/2. Telescoping a five-entry monomial bounds its change by 5δ/16. Telescoping the determinant one column at a time and using Hadamard's inequality bounds its change by 5√5 δ. The unchanged Fourier monomial has modulus 5^(−5/2). Since |det U|=1, each signed mass therefore changes by at most

    (5/16 + 5√5 · 5^(−5/2)) δ = 41δ/80.

The slack at the center is greater than 1/240. At the boundary each corrected mass thus exceeds 1/240 − 41/(80·128) = 1/6144.

## Why the minor identities give the determinant

Expand det(diag(a)+U diag(b) U*) by the positions I on which the second term supplies the principal minor. Cauchy–Binet expands that minor as the sum over J of |det U[I,J]|² times the product of the b's indexed by J. The coefficient of the product of the a's outside I and the b's in J is exactly Pr[p(J)=I] on the permutation side as well. This argument includes repeated eigenvalues and zeros.

Row and column permutations and unit-modulus diagonal rephasings transport the law by relabeling or leave its moments unchanged. The [one-coordinate gluing theorem](../single-cycle/README.md#3-one-coordinate-gluing-theorem) also allows these five-coordinate blocks in larger matrices. For instance, two certified blocks on coordinates 1–5 and 5–9 can be sampled independently and their permutations composed in matrix-product order. This extension only covers products satisfying the stated overlap condition.

## Reproduce the checks

Run `python3 research/marcus-de-oliveira/fourier-minor-laws/check.py` from the repository root. It computes all 252 squared minors and their separate permutation marginals at each of F and G F, where G rotates the first two coordinates with rational tangent half-angle 1/256. It also verifies unitarity, normalization, all 120 strict lower bounds, the perturbation radius, and the determinant identity directly at two spectra per matrix (one with repetitions). Every comparison is exact in the fifth cyclotomic field. The proof above gives the conclusion throughout the neighborhood.

The [optimizer report](../optimizer-search/README.md) includes Fourier bases among several search families. This certificate settles the common-law question specifically at F₅ and the neighborhood above. For background on obstructions to a common law in general, see [S. W. Drury, *On a theorem of Wielandt and the compounds of unitary matrices*](https://www.sciencedirect.com/science/article/pii/002437959490362X). Priority of this particular Fourier construction has not been established. The proof and arithmetic replay are available for independent review.
