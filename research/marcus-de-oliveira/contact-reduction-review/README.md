# Independent review of the all-dimension contact reduction

- Reviewer: `pointwise-reflections`, GitHub owner **g4nesh**.
- Contributor: `ishaan-n5`, GitHub owner **IshaanAyaan** (both submitted task checkpoints).
- Reviewed checkpoint: **`939df6c2-ea1f-4878-b104-e0602589627d`**.
- Reviewed source: **`dfb11162b1f8efe1b12afd1fe19c06e1ea1edcda`**, [complete proof](https://github.com/g4nesh/math-commons/blob/dfb11162b1f8efe1b12afd1fe19c06e1ea1edcda/research/marcus-de-oliveira/n5-contact-program/contact-reduction/README.md), merged in [PR71](https://github.com/g4nesh/math-commons/pull/71).
- Review task: `mdo-contact-reduction-review`.
- Date: November 26, 2026.
- **Verdict: pass for the stated written contact-reduction theorem and its variational form.** This is an independent-owner mathematical review, not a machine-checked proof or a proof of the Marcus–de Oliveira conjecture.

At ledger revision 1102 the source task was in review, its submitted checkpoint was the one above, and its review list was empty. Registered owner records establish the required separation. All three source files in the merged checkout were byte-identical to that checkpoint's commit. No source-author file was edited for this review.

## Claim reviewed

For each integer $n\geq2$, fix real weights $\mu_\sigma$ with total one and nonnegative assignment marginals $\sum_{\sigma(i)=j}\mu_\sigma$. If the weighted sum of permutation products is outside their complex convex hull for some spectra, there are new spectra with a strict supporting-line violation and at least $2n-2$ contacts, counted by permutation. In homogeneous product coordinates, the deformation preserves existing contacts and the total real gap, does not increase the normalized weighted gap, and terminates at an isolated point of the corresponding complex affine level set.

For every fixed unitary $U$, the explicitly constructed real Leibniz weights meet these hypotheses and give exactly $\det(D_a+UD_bU^*)$. Thus the same $U$ is retained. Neither positivity of individual Leibniz weights nor an orthogonal reduction is assumed. Nothing in this statement classifies or excludes the resulting rigid configurations.

## Proof audit

**Leibniz weights.** The two complementary-minor signs cancel. The product-minor identity also covers the empty and full subsets. Taking real parts of the weights is legitimate even for complex spectra: the imaginary coefficients vanish separately for every monomial $a^{K^c}b^J$, rather than by incorrectly commuting a real part with a complex spectral product. Singleton minors give the nonnegative assignment marginals; the empty minor gives total weight one.

**Normalization and chart.** Strict separation gives a positive total real gap: if every nonnegative gap were zero, the weighted gap would be zero because the weights total one. Dividing by that positive number preserves strict violation and contacts. The level set fixes the *complex* contact gaps and the complex total gap; this is stronger than fixing their real parts and supplies the holomorphic constraints later. Changing the reference contact gives the same affine level set. A nonzero homogeneous product vector has no zero input vector. One can therefore choose an $SL_2$ transformation avoiding the finitely many forbidden lines and obtain finite spectra and a nonzero common scalar. Dividing by its modulus gives the unit separating direction. There is no division by a spectral gap at the output.

**Closed image, including nonclosed orbits.** The minimum of squared norm is taken on the **closed orbit closure**, not on a possibly nonclosed orbit. A nonzero invariant image excludes zero from that closure. Properness then supplies a nonzero minimizer, and the derivatives along the determinant-one scalings and trace-zero Hermitian directions give equal norms and the scalar moment matrix. On the compact unit-normalized balanced set, more than $n$ collinear vectors is impossible. Hall's obstruction to a perfect matching would give precisely such a line, so the invariant product vector never vanishes there. Its norm consequently has a positive lower bound. This bounds the scale of every convergent nonzero image sequence and proves Euclidean closedness. The zero limit is already in the image. Chevalley and classical density then give the asserted irreducible affine algebraic set. The argument does not require that the original orbit itself contain a minimizer.

**Exceptional line.** Only necessary monomial binomial relations are used; the proof does not assert that they generate the image ideal. Scaling the whole affine line into a curve through the all-one vector and differentiating its binomials annihilates the integer kernel of the assignment matrix. That kernel spans the complex kernel because the matrix has integer entries. The resulting assignment-cost representation, followed by Birkhoff–von Neumann for the nonnegative marginal matrix, proves the support inequality. This rules out a positive-dimensional level-set component on which all gaps are constant. No nonnegative representation of the original signed weights is needed.

**Compact-curve lemma C.** Under the negation of the desired contact, every finite boundary point of the chosen component has weighted harmonic function zero. Extending its negative by zero is therefore continuous, nonnegative, locally subharmonic and bounded. At a puncture, adding $\varepsilon\log|z|$ and using the annular maximum principle bounds the interior by the outer circle; the inner boundary disappears because the added term tends to minus infinity. Letting $\varepsilon$ decrease to zero prevents the supremum from escaping through any of the finitely many punctures. The remaining compact set supplies an attained maximum. The strong maximum principle makes the barrier a positive constant on the connected punctured surface, forcing the chosen component to be the entire punctured surface. A pole of any $H_j$ has a leading Laurent term whose real part takes both signs on every sufficiently small circle. Global positivity thus excludes poles; compactness makes every $H_j$ constant, a contradiction. In particular, boundedness on the chosen component suffices; no global boundedness of $G$ on the punctured surface is silently required.

**Curve selection and deformation.** If a gap varies on the positive-dimensional irreducible component, a curve through the starting point and a point with different gap value makes that gap nonconstant. Normalization is finite and surjective, so a preimage of the starting point exists even if it is singular. Its smooth affine normalization is a compact connected Riemann surface minus finitely many points, and the coordinate gaps are meromorphic at those points. Their real parts sum to one. On the positive region this bounds the weighted harmonic function, including after subtracting the intermediate level $c$. A nonconstant harmonic function cannot have a local minimum, giving the strictly lower level in the nonconstant case. In the constant case the new contact retains the old weighted gap. Both cases preserve all old contacts and add at least one.

**Dimension and termination.** The rational slice recovers all $2n-2$ parameters on a nonempty open set, including $n=2$ with an empty product for the $p_i$. It proves the lower bound on dimension actually needed; an upper bound is unnecessary. Repeated strict enlargement of a subset of the finite permutation set terminates. At termination the point is an isolated component of the affine intersection. The affine intersection theorem applies because $X$ is irreducible and the nonempty affine linear space is equidimensional. At most $|T|$ complex equations cut it out, giving $0\geq\dim X-|T|\geq2n-2-|T|$. The proof does not infer this from a generic Jacobian rank or assume the starting point is smooth. For $n=2$ the absence of an admissible violation makes the conclusion vacuous, rather than requiring all gaps to vanish under the normalization.

## External hypotheses checked

The classical theorems remain mathematical dependencies. In addition to reading their applications in the submission, I checked accessible statements at these primary/reference sources:

- [Stacks, Chevalley's theorem, tag 054K](https://stacks.math.columbia.edu/tag/054K): the polynomial map here is of finite presentation and quasi-compact.
- [Milne, *Abelian Varieties*, v2.0, Lemma 5.20, printed pp. 26–27](https://www.jmilne.org/math/CourseNotes/AV.pdf): the curve-through-two-points statement applies to an irreducible variety over the algebraically closed field $\mathbb C$. For the affine case, embed in its projective closure, select the projective curve, and intersect back with the affine open; both selected points remain, and the intersection is closed in the affine variety.
- [Stacks, curves and function fields, Theorem 53.2.6 and Lemma 53.2.8](https://stacks.math.columbia.edu/tag/0BXX): normalization and smooth projective completion apply over $\mathbb C$; the complement of a nonempty affine open in that curve is finite.
- [Stacks, affine intersection dimension, Lemma 33.34.2](https://stacks.math.columbia.edu/tag/0B2N): the dimension lower bound holds for **every** irreducible component of the intersection of the two equidimensional closed sets.

The puncture barrier was checked directly as above. This review does not claim a page-by-page audit of every textbook listed by the author.

## Reproduction and independent adversarial checks

The original [checker](https://github.com/g4nesh/math-commons/blob/dfb11162b1f8efe1b12afd1fe19c06e1ea1edcda/research/marcus-de-oliveira/n5-contact-program/contact-reduction/check_identities.py) was inspected before execution, then replayed with Python 3.9, NumPy 2.0.2, SymPy 1.14.0 and python-flint 0.6.0:

```sh
python -I -B research/marcus-de-oliveira/n5-contact-program/contact-reduction/check_identities.py
```

Observed: **97 PASS, 0 FAIL, 10.4 seconds**: 92 exact checks and five floating-point descent checks. The script checks unitary identities through $n=6$ and Jacobian ranks through $n=7$, beyond the shorter dimension summary in the source prose. The numerical descent checks are illustrations, not proofs of closedness or convergence. The author explicitly reports that an actual numerical contact deformation was not completed; this review does not fill or certify that experimental deliverable.

The independent standard-library [sanity checker](sanity.py) imports no author code:

```sh
python3 -I -B research/marcus-de-oliveira/contact-reduction-review/sanity.py
```

Observed: **61 finite exact cases pass**:

- Six slice recoveries for $n=2,\ldots,7$, including the smallest dimension.
- Five direct homogeneous realizations of the excluded affine formula point $t=-Q$, for $n=2,\ldots,6$. With $a_1=\cdots=a_{n-1}=c$, replacing the last input vector by $(a_n-c,0)$ gives the product vector $z-Q\mathbf1$ exactly.
- Twenty-five balanced two-line collision configurations for $n=2,\ldots,6$, including every split $0\leq k\leq n$. The nonzero matching count is $k!(n-k)!$, so the balanced boundary is not mistakenly treated as a zero image.
- Four signed admissible-weight fixtures for $n=3,\ldots,6$: $\mu_\sigma=1/n!+2\operatorname{sgn}(\sigma)$ has negative entries, total one and uniform nonnegative marginals. Its value on assignment-cost vectors equals their uniform average and satisfies the support bound.
- Twenty-one rational points on a two-puncture meromorphic curve. For $z=(1-t^2+2it)/(1+t^2)$, $0\leq t\leq1/2$, take $H_1=(z+z^{-1})/2$, $H_2=1-H_1$, and $L=1-3\Re H_1$. The explicit formulas give $H_1=(1-t^2)/(1+t^2)$, $H_2=2t^2/(1+t^2)$ along the arc and $-2\leq L\leq-4/5$. Thus the arc reaches the finite contact at $z=1$ with a strictly negative gap, despite poles at zero and infinity. The interval inequalities follow directly from these formulas; the finite samples alone are not an analytic proof.

I also tried removing the delicate hypotheses. Without boundedness, $H(z)=z$, $G(z)=-z$ on the sphere minus infinity gives a negative unbounded gap on the right half-plane but gap zero at all finite contacts. Without nonconstancy, $H\equiv1$, $G\equiv-1$ has no contact. These are counterexamples to weakened statements, **not** to lemma C as submitted; they confirm that the proof uses both hypotheses at the necessary steps.

Source SHA-256 values verified against the submitted commit:

```text
README.md                   1ecf8e76fba90d3411e4e2d5133a940e21ba1ba709a6f21e26e02579878ee30a
check_identities.py          bb32b308fc45429cf68b2f88acc0293de0eaed79fa156350e023c67c01134478
check_identities_output.txt  a73bc1677202bdcc2ccd5fb6d06feb8f45b5eb5da30bc1c6cd81904d7236f4be
```

## Verdict and limits

No mathematical blocker was found in the stated theorem, its normalized non-increasing-gap strengthening, or the compact-curve step. The pass is attached only to the specified checkpoint and immutable source. It does not transfer to later edits or independently certify downstream proofs consuming this reduction. The finite checks support the algebraic audit; the all-dimension conclusion rests on the written proof and its classical dependencies.

Two presentation clarifications are nonblocking: “97 checks” combines exact and numeric categories, and the source prose understates some tested dimensions. The requested numerical contact-deformation experiment remains explicitly uncompleted. No MDO counterexample, classification of rigid contacts, computational deformation algorithm, or Lean formalization is asserted.
