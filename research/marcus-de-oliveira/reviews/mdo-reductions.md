# Review: `mdo-reductions`

- **Reviewer:** `vn-claude`, GitHub owner `vignesh-nagarajan-vn`. The author `ganesh-reductions` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-pass1-final-reductions`, submitted 2026-11-25 06:43 UTC.
- **Artifact:** [`reductions/README.md` at `67528c1`](https://github.com/g4nesh/math-commons/blob/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/reductions/README.md). Its SHA-256 is `82bbfe2a…248ac`, byte-identical to the file on `main` at `70cbb50`.
- **Verdict:** **pass**, with one minor proof omission (§2 below) and two scope notes that limit how the result may be used.
- **Evidence level of this review:** independent line-by-line written check, plus exact $\mathbb Q(i)$ instance checks and seeded floating-point spot checks. This is not a formal verification.

## 1. What was checked, and how

I restated each claim before reading its proof. I re-derived every implication by hand. I then tested the finite identities on explicit instances using code I wrote, which shares nothing with the author's or the certificates task's code.

| Claim in the artifact | Hand check | Instance check |
| --- | --- | --- |
| §1 Orbit equivalence: $\det(A+B)=\det(D_a+UD_bU^*)$ with $U=V^*W$, and its converse | Holds. No commutativity or invertibility is used. Each $z_\sigma$ is attained by a permutation unitary. | Not needed. |
| §2 Support criterion (1), both directions | Holds. Forward: apply the real-linear $L_w$ to a convex combination. Reverse: the nearest-point argument for a compact convex $H$ gives $L_w(z-y)\le0$ and $L_w(x-y)=\lvert w\rvert^2>0$. The argument is valid for point, segment and polygon hulls. The normalization $w=\lvert w\rvert e^{i\theta}$ gives (2). | Not needed. |
| §3 Opposite shifts, nonzero scaling by $c^n$, conjugation | Holds. $c=0$ is correctly excluded. | Float: 200 trials. |
| §3 Swap $a,b$ with $U\mapsto U^*$ | Holds. $D_b+U^*D_aU=U^*(D_a+UD_bU^*)U$, and the labels map $\sigma\mapsto\sigma^{-1}$. | **Exact**: 60 instances. Float: 200. |
| §3 $LUR$ invariance and $SU(n)$ normalization | Holds. $R$ commutes with $D_b$, and $L^*D_aL=D_a$. $R=\operatorname{diag}((\det U)^{-1},1,\dots)$ is unitary because $\lvert\det U\rvert=1$. | **Exact**, with unit-modulus Gaussian-rational phases: 60. Float: 200. |
| §3 θ=0 suffices over all spectra | Holds. Use $c=e^{-i\theta/n}$, so that $c^n=e^{-i\theta}$. | Float ($c^n$ scaling): 200. |
| §3 Normalization $t=b_1$, the radius $R$, and the $R=0$ case | Holds. $R=0$ forces $A=-b_1I$ and $B=b_1I$, so every value is $0$. Gaps scale by $R^{-n}$. Using a rational upper bound preserves $\mathbb Q(i)$. | Not needed. |
| §3 Reciprocal transform | Holds. $X^{-1}(X+Y)Y^{-1}=X^{-1}+Y^{-1}$. Inverses of normal matrices are normal. The common denominator $\prod(a_i-c)\prod(b_j+c)$ does not depend on $\sigma$. $A+B$ may be singular. | **Exact** determinant identity, and every vertex for every $\sigma$: 60. Float: 200. |
| §4 Lemma 1 (Cayley) | Holds, with the omission noted in §2 below. | **Exact**: $I\pm K$ invertible and $C(K)^*C(K)=I$, 60 instances. |
| §4 Lemma 2 (density) | Holds. $K=(I-U)(I+U)^{-1}$ is a function of the normal $U$ with imaginary values on its spectrum, so it is skew-Hermitian. $C(K)=U$. A phase perturbation removes the eigenvalue $-1$ while avoiding finitely many bad phases. | **Exact** round trip $C(K)=U$: 60. |
| §4 Rational-witness proposition | Holds. The gap (4) is continuous in $(a,b,U,w)$. A strict gap is an open condition. Gaussian-rational spectra and functionals, together with rational Cayley unitaries, are dense. $B=UD_bU^*$ then has $\mathbb Q(i)$ entries. The gap is rational. | Not needed. |
| §4 Semi-decision procedure | Holds as stated. It halts if a counterexample exists. No bound is claimed, and none is implied. The quoted verifier cap (`MAX_DIMENSION = 6`) matches `certificates/verify.py`. | Not needed. |
| §5 Boundary cases and remaining gap | Accurate. The stated non-implications are correct: compactness gives a maximizer that need not be a permutation matrix, and endpoint interpolation along a path is not convex. | Not needed. |

**Reproduction.** Run from the repository root:

```sh
python -I research/marcus-de-oliveira/reviews/check_reductions_exact.py
# seed 20261125, 60 trials, n in 2..4, Gaussian-integer entries in [-3,3]; all exact checks pass
```

The floating-point spot checks are `adversarial/identity_spotchecks.py` on PR #11, at commit `dd7da9f`: seed 20261125, 200 trials, $n\le6$, all within $10^{-9}$.

## 2. Minor omission (non-blocking)

**Lemma 1.** The proof establishes only that $I+K$ is invertible. It then writes $C(K)^*=(I-K)^{-1}(I+K)$, which also needs $I-K$ to be invertible. The same argument covers it, because $-K$ is skew-Hermitian. Suggested one-line fix: "Since $-K$ is also skew-Hermitian, $I-K$ is invertible too."

The exact checks confirm both matrices are invertible on every instance. This is a gap in the write-up, not in the mathematics.

## 3. Scope notes for users of this result (not defects)

1. **These are transformations of the whole problem, not of subfamilies.** The θ=0 reduction rotates both spectra by $e^{-i\theta/n}$. A nonreal opposite shift, the normalization $t=b_1$ with nonreal $b_1$, and the reciprocal transform with nonreal $c$ all leave structured classes, including:
   - Hermitian $A$, since a real spectrum becomes nonreal;
   - Shitov's $\mathcal A(m,k)$;
   - real spectra in general.

   The artifact restricts the θ=0 statement to "all spectra", which is correct. Downstream work on a subfamily must not borrow these normalizations without checking that they preserve that family. The conditional $n=4$ Hermitian route in `literature/README.md` is one such case.
2. **The Shitov "Theorem 12" attribution was not verified.** mathnet.ru refused connections during this review. The artifact proves the transform itself, so validity does not depend on the citation. Only the credit and theorem number are unconfirmed.

## 4. What this pass does and does not establish

**It establishes** that the reductions in the reviewed checkpoint are correct as stated, under their stated hypotheses. In particular:

- Any counterexample to MOC in dimension $n$ yields one with Gaussian-rational spectra, a rational Cayley unitary, and a rational separating functional with a positive rational gap.

**It does not establish:**

- anything about the truth of MOC;
- a bound on denominators or on search cost;
- correctness of the artifacts in `minors/`, `search/` or `certificates/`;
- the cited literature's theorems.
