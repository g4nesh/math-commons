# Review: `mdo-four-cycle-spectral-weights`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-search` and the supplement's author `ganesh-coordinator` both belong to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-four-cycle-spectral-final-20261125`, 2026-11-25 20:40:51 UTC, submitted as `mdo-four-cycle-spectral-submit-20261125` at 20:40:56 UTC. The task has no later work. The pinned commit is `fa8b4ea56e894f2287ce1a8f1c6a04059679de2b` (PR #21). The `four-cycle-spectral/` folder and its two dependencies are byte-identical on `main` at `1b98a07`.
- **Artifact** (SHA-256 of each file, recomputed from `git archive` of the pinned commit; all six match the folder's `SHA256SUMS`):
  - [`four-cycle-spectral/README.md`](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/README.md): `5f5df87987c3d5ae56a88135c3b1fae55044db6832bfe3e39f04f79116c93a6d`
  - [`four-cycle-spectral/certificate.json`](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/certificate.json): `73a5f6a51e1f699cfd114ad39922a7b89e243d3c329d305853eb28b8f0e95646`
  - [`four-cycle-spectral/fiber_certificate.json`](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/fiber_certificate.json): `d56447438ca64e58e7f776b998f6b384d49e7b9009e54e1308d588c0cdcc95c3`
  - [`four-cycle-spectral/fiber_certificate.py`](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/fiber_certificate.py): `91a1f30329478dbf23376e0b43fdae0a763cc2dd1038f8f93df6bc906e3fac17`
  - [`four-cycle-spectral/spectral.py`](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/spectral.py): `6c416aa47d5f44bd797b5ce1c6fb32e3d73ccb9f22721ca444104f2fa5f827b5`
  - [`four-cycle-spectral/test_spectral.py`](https://github.com/g4nesh/math-commons/blob/fa8b4ea56e894f2287ce1a8f1c6a04059679de2b/research/marcus-de-oliveira/four-cycle-spectral/test_spectral.py): `37e08aa89da0c75f9af24e8b55d11eedb29b96c4e958444e9fe9ecef785722a6`
  - Imported dependencies `single-cycle/results.json` (`7c04b6fe…173ce3`) and `certificates/verify.py` (`8cb1a728…bf78ff`) match the hashes recorded in `certificate.json`. `single-cycle/` is unchanged between `bbbf34f` (the commit pinned by the single-cycle task's final checkpoint, merged in PR #17) and this commit, as the README states.
- **Verdict:** **pass**. I found no defect. Two observations for users follow in §3; neither needs a change.
- **Evidence level of this review:** independent written re-derivation, plus an exact $\mathbb Q(i)$ recomputation with my own standard-library script, which imports no repository module and uses a different route (full symbolic Leibniz expansion in all four eigenvalues). I also re-ran the author's scripts and tests from a clean export. This is not a formal verification.

## 1. Claims restated, with their quantifiers

Throughout, $U$ (the $4\times4$ rational orthogonal matrix of README §1) and $b=(3-2i,2,1+2i,-i)$ are fixed, $A=D_a$, $B=UD_bU^*$, and $z_p(a)=\prod_j(a_{p(j)}+b_j)$ over the 24 column-to-row permutations $p$. This is the same hull as $\prod_i(a_i+b_{\sigma(i)})$ with $\sigma=p^{-1}$.

| README | Free | Fixed | Claimed certificate |
| --- | --- | --- | --- |
| §2 | real $t\in[-40,1]$, $a_0=-2+i(t-1)$ | $a_1=-1,\ a_2=i,\ a_3=1-i$ | three masses $w_j(t)=N_j(t)/K(t)$ on a fixed triangle, all $>0$ |
| §2, sign table | real $t$ | as above | this triangle's weights are $\ge0$ exactly on $[\alpha,\gamma]$, $\alpha\in(-41,-40)$, $\gamma\in(1,3/2)$ |
| §3 | $a_0\in\mathbb C$ | $a_1,a_2,a_3$ | five constant positive masses |
| §4a | $(a_0,a_1,a_2)\in\mathbb C^3$ | $a_3=1-i$ | twelve constant positive masses |
| §4b | $a\in\mathbb C^4$ | $U$, $b$ only | no single distribution works; seven-point real functional separates, gap $11201152/15625$ |

The checkpoint detail and the submit summary state the same quantifiers. Both say no general theorem and no pointwise counterexample is claimed. The README matches that.

## 2. What was checked, and how

| Claim | Hand check | Exact check (my code) |
| --- | --- | --- |
| §1 $U$ unitary; data equal the stored fixture; family normal and noncommuting | $A$ diagonal and $B$ unitarily diagonalizable, so both are normal. $[D_a,B]_{01}=(a_0-a_1)B_{01}$, and $\operatorname{Re}(a_0(t)-a_1)=-1$, so the entry is nonzero for every real $t$ once $B_{01}\ne0$. | $UU^*=U^*U=I$. $U,a,b$ and the triangle equal `fourCycleGreedyRuleFailure` in `results.json`. $BB^*=B^*B$, $B_{01}\ne0$. |
| §1 formula (1) and the vertex table | Only the $(0,0)$ entry moves, so $d$ is affine in $t$ with slope $i$ times the complementary principal minor. Each $z_p$ has one factor containing $a_0$. | Matches exactly. I specialized the full multi-affine expansion of $\det(D_a+B)$, not a cofactor. |
| §2 area and numerators | $[z_2-z_1,z_3-z_1]=500((t-2)^2+2)>0$. Cramer's rule in $\mathbb R^2$ gives the cyclic cross-product numerators. | Area $=(3000,-2000,500)$. $15625\cdot[z_{j+1}-d,z_{j+2}-d]$ equals $N_1,N_2,N_3$. $\sum N_j=K$. $\sum N_jz_j=Kd$ holds identically as cubics in $t$. |
| §2 positivity on the closed interval | The concavity, discriminant and monotonicity arguments are each valid. $K\ge 15625000$ on $\mathbb R$. | A different route: the exact minimum of each quadratic on $[-40,1]$, taken over the endpoints and the vertex, is $3403196$, $15269809024304/1550611$ and $5443030$, all positive. |
| §2 sign-change range $[\alpha,\gamma]$ | $N_1$ is concave with roots $\alpha<\beta$. $N_3$ is convex with roots $\gamma<\delta$. $N_2>0$ everywhere. The brackets are disjoint and ordered, so the nonnegative set is $[\alpha,\beta]\cap((-\infty,\gamma]\cup[\delta,\infty))=[\alpha,\gamma]$. | A sign change in each bracket. The discriminant of $N_2$ is $-244316944388864$. In floating point, $\alpha\approx-40.608$ and $\gamma\approx1.4918$. |
| §2 values and nonconstancy | | $w(0)$ equals the three stored `originalHullConvexWeights`, as uniqueness of barycentric coordinates requires. $w_1(0)=1100333/3906250$ and $w_1(1)=804497/2343750$. |
| §3 five masses for all $a_0\in\mathbb C$ | $d$ and every $z_p$ are holomorphic affine in $a_0$, so matching both complex coefficients proves the identity on $\mathbb C$. | The masses are positive and sum to $1$. $d(z)$ matches the README. The weighted constant and slope both equal those of $d$. |
| §4a twelve masses on $\mathbb C^3$ | $\det(D_a+B)$ and each $z_p$ have degree $\le1$ in each $a_i$ separately, so the corner-interpolation argument is valid. | The masses are positive and sum to $1$. All eight monomial coefficients in $a_0,a_1,a_2$ (at $a_3=1-i$) agree. This is a third route, besides the author's corners and principal minors. |
| §4b separator | Linearity of $\mathcal L$ gives $\mathcal L(d)=\sum_pw_p\mathcal L(z_p)\le0$ for any fixed distribution. The argument is valid. | $\max_p\mathcal L(z_p)=0$, attained by 10 of the 24 pairings. $\mathcal L(d)=11201152/15625$. |
| §4b pointwise inclusion at the seven diagonals | | Each $d(x^{(r)})$ lies in the convex hull of its own 24 values. I searched exact points, segments and triangles. |
| Stored `certificate.json` | | Every polynomial, mass and the six parameter checks equal my recomputation. |
| §5 reproduction | | All three README commands pass from a clean export (Python 3.11.7). The test suite reports `Ran 12 tests … OK`. |

**Robustness of my checker.** I ran it on three corrupted copies: two masses shifted by $\pm1/390625$, one separator coefficient changed, and one stored parameter-check weight changed. It rejected all three.

**Supplementary sanity check, not published.** I also ran a seeded numpy cross-check (seed 20261125). It used 2000 random points for each certificate, plus both endpoints of the interval. The largest relative error was below $10^{-14}$, and the smallest interval weight was about $2.47\times10^{-4}$, at $t=-40$. At the separator diagonal $a=(1,0,1,0)$, where $a_3=0$, the twelve masses miss the determinant by about $0.29$. This is consistent with §4b.

## 3. Observations (non-blocking, no change required)

1. **One of the seven separator diagonals is trivial.** At $x=(1,1,1,1)$, $A=I$, and all 24 pairing values equal $\det(I+B)=\prod_j(1+b_j)=48-24i$. Its inside-hull certificate is therefore automatic, and a triangle-only membership search will not find it. The other six diagonals have two-dimensional hulls (4 or 6 distinct pairing values each). For each, the determinant has strictly positive exact barycentric coordinates in some triangle of pairing values, so it lies strictly inside; its floating-point distance to the hull boundary is at least about $0.82$. The README's sentence on line 150 is true as written.
2. **The interval is tight at its left end.** The root $\alpha\approx-40.608$ of $N_1$ lies within $0.61$ of $-40$, so $w_1(-40)=3403196/13796875000\approx2.5\times10^{-4}$. The proof is exact, so this affects nothing. A user who perturbs $U$ or $b$ should not expect the left endpoint to survive.

## 4. Reproduction

Run from the repository root at the reviewed commit, with this review's script added:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/four-cycle-spectral/spectral.py --check research/marcus-de-oliveira/four-cycle-spectral/certificate.json
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/four-cycle-spectral/fiber_certificate.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/marcus-de-oliveira/four-cycle-spectral -p 'test_spectral.py' -v
python3 -I -B research/marcus-de-oliveira/reviews/check_four_cycle_spectral_weights_ishaan.py
# prints a JSON summary; exits nonzero on the first failed assertion (under 1 s)
```

The first two commands print `"exactParameterChecks": 6` and `"strictGap": "11201152/15625"`, respectively. My script has no randomness.

## 5. What this pass does and does not establish

**It establishes** that, for this fixed $U$ and fixed $b$, the checkpoint's claims are correct as stated:

- the closed-interval three-point certificate;
- the exact nonnegativity range $[\alpha,\gamma]$ of that triangle;
- the five-mass certificate for every complex $a_0$;
- the twelve-mass certificate for every $(a_0,a_1,a_2)\in\mathbb C^3$ at $a_3=1-i$;
- the exact impossibility of one fixed distribution over all of $\mathbb C^4$.

**It does not establish:**

- anything about MOC beyond these fixed-$U$, fixed-$b$ families;
- whether such certificates exist for other $a_3$, $U$ or $b$;
- the source fixture's no-common-minor property, which belongs to another task and was not re-reviewed here;
- novelty (none is claimed, and the artifact cites no external result).

My exact checks cover only the stated quantities and instances. The hand derivations are the argument for the continuous claims.
