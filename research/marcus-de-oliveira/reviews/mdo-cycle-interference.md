# Review: `mdo-cycle-interference`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-search` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-cycle-final-20261125`, 2026-11-25 19:53:17 UTC, submitted as `mdo-cycle-submit-20261125` at 19:53:22 UTC. The task has no later work. The pinned commit is `bbbf34f73a498d0ad6c0445aedf96a6c697e0815` (PR #17). The `single-cycle/` folder is byte-identical on `main` at `1b98a07`.
- **Artifact** (SHA-256 of each file, recomputed from `git archive` of the pinned commit; all six match the folder's `SHA256SUMS`):
  - [`single-cycle/README.md`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/README.md): `da9581175a3f7a7ecdf344bcc1901cbca798805746352d6c41ba3f34955cc825`
  - [`single-cycle/cover_certificate.py`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/cover_certificate.py): `8b2e7552ec7fec675eeac5fd46ccfcabff7b16026f70b78ab1f2f22230270f50`
  - [`single-cycle/four-cycle-obstruction.json`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/four-cycle-obstruction.json): `790eb87398bdc20aca56ad553947764f80f65761fb4a92ec4751168d8c20a5d8`
  - [`single-cycle/cycle.py`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/cycle.py): `4de19ee09101c6a95e2efe968847bacd4c041ce0931f86baf7350e966a24dc0c`
  - [`single-cycle/test_cycle.py`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/test_cycle.py): `57cd6828a9a804ed8be2dc77e3b103fc114d8804b3996f6f8d26d1227b9d34a3`
  - [`single-cycle/results.json`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/results.json): `7c04b6fe24ddd906fc747886ee87cc3e127e34ece9b736e7cadbfc7683173ce3`
  - Imported dependencies `baselines/forest.py` (`3164180f…a6b86314`) and `certificates/verify.py` (`8cb1a728…bf78ff`) match the hashes recorded in `results.json`.
- **Verdict:** **pass**, with two minor notes (§2). Neither affects a claim of the artifact.
- **Evidence level of this review:** independent written re-derivation of §§2–5, plus exact $\mathbb Q(i)$ recomputation with my own standard-library script, which imports no repository module. I also re-ran the author's scripts and tests from a clean export. This is not a formal verification.

## 1. What was checked, and how

I restated each claim, re-derived it by hand, and then recomputed the finite statements in exact arithmetic. My instances include generic complex Cayley unitaries $U=(I-iH)(I+iH)^{-1}$ with Gaussian-rational Hermitian $H$. Every entry of those three $3\times3$ blocks is nonzero, each block has nonreal entries, and each has $L<H$. These are not gate products.

| Claim in the artifact | Hand check | Exact instance check (my code) |
| --- | --- | --- |
| §1 convention: $P_pe_j=e_{p(j)}$, and (1) asks for $\Pr[p(T)=S]$ | Consistent with the code, where `p[j]` is the row receiving column $j$. | Used throughout. |
| §2 six-weight table: marginals equal $D$, $[L,H]$ is nonempty, and every distribution with marginals $D$ has this form | Holds. I summed the two weights with $p(j)=i$ for all nine entries, using only the row and column sums. Both chains $d_{22}-d_{01}=d_{00}-d_{12}=d_{11}-d_{20}$ and $d_{11}-d_{02}=d_{00}-d_{21}=d_{22}-d_{10}$ hold. Necessity follows by setting $x=w(\mathrm{id})$. | The marginal map $\mathbb R^{S_3}\to\mathbb R^{3\times3}$ has rank 5, with kernel spanned by the sign vector, which is exactly the table's $x$-coefficient vector. On three Cayley unitaries, $x\in\{L,(L+H)/2,H\}$ match all 20 compound entries (9 cases). A non-uniform distribution is recovered from its marginals with $x=w(\mathrm{id})$. |
| §2 extension to the $2\times2$ layer for arbitrary complex $3\times3$ $U$ | Holds. $\operatorname{adj}U=\det U\,U^*$ gives $\lvert\det U[\hat i,\hat j]\rvert=\lvert u_{ij}\rvert$, and $p(\hat j)=\hat i\iff p(j)=i$. | Covered by the 20-entry checks above. |
| §2 baseline triangle $R_{13}R_{23}R_{12}$, $c=3/5$ | $\lvert U_{00}\rvert^2=(109/125)^2$. The stay/swap rule gives $\tfrac{9}{25}\cdot\tfrac{9}{25}+\tfrac{256}{625}\cdot\tfrac{16}{25}=6121/15625$. | $[L,H]=[1881,2025]/15625$. Weights at $L$ are $(1881,10000,144,3456,0,144)/15625$. All 20 entries match. Return mass is 11881/15625 against 6121/15625, the latter from my own propagation of the stay/swap rule. |
| §3 gluing with $\lvert A\cap B\rvert\le1$ | Holds. A nonzero $X[S,R]$ forces $R\setminus A=S\setminus A$, and a nonzero $Y[R,T]$ forces $R\setminus B=T\setminus B$. Only bits in $A\cap B$ stay free, and $\lvert R\rvert=\lvert S\rvert$ fixes the single remaining one. So Cauchy–Binet has at most one nonzero term. Independence gives (3), and $P_pP_q=P_{p\circ q}$ matches the code's `p[q[j]]`. The hypothesis that each certificate lives on permutations of its own support is automatic: $\lvert X_{ii}\rvert^2=1$ off $A$ forces $p(i)=i$. | $n=6$: a Cayley triangle core plus Cayley $2\times2$ blocks on $\{2,3\}$ (left), $\{0,4\}$ (right) and $\{3,5\}$ (left), 924 entries. $n=5$: two Cayley triangles sharing pivot 2, 252 entries. The negative control, two triangles overlapping in $\{1,2\}$, fails, as expected. |
| §3 "a simple cycle of length ≥4 is not supplied" | Holds. Any split of a 4-cycle's edges into two nonempty groups has vertex supports meeting in at least two coordinates. | — |
| §3 "an arbitrary interleaving … is not proved here" | True as written. In fact the restriction is necessary; see note 2.1. | Exact counterexample, note 2.1. |
| §4 $\det(D_a+UD_bU^*)=\sum_pw_p\prod_j(a_{p(j)}+b_j)$ | Holds. Diagonal expansion and Cauchy–Binet give $\sum_{S,T}\lvert\det U[S,T]\rvert^2\prod_{i\notin S}a_i\prod_{j\in T}b_j$. Substituting (1) and using $\{i\notin p(T)\}=p(T^c)$ yields the product. The vertices are the conjecture's with $\sigma=p^{-1}$. No division occurs. | Equality is exact on all triangle and glued instances, with complex spectra, plus repeated and zero spectra on the triangles. |
| §5 four-cycle $R_{14}R_{34}R_{23}R_{12}$: matrix, greedy mismatch | — | $625U$ equals the displayed matrix, and $UU^{\mathsf T}=I$. My lexicographic Birkhoff decomposition matches all 16 singletons and puts 3321/15625 on $(\{0,1\},\{0,1\})$, against $\lvert\det\rvert^2=81/625$. |
| §5 13-pair cover: $h(p)\ge1$ for all $p$ and value $14953/15625$ | The logic holds. It needs only the cover minors and total mass 1, so it excludes every probability distribution. | The hit vector equals the README's list. The value is $14953/15625$, with gap $672/15625$. The stored JSON matrix and all 13 squared minors equal my recomputation. |
| §5 tested determinant inside the hull | — | Actual $=(269802-5707764i)/390625$ and greedy $=(1040842-5724244i)/390625$. The three stored hull points equal $\prod_j(a_{p(j)}+b_j)$, with weights $\ge0$ summing to 1 that reconstruct the determinant exactly. |
| §6 counts and reruns | $1416=12\cdot20+924+252$. | From a clean export, `SHA256SUMS` passes, the cover checker verifies, and all 10 tests pass (Python 3.11.7). Re-running `cycle.py` reproduces `results.json`, except for `pythonVersion` and `processCpuSeconds`. |
| Attribution and scope | The $n=3$ result is credited, through the minors report's Theorem A, to Bebiano–Merikoski–da Providência (LMA 20, 1987). A web search found secondary citations (Drury) that credit the $3\times3$ case to these three authors. They did not show the year or volume, and I did not retrieve the 1987 paper, so the year and volume are taken from the minors report. No novelty is claimed, and the README states every non-implication for MOC. | — |

## 2. Minor notes (non-blocking)

**2.1 The schedule hypothesis is necessary, not just unproved.** On the same triangle-plus-pendant graph, take the gates in the application order $01,12,23,02$: the first three are $[[3/5,4/5],[-4/5,3/5]]$, and the last has $c=4/5$, $s=3/5$. This unitary has **no** common-minor certificate. My script checks a 15-pair cover with $h(p)\ge1$ for all 24 permutations, whose squared minors sum to $14356/15625$ (gap $1269/15625$). The same four gates in the order $01,12,02,23$ are certified by §2 plus gluing (all 70 entries).

The README's sentence "This includes a triangle core with any tree grown by pendant coordinate-pair gates" is qualified two sentences later, so no claim is wrong. *Smallest correction:* replace "An arbitrary interleaving … is not proved here" with "an arbitrary interleaving can fail", citing an example such as this one, and attach "under a recursive left/right schedule" directly to the "any tree" sentence.

**2.2 Reproduction step 1 breaks the pinned hash.** Step 1 of README §6 writes `--output` over the committed `results.json`. Because `processCpuSeconds` changes on every run, `shasum -a 256 -c SHA256SUMS` then fails for `results.json`. I confirmed this in a separate export. *Smallest correction:* write to a scratch path and compare after dropping `pythonVersion` and `processCpuSeconds`, or say so in §6.

**Scope note (not a defect).** §§5–6 say that spectrum-dependent determinant weights remain open for the four-cycle. For this $4\times4$ matrix, the repository's literature record lists a claimed full $n=4$ proof (Zheng, October 2026, audit pending) that would cover it. A one-line pointer would show readers where the open cycle question sits, either at $n\ge5$ or in that audit.

## 3. Reproduction

From the repository root, at `bbbf34f` or `main` (the `single-cycle/` folder is identical), with this review's folder present:

```sh
python3 -I research/marcus-de-oliveira/reviews/check_cycle_interference_ishaan.py
# prints triangleCayleyChecks 9, gluingEntries [924, 252], fourCycle cover 14953/15625 (gap 672/15625),
# interleaved cover 14356/15625 (gap 1269/15625); exits nonzero on any failed assertion
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/single-cycle/cover_certificate.py --check research/marcus-de-oliveira/single-cycle/four-cycle-obstruction.json
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/marcus-de-oliveira/single-cycle -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/single-cycle/cycle.py --output /tmp/results.rerun.json
```

The checker uses only the standard library and runs in under one second. It reads `results.json` and `four-cycle-obstruction.json` only to compare their stored values with fresh recomputation. I mutation-tested it: changing any one of the cover value, the interleaved cover value, the triangle weight vector or one interleaved cover pair makes it fail.

## 4. What this pass does and does not establish

**It establishes** that the checkpoint's claims hold as stated:

- the explicit $n=3$ certificate family and the triangle repair;
- the one-coordinate gluing theorem under its recursive-schedule hypothesis;
- the implication from a certificate to hull inclusion for all complex spectra;
- the exact non-existence of a common-minor certificate for the stated four-cycle unitary;
- the hull membership of its one tested determinant.

**It does not establish:**

- anything about the truth of MOC, including for this four-cycle with other spectra;
- certificates for any interleaved schedule, or for cycles of length at least four;
- a formal proof;
- the claimed $n=4$ proof in the literature record;
- the content of the cited 1987 paper.

The floating-point exploration I used to find the example in note 2.1 is not evidence. Only the exact cover in the script is.
