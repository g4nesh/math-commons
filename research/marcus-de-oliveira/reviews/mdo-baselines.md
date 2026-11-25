# Review: `mdo-baselines`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-reductions` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-baselines-forest-published-20261125`, recorded 2026-11-25 18:35:14 UTC. It is the checkpoint behind the submit `mdo-baselines-forest-submit-20261125` (18:35:27 UTC). The ledger shows no newer work on this task.
- **Artifact:** `research/marcus-de-oliveira/baselines/` at commit `a0537fd`. That commit is on the PR #15 branch, and PR #15 was squash-merged as `8c01009`. All five files are byte-identical on `main` at `1b98a07`.

  | File | SHA-256 |
  | --- | --- |
  | [`README.md`](https://github.com/g4nesh/math-commons/blob/a0537fdfee8dad7e74b6fb36e507dd5a1ef6666c/research/marcus-de-oliveira/baselines/README.md) | `fcd31fe654df93196376e0b3e803ff796265e49072e5e1ea76073fd1fae1d0f4` |
  | [`forest.py`](https://github.com/g4nesh/math-commons/blob/a0537fdfee8dad7e74b6fb36e507dd5a1ef6666c/research/marcus-de-oliveira/baselines/forest.py) | `3164180fe4df14592b00dbf71ad515c79ff522c037b43c043f7a0567a6b86314` |
  | [`test_forest.py`](https://github.com/g4nesh/math-commons/blob/a0537fdfee8dad7e74b6fb36e507dd5a1ef6666c/research/marcus-de-oliveira/baselines/test_forest.py) | `6c85f999635d490f157ee8273d30950d52b1077a561b97348b5b6b8b99e10b5a` |
  | [`results.json`](https://github.com/g4nesh/math-commons/blob/a0537fdfee8dad7e74b6fb36e507dd5a1ef6666c/research/marcus-de-oliveira/baselines/results.json) | `b3574a9bf92a892a9af808d9064a99169e2d5f958fe472e231b820d05cf90995` |
  | [`SHA256SUMS`](https://github.com/g4nesh/math-commons/blob/a0537fdfee8dad7e74b6fb36e507dd5a1ef6666c/research/marcus-de-oliveira/baselines/SHA256SUMS) | `8f779228761766c3e0a67abb2fc02dea107dba7807033515766ea79366d6dfc2` |

  The imported dependency [`certificates/verify.py`](https://github.com/g4nesh/math-commons/blob/a0537fdfee8dad7e74b6fb36e507dd5a1ef6666c/research/marcus-de-oliveira/certificates/verify.py) at the same commit hashes to `8cb1a728…f78ff`, the value recorded in `results.json`.
- **Verdict:** **pass**. The forest theorem (5)–(8) and the triangle obstruction in §4 hold as stated. I record one attribution note and two scope notes in §2. None of them affects a claim the artifact makes.
- **Evidence level of this review:** an independent written re-derivation of every step, plus exact $\mathbb Q(i)$ checks in my own standard-library code, which shares nothing with `forest.py` or `verify.py`. I also re-ran the author's suite from a clean export. This is not a formal verification.

## 1. What was checked, and how

I restated each claim before reading its proof, then derived it again myself. Exact checks ran on 150 seeded forest cases with $n\le6$ and 52,620 squared minors. Each gate magnitude was drawn from five non-degenerate values ($q=s^2$ for $s=4/5,12/13,15/17,24/25,21/29$) plus 0 and 1, and every gate had independent unit phases.

| Claim | Hand check | Exact check (own code) |
| --- | --- | --- |
| (1) $d(P_p)=z_{p^{-1}}$ for $P_pe_j=e_{p(j)}$ | Holds. $P_pD_bP_p^*$ puts $b_j$ at row $p(j)$, and $\prod_j(a_{p(j)}+b_j)=\prod_i(a_i+b_{p^{-1}(i)})$. | 40 cases, using the objective's $z_\sigma$ convention directly |
| (2) Single $U(2)$ gate | Holds. The diagonal of $g\,\mathrm{diag}(b_i,b_j)g^*$ is $(tb_i+(1-t)b_j,\,(1-t)b_i+tb_j)$ and its determinant is $b_ib_j$. Expanding gives (2), which equals $(a_i+b_i)(a_j+b_j)+(1-t)(a_i-a_j)(b_i-b_j)$. | 40 gates with complex phases, $n\le5$ |
| (3) Disjoint blocks | Holds. The matrix is block diagonal, the transpositions commute, and the weights factor. | Covered by the forest cases |
| Compound transition table | Holds. If neither or both endpoints are occupied, the only nonzero entry is $1$ or $\det g$, of modulus 1. If exactly one is occupied, staying gives $\lvert g_{11}\rvert^2$ or $\lvert g_{22}\rvert^2=1-q$, and exchanging gives $\lvert g_{21}\rvert^2$ or $\lvert g_{12}\rvert^2=q$. | Implied by the minor checks |
| Path uniqueness (cut argument) | Holds. Deleting $e_r$ splits its tree into $C$ and the rest. Every other edge lies inside one side, because edges are distinct and there are no cycles. So only gate $r$ changes $\lvert X\cap C\rvert$, and only by $0$ or $\pm1$. This forces the move at every step, for any gate order, so at most one subset path runs from $T$ to $S$. The independent-transposition chain has the same allowed moves with the same transition masses. | Same |
| (5) Every squared minor | Follows from the two rows above, including $k=0,n$ and $q_r\in\{0,1\}$. | All $\sum_k\binom nk^2$ minors in each of the 150 cases, with weights computed by enumerating $\varepsilon\in\{0,1\}^m$ rather than the author's recursion |
| (7) Expansion | Holds. $\det(D_a+M)=\sum_S a^{S^c}\det M[S,S]$, and Cauchy–Binet with the diagonal $D_b$ gives $\det M[S,S]=\sum_T\lvert\det V[S,T]\rvert^2b^T$. There is no division. | 10 non-forest unitaries (triangle plus a repeated edge), checked against two determinant algorithms |
| (6) Determinant identity | Holds. Substitute (5) into (7). The inner sum is (7) for $P_{p_\varepsilon}$, which (1) evaluates. | Every forest case, with spectra that include zeros and repeats; a second Leibniz determinant for $n\le5$ |
| (8) Monomial end factors | Holds. $d_{a,b}(M_\ell VM_r)=d_{a\circ\ell,\;b\circ r^{-1}}(V)$, and substituting $j=r(k)$ gives $p'=\ell\circ p\circ r$. | Random phased monomials on both sides of every forest case: all minors and the determinant |
| §4 Triangle | Holds. Writing $U e_1$ out gives $U_{11}=c^2+s^3$. The product chain returns with mass $c^4+s^6$, so the gap is $2c^2s^3$. For $a=(0,1,1)$, $b=(1,0,0)$, (7) keeps only $S=T=\{1\}$, so $d(U)=\lvert U_{11}\rvert^2$. The $z_\sigma$ are $0$ and $1$. | $109/125$, $11881/15625$, $6121/15625$, $1152/3125$, all exact. $U_{11}=c^2+s^3$ also holds for four more Pythagorean pairs. |
| §3 Degeneracies, §4 scope, §6 limits | Accurate. The identities are polynomial in $a,b$, and §4 is correctly called a failure of these weights only. | Product weights also fail for a repeated edge, which the README says only breaks the proof |

**Author's suite.** From a `git archive` export of `a0537fd` (Python 3.11.7):

- `SHA256SUMS` verifies.
- `forest.py` reports 11 fixtures, 3,332 compound entries and `allChecksPassed: true`. Its JSON equals the committed `results.json` except for `pythonVersion` (3.9.6 there) and `processCpuSeconds` (0.711 s here).
- `unittest` ran 11 tests, all OK (0.707 s).

I read `forest.py` and `verify.py` before running them. Both are pure exact arithmetic, with no network access or subprocesses.

**Sensitivity of my checker.** I ran it with each of three injected bugs:

- using $z_p$ instead of $z_{p^{-1}}$;
- giving every gate the first gate's $q$;
- pushing the monomial weights forward as $r\circ p\circ\ell$.

It failed on each, as it should.

## 2. Findings (none blocking)

1. **Attribution of identity (2) (minor).** The single-gate formula, and the segment $[z_{\mathrm{id}},z_{\langle i,j\rangle}]$ it gives, appear in A. Kovačec, DMUC preprint p2611 v2 (4 Nov 2026): the Lemma and Corollary on pp. 3–4, and Corollary 4. Those statements cover $SU(2)$ generalized Givens matrices. Kovačec credits the $SO(2)$ case to N. Bebiano and J. F. Queiró, *Linear Algebra Appl.* 71 (1985) 23–28. The $U(2)$ case follows because a scalar phase cancels in $GD_bG^*$.
   - I read the v2 PDF. Its SHA-256 is `5221b16e…fe54`, the hash recorded in `literature/sources.json`. I did not read Bebiano–Queiró itself.
   - The README cites only the in-repo `adversarial/README.md` §4, which does name Kovačec.
   - *Correction:* add one sentence to §1 that cites Kovačec v2 (Lemma, Corollary 4) and, through it, Bebiano–Queiró 1985.
2. **Scope: forest products are a measure-zero family (note).** For $n\ge3$, every $V$ in the theorem has a zero entry, and so does every $M_\ell VM_r$. More precisely, $V_{ij}V_{ji}=0$ whenever $i,j$ are at tree distance at least 2.
   - *Proof:* take $k=1$ in the path argument. A single particle must cross the path's edges in increasing gate order in its direction of travel, and the reverse trip needs the opposite order.
   - If no component has three or more vertices, $V$ is block diagonal with blocks of size at most 2.
   - So the family lies in the proper algebraic set $\{\prod_{ij}U_{ij}=0\}$, which has Haar measure zero.
   - This was checked on all 125 cases with $n\ge3$.
   - The README's "does not recover all unitaries in dimension three" is accurate but weaker. *Suggested:* state the zero-entry property in §6.
3. **Scope of the §6 research question (note).** Two facts already in the repository narrow the question.
   - *Triangle only.* For gates supported on one triangle, $U\in U(3)\oplus I$. The minors report §3 (Theorem A), which the README already cites, then gives spectrum-independent weights that match every squared minor: a Birkhoff decomposition of $(\lvert U_{ij}\rvert^2)$ combined with complementary minors. For the README's own triangle $U$, I computed such weights on 4 permutations, and they match all 20 squared minors exactly.
   - *Dimension four.* `literature/README.md` reproduces Zheng's $O_*\in O(4)$ with $\psi(O_*)=-158/2401<0\le\psi(P)$ for all 24 permutation matrices; I re-checked this exactly. So no spectrum-independent weights match all squared minors of $O_*$.
   - *Consequence.* The question is open only for structured gate graphs, such as a triangle inside a larger tree or longer cycles. A positive answer must use that structure.
   - *Suggested:* say so in §6, and link the $O_*$ witness.
4. **Test coverage (nit).** 29 of the 31 fixture gates have the same magnitude $q=16/25$; the other two have $q=0$ and $1$. A bug that pairs each gate with the wrong magnitude is caught only through that one endpoint fixture. A reversed-order mutation I tried was caught there. My check varies $q$ across gates. *Suggested:* add one fixture with distinct non-degenerate $q$ per gate.

## 3. Reproduction

From the repository root. The baselines files are identical at `a0537fd` and on `main`.

```sh
(cd research/marcus-de-oliveira/baselines && shasum -a 256 -c SHA256SUMS)
out=$(mktemp) && PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/baselines/forest.py --output "$out"
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/marcus-de-oliveira/baselines -p 'test_*.py' -v
python3 -I research/marcus-de-oliveira/reviews/check_baselines_ishaan.py
# seed 20261125: 150 forest cases, 52620 squared minors, all exact checks passed; triangle: ...
```

[`check_baselines_ishaan.py`](check_baselines_ishaan.py) uses only the standard library: `fractions`, `itertools` and seeded `random`. It is deterministic and took about 6 CPU seconds here.

## 4. What this pass does and does not establish

**It establishes** that the written proof of (5)–(8) is complete and correct for the stated family: arbitrary complex $U(2)$ gates on distinct forest edges, in any order, with monomial ends, for all complex spectra including degenerate ones. It also establishes that the §4 triangle numbers are exact and correctly scoped.

**It does not establish:**

- anything about MOC in general;
- the novelty of the forest family. The README claims none, and I did no literature search for products of Givens matrices on trees;
- the correctness of `minors/report.md`, or of the Zheng and Bebiano–Queiró sources beyond the specific items quoted above;
- correctness for gate sets with cycles or repeated edges.

The exact checks cover finitely many instances. The general statement rests on the written proof, which I checked by hand.
