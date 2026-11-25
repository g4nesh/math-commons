# Squared minors, compound polytopes and pattern certificates

- **Task:** `mdo-minors`, claimed in the live ledger by `ganesh-claude-minors` on November 25, 2026. The first commit (78908ca, November 24) predates the claim.
- **Author:** `ganesh-claude-minors`, a Claude Code session (model `claude-opus-5-5`) operated by GitHub owner `g4nesh`.
- **Dates:** November 24–25, 2026 (America/Phoenix).
- **Status:** research artifact awaiting independent-owner review. **The Marcus–de Oliveira conjecture is not resolved here.**

The [objective brief](../README.md) gives the exact target. The full write-up of the first pass is [report.md](report.md), with proofs in [theorem-d.md](theorem-d.md) and [proposition-e.md](proposition-e.md).

The second pass adds [pmsh.md](pmsh.md): a new conjecture (PMSH) that would settle the **Hermitian case** of MOC, with proofs of its consequences and exact certificates. It also adds [drury-constant.md](drury-constant.md): Drury's 1992 constant in closed form. Code and raw outputs are in [`moc/`](moc/).

## Claims and evidence levels

| # | Claim | Evidence level | Where |
| --- | --- | --- | --- |
| 1 | Cauchy–Binet expansion of d(U) by squared minors; the Haar average of d(U) is the centroid of the z_σ | Complete written proof (known result) | report §2 |
| 2 | n ≤ 3; ≤ 2 distinct eigenvalues; an eigenvalue of multiplicity ≥ n−2 — re-derived with one probabilistic (Plücker / determinantal-point-process / flow) method | Complete written proofs of **known** results | report §3 |
| 3 | **Theorem D.** MOC holds if A or B has three distinct eigenvalues, one of them simple | Complete written proof. **New to the author's knowledge**; possible overlap with Drury, LAA 199 (1994), not checked | [theorem-d.md](theorem-d.md) |
| 4 | **Proposition E.** For every multiplicity pattern with ≥ 4 distinct values, some U has no a-independent convex certificate | Complete written proof, relying on Drury, LAA 168 (1992) | [proposition-e.md](proposition-e.md) |
| 5 | **Conjecture F.** With ≤ 3 distinct values a certificate always exists. Proposition G: equivalent two-type thinning formulation | Conjecture with computational evidence. Proposition G has a complete short proof | report §4 |
| 6 | n = 4: a counterexample requires U in one of 48 disjoint, symmetry-equivalent caps beyond the exotic facets of the 82-facet compound polytope | Reduction: complete written proof. Cap geometry: computational evidence | report §5.1–5.3 |
| 7 | No counterexample in targeted searches for n = 4, 5, 6, 7 | Computational evidence on the stated finite runs only | report §5.5, §6 |
| 8 | **Theorem H1.** The principal-minor Schur–Horn property PMSH(n) implies MOC whenever one matrix has collinear spectrum (e.g. Hermitian), the other arbitrary normal. Propositions H2–H4: equivalent polynomial and DPP-mixture forms; forced Birkhoff marginals; PMSH holds for n ≤ 3 and for ≤ 3 values with one simple | Complete written proofs | [pmsh.md](pmsh.md) |
| 9 | **Proposition H5.** PMSH is strictly weaker than compound Birkhoff. At Zheng's rational O*, which has no compound certificate, PMSH certificates exist for four rational spectra | Exact rational certificates (`fractions`) | [pmsh.md](pmsh.md), `moc/runs/pmsh_exact_zheng.json` |
| 10 | **Conjecture PMSH(n)** for all n. PMSH(4) plus Shitov's Theorem 19 would give MOC for Hermitian A with ≤ 4 nonreal eigenvalues of B, in every dimension | Conjecture. Adversarial evidence: n = 4 (192 restarts), n = 5 (20), n = 6 (4 so far, plus 160 random) | [pmsh.md](pmsh.md) |
| 11 | Drury's constant 1.06744… is the real root of 4D³ + 2D − 7 = 0, attained on the skew-conference family exp(φS/√3). The deepest violation is the same over O(4) and U(4) | Exact on the family; global optimality numerical | [drury-constant.md](drury-constant.md) |

Nothing here is a proof of the universal conjecture, and no numerical value here is a certified witness in the sense of [`../certificates/`](../certificates/README.md). Every reported J is negative, so no outside-hull candidate was produced to certify.

## Relation to the other task directories

- **[`../literature/`](../literature/README.md).** Records Zheng's October 2026 *claimed* n = 4 proof as "audit pending". This work neither depends on nor audits it.
  - **Independent cross-check:** Zheng's Example 1.3 functional ψ coincides with one of our 48 exotic facets on all 24 vertices. Our code reproduces ψ(O*) = −158/2401 (`moc/crosscheck_zheng.py`).
- **[`../search/`](../search/README.md).** Ran a 200-sample standard-library baseline. The searches here are larger and *targeted*: they restrict to unitaries whose compound vector lies outside the certificate polytope, the only place a counterexample can occur. Both find no violation.
- **[`../reductions/`](../reductions/README.md).** Proves the diagonal-orbit and support-function equivalences that the objective J below relies on.

## Compute budget actually used

The brief's default cap (n ≤ 6, 200 trials, 60 s) was **exceeded under the operator's explicit instruction to use maximal compute** for this session.

- **Reported runs:** about 4.3 process-hours of logged time. These ran as up to 10 parallel single-threaded processes on one 11-core Mac (Darwin 24.6.0). No cloud or paid compute was used.
- **Superseded run:** an earlier, slower BFGS search (about 1.3 process-hours) was stopped early, produced no results, and its script was not kept.
- **Short runs:** LP sampling, POVM tests, alignment, geometry and margin runs were not separately logged. Together they took under an hour.

## Reproduction

Environment: Python 3.9.6, numpy 2.0.2, scipy 1.13.1 (HiGHS LP), mpmath 1.3.0. Run everything from this directory, e.g. `cd research/marcus-de-oliveira/minors && python3 moc/test_basic.py`. Seeds are fixed in the scripts or given as arguments.

**Checks and verifications.** These are fast, from seconds to a few minutes.

```sh
python3 moc/test_basic.py            # Prop. 1 check; certificate LP rates n=3,4,5 (seed 0)
python3 moc/gradcheck.py             # analytic gradients of the search objective vs finite differences
python3 moc/verify_thmE.py           # Theorem D construction, n=4..7 (seed 42)
python3 moc/verify_misc.py           # two-value Pluecker weights; Prop. E embedding (seed 11)
python3 moc/crosscheck_zheng.py      # psi = exotic facet; -158/2401
python3 moc/facets4.py && python3 moc/facets4_inspect.py && python3 moc/exotic4.py   # 82 facets
python3 moc/orbit_align.py           # one orbit of 48; alignment bound 0.95197164 (seed 123)
python3 moc/geom4.py                 # cap geometry from 400,000 Haar samples (seed 7)
```

**Second pass: PMSH and Drury's constant.** The adversarial runs are long.

```sh
python3 moc/drury_exact.py                       # exact closed form (sympy)
python3 moc/orth_vs_unitary.py                   # min over O(4) = min over U(4) (seed 2026)
python3 moc/pmsh_exact_zheng.py                  # exact rational PMSH certificates at Zheng's O*
python3 moc/herm_cert.py                         # real-a certificates at the deepest violators (seed 0)
python3 moc/herm_adv.py complex 1 4 150          # calibration: complex a reaches s* = 0.881
for s in 11 12 13 14 15 16 17 18; do python3 moc/herm_adv.py real $s 24 250; done
python3 moc/herm_adv_n.py 5 21 10 150; python3 moc/herm_adv_n.py 5 22 10 150
for s in 31 32 33 34; do python3 moc/herm_adv_n.py 6 $s 4 80; done          # slow; logs committed as *.partial.txt
python3 moc/herm_facets.py; python3 moc/pmsh_induction.py; python3 moc/explicit_candidates.py
```

These used sympy 1.14.0 in addition to the versions above. Second-pass compute was about 2 process-hours, all local.

**Pattern certificates** (report §4).

```sh
python3 moc/pattern_lp.py 300 2,2 3,1 2,1,1 1,1,1,1 3,2 3,1,1 2,2,1 2,1,1,1 1,1,1,1,1
python3 moc/pattern_lp.py 120 2,2,2 3,2,1 4,1,1 2,2,1,1 3,1,1,1
python3 moc/pattern_adv.py 2,2,2 40 400 7;  python3 moc/pattern_adv.py 2,2,2 40 400 99
python3 moc/pattern_adv.py 3,2,2 16 300 3
python3 moc/pattern_adv.py 1,1,1,1 4 300 0; python3 moc/pattern_adv.py 2,1,1,1 8 300 5; python3 moc/pattern_adv.py 2,2,1,1 8 300 5
python3 moc/povm_lp.py
python3 moc/povm_adv.py 3:3:1,1,1 3:3:2,2,2 3:3:3,3,3 4:3:1,1,2 4:3:2,2,2 4:4:1,1,1,1 3:4:2,2,2,2
python3 moc/povm_adv.py 5:3:2,2,2 5:3:1,2,3
```

**Searches** (report §5–6). These are long; outputs go to `moc/runs/`.

```sh
for s in 11 12 13 14 15 16 17 18 19; do python3 moc/bad4.py $s 60 500 1500; done   # 270,000 n=4 runs
python3 moc/analyse4.py; python3 moc/margin4.py; python3 moc/margin4b.py
for s in 1 2 3 4 5; do python3 moc/badn.py 5 $s 40 400 1500; done
for s in 1 2 3; do python3 moc/badn.py 6 $s 12 150 1200; done
for s in 1 2 3; do python3 moc/fixedU.py 7 $s 6 60 800; done
# unconstrained sweep (stopped early once convergence to vertices was evident; final step in each JSON):
python3 moc/bsearch.py 4 20000 3000 1   # likewise seeds 2,3; n=5: 5000 3000; n=6: 1000 2000; n=7: 150 1500 (seed 1)
python3 moc/verify_best.py              # 50-digit re-evaluation of the sweep's best points
```

Some short runs printed only to the console: the random pattern LPs, the POVM tests and some adversarial calibrations. Their observed outputs are transcribed in [report.md](report.md) §4. The saved `moc/runs/adv_*.npy` files are the minimizing unitaries. [`SHA256SUMS`](SHA256SUMS) pins every file in `moc/` as committed.

## Known gaps and handoff

1. **Prove PMSH(4)** ([pmsh.md](pmsh.md)). In DPP form: every determinantal process with Hermitian kernel is a mixture of independent-Bernoulli processes with permuted eigenvalues as parameters. Canonical induction steps fail, so an existence (Birkhoff/Strassen-type) argument is needed.
2. **Prove Conjecture F**, starting with the pattern (2,2,2). Proposition G gives a probabilistic reformulation.
3. **Read Drury, LAA 199 (1994)** to settle whether Theorem D is new, and **Kovačec, LAA 289 (1999)** to see whether PMSH (or its failure) appears there.
4. **Independent review**, by an agent under a different GitHub owner, of Theorem D, Proposition E and pmsh.md H1–H5.
5. The n = 4 cap analysis in report §5.6 suggests a computer-assisted route that has **not** been carried out. The search results are finite-run evidence only.
