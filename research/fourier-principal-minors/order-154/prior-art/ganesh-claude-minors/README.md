# Is P(154) already known? Prior-art report

Task `fourier-154-prior-art` · agent `ganesh-claude-minors` (GitHub owner g4nesh) · search cutoff **2026-11-27, 20:15 UTC**

**P(154)** is the statement that every principal minor of the unpermuted complex Fourier matrix F₁₅₄ = (ζ^(ij)), 0 ≤ i, j ≤ 153, with ζ a primitive 154th root of unity, is nonzero.

## Verdict

**No prior proof located within the documented search.**

- No source I found states P(154), lists 154 in a table or certificate, or proves a family that contains 154 or a square-free multiple of 154.
- Prior theorems do reduce P(154) to a finite computation. The lifting theorem of Caragea, Lee, Malikiosis and Pfander ([Theorem 1.4 of v3](https://arxiv.org/abs/2505.24326v3); it is Theorem 1.2 of [v1](https://arxiv.org/abs/2505.24326v1), 30 May 2025) gives P(154) as soon as every principal minor of F₁₄ is nonzero in characteristic 11, or every principal minor of F₂₂ is nonzero in characteristic 7. I found no source that reports either computation.
- A proof of P(154) by that route would therefore be a new finite verification fed into a known theorem. Gu, Zhou and Wang used exactly this pattern for orders 70 and 143 in August 2026. The method is prior art; the order-154 finite certificate is not in any source I located.

An absence of search results does not establish priority. The verdict is bounded by the search in [search-log.md](search-log.md) and by the limits listed under [Residual uncertainty](#not-accessed-and-residual-uncertainty).

**How sure.** I read the four papers named in the task page in full, and none proves P(154): that part is certain up to my reading of the arXiv HTML renderings. For the other located sources I read every numbered statement, and their stated hypotheses exclude 154 or leave a finite hypothesis open. Whether some source outside this search proves P(154) cannot be settled by searching.

## Substituting 154 into every candidate theorem

154 = 2 · 7 · 11, so the ways to write it as a prime times a cofactor are 2 × 77, 7 × 22 and 11 × 14. The arithmetic used below is reproduced by [arithmetic_checks.py](arithmetic_checks.py).

### Unconditional results

| # | Source and statement | Substituting 154 | Covers P(154)? |
|---|---|---|---|
| 1 | Caragea–Lee, [v4](https://arxiv.org/abs/2409.09793v4) Thm 1.1: for square-free N ≥ 4, every r × r principal minor is nonzero for r ∈ {2, 3, N−3, N−2} | 154 is square-free | Ranks 2, 3, 151, 152 only |
| 2 | CLMP, [v3](https://arxiv.org/abs/2505.24326v3) Thm 1.3 (from Cor 5.2 and Prop 5.4): N square-free of the form 2p, 3p, 5p, 6p or 7p with p prime | 154 = 2 · 77 = 7 · 22, and 77, 22 are not prime; 3, 5, 6 do not divide 154 | No |
| 3 | CLMP v3 Thm 1.5: N = p₁⋯p_k with p₁ < ⋯ < p_k and p_(j+1) > (P_j/2)^(P_j φ(P_j)/4), P_j = p₁⋯p_j | (2, 7, 11): j = 1 needs 7 > 1, true; j = 2 needs 11 > 7²¹ = 558 545 864 083 284 007, false | No |
| 4 | Two-prime and prime-size results: CLMP v3 Cor 5.1–5.3; Loukaki [v5](https://arxiv.org/abs/2412.08600v5) Thm A, Thm B, Cor 1, Cor 2; Zhou [v3](https://arxiv.org/abs/2505.01189v3) Cor 1; Emmrich–Kunis [v2](https://arxiv.org/abs/2506.02947v2) Thm 16, Cor 17 | 154 has three prime factors | No |
| 5 | Zhou v3 Cor 2: F₂^(⊗k) ⊗ F_q for odd prime q and k ≤ 3, after a permutation (identity for k = 1) | 154 is not of the form 2^k q | No |
| 6 | CLMP v3 Thm 5.7: every "arithmetic N′-principal" minor is nonzero, for all square-free N = pN′ | Applies to 154, but only to that subclass of minors (Definition 5.6) | Partial class only |
| 7 | Gu–Zhou–Wang, [v1](https://arxiv.org/abs/2608.17746v1) Thm 1.1: N ∈ {70, 143} | 154 is neither | No |
| 8 | Wang–Zhang, [v1](https://arxiv.org/abs/2608.22119v1) Thm 1.1: P(10p) for every prime p ∉ {2, 5, 11} | 10 does not divide 154 | No |
| 9 | Wang–Zhang v1 Cor 1.2: P(10pq) for p ∉ {2, 5, 11} and q > (5p)^(10p(p−1)) | 10 does not divide 154 | No |
| 10 | EPOCH, [v1](https://arxiv.org/abs/2610.06986v1) Supplement K.5: lifting conclusion for orders 11p with prime p outside {3, 11, 23, 67, 199, 419}, and 10p outside {2, 3, 5, 11, 31} | 154 = 11 · 14 and 14 is not prime | No |

### Conditional lifting theorems

Each theorem below has a finite hypothesis about the cofactor. I substituted all three factorizations into each.

**CLMP v3 Theorem 1.4.** Let N = pN′ be square-free with p prime, and suppose all principal minors of F_N′ are nonzero in characteristic p. Then every N′-principal minor of F_N, hence every principal minor, is nonzero. "Nonzero in characteristic p" means p does not divide the norm N_(ℚ(ζ_N′)/ℚ)(det F_N′[S, S]) for any S (v3 Section 3.2; the proof uses exactly p ∤ N(det Z_k) at equation (5.11)). So the minors must be nonzero at every prime ideal above p.

| (p, N′) | Structural hypotheses | Finite hypothesis | Status in the located literature |
|---|---|---|---|
| (11, 14) | hold | the 2¹⁴ principal minors of F₁₄ are nonzero at both primes above 11 (ord₁₄ 11 = 3, residue field 𝔽₁₃₃₁) | Not reported anywhere. CLMP Prop 5.5 gives it only for p > 7²¹ |
| (7, 22) | hold | the 2²² principal minors of F₂₂ are nonzero at the single prime above 7 (ord₂₂ 7 = 10, residue field of 7¹⁰ elements) | Not reported anywhere |
| (2, 77) | hold | the principal minors of F₇₇ are nonzero at the primes above 2 | **Fails.** CLMP Prop 4.5 computes the norm 2³ · 7³ for the principal minor {0, 1, 3} of F₇. The index set {0, 11, 33} gives a Galois conjugate of that minor inside F₇₇ (the computation in Wang–Zhang Prop 7.1), so a principal-minor norm of F₇₇ is even |

**Zhou v3 Theorem 1** (the same statement is in [v1](https://arxiv.org/abs/2505.01189v1), 2 May 2025). F_m ⊗ F_n is principally nonsingular if (i) the ideal ⟨1 − ζ_n⟩ is prime in ℤ[ζ_mn] with residue characteristic 0 or at least n, and (ii) F_m is principally nonsingular over that residue field (his statement allows any matrix F of powers of ζ_m; I take F = F_m). For prime n coprime to m the quotient is ℤ[ζ_m]/(n), so (i) holds exactly when n is inert in ℚ(ζ_m), that is ord_m(n) = φ(m). That reformulation is mine.

| (m, n) | Condition (i) | Condition (ii) |
|---|---|---|
| (14, 11) | fails: ord₁₄ 11 = 3, not 6 | not reached |
| (22, 7) | holds: ord₂₂ 7 = 10 = φ(22), residue characteristic 7 ≥ 7 | F₂₂ principally nonsingular over the field of 7¹⁰ elements: not reported anywhere |
| (77, 2) | fails: ord₇₇ 2 = 30, not 60 | not reached |

A result for F₂₂ ⊗ F₇ transfers to the unpermuted F₁₅₄ by CLMP v3 Lemma 3.4: when gcd(m, n) = 1, the principal minors of F_m ⊗ F_n are Galois conjugates of those of F_mn up to sign.

**Wang–Zhang v1 Theorem 4.3** (active-rank norm budget). If M is square-free, P(M) holds, q ∤ M is prime and β(M, q) < φ(M), then P(Mq). P(14), P(22) and P(77) hold by CLMP Cor 5.2. The budgets β(14, 11) and β(22, 7) need the q-adic valuations of the principal-minor norms of F₁₄ and F₂₂; Wang–Zhang compute such a spectrum only for base 10 (Prop 5.1, Table 2). At (77, 2) the hypothesis fails: the F₇ minor above contributes valuation at least 30 at ranks 3 and 74, so β(77, 2) ≥ 60 = φ(77). That last inequality is my derivation from the cited norm. Their Lemma 4.1 (one safe prime ideal per occupancy flag) is the only located statement that could still use 2 as the multiplier, and it needs prime-ideal data for F₇₇ that no source reports.

### Indirect routes

- **Upward inheritance.** Wang–Zhang Prop 7.1: a zero principal minor of F_d gives one in F_dm when gcd(d, m) = 1. So P(154m) for one coprime m would imply P(154). No located theorem proves P(N) for a square-free multiple N of 154. Rows 2, 4, 5, 7, 8 and 10 contain no multiple of 154: their orders have at most two distinct prime factors, or have the form 6p or 10p, or are 70 and 143. Row 9 would need 770 = 10 · 7 · 11 with 11 > 35⁴²⁰, or p = 11, which is excluded. A Theorem 1.5 tower through 2, 7 and 11 has 11 directly after 7 with P_j ≥ 14, which again needs at least 11 > 7²¹.
- **All but finitely many multipliers.** CLMP (Section 1, after Thm 1.4) and Wang–Zhang (Section 7.1) note that P(Mq) holds for every prime q outside the support of the norm product of the principal minors of F_M. For M = 14 that support is not computed in any located source, so the remark does not decide q = 11.
- **Support-size uncertainty principles.** F_N[K, K] is singular exactly when some nonzero f is supported in K with Fourier transform supported in the complement of K. A bound on support sizes alone (Donoho–Stark, Tao, Meshulam) could exclude this only by forcing |supp f| + |supp f̂| > N, and for composite N the indicator of a subgroup of order d has |supp f| + |supp f̂| = d + N/d ≤ N. I therefore did not treat support-size theorems as candidates. Zhou ([2608.06335v2](https://arxiv.org/abs/2608.06335v2)) discusses this disjoint-support formulation and keeps the square-free statement as a conjecture (Conjectures 2 and 3).
- **Direct numerics.** CLMP v3 (Section 1) says brute-force checks reach about N = 45; the George Mason poster reaches N = 30. F₁₅₄ has 2¹⁵⁴ principal minors.

## Sources inspected

All reading was done on 2026-11-27 through arXiv's HTML renderings unless noted. "Full" means every section; "statements" means the abstract and every numbered theorem, corollary, proposition and conjecture, without the proofs.

| Source | Version read (date) | Earlier versions | Depth | Bearing on 154 |
|---|---|---|---|---|
| Caragea, Lee, *On the principal minors of Fourier matrices*, [arXiv:2409.09793](https://arxiv.org/abs/2409.09793v4) | v4 (20 Apr 2025) | v1 15 Sep 2024 | Full | Conjecture 1.3(b); ranks 2, 3 and complements |
| Caragea, Lee, Malikiosis, Pfander, *Principal minors of Fourier matrices of square-free order*, [arXiv:2505.24326](https://arxiv.org/abs/2505.24326v3) | v3 (7 Dec 2025) | v1 30 May 2025 and v2 30 Jul 2025 (statements read; neither mentions 154 or 14p) | Full | Lifting theorem; families 2p, 3p, 5p, 6p, 7p; says 70 is the smallest unresolved three-prime order |
| Gu, Zhou, Wang, *Principal nonsingularity of the Fourier matrices of orders 70 and 143*, [arXiv:2608.17746](https://arxiv.org/abs/2608.17746v1) | v1 (18 Aug 2026) | none | Full | Orders 70 and 143 only |
| Wang, Zhang, *Prime multipliers of order ten*, [arXiv:2608.22119](https://arxiv.org/abs/2608.22119v1) | v1 (22 Aug 2026) | none | Full | Base 10 only; P(110) left open |
| Guo, Mo, Li, Deng, *EPOCH: Reliable Discovery through Evidence-Governed Search*, [arXiv:2610.06986](https://arxiv.org/abs/2610.06986v1) | v1 (4 Oct 2026) | none | Supplement K.5 only, the one passage a text search for Fourier, principal minor, Chebotarev or 154 finds | Bases 10 and 11 |
| Loukaki, *Chebotarev's theorem for groups of order pq and an uncertainty principle*, [arXiv:2412.08600](https://arxiv.org/abs/2412.08600v5) | v5 (29 Aug 2025) | v1 11 Dec 2024 | Statements | n = 2p and n = pr |
| Zhou, *Principal non-singularity of Fourier matrices on ℤ_p × ℤ_q and ℤ₂^k × ℤ_q*, [arXiv:2505.01189](https://arxiv.org/abs/2505.01189v3) | v3 (27 Feb 2026) | v1 2 May 2025 (statements read) | Statements | Conditional Theorem 1; two-prime corollaries |
| Emmrich, Kunis, *Real and finite field versions of Chebotarev's theorem*, [arXiv:2506.02947](https://arxiv.org/abs/2506.02947v2) | v2 (30 Sep 2025) | v1 3 Jun 2025 | Statements | Prime size only |
| Zhou, *On number of cyclic n-roots and disjointness of Fourier supports*, [arXiv:2608.06335](https://arxiv.org/abs/2608.06335v2) | v2 (12 Aug 2026) | v1 6 Aug 2026 | Statements | States the conjecture; no composite order proved |
| Cabrelli, Molter, Negreira, *Weaving Riesz bases*, [arXiv:2404.02229](https://arxiv.org/abs/2404.02229) | abstract (v1 2 Apr 2024) | | Abstract; scope as described in Caragea–Lee v4 | First explicit statement of the conjecture; 2 × 2 minors |
| Brennan, Tuttle, Wigginton (dir. Walnut), *Principal Minors of the Fourier Matrix*, George Mason MEGL poster, Fall 2023 | PDF | | Full (text extracted from the PDF) | Numerics to N = 30; N = pq stated as a conjecture |
| Caragea, Pfander, Lee, Malikiosis, *Towards a conjecture concerning minors of Fourier matrices*, SampTA 2025, doi:10.1109/SampTA64769.2025.11133507 | IEEE Xplore abstract | | Abstract only | Permutation conjecture, false at N = 16; "partial results and numerics" |
| Malavolta, Rosen, [arXiv:2606.08662](https://arxiv.org/abs/2606.08662v1); Narayanan, [arXiv:2607.15848](https://arxiv.org/abs/2607.15848v1) | v1 | | Abstract | Citers of CLMP; prime order only |

Sources that searches surfaced and I judged out of scope from the title or abstract alone: arXiv:2511.08504 (OTFS channel estimation, a citer), arXiv:2310.09992, arXiv:1807.07648, arXiv:math/0611493, arXiv:1910.05652, arXiv:2610.20653, arXiv:2610.24518, arXiv:2510.27670, and the support-size uncertainty papers returned by web query 18.

## Dates and status

| Item | Earliest public date found | Latest version checked | Status |
|---|---|---|---|
| The conjecture for square-free N | 2 Apr 2024 (Cabrelli–Molter–Negreira v1, per Caragea–Lee v4); 15 Sep 2024 (Caragea–Lee v1) | Caragea–Lee v4, 20 Apr 2025 | Cabrelli–Molter–Negreira is published (J. Fourier Anal. Appl. 31 (2025), per the later papers' reference lists). Caragea–Lee is cited as Commun. Korean Math. Soc. (2026), online first, by Gu–Zhou–Wang. I opened neither journal version |
| Lifting theorem, inert case | 2 May 2025 (Zhou v1, Thm 1) | Zhou v3, 27 Feb 2026 | Published as Demonstr. Math. 59 (2026), doi:10.1515/dema-2025-0242, per Google Scholar and Semantic Scholar. Journal version not opened |
| Lifting theorem, general | 30 May 2025 (CLMP v1, Thm 1.2) | CLMP v3, 7 Dec 2025 | arXiv preprint. The arXiv record showed no journal reference on the access date |
| First finite certificate fed into the lifting theorem for a three-prime order | 18 Aug 2026 (Gu–Zhou–Wang v1, order 70) | v1 | arXiv preprint, AI-assisted by its own acknowledgment; no peer-reviewed version located. The verifier scripts are said to ship with the TeX source; arXiv lists no ancillary files and I did not open the source bundle |
| Fixed-base norm spectra | 22 Aug 2026 (Wang–Zhang v1, base 10); 4 Oct 2026 (EPOCH v1, bases 10 and 11) | v1 of each | arXiv preprints, AI-assisted (Wang–Zhang, by their acknowledgment) or reported from an AI research agent (EPOCH); unreviewed |
| Two-prime results | 11 Dec 2024 (Loukaki v1) | v5, 29 Aug 2025 | Published: Bull. Lond. Math. Soc. 57 (2025), doi:10.1112/blms.70192, per the arXiv record |
| **P(154)** | none found | | No statement, table entry or certificate located up to the cutoff |

## Not accessed and residual uncertainty

- **Abstract only.** The SampTA 2025 paper: the OpenReview copy sits behind a browser verification page, and on IEEE Xplore I read the abstract only. Its abstract concerns the permutation conjecture. The same authors' later v3 (Dec 2025) lists only 2p, 3p, 5p, 6p, 7p and calls 70 unresolved, so I infer the conference paper does not prove P(154); I did not confirm that from its text.
- **Not opened.** The journal versions named above; the body of the EPOCH paper outside Supplement K.5; the TeX-source bundles of Gu–Zhou–Wang and Wang–Zhang; Caragea's 2024 thesis and the Lee–Pfander–Walnut paper (both cited for conjectures only). CLMP announce a second manuscript on non-square-free orders as "in preparation"; I found no posted version.
- **Blocked services.** zbMATH Open returned a bot-verification page. MathSciNet was not tried (subscription). The doi.org resolver could not be opened from this session. Semantic Scholar rate-limited the citation list of arXiv:2608.17746.
- **Index lag.** The EPOCH paper appears in neither citation index and surfaced only through web search. Recent or oddly worded papers can be missed the same way. arXiv's own search covers metadata, not full text.
- **Outside any index.** Talks, theses, private computations, code repositories and unposted manuscripts. A web search for a public order-154 verifier or certificate found none.
- **Depth of reading.** For Loukaki, Zhou and Emmrich–Kunis I read the statements, not the proofs. I relied on arXiv HTML renderings rather than the PDFs; where the text extraction of the Theorem 1.5 bound was garbled I confirmed it from the proof and from the LaTeX of v1 and v2.
- **Soundness of the cited papers is not assessed here.** This report asks what the literature claims, not whether those claims are correct. The August 2026 preprints in particular are unreviewed.

## Leads for the 14p audit

1. **What the literature already gives for 14p.** P(42) by CLMP Prop 5.4 (42 = 6 · 7). P(70) by Gu–Zhou–Wang Thm 1.1 and by Wang–Zhang Thm 1.1 (70 = 10 · 7). P(14p) for every prime p > 7²¹ by CLMP Thm 1.5 with (2, 7, p). 28 and 98 are not square-free. I located nothing for primes 11 ≤ p < 7²¹.
2. **The template is published.** Wang–Zhang carry out the full programme for base 10: norm spectrum of the base, ordinary lifting for primes outside its support, prime-ideal analysis for the exceptional ones (Lemma 4.1, Thm 4.3). A base-14 classification would follow that template, so its novelty would lie in the base-14 spectrum and the treatment of its exceptional primes, not in the method.
3. **Small multipliers can be exceptional.** For base 10 the support is {2, 3, 5, 11, 31}; 31 was rescued by a safe prime ideal and P(110) was left open (Wang–Zhang Props 6.1, 6.2). Expect the same three outcomes for base 14: safe, rescued, undecided.
4. **Other presentations of 14p.** Multiplier 7 over base 2p needs F_2p in characteristic 7. Multiplier 2 over base 7p fails for ordinary lifting for every p, by the F₇ minor with norm 2³ · 7³. Zhou's Theorem 1 with multiplier p applies only when p is inert in ℚ(ζ₁₄), that is p ≡ 3 or 5 (mod 7). Wang–Zhang (Section 7.2) note that reordering the factors can remove an obstruction.
5. **Upward inheritance cuts both ways.** A zero principal minor at any order d dividing 14p (coprime cofactor) would propagate to 14p, and any future theorem for a square-free multiple of 14p would give P(14p) (Wang–Zhang Prop 7.1).
6. **Priority risk is high and current.** Three AI-assisted preprints that apply the lifting theorem to new orders or bases appeared between 18 Aug and 4 Oct 2026. Recheck arXiv (math.NT, math.FA, cs.AI) and the citers of arXiv:2505.24326 immediately before any 14p claim.

## Independence

I did not open `order-154/seed.md`, `order-154/verify_seed.py`, `order-154/replay.md`, anything under `order-154/replays/` or `family-14p/` before writing this verdict (I saw directory and file names in a listing, not contents). Two things reached me through the commands the task page requires: the title of commit `82e1861` on `main` in the git log, and the one-line abridged finding of `fourier-154-replay` printed by `context --objective`. Both point to CLMP Theorem 1.4 with (p, N′) = (11, 14), and the second says the replay found that finite hypothesis to hold over 𝔽₁₃₃₁. I did not use them as evidence. The substitution tables above cover all three factorizations from the papers alone, and the verdict concerns what the literature contains, which that exposure does not change.

The inferences marked as mine (the inertness reformulation of Zhou's condition, the failure at (2, 77), the bound β(77, 2) ≥ 60, the remark on support-size principles) are short derivations from cited statements, not literature claims, and deserve the reviewer's check.
