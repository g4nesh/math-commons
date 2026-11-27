# Independent review: fourier-154-prior-art

- Reviewer agent / GitHub owner: vn-claude (Claude Code worker, model Claude Sonnet 5.5) / vignesh-nagarajan-vn
- Contributor owner(s): g4nesh (author agent ganesh-claude-minors)
- Exact checkpoint ID: `7fb0f261-bff0-4fbf-99c0-70be751a944a`
- Artifact commit(s): `092a7cecb7566d3a4f21101e7871631461aa3305` (adds only `order-154/prior-art/ganesh-claude-minors/README.md`, `search-log.md`, `arithmetic_checks.py`; the folder is unchanged on the working tree I read)
- Verdict: pass (with optional wording corrections listed under Findings)
- Date: November 27, 2026

## Claim and evidence level

The checkpoint claims: no prior proof of P(154) (every principal minor of the unpermuted complex F_154 is nonzero) was located within a documented search with cutoff 2026-11-27 20:25 UTC. It further claims that P(154) is implied only conditionally, through CLMP (arXiv 2505.24326v3) Thm 1.4 at (p, N') = (11, 14) or (7, 22) and Zhou (arXiv 2505.01189v3) Thm 1 at (m, n) = (22, 7), each needing a finite-characteristic check of F_14 or F_22 that no located source reports, while (2, 77) fails because the F_7 minor {0,1,3} has norm 2^3 * 7^3.

Evidence level: a bounded literature search (my own plus the author's log) and exact integer arithmetic. It is not a proof that no source exists, and the report says so. The claim as worded ("no prior proof located within the documented search", category 3 of 4) is correctly bounded.

## Checks performed

All access was on 2026-11-27 through the WebFetch and WebSearch tools of the Claude Code session. The arXiv HTML pages were summarized by the fetch tool's small model, so my theorem quotes are tool-extracted text, not my own reading of the PDFs. Where that matters I say so.

### Own search log

| # | Service | Query or page | Result |
|---|---|---|---|
| 1 | Web search (standard) | principal minors Fourier matrix square-free order 154 nonzero | Caragea-Lee, CLMP, Gu-Zhou-Wang, Wang-Zhang, GMU poster. Nothing for 154 |
| 2 | Web search (standard) | Chebotarev theorem composite order Fourier matrix principal minors three primes 2026 | Same papers plus Zhou, Emmrich-Kunis, Loukaki, an uncertainty-principle paper. No three-prime result for 154 |
| 3 | Web search (extended) | "Fourier matrix" "principal minors" order 154 OR "14p" OR "2·7·11" arXiv November 2026 | Same set; nothing on 154 or 14p. The tool reported it had simplified the operators |
| 4 | Web search (extended) | arXiv 2026 principal minors Fourier matrices 14 11 characteristic lifting theorem new orders 154 105 130 110 | Same set plus EPOCH 2610.06986 and one new hit, 2611.25715 (row 8). Nothing on 154 |
| 5 | arXiv listing search | "principal minors" Fourier, newest first | 6 hits: 2610.20653 (unrelated), 2608.22119, 2608.17746, 2505.24326, 2505.01189, 2409.09793 |
| 6 | arXiv listing search | Chebotarev Fourier minors | 3 hits: 2506.02947, 2505.24326, 2505.01189 |
| 7 | Google Scholar | "principal minors" "Fourier matrix" 154 | Redirected to a sync page; not completed (blocked) |
| 8 | arXiv abstract | 2611.25715 (Deng, Zhang, 19 Nov 2026, "Real Principal-Minor Detectors") | Surfaced by query 4. Concerns symmetric matrices over F_2 and real detectors; not Fourier matrices, not 154. Out of scope. Not in the author's log |
| 9 | arXiv abstract | 2505.24326 | v1 30 May 2025, v2 30 Jul 2025, v3 7 Dec 2025; no journal reference listed |

Not done: zbMATH, MathSciNet, Semantic Scholar citation API, Google Scholar cited-by lists, IEEE/SampTA full text, journal versions, TeX source bundles. Nothing I ran contradicts the author's search log. The author ran more queries than I did (19 web, 9 arXiv, 7 Scholar, 5 Semantic Scholar); I treat that log as additional bounded coverage that I did not re-execute.

### Cited theorems: number and hypotheses

| Cited statement | What I opened | My verification | Covers 154? |
|---|---|---|---|
| CLMP v3 Thm 1.3 (2p, 3p, 5p, 6p, 7p with p prime) | arxiv.org/html/2505.24326v3 | Number and text confirmed. 154 = 2*77 = 7*22 with 77 and 22 not prime; 3, 5, 6 do not divide 154 | No |
| CLMP v3 Thm 1.4 (lifting) | same | Number confirmed. Hypotheses: N = pN' square-free, p prime, all principal minors of F_N' nonzero in characteristic p. Conclusion: every N'-principal (hence every principal) minor of F_N nonzero. The replay's raw-LaTeX check (`replays/claude-opus-5-5-fourier-154-replay/source-check.md`, section d) shows the proof needs p not dividing the norm of each det Z_k, so all primes above p. All three factorizations of 154 satisfy the structural hypotheses | Only conditionally |
| CLMP v3 Thm 1.5 (separated primes) | same | Number and bound confirmed: p_(j+1) > (P_j/2)^(P_j phi(P_j)/4). For (2, 7, 11), j = 2 needs 11 > 7^21 = 558545864083284007 (I re-ran the arithmetic): false | No |
| CLMP v3 Cor 5.2 (N = pq, p in {2,3,5,7}) | same | Confirmed. Gives P(14), P(22), P(77) over C. Not a three-prime statement | No |
| CLMP v3 Prop 5.4 (6p) | same | Confirmed; 154 is not 6p | No |
| CLMP v3 Prop 4.5 and the {0,1,3} minor of F_7 | same | Fetch tool reported norm 2^3 * 7^3, divisible by 2 and 7 only. I re-ran the exact computation in `arithmetic_checks.py`: 2744 = 2^3 * 7^3, true | n/a |
| Caragea-Lee v4 Thm 1.1 (ranks 2, 3, N-3, N-2) | Not opened by me; only consistent with search summaries and CLMP's statement of Conjecture 1.2 | Not independently verified from the paper | Ranks only |
| Zhou v3 Thm 1, Cor 1, Cor 2 | arxiv.org/html/2505.01189v3 | Thm 1 confirmed: F tensor F_n principally nonsingular if (i) the ideal <1 - zeta_n> is prime in Z[zeta_mn] with residue characteristic 0 or at least n, and (ii) f(F) is principally nonsingular in that quotient. Cor 1 needs distinct odd primes p, q with q generating Z_p^* and q > Gamma_p; Cor 2 is F_2^(tensor k) tensor F_q, k in {1,2,3}, after a permutation. 154 matches neither. The author's inertness reformulation (condition (i) iff ord_m(n) = phi(m), n prime coprime to m) is standard and I accept it; the order values ord_14(11)=3, ord_22(7)=10=phi(22), ord_77(2)=30 re-ran correctly | Only conditionally, at (22, 7) |
| Loukaki v5 Thm A, B, Cor 1, 2 | arxiv.org/html/2412.08600v5 | Confirmed: Thm A is n = 2p; Thm B is n = pr with p primitive mod r and p > Gamma_r; Cor 1 and 2 are the principal-submatrix versions. Loukaki remarks that the inductive proof extends to more than two primes, but her application stays at two (per the fetch tool); the report does not record that remark, see Findings | No |
| Gu-Zhou-Wang v1 Thm 1.1 (N in {70, 143}) | arxiv.org/html/2608.17746v1 | Confirmed. Prop 2.1 there restates CLMP Thm 1.4 via the norm condition. Neither 154, 14, 22, 77 nor 11p appears | No |
| Wang-Zhang v1 Thm 1.1 (10p, p not in {2,5,11}), Cor 1.2, Thm 4.3, Prop 7.1 | arxiv.org/html/2608.22119v1 | Thm 1.1, Thm 4.3 and Prop 7.1 are as the report says. Cor 1.2 reads "all but finitely many primes q" in my fetch; I could not confirm the explicit bound q > (5p)^(10p(p-1)). No mention of 154 or base 14 | No |
| EPOCH 2610.06986v1 Supplement K.5 | arxiv.org/html/2610.06986v1 | Title "Fourier principal nonsingularity for orders 10p and 11p"; exception supports {3,11,23,67,199,419} for 11p and {2,3,5,11,31} for 10p; orders 187 and 130 named as included. 14p and 154 not mentioned | No |

### Factorization logic

154 = 2 * 7 * 11 has exactly three prime-times-cofactor splittings: 11 x 14, 7 x 22, 2 x 77. I re-derived each.

- (p, N') = (11, 14): the structural hypotheses of Thm 1.4 hold. The finite hypothesis is that the principal minors of F_14 are nonzero at both primes above 11 (ord_14(11) = 3, residue field of 1331 elements). Not reported in the literature, apart from CLMP's size-based Prop 5.5 bound.
- (7, 22): a single prime above 7 (ord_22(7) = 10). Not reported in the literature.
- (2, 77): the F_7 minor {0,1,3} has norm 2^3 * 7^3. In Z_77 the index set {0, 11, 33} has entries zeta_77^(121 ab) = zeta_7^(4ab) for a, b in {0,1,3}, a Galois conjugate of that F_7 minor, so F_77 has a principal minor with even norm and the hypothesis of Thm 1.4 fails at p = 2. I checked this by hand. It closes this route; it does not say P(154) is false.
- Upward inheritance (Wang-Zhang Prop 7.1, confirmed): a zero minor at order d propagates to dm when gcd(d, m) = 1. The divisors 2, 7, 11, 14, 22, 77 all satisfy P over C (CLMP Cor 5.2 and Chebotarev), so nothing is inherited downward onto 154. No located theorem proves P(N) for a square-free multiple of 154, so no implication runs the other way either.

Within what I found, no known reduction from a covered smaller order implies 154.

### Reproduction

I read `arithmetic_checks.py` (standard library only, short) and ran it. Observed: ord and phi table for (14, 11), (22, 7), (77, 2); Thm 1.5 bound 558545864083284007, so "11 > bound" is False; norm of the F_7 {0,1,3} minor 2744 = 2^3 * 7^3, True; beta(77, 2) >= 60 = phi(77). All as the report states.

I also read order-154 `seed.md`, `replay.md` and the replay's `source-check.md` as background. I did not rely on them for the verdict, and the replay's finite-field result is outside this review's scope.

## Independent adversarial attempt

I looked for a source that proves P(154) directly or through a hypothesis that 154 satisfies: the 2p, 3p, 5p, 6p, 7p families; N = pq results; the Thm 1.5 separated-primes tower; Zhou's inert-prime theorem and its 2^k * q corollary; orders 70 and 143; the base-ten 10p and 10pq family; EPOCH 10p and 11p; arithmetic N'-principal minors; and general uncertainty-principle results. I found none that covers 154. The closest implicit routes are the conditional ones above, and each needs a computation (F_14 mod 11 or F_22 mod 7) that no source I opened reports.

Limits: I did not search zbMATH or MathSciNet, did not read SampTA 2025 full text, journal versions, TeX sources or the D. G. Lee preprints, and my Google Scholar attempt failed. A paper posted after about 20:25 UTC on the cutoff day, or any unindexed talk, thesis or repository, would not be seen. A literature search can only ever be bounded.

## Findings

No fatal gap. Items for the author or synthesis, none blocking.

1. **Repairable, wording (Prop 5.5 attribution).** The README says CLMP Prop 5.5 gives the (11, 14) hypothesis "only for p > 7^21". Prop 5.5 as extracted says p > m^(m phi(N)/2) where m is the minor size, so for N = 14 the bound is 7^21 at m = 7 but 14^42 at m = 14, and no single threshold is stated there. The 7^21 figure matches the Thm 1.5 bound at P_2 = 14, which I confirmed, and presumably uses complementary-minor symmetry to restrict to m <= 7. The report should say that or cite Thm 1.5. The conclusion that 11 lies below every such threshold is unaffected.
2. **Repairable, wording (Prop 7.1 attribution).** The report says the index-set {0, 11, 33} computation is in Wang-Zhang Prop 7.1. That proposition concerns a zero complex minor propagating from F_d to F_dm; it does not mention F_7 mod 2 or order 77. The embedding argument is the same and the conclusion is correct, but it is the author's own derivation, as the report's closing paragraph partly says.
3. **Optional, completeness.** Loukaki's remark that her inductive proof extends beyond two primes is not recorded. A sentence noting that it still yields no result for 154 would help synthesis.
4. **Optional, completeness.** The report cites CLMP v1 as Thm 1.2 for the lifting statement. I did not open v1 and could not confirm that numbering. The v3 numbering (Thm 1.4) is confirmed.
5. **Optional, log.** arXiv 2611.25715 (19 Nov 2026) appears in my search but not the author's log. It is unrelated, so this is a bookkeeping gap, not a missed prior-art item.
6. **Not independently checked.** Journal-status claims (Commun. Korean Math. Soc., Demonstr. Math. 59, Bull. Lond. Math. Soc. 57, J. Fourier Anal. Appl. 31); the GMU poster contents; the SampTA abstract; the EPOCH body beyond K.5; Wang-Zhang's explicit bound in Cor 1.2. CLMP's "approximately N = 45" brute-force sentence is confirmed.
7. The report discloses its exposure to the replay's one-line summary and the commit title 82e1861. That is handled adequately, since the verdict concerns literature content, which that exposure does not change.

## Verdict rationale

Pass. The checkpoint's claim is correctly bounded: "no prior proof located within the documented search", with cutoff, services, queries and limits stated, and an explicit statement that absence of results does not establish priority. Every theorem I could open carries the number and hypotheses the report gives, apart from the two attribution slips above. For each of the cofactors 14, 22, 77 the substitution is correct, and I re-derived the one definitive negative (route (2, 77)). My independent search found no source covering 154 directly or through a theorem whose hypotheses 154 satisfies, apart from the conditional CLMP Thm 1.4 and Zhou Thm 1 routes the report already identifies.

What this pass establishes: the report's conclusion is consistent with an independent bounded search and with the cited statements. What it does not establish: that no proof of P(154) exists anywhere; that the cited papers are correct (several are unreviewed AI-assisted preprints); or the truth of P(154), which belongs to the replay and synthesis tasks. I recommend that the wording fixes in findings 1 and 2 be made in synthesis or a follow-up commit; they do not change the verdict.

Reviewer: agent vn-claude (Claude Code worker, model Claude Sonnet 5.5), owner vignesh-nagarajan-vn, November 27, 2026.
