# Search log for the P(154) prior-art report

All searches and readings were made on **2026-11-27** between about 20:00 and 20:15 UTC. The verdict and the reading of each source are in [README.md](README.md).

No query below returned a source that states P(154), tabulates 154, or reports a finite-characteristic computation for F₁₄, F₂₂ or F₇₇. An absence of results does not establish priority.

## Services

| Tag | Service | What it covers | Limits |
|---|---|---|---|
| W | Web search tool of the Claude Code session ("standard" or "extended" mode) | General web index | US index; returns about ten links and a machine summary per query. I used the links and read the primary pages myself |
| A | arXiv listing search, `arxiv.org/search`, all fields, newest first | Titles, abstracts, authors, comments | Metadata only, no full text |
| G | Google Scholar | Full text of indexed papers; cited-by lists | Coverage and ranking are opaque; recent papers lag |
| S | Semantic Scholar Graph API, `/paper/arXiv:ID/citations` | Citing papers | Rate-limited without a key; recent papers lag |
| R | Direct reading in a browser | arXiv abstract and HTML pages, IEEE Xplore abstract page, one PDF | HTML renderings, not PDFs |

## Web search (W)

| # | Query | Mode | Relevant results |
|---|---|---|---|
| 1 | Caragea Lee conjecture principal minors Fourier matrix square-free order nonzero | standard | The four papers named in the task page |
| 2 | "principal minors" "Fourier matrix" order 154 | standard | Same papers and the George Mason poster; nothing on 154 |
| 3 | "principal nonsingularity" Fourier matrix square-free order | standard | Adds Zhou arXiv:2505.01189 and Loukaki arXiv:2412.08600 |
| 4 | Fourier matrix principal minors nonzero "14p" prime multiplier square-free lifting | extended | Adds EPOCH arXiv:2610.06986 and the SampTA 2025 abstract on OpenReview; nothing on 14p |
| 5 | "Principal nonsingularity of the Fourier matrices of orders" 70 143 cited follow-up 2026 | extended | Wang–Zhang, EPOCH, a third-party summary page; no further follow-up |
| 6 | Fourier matrix principal minors square-free conjecture 2026 arXiv "110" OR "154" OR "105" OR "130" three primes proved | extended | Nothing on 105, 110, 130 or 154 beyond Wang–Zhang leaving 110 open |
| 7 | "principal submatrices" Fourier matrix invertible square-free N Chebotarev composite order uncertainty principle cyclic group | standard | Adds Zhou arXiv:2608.06335 |
| 8 | "Fourier matrix" principal minors "2·7·11" OR "2 \cdot 7 \cdot 11" OR "14 \cdot 11" OR "N=154" Chebotarev | extended | No page with those strings; adds arXiv:1807.07648 (out of scope) |
| 9 | Caragea Lee Malikiosis Pfander "Principal minors of Fourier matrices of square-free order" journal published | standard | arXiv record only; no journal version found |
| 10 | mathoverflow OR math.stackexchange "DFT matrix" OR "Fourier matrix" "principal minors" nonzero square-free | standard | No forum thread; adds arXiv:math/0611493 (out of scope) |
| 11 | github Fourier matrix principal minors verifier finite field certificate "Verify_F70_Principal_Minors" OR "F154" OR "order 154" OR "P(154)" | extended | No repository or certificate for 154 |
| 12 | "Fourier matrix" "principal minors" "characteristic" norm spectrum "order 14" OR "F_{14}" OR "base 14" OR "14p" cyclotomic | extended | Nothing for base 14 |
| 13 | Chebotarëv theorem composite order "square-free" Fourier minors nonvanishing new result October OR November 2026 | extended | Nothing dated after August 2026 |
| 14 | Caragea Lee "On the principal minors of Fourier matrices" Communications of the Korean Mathematical Society CKMS c250267 | standard | Journal page not found |
| 15 | Zhou "Principal non-singularity of Fourier matrices" Demonstratio Mathematica 2026 dema-2025-0242 | standard | Journal page not found; arXiv v3 only |
| 16 | "Fourier matrix" OR "DFT matrix" "all principal minors" nonzero "square-free" OR "squarefree" -arxiv 154 OR "three primes" uncertainty "support" disjoint | standard | Nothing new; adds arXiv:1910.05652 (out of scope) |
| 17 | "Fourier matrix" principal minors "22 \times 7" OR "22×7" OR "77×2" OR "order 22" OR "F_{22}" "characteristic 7" square-free lifting | standard | Nothing for F₂₂ in characteristic 7 or F₇₇ in characteristic 2 |
| 18 | discrete uncertainty principle cyclic group square-free order function and Fourier transform disjoint supports "principal minors" theorem three distinct primes | standard | Uncertainty-principle papers on support sizes; no three-prime theorem |
| 19 | Chebotarev cyclic group Fourier matrix "product of three primes" OR "three distinct primes" principal minors nonzero proved Fourier "2pq" OR "pqr" | extended | No three-prime result beyond 6p, 70 and the 10p family |

## arXiv listing search (A)

| Query | Results | Records |
|---|---|---|
| "principal minors" Fourier | 6 | 2409.09793, 2505.01189, 2505.24326, 2608.17746, 2608.22119, and 2610.20653 (unrelated) |
| "principal nonsingularity" OR "principal non-singularity" OR "principally non-singular" | 2 | 2505.01189, 2608.17746 |
| Chebotarev "square-free" Fourier minors | 1 | 2505.24326 |
| "Fourier matrix" minors | 8 | The above plus 2506.02947, 2310.09992 (unrelated), 1807.07648 (unrelated) |
| "Fourier matrices" "square-free" | 2 | 2409.09793, 2505.24326 |
| Chebotarev "roots of unity" minors | 3 | 2506.02947, 2505.24326, math/0312398 |
| Fourier "disjoint supports" cyclic | 1 | 2608.06335 |
| Cabrelli Molter Negreira weaving Riesz bases | 1 | 2404.02229 |
| "principal submatrices" Fourier OR "principal submatrix" Fourier | none returned | The page showed no result list; not counted as a search |

## Google Scholar (G)

| Query | Results |
|---|---|
| "principal minors" "Fourier matrix" "154", since 2024 | 2, both unrelated (arXiv:2610.24518, arXiv:2510.27670) |
| "principal minors" "Fourier matrices" "square-free", since 2024 | 14: 2409.09793, 2505.24326 (and its v1 title "Chebotarev type results for composite group order Part 1"), 2608.17746, 2505.01189 and its journal version, SampTA 2025, 2412.08600, 2608.22119, 2608.06335, 2506.02947, 2606.08662, 2607.15848, 2511.08504 |
| Cited-by list of "Principal minors of Fourier matrices of square-free order" | 6: 2511.08504, Loukaki (Bull. Lond. Math. Soc.), 2608.17746, 2608.06335, 2606.08662, 2608.22119 |
| "Principal nonsingularity of the Fourier matrices of orders" OR "Prime multipliers of order ten" | 1 record (2608.17746), cited by 1; I did not open that citing record |
| "principal nonsingularity" OR "principal non-singularity" Fourier, since 2024 | About 52; the first page has the four known papers and unrelated Demonstratio Mathematica articles. Later pages not read |
| "Fourier" "principal minors" "lifting theorem" "characteristic" "square-free" 14, since 2025 | 2: 2608.17746, 2608.22119 |

## Semantic Scholar citations (S)

| Paper | Citing papers returned |
|---|---|
| arXiv:2505.24326 | 7: 2608.22119, 2608.17746, 2608.06335, 2607.15848, 2606.08662, 2511.08504, 2412.08600 |
| arXiv:2409.09793 | 8: 2608.22119, 2608.17746, 2608.06335, Zhou (Demonstr. Math., doi:10.1515/dema-2025-0242), SampTA 2025 (doi:10.1109/SampTA64769.2025.11133507), 2506.02947, 2505.01189, 2412.08600 |
| arXiv:2505.01189 | 1: its own journal version |
| arXiv:2608.22119 | 0 |
| arXiv:2608.17746 | No data: rate-limited on three attempts |

## Direct reading (R)

The version and depth of reading for each source are in the "Sources inspected" table of [README.md](README.md). Pages opened: arXiv abstract pages for every arXiv identifier listed there (version histories and dates come from those pages); arXiv HTML renderings of 2409.09793v4, 2505.24326v3, 2505.24326v1, 2608.17746v1, 2608.22119v1, 2610.06986v1, 2412.08600v5, 2505.01189v3, 2505.01189v1, 2506.02947v2 and 2608.06335v2; the ancillary-file pages of 2608.17746v1 and 2608.22119v1 (both report no ancillary files); the IEEE Xplore abstract page of document 11133507; the poster PDF on `megl.science.gmu.edu`.

## Could not be opened

| Target | What happened |
|---|---|
| `openreview.net` copy of the SampTA 2025 paper (id `becyRjuIVt`) | Browser verification page, in the browser and through a fetch tool. Not bypassed |
| `zbmath.org` search | Bot-verification page. Not bypassed |
| `doi.org` links for the Commun. Korean Math. Soc. and Demonstr. Math. versions | The resolver could not be opened from this session |
| MathSciNet | Not tried; subscription service |
