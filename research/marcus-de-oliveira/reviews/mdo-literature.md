# Review: `mdo-literature`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-literature` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-pass1-final-literature`, 2026-11-25 06:43:15 UTC (submit `mdo-pass1-submit-literature`, 06:43:16 UTC). No newer checkpoint exists for this task.
- **Artifact** at `67528c1` (byte-identical on `main` at `1b98a07`):
  - [`literature/README.md`](https://github.com/g4nesh/math-commons/blob/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/literature/README.md), SHA-256 `85d06dd8…acbf2b`
  - [`literature/sources.json`](https://github.com/g4nesh/math-commons/blob/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/literature/sources.json), SHA-256 `6113e91d…18b98b4`
  - [`literature/check_source_examples.py`](https://github.com/g4nesh/math-commons/blob/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/literature/check_source_examples.py), SHA-256 `1d757137…a7213`
  - [`literature/source-example-results.json`](https://github.com/g4nesh/math-commons/blob/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/literature/source-example-results.json), SHA-256 `7173175b…68ba2b`
- **Verdict:** **changes requested.** The three exact checks reproduce and every source statement I could reach is accurate. However, the lane's central status claim is out of date. Zheng's supplements S0, S1 and S2, which the README presents as "the concrete audit blocker", have been public in the author's GitHub repository since 2026-11-03. The lane also does not cite the established cases to primary sources. The corrections are small; see §2.
- **Evidence level of this review:** independent review. I checked the source claims against the primary sources I could retrieve, re-ran the checker, and recomputed every example with my own exact code. This is not an audit of Zheng's proof, and it does not establish any theorem.

## 1. What was checked, and how

| Claim in the artifact | How I checked it | Result |
| --- | --- | --- |
| Kovačec DMUC p2611 is v2 of 2026-11-04, 10 pages, 927388 bytes, with the stated SHA-256 | Downloaded the PDF from `mat.uc.pt` and hashed it. Read the `Last-Modified` header and the cover line. | **Holds**: same bytes and hash; `Last-Modified` 2026-11-04; the cover carries "DMUC Preprint 26-11, 07 Apr 2026" and "[v2] Wed, 04 Nov 2026". |
| Kovačec's unnumbered Theorem assumes simple connectedness and a Jordan polygon, and concludes an inner inclusion. The former Conjecture 7 has counterexamples. | Extracted the PDF text and read the Theorem and abstract. | **Holds**: the conclusion is cl$(J(A,B))\subseteq\Delta(A,B)$ provided $\Delta(A,B)$ is simply connected. The abstract says Conjecture 7 had counterexamples. |
| Corollary 3 prints the outside factor with $b_\nu$, where it should be $b_{\pi(\nu)}$ | Text extraction shows $b$ with a single subscript in the outside product, but $b_{\pi(i)}$ inside the bracket. Recomputed exactly, with $G=I$ and also with the rotation $c=3/5,\ s=4/5$. | **Holds**. For $G=I$ the actual value is 343, against 441 with the printed factor. For the rotation I derived the Givens bracket myself; the actual value is 8463/25, against 10881/25 with the unpermuted outside factor. |
| Goyal witness: $a,b,U$ as stated, $\det=2/3+2i/9$, inside the original hull | Data matched against Kovačec v2 §4. Recomputed with Gaussian elimination over $\mathbb Q(i)$. Found every $\sigma$-point triangle that contains the determinant, using exact barycentric coordinates. Then checked the author's weights. | **Holds**: $U$ is orthogonal, the determinant is as stated, and it lies in 4 triangles. The weights $(0,0,1/3,13/27,5/27,0)$ are exact. |
| $O_*$, $\psi$, $\psi(O_*)=-158/2401$, and the histogram $\{0{:}13,1{:}10,2{:}1\}$. The separating argument rules out common permutation weights, not MOC. | Matched $O_*$ and $\psi$ against Zheng's paper, Example 1.3. Recomputed both independently. Checked the argument by hand: $\psi$ is a signed sum of squared minors, so it is linear in $r(V)$ and nonnegative at every $P_\sigma$. | **Holds**. The scope statement is correct. Attribution: see §2.3. |
| Zheng preprint: DOI, 66 pages, Theorems 1.1 and 1.2, S0, S1, S2, Appendices F and G, 82, 6846, 214, software versions | DataCite record of the DOI. Text of `paper/moc4.pdf` from the author's repository: 66 page objects, SHA-256 `1f3dafad5021bbeff8b4a9cbc79de45846c47c7bfb1d542167e8a8c29aa1c138`. The ResearchGate PDF URL returned 403 for me as well. | **Holds** as a description. "Supplements not located" is **out of date** (§2.1). |
| Shitov 2025, Theorem 1 scope and date | OpenAlex abstract and Crossref (Izv. Math. 89(4), 862–869; 2025-08-25). mathnet.ru returned 403 and then 503. | Theorem 1 scope **holds**. Definition 5 and Theorem 19 were **not checked** (§3). |
| Drury 2006; Rodtes; the 2025 survey; Mulero-Martínez 2026 | ELA landing page and abstract, the arXiv abstract page, and Crossref for each DOI. | Scope descriptions **hold**. Rodtes has a published version that is not cited (§2.4). |
| Lane checker reproduces its recorded output | Ran it from a `git archive` export of `67528c1`, under `python3 -I` (3.11.7). | **Holds**: exit code 0. The output equals `source-example-results.json` exactly. |

## 2. Findings, most severe first

### 2.1 The n = 4 audit blocker is not current (blocks a pass)

The README's section "Four-dimensional proof: the concrete audit blocker" and handoff item 1 say S0, S1 and S2 were not located. `sources.json` says `"supplements_status": "... not located in bounded search ..."` and `"sha256": null`.

The author's public repository [`DongzheZheng/The-Marcus-de-Oliveira-conjecture-for-n-4`](https://github.com/DongzheZheng/The-Marcus-de-Oliveira-conjecture-for-n-4) was created on 2026-11-03, about three weeks before this checkpoint; its last push was also on 2026-11-03. Its README maps the material as follows:

- S0 is `docs/notation.md`.
- S1 is `cases/` (T01–T90 plus `index.tsv`).
- S2 is `docs/verification-scope.md`, with programs in `checks/` and data in `data/`.

At head commit `53ec16567d2307266ce52a8e9678410b2457cc27` the repository has 1472 files. These include the 66-page `paper/moc4.pdf` with the hash given above. `docs/provenance.md` lists that same hash for the paper version of 24 October 2026. The repository's own scope statement says the following:

- five portable exact checks, "not … a new independent verification of all ninety spectral inequalities";
- an incomplete Lean formalization;
- a note of 29 October 2026 by Ma and Zhang, *as reported by the author*. It checked the contact reduction and the finite classification, and left the ninety inequalities unchecked. I did not locate the note itself.

The README did label its search as bounded, so it is not dishonest. But the "blocker" framing and the priority handoff now send effort to a solved retrieval problem.

**Smallest correction:**

1. In `sources.json` `zheng-n4`, add the repository URL, commit `53ec165` and the PDF SHA-256. Change `supplements_status` to "public; not yet replayed or audited".
2. Rewrite the blocker section and handoff item 1. Point them at replaying `checks/` from a pinned export, and at auditing the ninety case proofs. The other three proposed audit tasks remain sound as written.

### 2.2 Established cases are not cited to primary sources

The task asks for "established cases". The lane cites three primary case results (Shitov, Drury 2006, Rodtes) and otherwise defers to a survey ("established special-case context"). It omits the basic cases. Zheng's §1 and §1.3 ("Related work") attribute the following. I confirmed each title and venue on Crossref; the case descriptions are Zheng's.

| Case | Primary source |
| --- | --- |
| n = 3 | Bebiano, Merikoski, da Providência, *LMA* 20(2) (1987) 167–170, doi:10.1080/03081088708817750 |
| Both matrices Hermitian | Fiedler, *PAMS* 30 (1971) 27–31, doi:10.1090/S0002-9939-1971-0286814-1 |
| One matrix with at least n − 2 equal eigenvalues | Kovačec, *LAA* 201 (1994) 91–97, doi:10.1016/0024-3795(94)90107-4 |

§1.3 also credits the following, which I did not open:

- de Oliveira 1982;
- Bebiano–Queiró 1985;
- da Providência–Bebiano 1986;
- Bebiano–da Providência 1988 (a pair of unitaries);
- Bebiano–Kovačec–da Providência 1994;
- Drury 2007, which Zheng says the concyclic-spectra extension of Drury 2006 also uses.

It credits Marcus 1973 for n = 2.

**Smallest correction:** add a short known-cases table with primary DOIs. Take each case's scope from the paper itself, not from a survey.

### 2.3 The common-weight obstruction lacks its primary attribution (minor)

The README calls the obstruction "known" but cites only Zheng's Example 1.3. Zheng's own text attributes the failure of common weights for n = 4 to Drury, *LAA* 168 (1992) 251–257, doi:10.1016/0024-3795(92)90296-M (title confirmed on Crossref). His footnote 2 attributes an equivalent form of $O_*$ to Holden's unpublished MI-05 report (13 September 2026).

**Smallest correction:** add these two attributions in the section heading's paragraph.

### 2.4 Smaller citation notes (minor, non-blocking)

- **Rodtes.** The work is published as *Oper. Matrices* 15(1) (2021) 127–130, doi:10.7153/oam-2021-15-09. Only arXiv v1 is cited.
- **Kovačec v2 on Zheng.** Kovačec does more than acknowledge the claim "tentatively":
  - the introduction says the paper "quite probably proves MOC for n = 4";
  - the introduction still states that MOC "is open even in the case n = 4";
  - §4 reports that Zhang and Ma are formally checking parts of it.

  The name appears as "Denzel Zhang" in Kovačec's introduction, as "Denzel Zheng" in the paper and repository, and as "Dongzhe Zheng" in DataCite. The manifest should record these variants.
- **Missing manifest entries.** The Goyal [G] and Zhang–Ma [ZM] preprints that Kovačec cites are not in the manifest, although the Goyal witness is used.
- **Mulero-Martínez entry.** It should carry the title, *LAA* 740 (2026) 19–38.

## 3. Reproduction

From the repository root at `67528c1` (with this review file added):

```sh
python3 -I research/marcus-de-oliveira/literature/check_source_examples.py   # lane checker; output equals source-example-results.json
python3 -I research/marcus-de-oliveira/reviews/check_literature_ishaan.py    # independent exact re-check; exit 0
curl -sL https://www.mat.uc.pt/preprints/ps/p2611.pdf | shasum -a 256         # 5221b16e…3fe54 on 2026-11-25
curl -sL https://raw.githubusercontent.com/DongzheZheng/The-Marcus-de-Oliveira-conjecture-for-n-4/53ec16567d2307266ce52a8e9678410b2457cc27/paper/moc4.pdf | shasum -a 256   # 1f3dafad…c138
```

`check_literature_ishaan.py` uses the standard library only and is deterministic. Its SHA-256 at the time of writing is `ddb9a39a…50fecbe`. If the recorded results file is present, the script also cross-checks it.

## 4. Limits of this review

- **Zheng's proof.** I did not audit it. I did not read the ninety case proofs, run any program from Zheng's repository, or verify Theorem 1.2. I did not compare the ResearchGate bytes with the GitHub PDF, because ResearchGate returned 403. Their identity rests on the author's provenance file and the matching page count and content.
- **Shitov.** I could not reach Definition 5 and Theorem 19 (mathnet.ru returned 403 and 503). The lane's conditional inference ($m=k=4$) is correct *if* the theorem is as the lane quotes it. I did not confirm that the hypotheses admit $m=k=4$.
- **Text extraction.** The PDF text extraction loses Greek glyphs. My reading of Corollary 3 rests on the subscript structure, $b$ with one index outside the bracket against $b$ with an argument $(i)$ inside it. I did not check the page numbers the lane cites.
- **Freshness of the search.** I made no systematic search for 2026 work beyond what the sources above cite. "Current" here means checked against those sources and the author's repository on 2026-11-25.
- **Scope of a pass.** A later pass on this lane would establish only that the inventory is accurate and current. It would say nothing about MOC or the n = 4 claim.
