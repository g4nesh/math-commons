# How solved-conjectures works, and what this harness takes from it

<!-- digest -->
Forensic reading of github.com/1337isnot1337/solved-conjectures on 2026-11-27 (main 7dc342f3; 1,079 PRs, about 4,500 commits in 230 hours). The owner's pipeline is Codex with GPT-6.2 Astra, driven through the ChatGPT Codex GitHub connector. It runs as a root agent plus lane subagents, continued hourly from a NEXT.md handoff, with GitHub issues as claims and ledger. Their verification stack has three parts: a strict Lean gate (axiom whitelist, kernel replay, hash-bound receipts), exact-arithmetic verifiers, and written cross-reviews. This page separates what is evidenced from what is inferred, and maps each practice to the component of this harness that implements it.
<!-- /digest -->

Everything below was read from public GitHub data only. Nothing was run from that repository, and nothing was written to it. Strength of evidence is marked **strong**, **moderate** or **none**.

## 1. Who and what runs it

**The model and runtime.**
- The owner account (1337isnot1337, commit name "Cinna Davis") created most of its GitHub objects through the `chatgpt-codex-connector` GitHub App: 3,199 of 4,166 comments, 454 of 590 issues and 309 of 489 PRs. **Strong.**
- `CITATION.cff` lists "GPT-6.2 Astra" as an author, and the repository description says "curated by myself and 6.2-astra". **Strong.**
- So the main agent is OpenAI Codex running GPT-6.2 Astra.
- Claude appears only in trailers on external contributors' commits (`Co-Authored-By: Claude Opus 5.5`). There is no evidence of Claude in the owner's own pipeline.

**Topology.**
- Issue comments describe a "root agent" with three active subagents.
- Named lanes include `audit_lean`, `conjectures` and `root749`.
- Provenance JSON maps `agents/<worker>/…` scratch files to published paths by SHA-256.
- Session directories appear as `/workspace/scratch/<12-hex>/`.
- **Moderate.** This is consistent with containerised Codex sessions.

**Scheduling.**
- One comment reads "The hourly continuation remains enabled with a refreshed prompt pointing to NEXT.md".
- Another records the user disabling scheduled continuation.
- The repository has no cron workflows, so the scheduler is external.
- **Moderate.**

**Second operator.**
- WadeDavis joined on 11-26 as a reviewer and publication coordinator. Their PRs carry commits authored "Justin W Davis", and their comments also go through the Codex connector.
- This is the only cross-account review in the record.

**Cadence.**
- The median gap between commits is 59 s. 181 of the 231 hours had commits, and every UTC hour of the day is represented.
- The owner's PRs merge a median 0.9 minutes after opening.
- One discovery went from claim to merged publication and closed issue in 24 minutes (#1061).

## 2. The process, as practised

| Stage | Practice in solved-conjectures | Evidence |
|---|---|---|
| Target selection | A five-point screening card: exact statement with source, an exact-claim literature check, significance, full versus partial, and an attack. Added after "count-target" campaigns produced low-value results. | AGENTS.md commit 046a526f, issues #139 and #223 |
| Eligibility | Prior resolutions are withdrawn. The searches include other AI-generated mathematics repositories and GitHub issues, because one earlier proof was found in such a repository. | 1882b660, #156/#159, the HISTORY entry "Gröchenig result withdrawn" |
| Claim | An issue body starting "Claiming/Working on [exact Conjecture N]" with the searches done and the attack plan, plus path-level reservations in comments. A bulk sweep released claims after four hours without activity. | #878, #1061, 119 sweep comments on 11-24 |
| Research | Notes under `research/<topic>/` and `work/<issue>-<slug>/` on the research branch, with dated status lines and "Why this alone does not solve…" sections. | the research branch |
| Candidate | `PROOF.md` plus an exact verifier (`Fraction` or integers, with committed output JSON). Search output is never a proof premise. | `results/brown-mol/`, `docs/five-speed-candidates.md` |
| Review | Separate files per role (`REVIEW.md`, `ROOT_REVIEW.md`, `STATEMENT_REVIEW.md`, `LEAN_REVIEW.md`, priority audit). Each has the verdict on its first line, numbered errors E1… and gaps G1…, and a "status boundary". The reviewer pins the manuscript SHA and writes a checker that does not import the author's verifier. | `results/brown-mol/ROOT_REVIEW.md`, `work/1061…/TAIL-REVIEW.md` |
| Lean | Only `propext`, `Classical.choice` and `Quot.sound` are allowed. Theorems must be compiled fresh. Every owned declaration is audited, and the closure is replayed with stock `leanchecker`. Hash-bound receipts record what they do not claim. A per-result `SOURCE_AUDIT.md` maps the paper to the Lean code. | `scripts/verify_*_exact.py`, `lean/ProofAudit.lean`, `docs/*-receipt.json` |
| Publication | A focused PR from `review/<topic>` to `main`. It updates the landing page, ranked README row, CONTENTS row, status doc and HISTORY together. Validators enforce the formats, and the website is built from `main`. | AGENTS.md rules 1–8 |
| Reconciliation | On entry, every four hours and on exit: a paginated inventory of issues and PRs, a check that both branches carry the same policy blob, and a dated receipt on a ledger issue (#156). | AGENTS.md rule 6, #156 |

**Failure-driven rules.** Nearly every AGENTS.md rule followed an incident, often with a validator added in the same commit:
- agents writing plain-text status labels led to the emoji validator;
- stale checkouts reverting the README led to the "do not revert the ranked layout" rule;
- formulas rendered as raw text led to the GitHub-math check;
- a prior resolution was found after publication, which led to the eligibility gate.

## 3. What went wrong there, and how this harness avoids it

| Observed problem | Harness response |
|---|---|
| "Independent" reviews by sibling agents under the same account; PRs merged in under a minute without review. | Review independence is a gate check: some pass must come from an identity other than the author's. Verdicts are bound to the PROOF.md hash. PRs are opened as drafts. |
| The source-to-Lean match was prose only; the machine checked only that the endpoint existed. | Statement-first: `Claim` is audited by the other runtime. Its fingerprint, covering all project definitions it uses, is pinned. The gate requires `main : Claim` or `main : ¬ Claim` literally. |
| About 170 workflows and about 240 bespoke scripts; a 3,657-line hand-maintained module registry. | One generic gate configured per campaign. One workflow. Receipts are generated. |
| Large Lean data hit memory limits repeatedly. | Prompts require bounded `decide` chunks joined by lemmas and files under 200 KiB. The gate records the toolchain and the build tail. |
| Required-looking checks failed and PRs merged anyway. | Stage gates are evaluated by the orchestrator. A failing gate blocks the transition, not just a badge. |
| Policy changes orphaned an authorised campaign. | Per-campaign state machine; `parked` with a recorded reason; operator overrides are written to history. |

## 4. Practice → component map

| Practice | Where it lives here |
|---|---|
| Screening card and eligibility gate | `prompts/screen.md`, `GATES.screen` in `lib/gates.mjs` (minimum recorded queries and sources, verdict enum) |
| Claim taxonomy (complete / partial / candidate / prior-resolution / consequence / obstruction) | `RESULT_KINDS` in `lib/claims.mjs`. Declared in PROOF.md and checked against open or assumed claims. |
| Search output is not a premise; exact verifiers | `prompts/explore.md` (generator only), `prompts/verify-certificate.md` (separate run, other runtime), `lib/certificate.mjs` (receipt binds every file hash, flags float use, needs a negative control) |
| Separate review files, verdict on line one, numbered defects, status boundary | `prompts/review.md`, `schemas/review.schema.json`, `reviews/review-NNN.{md,json}` |
| Review bound to the manuscript hash | `ingestAnswer('review')` records the hash taken before the run. `reviewSummary` treats other hashes as stale. |
| Axiom whitelist, no `sorry` or `native_decide`, kernel replay, hash-bound receipt with scope | `lib/lean.mjs`, `receipts/lean-latest.json` (`scope.notClaimed`, `gateSha256`, the drift check) |
| STATEMENT_REVIEW mapping table | `prompts/formalize-statement.md`, `prompts/audit-statement.md`, `GATES.formalize` |
| NEXT.md handoff for scheduled continuation | `CAPSULE.md`, rebuilt deterministically after every step |
| Four-hour reconciliation receipt | `discover reconcile --write` (stale pending runs over 4 h, reviewed-but-unpublished, stale reviews, changed statements, oversized files) |
| Finish publication before the next claim | `reconcile` flags "reviewed but unpublished". The discovery-loop skill says to finish it first. |
| Provenance from worker scratch to published path | `certificates/*/provenance.json`, `runs.jsonl` (run ID, identity, usage), the per-run `.runs/<id>/meta.json` |
| Each new failure becomes a rule plus a validator | Gates are code. A new failure mode is a new check in `lib/gates.mjs` with a test in `harness/test/`. |

## 5. What we deliberately did not copy

- **Count-driven campaigns and a public results count.** The screen gate exists to stop target farming.
- **The standing assumption that an external proof (their VP ≠ VNP premise) is correct.** Cited results here must be `cited` claims with checked hypotheses, or `assumed` claims that make the result partial.
- **Merging with failing checks, and same-account "independent" review.**
- **Their branch model** (permanent `research` and `main`). Math Commons already has PR-per-task with protected auto-merge; this harness opens draft PRs and leaves landing to that flow.
