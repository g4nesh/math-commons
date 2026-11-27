---
name: math-review
description: Independent adversarial review of a mathematical proof or certificate bound to an exact file hash - target match, claim-by-claim re-derivation, hypothesis checks for citations, re-running verifiers, numbered defects with severity, status boundary. Use when asked to referee, review, check or audit a proof, a PROOF.md, a certificate, or another agent's result.
---

# Independent mathematical review

Follow `harness/prompts/review.md`; the essentials:

- Record the sha256 of the file you review (`shasum -a 256 PROOF.md`). Your verdict is about that exact text.
- Do not read the author's attempts, notes or chat. Independence is the point.
- **Target match first**: compare the proved statement with STATEMENT.md#exact-statement word by word. A weaker statement cannot be a complete result.
- Re-derive each claim in dependency order. For citations: exact source statement, every hypothesis checked here. Re-run every verifier; write a small exact check of your own for any key computation.
- Hunt the usual failures: dropped cases, WLOG that loses generality, wrong inequality direction, unjustified limit/sum exchange, finite check presented as general, a cited theorem outside its hypotheses.
- Defects D1, D2, … with severity blocking/major/minor, location, description, repair (give the closing argument if you can).
- Verdict on the first line: pass only with every claim checked and nothing blocking or major. End with a status boundary: what this review did not establish.

In a campaign, record your review through the harness (`discover step SLUG --runtime manual --as <your-identity>` at the review stage, then `discover record`) so it is bound to the hash and counted.
