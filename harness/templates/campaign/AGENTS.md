# Campaign {{SLUG}}: rules for agents working in this directory

This directory is one discovery campaign run by `harness/bin/discover.mjs`. Start from `CAPSULE.md`; it is regenerated after every step and points to everything else by locator.

- Do not edit `campaign.json`, `CAPSULE.md`, `screen.json`, `directions.json`, `ledger.jsonl`, `runs.jsonl`, `receipts/`, `reviews/review-*`, or any `receipt.json`. The orchestrator owns them.
- `STATEMENT.md` is frozen after screening. Changing it sends the campaign back to screening.
- `PROOF.md` edits void every review bound to its previous hash.
- Lean: the statement (`formal/**/Statement.lean` and every definition it uses) is pinned by fingerprint after the statement audit. Prove exactly the pinned `Claim` (or `¬ Claim`).
- Run a step through the harness (`node harness/bin/discover.mjs step {{SLUG}}`), or, in an interactive session, take the prompt from `discover step {{SLUG}} --runtime manual` and finish with `discover record`.
- Evidence levels, scope discipline and file-writing rules: `harness/prompts/_preamble.md`.
