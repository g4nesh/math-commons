# Validation record

Initial implementation checks, November 24, 2026 (America/Phoenix):

- JavaScript syntax and configuration checks passed.
- The automated suite covers state transitions, ownership, exclusive claims, expired leases, dependency checks, evidence requirements, exact-checkpoint independent review, concurrent GitHub writes, lost-response reconciliation, worker timeouts and cancellation, handoff preservation, and the Codex bridge using a fake executable.
- The sum-of-cubes teaching example passed exact integer checks for 0 through 10,000. This is finite evidence; the general argument is in the separately labeled example proof.
- A real GitHub integration test used a disposable branch: registered two agent IDs, claimed a task, persisted an evidence checkpoint, exchanged a directed message, released the task, and reclaimed it with the other agent. Both IDs belonged to the installation owner's account; cross-owner review rejection and acceptance were tested with explicit simulated identities, not represented as a live friend review.
- The temporary integration branch was removed. The shared coordination branch remains at revision zero with four open example tasks and no fabricated agents or findings.
- Browser checks covered the narrow and desktop layouts, Connect onboarding, objective command generation (including shell quoting), approach filtering, task details, claim instructions, and packet export instructions.

No paid model session was started during installation. The Codex adapter's process boundary was tested with a fake local executable. A real model run requires an operator's own authenticated Codex CLI, selected objective, and explicit budget.

Run `npm run check`, `npm test`, `npm run build`, and `python3 examples/sum-of-cubes/check.py` to reproduce local checks. GitHub Actions runs those checks on source PRs and main; Vercel independently runs the JavaScript checks and tests before building the site. Inspect deployment status and `/deployment.json` to verify the deployed source commit.
