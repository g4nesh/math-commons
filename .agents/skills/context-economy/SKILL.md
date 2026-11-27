---
name: context-economy
description: Keep long mathematical research sessions cheap and resumable - read the capsule not the history, expand locators on demand, write digests, return locators instead of content. Use when starting or resuming campaign work, when context is getting long, or before handing off.
---

# Context economy

Tokens spent re-reading are tokens not spent on mathematics. The harness is built so that every step can start fresh.

## Reading

- Start from `campaigns/SLUG/CAPSULE.md` (≤ ~1.5k tokens, regenerated deterministically after every step). It lists claims, certificates, reviews, directions, killed approaches and the next action, each with a locator.
- Expand one thing at a time: `node harness/bin/discover.mjs show campaigns/SLUG/PROOF.md#c4` or `…:L120-L160`. List a file's sections with their token cost first: `discover outline FILE`.
- See what changed since you last looked: `discover changed SLUG --since 2026-11-27T12:00` (or a git ref).
- Never cat a large JSON or log; select the field (`jq '.failures'`) or the tail.

## Writing

- Put a `<!-- digest -->…<!-- /digest -->` block (≤ 10 lines: result, evidence level, what is not claimed) at the top of any artifact over a page. Context packs send the digest instead of the file when budgets are tight.
- First line of every attempt/review: outcome and status, so a reader can stop after one line.
- Long proofs: PROOF.md under 200 KiB; case analyses in `proof/<topic>.md` linked from the claim.
- Final answers return locators and a ≤ 3-sentence summary, never file contents.

## Handoff

Before stopping, make sure the capsule's "Next action" is true (run `discover capsule SLUG`), commit artifacts, and record failed approaches in the ledger (they become the capsule's "Do not retry" list).

## Check a prompt's cost

`discover pack SLUG --role review` shows exactly what a role would receive, item by item with token estimates, and what was omitted.
