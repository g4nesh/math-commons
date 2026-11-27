# Context economy: fewer input tokens, more useful output

<!-- digest -->
Every model run is a fresh session that starts from about 2–25k tokens of deliberately chosen context instead of a growing transcript. Four mechanisms do the work:
1. A deterministic capsule (≤ `context.capsuleTokens`, default 1,500) that replaces history.
2. Role-specific context packs with token budgets, priorities and graceful degradation (full → producer digest → outline → locator).
3. Prompts ordered for prefix caching.
4. Outputs written to files with only locators returned.
Agents fetch anything else on demand with `discover show|outline|changed`.
<!-- /digest -->

## Why fresh sessions

Long research sessions get worse as they grow: earlier dead ends stay in context, the window fills with tool output, and resumption after a crash or limit hit loses everything. solved-conjectures handled this with an hourly continuation prompt pointing at a NEXT.md file. This harness makes that the only mode:
- every role run uses `--no-session-persistence` (Claude) or `--ephemeral` (Codex);
- continuity lives in files the orchestrator maintains.

A run that dies costs only that run. `discover step` resumes from the record.

## 1. The capsule (`CAPSULE.md`)

`lib/capsule.mjs` builds it after every step **from files, without a model**, so it is free and cannot drift from the record. It holds:
- stage, status, polarity and the computed label;
- the statement (truncated, with its hash and frozen state);
- claims with status and dependencies;
- certificate and review status;
- each direction with its last attempt;
- **Do not retry**: the last eight failed or killed approaches, with their failure points;
- next action: the stage's standing instruction plus the first failing gate checks with their fixes;
- the last transition and the budget.

Over budget, it drops sections in a fixed order (Do not retry, Directions, Certificates, Reviews, Claims) and leaves a pointer. A human or an interactive agent resumes by reading this one file.

## 2. Context packs

Each prompt's frontmatter lists `inputs` with priorities and a token `budget`. `buildPack` in `lib/prompt.mjs` resolves them:

| Input | Resolves to |
|---|---|
| `file:PATH`, `section:PATH#anchor` | the file or one section |
| `digest:PATH`, `outline:PATH` | the producer's digest block, or headings with line ranges and token costs |
| `direction` | this worker's direction JSON only |
| `ledger:direction` | the last 10 attempts on this direction, one line each |
| `promising-attempts` | full text of up to 4 live promising attempts (shrinks to one-liners) |
| `certificates:index` | one line per certificate: receipt state, claim, not-claimed, verify command |
| `certificate` | one certificate's manifest and file list (for the verifier, never the generator's reasoning) |
| `reviews:defects` | open blocking/major defects with suggested repairs, then minors |
| `lean:statement`, `lean:receipt`, `lean:project` | the statement files with the pin and closure; the last gate failures and build tail; the project layout |

Items are added by priority. An item that does not fit is replaced by its producer digest, then its outline, then left out and listed under "Available on demand" with its token cost. `discover pack SLUG --role R` shows the exact manifest before you spend anything.

**Independence by construction.** The pack is also the isolation boundary:
- the reviewer's inputs exclude attempts, the ledger and the map;
- the verifier sees only the certificate and the statement.

So the cheaper prompt is also the more independent one.

## 3. Prompt layout for caching

`assemblePrompt` orders every prompt as follows:

1. The shared preamble (static).
2. The role prompt (static per role).
3. Campaign constants: the exact statement (changes only at screening).
4. The capsule (changes per step).
5. The context pack.
6. The step header (run ID, role, write scope, task).

Anthropic and OpenAI both cache by prefix. Repeated runs of a role on a campaign therefore share everything up to the capsule, and runs of the same role across campaigns share the first two blocks. A unit test asserts that two runs share the prefix byte for byte.

## 4. Output discipline

- Roles write long content to files in their scope and return a small JSON object (strict schema). Nothing large returns through the orchestrator, so the user's own session stays small too.
- Artifacts over a page start with a `<!-- digest -->` block (≤ 10 lines). The draft gate requires one in PROOF.md. Packs prefer digests when budgets are tight.
- Attempts and reviews start with a one-line outcome/verdict, so a reader can stop after one line.
- Files stay under `context.maxArtifactBytes` (200 KiB, the GitHub preview limit); `reconcile` flags violations. Split long case analyses into `proof/<topic>.md`.

## 5. On-demand retrieval

| Command | Use |
|---|---|
| `discover show campaigns/X/PROOF.md#c4` | one section by heading anchor |
| `discover show campaigns/X/attempts/D1/002.md:L40-L80` | a line range |
| `discover outline FILE` | headings with line ranges and token estimates, before reading anything long |
| `discover changed X --since 2026-11-27T10:00` or `--since <git-ref>` | what changed since you last looked |
| `discover capsule X` | rebuild and print the capsule |

## 6. Maximising useful output per token

- **Effort routing.** Configure `runtimes.*.effort` globally and `effort:` per role. Bookkeeping (capsule, gates, receipts, publication) costs no model tokens at all.
- **Small, decisive first tests.** Every direction carries a falsification test that runs first, so bad directions die cheaply. The "Do not retry" list stops re-exploration.
- **One complete pass beats several partial ones.** The repair prompt asks for every defect to be fixed in one edit, because each PROOF.md edit voids all reviews.
- **Budgets.** `budget.tokens` per campaign; `step` refuses model runs past it. `runs.jsonl` holds runtime-reported usage (input including cache, output, cached input) from the existing Math Commons collectors.
