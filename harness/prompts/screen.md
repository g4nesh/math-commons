---
role: screen
schema: screen
budget: 8000
effort: high
maxTurns: 80
writes: ["SCREEN.md", "STATEMENT.md", "sources/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}]
---
# Step: screen the target before any research

Decide whether this campaign may proceed. This gate exists because the most expensive failure of an autonomous research program is spending days "discovering" something already in the literature, or polishing a trivial target. Be skeptical; a rejection is a successful screen.

## 1. Pin the exact statement

Open the primary source. Copy the statement verbatim into STATEMENT.md#exact-statement (you may fix only formatting) with: source version (e.g. arXiv v2), statement number, printed page or section, and every definition and convention it depends on (normalisations, whether graphs are simple, whether 0 ∈ ℕ, …). Under `## Scope` say what a complete resolution must cover and which weaker results would be partial. If the source is ambiguous, list the readings and choose none: return verdict `needs-clarification`.

## 2. Prior-resolution search (bounded but real)

Search, at minimum:
- the statement number and the authors' names on arXiv full text and Google Scholar / Semantic Scholar, including papers citing the source;
- alternative names and equivalent formulations of the statement, and its direct consequences of established theorems;
- MathOverflow / Math StackExchange, the relevant problem databases and surveys;
- formal libraries (Mathlib) and other AI-generated mathematics repositories on GitHub (search their issues and READMEs), which have already resolved many "open" problems;
- this repository: existing campaigns and `results/`.

Record every query and where you ran it, every source you actually read with its URL and what it says about the exact claim, and every search that failed (timeouts and rate limits are failures, not negative evidence).

## 3. Decide

- `prior-resolution`: the exact statement (or something that implies it) is settled. Cite it precisely.
- `low-impact`: open but not worth a campaign under the significance policy (isolated finite case, record chasing, routine corollary, artificially narrow variant). Say why.
- `needs-clarification`: the statement cannot be pinned.
- `eligible`: open as far as a careful search shows, significant, and you can name an attack.

For `eligible`, argue significance concretely (what a full proof or disproof would change, who has worked on it, what it connects to), separate full from partial resolution, and give 1–3 attacks each with a cheap first test.

## Files

- `SCREEN.md`: first line `**Decision:** <verdict> (<date>)`, then sections `## Exact target`, `## Prior-resolution search` (table: query | where | result), `## Sources read`, `## Search failures`, `## Significance`, `## Full vs partial`, `## Attacks`.
- `STATEMENT.md`: filled exact statement, source and scope.
- `sources/`: optional short excerpts (fair-use length) with citations; never whole papers.
