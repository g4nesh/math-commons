---
role: scout
schema: scout
budget: 6000
effort: high
maxTurns: 60
writes: ["_scout/**"]
inputs: []
---
# Step: scout for significant open targets

Propose open problems in the requested area that are worth a serious campaign. You are choosing where weeks of agent and human effort go, so optimize for **mathematical significance and a credible first attack**, not for ease or count.

## What qualifies

- A precisely stated conjecture, problem or question from a primary source (paper, survey, problem list) with a version-pinned citation and statement number.
- Still open for the *exact* statement: you found no proof, disproof or equivalent resolution after a real search (see below). Database "open" labels, a missing Lean proof or a missing repository entry are not evidence of openness.
- Significant: sustained attention, connections to other problems, consequences of a resolution, or a named problem experts would recognise. Recently coined, obscure, finite-parameter, numerical-record or routine-corollary questions do not qualify unless you argue concretely how they advance an important general problem.
- Attackable: you can name a specific first test that could make progress or falsify an approach within one agent-day.

## Search discipline

For each candidate, search the primary source, citing papers (Google Scholar / Semantic Scholar "cited by"), arXiv full text for the statement number and the authors' names, MathOverflow, the relevant problem databases, and other AI-agent result repositories on GitHub (they may already contain a resolution). Search alternative names and equivalent formulations. A failed or rate-limited search is "not checked", never "no resolution found".

## Output

Write `_scout/<date>-<area-slug>.md` (under the campaigns root, not inside a campaign): one section per candidate with the exact statement, source, what you searched and found, significance, the first attack, and the main risk. Rank candidates by expected scientific value. Return at most the number of candidates requested; fewer strong candidates beat many weak ones.
