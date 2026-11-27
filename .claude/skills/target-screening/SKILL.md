---
name: target-screening
description: Screen a proposed open problem before investing in it - pin the exact statement and source, run a real prior-resolution search (including other AI-math repositories), argue significance, separate full from partial resolution, name attacks. Use when choosing a research target, scouting conjectures, or checking whether a "new" result is already known.
---

# Target screening

Follow `harness/prompts/screen.md`. The decision is one of `eligible`, `prior-resolution`, `low-impact`, `needs-clarification`.

- **Exact statement**: verbatim from a version-pinned primary source with statement number and page; all conventions.
- **Prior-resolution search**: statement number and author names on arXiv full text and Scholar (citing papers), equivalent formulations and consequences of known theorems, MathOverflow, problem databases, Mathlib, and other AI-generated mathematics repositories (README and issues). Record every query; a failed search is "not checked", not "nothing found".
- **Significance**: what a resolution changes, who worked on it, what it connects to. Reject record-chasing, isolated finite cases, routine corollaries and artificial variants.
- **Full vs partial**: what must be covered for a complete resolution.
- **Attacks**: 1-3, each with a cheap first test.

Repeat the duplicate check before publication: a concurrent agent may have resolved it.
