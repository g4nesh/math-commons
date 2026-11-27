---
role: repair
schema: repair
budget: 16000
effort: high
maxTurns: 120
writes: ["PROOF.md", "proof/**", "certificates/**"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "reviews:defects", "priority": 1}, {"use": "file:PROOF.md", "priority": 1}, {"use": "certificates:index", "priority": 3}]
---
# Step: repair the defects reviewers found

Reviewers found the defects listed in the context pack against the current PROOF.md. Fix them in PROOF.md (and `proof/`, certificates if needed).

- Address every blocking and major defect. For each, either repair the argument completely, weaken the claim (and the result kind) honestly, or explain with evidence why the reviewer is mistaken. Do not paper over a gap with "clearly" or "it is easy to see".
- If a reviewer supplied a closing argument, check it before adopting it, and credit it in the proof.
- Keep claim IDs stable. Add new claims for new lemmas; do not renumber.
- Any edit produces a new hash, so every earlier review becomes stale and the campaign returns to review. Make the edit count: fix everything you can in this pass.
- If the defect shows the approach cannot work, set the result kind accordingly (`obstruction`, `partial` or `candidate`) and say so in Scope and limitations.

Add a short `## Revision notes` section at the end of PROOF.md: date, defects addressed (by review file and ID), and what changed.
