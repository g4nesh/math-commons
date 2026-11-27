---
role: verify-certificate
schema: verifier
budget: 6000
effort: high
maxTurns: 60
writes: ["certificates/{{CERTIFICATE}}/verify.py", "certificates/{{CERTIFICATE}}/manifest.json", "certificates/{{CERTIFICATE}}/VERIFIER.md"]
inputs: [{"use": "file:STATEMENT.md", "priority": 1}, {"use": "certificate", "priority": 1}]
---
# Step: write an independent verifier for one certificate

Another run produced `certificates/{{CERTIFICATE}}/`. Write a verifier that would catch a mistake in that run. You have deliberately **not** been shown the generator's reasoning; do not open the generator script except to learn the file format if `manifest.json` does not describe it.

## Requirements

1. Derive what must be checked from the exact statement and the manifest's `claim`, not from the generator. Write it down first in `VERIFIER.md` as a checklist: each condition of the statement the certificate must satisfy.
2. `verify.py`: Python 3.9+ standard library only; integers and `fractions.Fraction`; no floats, no `math` module, no numpy. Parse `certificate.json`, check every condition, print one line per check with its result, print `VERIFIED` and exit 0 only if every check passes; otherwise print the failure and exit 1. Keep it short and readable: a referee will read every line.
3. Check the conventions (indexing, simple vs multi-graphs, closed vs open intervals) against STATEMENT.md, not against the generator.
4. Include at least one negative control: a deliberately corrupted copy of the data (in memory) that the verifier must reject, and assert that it does.
5. Set `verifier.path` to `verify.py` and `verifier.command` to `["python3", "-I", "-B", "verify.py"]` in `manifest.json`; set `verifier.runId` to the run ID of this step. Do not modify the certificate data.

Then run it once yourself and report the outcome. The orchestrator re-runs it and writes the receipt; a verifier that passes only on your machine fails the gate.
