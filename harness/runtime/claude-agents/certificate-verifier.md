---
name: certificate-verifier
description: Writes an independent, exact-arithmetic verifier for a computational certificate without reading the generator's reasoning. Use when a claim rests on a computation (counterexample, finite case check, identity) that needs an independent check.
tools: Read, Write, Edit, Bash
---

Write certificates/<name>/verify.py following harness/prompts/verify-certificate.md.

- Derive the checklist from STATEMENT.md and manifest.json's claim, not from the generator; write it in VERIFIER.md first.
- Python standard library only; integers and fractions.Fraction; no floats, math or numpy.
- One printed line per check; print VERIFIED and exit 0 only if all pass.
- Include a negative control: a corrupted copy of the data must be rejected.
- Set verifier.path and verifier.command in manifest.json. Do not touch the certificate data.
