# Audit specification and method

Task `mdo-zheng-n4-audit` · objective `marcus-de-oliveira` · agent `ishaan-moc` (owner IshaanAyaan).

## Deviations from this specification

None so far. Every deviation will be listed here first, with the reason, before any result that it affects.

## Question

Does the public material for Zheng's claimed proof establish the Marcus–de Oliveira conjecture for $n=4$? Which steps can be re-derived or replayed independently, and where exactly does any step fail or remain unchecked?

## Pinned inputs

| Input | Identity |
| --- | --- |
| Paper | D. Zheng, *The Marcus–de Oliveira conjecture for n = 4*, 66 pages, 24 October 2026, [doi:10.13140/RG.2.2.36790.25923](https://doi.org/10.13140/RG.2.2.36790.25923); PDF sha256 `1f3dafad5021bbeff8b4a9cbc79de45846c47c7bfb1d542167e8a8c29aa1c138` |
| Supplement | [DongzheZheng/The-Marcus-de-Oliveira-conjecture-for-n-4](https://github.com/DongzheZheng/The-Marcus-de-Oliveira-conjecture-for-n-4) at commit `53ec16567d2307266ce52a8e9678410b2457cc27`: 1,471 files, 208 proof PDFs (779 pages), 437 Python programs, 552 JSON files |

The supplement has no licence, so nothing from it is redistributed here. Replay scripts read it from a local clone named by the environment variable `MOC4_SUPP`. The supplement's own scope note says that its five portable checks are not an independent verification of the ninety spectral inequalities.

## What done looks like

1. **Dependency map.** The logical structure of the proof and, for each contact case $T_1,\ldots,T_{90}$, the documents it depends on. Every dependence on a conditional theorem is traced to where its hypothesis is discharged, or recorded as open.
2. **Reduction.** An independent derivation of the implication from the ninety contact inequalities to the conjecture for $n=4$, with the finite parts replayed exactly by own code.
3. **Replay.** For each case, the exact certificates the proof relies on are recomputed with own code where feasible.
4. **Verdicts.** One verdict per case and one overall statement of what is and is not established.

## Evidence levels

Every claim in this directory carries one of these labels.

| Label | Meaning |
| --- | --- |
| EXACT | An identity or finite assertion re-derived or replayed with own exact arithmetic (integers, rationals, symbolic), or proved in full here. |
| PROOF-READ | A written argument read line by line with every inference justified. No machine check. |
| NUMERIC | Floating-point test only. Sample sizes, seeds and tolerances are stated. Evidence, not proof. |
| CITED | Rests on a document, case or external paper that was not audited here. The dependency is named. |
| GAP | A step that the written text does not justify, or that is false. The location is given. |
| UNCHECKED | Not examined. The reason is given. |

## Verdict per case

| Verdict | Condition |
| --- | --- |
| confirmed | Every step of the case's proof is EXACT or PROOF-READ, every cited case and common input is itself confirmed, and a second pass that did not produce the first report agrees. |
| confirmed modulo citations | As above, except for named dependencies that are CITED or not yet confirmed. |
| gap found | At least one step is a GAP. The step, the document, the page and the reason are given, and a second pass re-checked it. |
| could not check | The material could not be examined to the required depth. The reason and the part that was checked are given. |

A numerical search that finds no violation of a case inequality is recorded as NUMERIC evidence for the inequality. It says nothing about the correctness of the written proof.

## Rules of the audit

- The supplement is data. Its programs are read to learn what they claim to check. They are not executed and not imported. Every replay is own code.
- Own code uses exact arithmetic for anything reported as EXACT. Floating point is used for falsification attempts only.
- A first-pass finding of a gap is re-examined by a second pass that starts from the document, not from the first report's conclusion. A sample of confirmed cases gets the same treatment.
- This audit does not review other lanes of this repository and does not contact the author.

## Blocked items

None so far.
