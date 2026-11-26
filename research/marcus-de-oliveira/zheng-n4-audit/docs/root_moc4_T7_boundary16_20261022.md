# root_moc4_T7_boundary16_20261022 (seven-contact face T_8 ∪ {16})

Evidence level: EXACT for every displayed identity, including the two interval derivations (10 checks, 0 failures). PROOF-READ for the sign conclusion. CITED for the SOC atom and the collision/circle layers.

PDF: proofs/manin_transfer_20261021/root_moc4_T7_boundary16_20261022.pdf (1 page). Replay: `replay/root_moc4_T7_boundary16_20261022__check.py` and its output.

## Deviations and limits

- The replay rebuilds the T_8 source curve, using the document's sign convention for the bracket with ∞. It then specializes to s = 1 exactly.
- The document does not display the derivations of LV and LG. I derived them myself and verified exactly: g3 >= d is equivalent to 2(1+X)Re w <= LG, and g10 >= 0 is equivalent to 2(1+X)Re w >= LV, using the cross-multiplied (6), which also holds at RX = 1.

## A. Statements

- **Claim (p. 1).** Let x be in X with g >= 0 and contacts ⊇ {0,2,5,6,7,8,16}. Generators: {0,2,3,5,6,39,45,46}. Then Phi >= 0 for all O. In the nonminor case d > 0 the canonical coefficients satisfy ab >= d^2, so the SOC atom on (f2, f39; f46) applies. Collision and circle strata are cited as separately safe.

## B. Dependencies

- T_8 chart identities (root_moc4_T7_uniform): audited in this batch, EXACT.
- SOC atom (F03): CITED. Its orbit relation was checked in the T_8 replay.
- Collision and circle strata (paper B.1, B.2; three-A/three-C circle theorem): CITED.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | g16 = (1-s)g12 and g12 > 0 (from d = s g12 - g18 > 0) give s = 1 | EXACT + PROOF-READ |
| 2 | the s = 1 forms of L, J, Z; d = tR^2/|w|^2; g3 = (1+X)LR^2/|(R+1)w + R^2(1+X)|^2 | EXACT |
| 3 | LV <= 2(1+X)Re w <= LG | EXACT (auditor's derivation) |
| 4 | LG - LV = (X-R)(1+d)L/((R+1)d), so R <= X | EXACT + PROOF-READ |
| 5 | g23 = R - d >= 0, so d <= R <= X, and a = Xd >= d^2 = d^2 b (b = 1) | EXACT + PROOF-READ |

## D. Falsification

All claims are identities and were replayed exactly. Nothing was left that numerical testing could falsify.

## E. Verdict

**confirmed-modulo-citations.** The citations are the SOC atom and the separately cited collision/circle layers.
