# Audit: root_inventory50_bounded_native_conditional_20261023

Evidence level: EXACT for every algebraic identity and combinatorial claim (`replay/root_inventory50_bounded_native_conditional_20261023__check.py`, 31 checks, all pass); PROOF-READ for the analytic curve argument and the wall restoration (puncture set and Rouche quantifiers made explicit by the auditor); CITED for the boundary cases, the seven-element hypothesis D, the bounded-negative-part lemma and Appendices B and C; NUMERIC for the falsification run.

PDF: `proofs/manin_transfer_20261021/root_inventory50_bounded_native_conditional_20261023.pdf`, 3 pages. Role: conditional reduction for paper case T_51 (also listed by the catalogue for T_87).

## Deviations and limits

- The supplement's checker ("148 exact groups") was not run or read; all identities re-derived from the PDF formulas and the definitions.
- The Laurent analysis at W = 0 was done by exact power series of the branch to first order and of numerator and denominator of H14 to order W^1. "H4 + H14 holomorphic" is deduced from: both poles simple, W a local parameter on both branches (simple roots), residues r and -r. The direct expansion of the sum was dropped because it exceeded the 30 s light-job limit while the heavy pool had no free slot for 15 minutes.
- Numeric falsification is small: 4 stratum points (seeds 0-3), 12 BFGS restarts each.
- The document does not name the six-element pullbacks for ties 0, 7, 16; resolved by the auditor.

## A. Statements

**Theorem 1 (p. 1).** Assume (H_D): (*) holds for every x in X with P c(x) >= 0 whose contacts contain some image of D = {0,1,2,7,9,16,23} (Table 5 T_12^(7)). Then (*) holds for T = {1,3,5,10,11,12} (= paper T_51 literally): for every x in X with P c(x) >= 0 and contacts containing T (non-strict), Phi_{c(x)}(O) >= 0 for all O in SO(4). Further inputs: "the 60 preceding" six-element cases, B = T_8^(7), C = T_11^(7), Appendices B, C ("colliding spectra satisfy the support inequality by the global theorems"). With those all spectra are covered.

Sub-claims: (1) unique nonnegative expansion in f0,f1,f2,f8,f10,f12,f30,f59 with d = g4 the only non-minor coefficient, L_O >= d f30(O); ties; inverse (2); H levels; (3) H4, H0; the 105-term curve num(H12 - h) of degree (2 in Z, 5 in W); W = 0 local models (finite branch Z(0) = -1/(r+s) with Z'(0) as displayed; second branch V = WZ -> (r+s)/(r+s-u)); H14 has a simple pole of residue -r on both branches, so H4 + H14 is holomorphic; (4) 0 <= d <= d + g14 = -eta Im(H4 + H14); Z^2-coefficient of the curve; Rouche restoration.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| D = {0,1,2,7,9,16,23} | seven-contact (hypothesis) | T_12^(7); tie 20 (a group element maps D onto T+{20}) | EXACT action |
| B = {0,1,2,7,9,10,16} | seven-contact | T_8^(7); ties 6, 13 | EXACT action, CITED |
| C = {0,1,2,7,9,16,17} | seven-contact | T_11^(7); tie 14 | EXACT action, CITED |
| ties 0, 16 | paper cases (unnamed: "strictly prior closed pullbacks") | T_9, T_11, T_17, T_38, T_65, T_76 map into T+{j} (stages 3-6 < 7) | EXACT existence, CITED |
| tie 7 | paper cases (unnamed) | T_42, T_65, T_71 | EXACT existence, CITED |
| bounded-negative-part lemma | common lemma | final_referee_bounded_bad_harmonic_lemma (its hypotheses 1-3 match the use, read) | CITED |
| collision theorems | common lemmas | paper App. B, App. C | CITED |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = paper T_51 | identity | EXACT | |
| 2 | contact rank 6; kernel basis f0,f1,f2,f8,f10,f12,f30,f59; rows 20,16,13,14,7,6,4,0 form I8 | identity | EXACT | equals lib generators(T); so each coefficient is a gap and is >= 0 |
| 3 | only f30 non-minor, coefficient d = g4; L_O >= d f30(O) | identity + inequality | EXACT + PROOF-READ | 30 exact random instances |
| 4 | ties outside {0,6,7,13,14,16,20} force d = 0; 6,13 -> B, 14 -> C, 20 -> D, 0,7,16 -> earlier cases | case split | EXACT | |
| 5 | raw ratios, inverse (2), denominators sum to 1 + W | identity | EXACT / PROOF-READ | |
| 6 | contact levels; H4 = r(W+1)/W, H0 = (W-r)(Z-1) | identity | EXACT | from raw sigma-point products |
| 7 | curve numerator: 105 terms, degree 5 in W and 2 in Z; last-coordinate denominator of degree 1 in Z | identity | EXACT | |
| 8 | punctures with W != 0 and W = infinity have bounded d | analytic | PROOF-READ | H4 |
| 9 | F(0,Z) = r^3((r+s)Z+1); Z(0), Z'(0) | identity | EXACT | |
| 10 | finite branch: H14 simple pole, residue -r | identity | EXACT | power series |
| 11 | every monomial has deg_W - deg_Z >= -1; lowest weighted polynomial -r^3 V((r+s-u)V - r - s); V(0) simple | identity | EXACT | |
| 12 | V-branch: H14 simple pole, residue -r | identity | EXACT | power series in W with V(W) |
| 13 | H4 + H14 holomorphic on both branches; these are the only branches over W = 0 (degree-2 Z projection, one finite root, leading coefficient vanishing simply at W = 0) | analytic | PROOF-READ from 9-12 | |
| 14 | (4): 0 <= d <= g4 + g14 = -eta Im(H4 + H14) bounded at W = 0 | inequality | PROOF-READ | g14 >= 0 on the support region |
| 15 | generic conditions proper together, witness (2,3,4,5) | identity | EXACT | |
| 16 | nonconstant coordinate (H4 or H0) | analytic | PROOF-READ | |
| 17 | harmonic exclusion | analytic / citation | PROOF-READ + CITED | puncture set as in the T_50 conditional report, note (a): all non-regular points go into P; d is bounded near each |
| 18 | Z^2 coefficient -W(r-W)^2((u-1)W+r)(r+s-u-W) nonzero for r != 0, so the last germ is nonconstant; Rouche restores generic parameters | identity + analytic | EXACT + PROOF-READ (argument as in the T_50 report, note (b)) | |

## D. Falsification

- EXACT instances of L_O >= d f30(O), d = g4 (30 random rational points of N_T times rational SO(4)).
- NUMERIC: 4 stratum points (seeds 0-3, pi* = pi_1), 12-restart BFGS over SO(4): worst min Phi = +1.3e-15. No violation.

## E. Verdict

**confirmed-modulo-citations.** Citations: hypothesis D; B, C; the six-element cases for ties 0, 7, 16 (earlier stages, resolved by the auditor); the bounded-negative-part lemma; Appendices B, C. Presentation remarks: pullbacks for ties 0, 7, 16 unnamed; puncture set described loosely.
