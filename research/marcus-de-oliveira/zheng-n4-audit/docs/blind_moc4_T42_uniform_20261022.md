# Audit of blind_moc4_T42_uniform_20261022

Evidence level: EXACT for the reduction, cubic restriction and Q=0 guard; NUMERIC (exact-at-points) for the radial factor and its lower bound; CITED for the capacities and the 60-table certificate.

Auditor batch: C18-T43-T45. Source: `proofs/manin_transfer_20261021/blind_moc4_T42_uniform_20261022.pdf` (3 pages). Replay: `replay/blind_moc4_T42_uniform_20261022__check.py` (output in `replay/blind_moc4_T42_uniform_20261022__check_output.txt`).

## Deviations and limits

- Batch budget about one hour for five documents; depth was put on the three primary proofs.
- Nothing from the supplement was executed or imported; no supplement JSON was loaded. All polynomials were rebuilt from the paper's definitions with the shared library (cubic (A.2), P, f_j of Appendix D).
- The large stored certificates (1000+-term radial factors, Bernstein tables) were NOT re-expanded symbolically. Instead the radial factor was reconstructed EXACTLY at sample points: G^2 = -Q^2 det Hess C(g - (C/Q) kappa)/4 by exact 14x14 rational determinants, interpolated along lines and square-rooted exactly; sign fixed by the d^9 coefficient = +1. Polynomial identities checked this way at finitely many random rational points are evidence (exact arithmetic, probabilistic), not proofs; they are labelled NUMERIC (exact-at-points).
- Library value check: the library cubic gives C, Q and the base normalisation G(d=1, rest 0) = 1 exactly as the documents state (constant factor 1); no disagreement with the library was found.

## A. Statements

- **Theorem 1** (p. 1) - paper case(s) T_43: For every x in X and pi* with P c(x) >= 0 and (P c(x))_pi = 0 for pi in {1,2,4,6,7,8} (non-strict: extra contacts allowed), Phi_{c(x)}(O) >= 0 for all O in SO(4). No restriction on the spectra: the proof only uses that Re x (height restored) satisfies C = 0 and det Hess C = 0, which holds for every x in X, including repeated eigenvalues and a_i+b_j=0.

The contact set equals the paper's Appendix E set literally (no group element needed), and the generator list equals lib `generators(T)`. The statement is non-strict, covers all O in SO(4) (O(4) equivalently) and all spectra; it is at least as strong as (*).

## B. Dependencies

- [core-doc] old-H capacities (3,5,54),(23,5,54); chains centred at f5 with ends (f0,f3),(f23,f6); SOC tuple (5,2,54) (coordinate images of earlier proved capacities) - used for: as>=1 safe; a>=J safe. Statement match checked: no (CITED).
- [common-lemma] H/chain resource-allocation lemma (proof sketched in the doc; capacity forms not stated) - used for: definition of J and safety for a >= J. Statement match checked: no (CITED).
- [common-lemma] five-chart coverage of the resource lemma ('already proved') - used for: Bernstein tables are computed per KKT chart. Statement match checked: no (CITED).
- [other-doc] supplement uniform/full-radial verifiers: 60 exact tables (L0,B1,B2 coefficients in Y nonnegative) - used for: L0 >= (1-a/J)^2, B1,B2 >= 0. Statement match checked: no (CITED).
- [common-lemma] radial necessity: Re x satisfies C=0 and det Hess C=0 (secant rank <= 12, cf. HMSV) - used for: hard cost must have G=0. Statement match checked: no (CITED).
- [core-doc] symmetry group G (paper Sect. 4.2) - used for: swap f22<->f54. Statement match checked: yes.

## C. Steps

| step | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | generator inventory {0,1,2,3,4,5,6,21,22,23,53,54}, nonminor 22,54 | identity | EXACT | lib generators(T) |
| S2 | four linear identities among f's | identity | EXACT | profiles compared on all 24 permutations |
| S3 | g in G fixing T and swapping f22<->f54 (row/col perms as stated) | identity | EXACT | found ((1,0,3,2),(0,2,1,3),0) in lib convention; 2 such elements |
| S4 | cancellation/absorption gives (1) or a minor sum | case split | PROOF-READ | coefficient bookkeeping with S2 |
| S5 | basis of (1) has rank 8 and common zero set exactly T | identity | EXACT | rank and zero sets |
| S6 | C(g+z kappa)=C+zQ with the displayed C, Q | identity | EXACT | sympy from library cubic; constant 1 |
| S7 | hard region a<J: -MC identity, MT0-bc >= 0 incl. b,c terms; M=0 case; hence C<0, Q!=0 | identity/inequality | EXACT | sympy; J <= 1+1/M is immediate (delta=0) |
| S8 | Q^4 det Hess = -4 Q^2 G^2, G of degree 9 with d^9 coefficient 1 | identity | NUMERIC | exact-at-points: perfect square along 12 lines, base value 1 |
| S9 | five KKT charts and J formulas; chart coverage | citation | CITED | 'already proved five-chart coverage' |
| S10 | 60 tables: coefficients of L0,B1,B2 in Y nonnegative; L0(a=0) >= 1 | certificate | CITED | supplement verifier; not replayed; conclusion sampled in S12 |
| S11 | (1+aS)^2 G(a,S/(1+aS)) = L0+S B1+S^2 B2 and hence G >= (1-as)^2(1-a/J)^2 | identity | PROOF-READ | two-line algebra, re-derived by auditor |
| S12 | G >= (1-as)^2(1-a/J)^2 > 0 on the hard region | inequality | NUMERIC | 12 exact rational hard points (J by 1-D convex minimisation in floats); min ratio 571 |
| S13 | zero-resource and boundary limits (J_eps increases to J; closedness) | analytic | PROOF-READ | argument read; J_eps <= J monotone because adding resource lowers each energy |
| S14 | as>=1 or a>=J safe | citation | CITED | capacities |

## D. Falsification

- NUMERIC: G^2 = -Q^2 detH/4 is a perfect square along 12 exact lines (degree <= 18)
- NUMERIC: G >= (1-as)^2(1-a/J)^2 at 12 exact hard points (seed 20261125), min ratio 571.7
- NUMERIC: rank Hess C(Re x) = 12 at 5 random x in X (supports det=0 necessity)
- NUMERIC: 12 realizable T_43 points (lib stratum sampler, seeds 0-11, margin 1e-3): none hard; min Phi over SO(4) ~ 7e-12 (6 restarts)

Exact replays:
- EXACT: generators and nonminor set of T_43
- EXACT: four generator identities
- EXACT: swap symmetry in G
- EXACT: rank 8 / common zeros of basis (1)
- EXACT: C and Q of the cubic restriction (constant 1)
- EXACT: -MC identity and nonnegativity of MT0-bc at d=1
- EXACT: M=0: C=-(1+p)(1+q)(1+r)
- EXACT: G(base)=1 normalisation

## E. Verdict

**confirmed-modulo-citations**. Primary proof of paper T_43={1,2,4,6,7,8} (literal match): reduction to an 8-generator cone, cubic guard and C,Q formulas replay exactly; the radial-factor bound G >= (1-as)^2(1-a/J)^2 holds at all 12 exact sample points; rests on cited SOC/old-H/chain capacities, the resource lemma and the unreplayed 60 Bernstein tables.

No gap found. The only links not checked by me are the cited capacities, the resource lemma's chart coverage and the stored table certificate; the quantitative bound those tables imply was tested directly with exact radial values and holds with a large margin at every sample.


Unchecked:
- full symbolic replay of G (1398 terms) and of the 60 Bernstein tables (time)
- the KKT chart coverage of the resource lemma
- the forms of the old-H/chain/SOC capacities (other documents)
