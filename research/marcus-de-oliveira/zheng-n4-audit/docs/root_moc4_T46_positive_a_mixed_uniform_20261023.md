# Audit of `root_moc4_T46_positive_a_mixed_uniform_20261023`

Evidence level: EXACT for all twelve chart sign tables (own replay, ten with term counts equal to the document's; the two square tables re-proved by an own AM-GM certificate), for G (1556 terms) and C0, and for the capacity algebra; CITED for radial necessity, the six-ray coherent capacity and the (2) capacity; NUMERIC not run for T46.

PDF: `proofs/manin_transfer_20261021/root_moc4_T46_positive_a_mixed_uniform_20261023.pdf` (3 pages). Auditor label C24-T47.

## Deviations and limits

- The two 'above e0, b > q' tables were cleared by my own clearing (1,869,588 and 50,618 terms vs the document's 1,740,714-term residual and 45,903 terms); I did not use the document's 1718/56 squares but built my own AM-GM pair certificate (1829 and 58 positive terms absorbed).
- The supplement name T46 is the paper's T_47 (contact set {1,3,5,6,7,8}, identity element). N_T has 10 generators (f0,f1,f2,f3,f4,f5,f6,f27,f28,f59); the 8-ray 'canonical source' omits f0,f3, whose absorption via f28+f0 = f1+f2+f5+f27 and f28+f3 = f4+f5+f27 (both replayed EXACT) is done in the T46 assembly note, not here.
- G (the radial polynomial) was rebuilt from the paper's definitions only: w(g,h) = sum of the source coefficients times the Appendix D vectors f_j plus h*kappa, C(w) = C0 + hQ, G^2 := -Q^2 det Hess C(w(g,-C0/Q))/4 (d = 1). The 14x14 determinant was computed by Kronecker substitution (degree bounds found by univariate interpolation) and every G was cross-checked against an independent exact Fraction determinant at random rational points. G is only defined up to a global sign; the sign was fixed once (positive c^2-coefficient at the all-ones point) and the document's sign claims were then reproduced.
- The chart clearings are my own (each variable cleared by its own degree, or by a common degree where the document says so). Where my term counts equal the document's this is noted; where they differ, the sign conclusion was re-proved for my clearing.
- Nothing in the supplement was executed; its JSON certificates were not used.
- Cited inputs (radial necessity G = 0 for spectrally realizable sources, the negative-Phi exclusion, and the pointwise capacities) were not audited here and are CITED.
- Falsification: the exact certificates prove the radial claims outright. A floating-point realizability test was run only on a few points (budget): on 2 random realizable T_65-stratum points, |G| relative to the absolute-term sum was <= 6e-12 (consistent with G = 0 being necessary), Phi > 0 and H6 > 0. No random realizable point fell in any hard region (too few samples to mean much).

## A. Statements

**Theorem 1** (contact set [[1, 3, 5, 6, 7, 8]]). For x in X with c(x) in N_{T_47} written as the canonical source g = q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4 (all >= 0), if qb <= d^2 and max(q,b) > d, then Phi_{c(x)}(O) >= 0 for all O in SO(4), for every a >= 0. Proof covers d>0 and q,p,j,c,e,b>0 itself; zero coefficients and d=0 are delegated (transfers, minor sum). Spectra: all realizable ones (no genericity assumed), via the necessity G=0.

## B. Dependencies

| kind | reference | used for | match checked |
|---|---|---|---|
| core-doc | T46 radial necessity (det Hess = -4 G9^2/Q^2 => G9 = 0 at realizable sources, incl. collision loci) | contradiction in the hard region | False |
| common-lemma | final_referee_T46_six_ray_coherent_capacity_20261023 (six-ray coherent coefficient bound with beta=-r, v=-rt/q) | capacity (1) p >= P* | False |
| other-doc | root_moc4_T46_old_four_capacity (capacity (2) ab(e-1) >= a+b+e+1) | capacity (2) | False |
| other-doc | T46 assembly: zero-coefficient transfers (to T28,T35,T0,T38,T41,T6) and the a=0 face T u {14} | zero coefficients | False |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Contact set {1,3,5,6,7,8} = paper T_47; source rays are 8 of the 10 generators | identity | EXACT | lib |
| 2 | C0 = -(q+1)L - jM, M=(c+1)(a+1+e), L=M+b(a+1); Q = jp - j(c+e+1) - (c+1)(e+1) | identity | EXACT | own computation |
| 3 | G9 has 1556 terms, degree 9 | identity | EXACT | own Kronecker determinant + point checks |
| 4 | Realizable => G9 = 0 | citation | CITED | not audited |
| 5 | Capacity (1) follows from the cited six-ray bound with t = eta j/(D0+eta^2 j): thresholds become P* and rho^2 eta^2 j/N < c | identity | EXACT | sympy simplification of the three substitutions and of N - (rho^2/c)eta^2 j |
| 6 | Cited six-ray bound itself (both mixed orders, r=1 limit) | citation | CITED | - |
| 7 | Capacity (2) | citation | CITED | - |
| 8 | Charts (3)-(4) with S>0 are onto the complement of (1),(2) in the sector (12 charts) | case split | PROOF-READ | b(e-1)-1 > 0 iff e > e0; (q,b)=(1+A, r^2/(1+A)) or swapped covers qb=r^2<=1, max>1; p = P* S/(1+S) |
| 9 | Ten cleared tables have only negative coefficients and an H-free monomial | inequality | EXACT | own clearing; term counts 970413, 347708, 1861444, 986482, 327932, 19309, 6782, 34124, 23910, 8054 equal the document's |
| 10 | Two 'above, b>q' tables are negative via squares | inequality | EXACT | own AM-GM certificate on my clearing (different from the document's squares) |
| 11 | All chart denominators positive | inequality | EXACT | positive-coefficient polynomials in chart variables |

## D. Falsification

- none for T46 within budget (the exact tables settle the radial claim)

## E. Verdict

**confirmed-modulo-citations**. Proves the qb<=1, max(q,b)>1 sector of supplement T46 (= paper T_47) by excluding realizable sources below two cited capacities; all 12 chart sign tables replayed exactly (10 with identical term counts, 2 by own AM-GM certificates); rests on cited radial necessity and capacities.


Unchecked:

- cited six-ray coherent capacity and capacity (2)
- radial necessity
- zero-coefficient transfers and a=0 face (assembly)

Replay: `replay/root_moc4_T46_positive_a_mixed_uniform_20261023__check.py`, run from `replay/` as `python -I -B <script> ../lib "=1"` and `python -I -B <script> ../lib "<1"` (two runs to fit the 600 s job limit; both are in the output file).

