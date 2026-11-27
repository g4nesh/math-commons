# Order-154 replay: reduction via CLMP Theorem 1.4

Sources read (arXiv HTML, via a page-summarising fetch tool, so statements below are faithful paraphrases, not verbatim quotes; see open points): Caragea-Lee-Malikiosis-Pfander (CLMP), arXiv 2505.24326v3, Sections 1, 3, 4 (summary), 5 (summary of proof); and the abstract of 2409.09793v4. No repository material under research/fourier-principal-minors/ other than this folder was read.

## 1. Theorems used

Notation: F_M = (omega^(ij)), 0<=i,j<M, omega a primitive M-th root of unity (CLMP use e^(-2 pi i/M); any primitive root gives Galois-conjugate minors, so nonvanishing of all principal minors is independent of the choice).

- **Theorem 1.4 (CLMP).** Let N = p N' be square-free, p prime (so p does not divide N'), and assume all principal minors of F_{N'} are nonzero in characteristic p. Then every N'-principal minor of F_N is nonzero; in particular every principal minor of F_N is nonzero.
  - N'-principal minor of F_N (Def. 3.2): submatrix F_N[K,L], |K|=|L|, such that for every residue a mod N' the number of elements of K congruent to a equals that of L. A principal submatrix (K=L) is N'-principal.
  - "Nonzero in characteristic p" (Sec. 3.2, Def. 3.6): reduce Z[omega_{N'}] modulo a prime ideal above p (p does not divide N', so unramified; residue field F_{p^f}, f = ord_{N'}(p)); equivalently, fix a primitive N'-th root zeta in F_{p^f} and use (zeta^(ij)). The condition used in the proof (Sec. 5, step 8) is that p does not divide the norm N_{Q(omega)/Q}(det F_{N'}[K,K]). That is: the minor is nonzero modulo EVERY prime above p, equivalently nonzero for every primitive N'-th root zeta in F_{p^f} (all Galois conjugates).
- **Proof skeleton (Sec. 5), for the record.** Lemma 3.5 swaps F_N for F_{N'} (x) F_p (Galois-conjugate minors). Single-residue-class case is a minor of F_p, nonzero by Chebotarev (Thm 1.1). Multi-class case: differentiate the Vandermonde-type factorisation, evaluate at powers of zeta = omega^p, apply the block-determinant Lemma 3.1; get p*unit dividing a quantity whose norm is a product of N(det Z_k) (Z_k leading principal minors of a row/column-permuted F_{N'}), N(zeta)=+-1, N(zeta^a - zeta^b) | power of N' (prime to p), factorial and Vandermonde-difference factors prime to p. Contradiction.
- **Prop. 3.3** (unitary complement): A[K,L] invertible iff A[K^c,L^c] invertible; remark extends it to reductions mod q when det A is nonzero mod q. For F_{14} mod 11, det^2 = +-14^14, a unit mod 11.
- **Prop. 3.7**: norm of D_{K,L} is invariant under x -> ax+b (a a unit) applied to K or L.
- **Thm 1.1 (Chebotarev)**: all minors of F_p nonzero.
- (Context: Prop. 5.4 treats N=6p exactly this way: Theorem 1.4 reduces to F_6 in char p, minors of size <=3.)

## 2. Application to N = 154 = 2*7*11

Three factorisations N = p N': (p,N') = (11,14), (7,22), (2,77). Choose **p = 11, N' = 14** (smallest check: 2^14 subsets; 2^22 and 2^77 are worse).

- p=11 does not divide 14. ord_14(11): 11^2=121=9, 11^3=99=1 mod 14, so f=3. Residue field **F_{11^3} = F_1331**; 1330 = 14*95, so primitive 14th roots exist. Construct F_1331 = F_11[x]/(monic irreducible cubic), take a generator g of the multiplicative group (order 1330), **zeta = g^95**.
- Q(omega_14) has degree phi(14)=6; Gal = Z_14^* = {1,3,5,9,11,13}; Frobenius subgroup <11> = {1,9,11}. So 11 splits into r = 6/3 = 2 primes, each with residue field F_1331. They correspond to the two cosets of <11>: choose zeta and zeta' = zeta^(-1) (13 = -1 is not in <11>). **Both must be checked** (the norm condition covers all primes above 11).
- Matrices: M_zeta = (zeta^(ij))_{0<=i,j<14} over F_1331 (no normalisation; unnormalised entries are algebraic integers, no factor 1/sqrt(14) is used; coefficient ring Z[omega_14] -> F_1331). M_{zeta^-1} likewise.
- Statement to certify (**Hyp(14,11)**): for every S subset of Z_14 and each zeta in {zeta, zeta^-1}, det M[S,S] != 0 in F_1331 (empty determinant = 1).

## 3. Exact finite check set and validity of symmetry reductions

Full set: 2^14 = 16384 subsets S, times 2 roots = 32768 determinants (sizes <=14). This is exhaustive with no symmetry argument and is the recommended primary run (trivial cost).

Optional quotients (derived here, use only as cross-check or speedup):
1. **Translation** S -> S+b: zeta^((i+b)(j+b)) = zeta^(ij) zeta^(bi) zeta^(bj) zeta^(b^2); diagonal scalings by nonzero elements, so det changes by a nonzero factor. Valid for the same zeta.
2. **Unit multiplication** S -> aS, a in Z_14^*: entries become zeta^(a^2 ij), and a^2 in {1,9,11} = <11> (squares of units mod 14 are 1,9,11,11,9,1). zeta^(11) is the Frobenius image of zeta, so det_{aS}(zeta) = Frob^k(det_S(zeta)); zero iff zero. So for a FIXED zeta, all of Aff(Z_14) (order 14*6 = 84, including S -> -S) preserves vanishing. (Prop. 3.7 gives the same at the norm level.) Number of orbits: not computed by hand (Burnside count by program).
3. **Complement** (Prop. 3.3 / Jacobi): F^{-1} = (1/14) M_{zeta^-1}, so det M_zeta[S^c,S^c] = +- det M * 14^(-|S|) * det M_{zeta^-1}[S,S]; complement swaps zeta <-> zeta^-1 (the two primes), and det M != 0 mod 11. Hence the check on {all sizes, both zeta} equals {sizes 0..7, both zeta}: sum_{k=0..7} C(14,k) = (16384 + 3432)/2 = 9908 subsets, times 2 roots = 19816 determinants. Note: complement does NOT reduce the zeta count unless one checks all sizes for one zeta; do not combine complement with "only zeta".
4. Sizes 0 and 1 are trivial (1 and 1); size 2 and arithmetic progressions are handled by CLMP Prop. 4.3 but need not be assumed.

## 4. Implication chain

(A) Hyp(14,11) [all principal minors of F_14 nonzero in char 11, both primes above 11, via norm criterion Def. 3.6 (3)<=>(4)] 
=> (B) Theorem 1.4 with N = 11*14, p=11, N'=14: every 14-principal minor of F_154 is nonzero over C 
=> (C) every principal minor of F_154 (principal => 14-principal, Def. 3.2) is nonzero; the empty minor is 1. Over any primitive 154th root: nonvanishing is Galois-invariant, so independent of zeta.

Failed check: if some S has det M_zeta[S,S] = 0 in F_1331 (for either root), then the hypothesis of Thm 1.4 fails for (p,N')=(11,14). This only kills this lifting route (and means 11 divides the norm of an algebraic integer); it says nothing about det F_154[S',S'] over C, nor about F_14 over C (Theorem 1.3 gives 2p, i.e. F_14 itself fine over C). Fallbacks: (7,22) (check 2^22 subsets in F_{7^10}; ord_22(7)=10=phi(22), one prime) or (2,77) (infeasible). Failure of all routes would still not disprove the claim. A successful run does not need the other routes.

## 5. Algorithm specification (exact arithmetic only)

1. Find monic irreducible cubic c(x) over F_11 (no roots mod 11; test all 11 values). Represent elements as triples in F_11; multiply mod c.
2. Find g of order 1330: test g^(1330/2), g^(1330/5), g^(1330/7), g^(1330/19) != 1 (1330 = 2*5*7*19). Set zeta = g^95; assert zeta has order exactly 14 (zeta^7 != 1, zeta^2 != 1, zeta^14 = 1). zeta_inv = zeta^13.
3. Precompute powers z[k] = zeta^(k mod 14); matrix entry (i,j) = z[(i*j) mod 14].
4. For each root in {zeta, zeta_inv}, each bitmask S in [0,2^14): build M[S,S] (size s), compute det by Gaussian elimination over F_1331 (needs only field inverse via a^(1329)); det of empty = 1. Record any S with det=0.
5. Output: count checked, list of zero cases. Pass criterion: zero cases empty for both roots.
6. Self-tests:
   - Field: assert zeta has order 14, sum_{j} zeta^(ij) = 0 for i not 0 mod 14, and det of full M_zeta nonzero (also det^2 = +-14^14 mod 11 up to sign convention).
   - Empty set gives 1; singletons give 1; full set equals det M (check M * M_{zeta^-1} = 14 I).
   - Complement identity: check det M_zeta[S^c,S^c] = eps * det M * 14^(-|S|) * det M_{zeta^-1}[S,S] with eps = +-1 for all S (cross-validates the determinant routine and the two-prime bookkeeping).
   - Symmetry: dets for S, S+b, aS agree up to nonzero (translation, unit) and vanish together.
   - Zero detection: rerun the same determinant routine with N not square-free, e.g. over F_p with a primitive 4th root, S={0,2}: matrix [[1,1],[1,1]] det 0; N=9, S={0,3,6}: all-ones det 0 (N=4, N=9 are cases where principal minors vanish, 2409.09793). Also a char-p failure: per CLMP Thm 4.6 (as summarised), F_11 fails the 3-Chebotarev property (3 divides (11+1)/4); an F_11 minor vanishing mod 3 should be detected by the routine over F_{3^5} (ord_11(3)=5). (Unverified as to which minor; use only as a smoke test.)
   - Cross-check of the reduction on a small case: N=14=7*2, p=7, N'=2 (F_2 principal minors 1,1,-1,-2 nonzero mod 7); N=22 (p=11,N'=2), N=6*p, N=42 (p=7,N'=6, F_6 in char 7, 7=1 mod 6 so the field is F_7) with minors up to size 3. Compare, for N=14 and N=22, against a direct complex check done exactly (e.g. determinant over F_q with q = 1 mod N large prime: nonzero there proves nonzero over C; or exact Z[omega] norms) for all 2^14 / 2^22 subsets, and check that every principal minor passes in both, and the char-p check implies it. Never use floating point.

## 6. Open points (not confirmed, flagged)

1. My reading of the paper came through a summarising fetcher, not the raw text: exact numbering and phrasing of Def. 3.2, Def. 3.6, Lemma 3.5 and the Section 5 proof are second-hand. The last 8.6k characters were read only as bibliography. Someone should read Section 5 directly, especially (5.12) and the claim that the Z_k are principal minors (not Galois-twisted ones) of F_{N'} with root omega^p; if Z_k could be Galois-twisted principal minors, the all-primes-above-p reading (checking both zeta and zeta^-1) covers it anyway, which is why I check both.
2. "Nonzero in characteristic p" is read as the norm condition (all primes above p). Single-prime checking would be insufficient unless complement/units symmetry closes it; here it does not (-1 not in <11>).
3. Whether Theorem 1.4 requires p < N' or other ordering: none stated in the summaries; hypotheses as stated above only.
4. Theorem 4.6 smoke test details (which minor, which field) unverified.
5. The Burnside orbit count under Aff(Z_14) (with or without complement) is not given; only the unreduced counts 16384/9908 per root are exact.
6. Orbit/quotient claims in section 3 are my own derivations, not from the paper (except Prop. 3.3 and 3.7).
