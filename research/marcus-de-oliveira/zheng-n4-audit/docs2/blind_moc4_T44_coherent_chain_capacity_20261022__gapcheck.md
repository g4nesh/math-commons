Evidence level: EXACT. The Gram identity of Proposition 1 and the b = 0 Gram of Proposition 2 were replayed by the auditor's own exact code (python-flint, sympy) from the paper's definitions. The node functions came from the supplement's JSON data and served only as a candidate certificate. A floating-point sanity check (NUMERIC) is reported but not relied on.

# Gap check: blind_moc4_T44_coherent_chain_capacity_20261022 (label C18-T43-T45-gapcheck)

PDF: proofs/manin_transfer_20261021/blind_moc4_T44_coherent_chain_capacity_20261022.pdf (2 pages).

## Deviations and limits

- The document never states the eight node functions ("diagonal-coefficient Gram coordinates") or the "ten relations". To test the step I took a candidate set of node functions from the supplement's data files:
  - `blind_moc4_T44_cubic_reduction_20261022.json`: the keys `features` (two blocks of eight cubic polynomials in quaternion coordinates p0..p3, q0..q3) and `left_four` (four signed variable permutations).
  - `blind_moc4_T44_bad18_budget_PSD8_20261022.json`: the key `T` (an 8x8 rational change of basis).

  Both files were read with `json.load`, and the polynomial strings were parsed by my own tokenizer, with no `eval`. No supplement script was executed or imported. I read the companion `.py` only to learn which data keys it uses. Its own assertions were not used.
- I did not use the supplement's `gram_columns`, `cost_columns` or `R` data. Every Leibniz sum was rebuilt from the paper's f_j (moc4lib, Appendix D) and the definition (3.4).
- The quaternion convention (O e_j = p e_j q) was found by trying 8 conventions numerically. Everything after that choice is exact.
- Only the disputed step is examined: Prop. 1 and the b = 0 Gram of Prop. 2. The optimisation in Prop. 2 and the scalar test (1) are outside this check.

## The gap as reported

p.1, proof of Prop. 1 (and the b = 0 Gram in Prop. 2). The document asserts that the Gram matrix G = V V^t of the eight listed vectors in R^3 "satisfies all ten original coefficient relations" in previously verified Gram coordinates. It also asserts that G therefore represents the Leibniz sum of the cost

  g* = p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15 + q f2,

with d = 1, j = r^2(1+v^2), c = k^2(1+v^2), a = e = v^2, b = w^2, q = 0, p = x^2 + y^2. Neither the coordinates nor the relations are given. Only the companion program is cited.

## Classification: FILLED (severity: presentation)

The statement is true and admits a complete, checkable sum-of-squares proof. The missing data (the node functions) can be written down explicitly, and the identity was replayed exactly. The argument follows.

### Setting

For unit quaternions p, q, let O(p,q) be the matrix of x -> p x q on H = R^4 (basis 1, i, j, k). Its entries are bilinear in (p, q). The code checks two identities exactly: O(p,q)^t O(p,q) = |p|^2|q|^2 I and det O(p,q) = (|p|^2|q|^2)^2. So O(p,q) is in SO(4) when |p| = |q| = 1. Every element of SO(4) has this form. This is the standard double cover S^3 x S^3 -> SO(4), x -> p x q', composed with the inversion q' -> q-bar on S^3.

The eight sources and their nodes are:

| node | J | C | E | A | B | Q | P | D |
|---|---|---|---|---|---|---|---|---|
| source | f53 = n13,13 | f12 = n14,24 | f8 = n12,23 | f15 = n12,14 | f7 = n3,1 | f2 = n13,34 | f1 = n12,34 | f18 = h3412,- |

Every source except f18 is a minor vector. By (4.2), each such source has Phi = det(minor)^2 >= 0 on SO(4). This is checked against moc4lib.

### Node functions

Let phi^(beta,gamma) = (phi_1, ..., phi_8) be the eight cubic features of block beta (two blocks), evaluated at the gamma-th signed variable permutation of `left_four` (four copies). That gives eight copies in all. Define the R^8-valued node functions

  m_K = ( sum_l (T^{-1})_{K,l} phi_l^(beta,gamma) )_{(beta,gamma)},   K = J, C, E, A, B, Q, P, D,

and the pairing M_KL = <m_K, m_L> (a sum over the eight copies). The matrix M(p,q) is a Gram matrix, so it is PSD at every point.

### Step 1 (EXACT)

I solved symbolically (sympy) for the conditions on a symmetric matrix M under which sum_{K,L} G_KL M_KL equals the Leibniz-sum target for all r, k, v (r != k) and w != 0. The target is

  r^2(1+v^2) Phi_J + k^2(1+v^2) Phi_C + v^2 Phi_E + v^2 Phi_A + w^2 Phi_B + (x^2+y^2) Phi_P + Phi_18.

The general solution consists of eleven linear relations:

- R1-R4: M_BB, M_CC, M_JJ, M_PP equal Phi_B, Phi_C, Phi_J, Phi_P.
- R5: M_JC = 0.
- R6: 2 M_JP - M_DD = -Phi_18.
- R7: M_CP - M_BP + M_JP = 0.
- R8: M_JE + M_JA + M_DD = Phi_18.
- R9: 2(M_CE + M_CA + M_BP) - M_DD = -Phi_18.
- R10: M_EP + M_AP - M_PD = 0.
- R11: M_EE + M_AA + M_DD + 2M_EA - 2M_ED - 2M_AD = Phi_A + Phi_E.

Substituting R1-R11 back makes the identity hold identically in r, k, v, w. This is the document's "ten relations", up to how they are counted.

M_AA is not forced to equal Phi_A, and M_DD is not forced to equal Phi_18. This is consistent with the document's remark that the D diagonal 1 + v^2 is not the coefficient d = 1. Phi_18 takes negative values, so a naive reading with M_DD = Phi_18 could not work.

### Step 2 (EXACT)

For the node functions above, each relation "M <-> Phi" was checked as the polynomial congruence

  |p|^2 * M(p,q) - 2 * Phi(O(p,q)) == 0   mod (|q|^2 - |p|^2).

Here Phi(O(p,q)) is the Leibniz sum, of degree (4,4), of the unnormalised matrix. The congruence is tested by division with remainder in lex order with q0 first, where the single generator is a Groebner basis. On |p| = |q| = 1 it says M = 2 Phi. All eleven relations pass.

### Step 3 (EXACT, independent of steps 1-2)

At three rational parameter points I checked the full identity

  |p|^2 * sum_{copies} sum_{s=1..3} ( sum_K V_{Ks} m_K )^2 == 2 * Phi_{g*}(O(p,q))   mod (|q|^2 - |p|^2).

The points (r, k, v, w) were (1, -1, 0, 1), (3/2, -2/5, 1/3, 7/4) and (-2, 5/3, -3/7, 1/2).

### Conclusion

For every O in SO(4), write O = O(p,q) with |p| = |q| = 1. Then

  2 Phi_{g*}(O) = sum over 8 copies and s = 1..3 of ( sum_K V_{Ks} m_K(p,q) )^2 >= 0.

Raising p, j, c, e, a, b, q above their threshold values adds nonnegative multiples of Phi of minor vectors, each equal to det(minor)^2 >= 0. This proves Prop. 1 for all real r, k, v with r != k and all w > 0.

The stated diagonal (r^2(1+v^2), k^2(1+v^2), v^2, v^2, w^2, 0, x^2+y^2, 1+v^2) of V V^t was also checked symbolically.

### b = 0 Gram of Prop. 2 (EXACT)

The vectors are J = (0, s, 0), P = (0, -1/s, 0), D = (1, 0, 0) with s = sqrt(j), and all other nodes zero. Then

  sum G_KL M_KL = s^2 M_JJ + M_PP/s^2 - 2 M_JP + M_DD.

By R3, R4 and R6 this equals j Phi_J + Phi_P/j + Phi_18, the Leibniz sum of j f53 + (1/j) f1 + f18. That cost is therefore a sum of squares on SO(4), with no division by w. Checked symbolically.

### NUMERIC sanity check (not relied on)

Phi_{g*} was minimised over SO(4) (expm of skew matrices, BFGS) for 150 random (r, k, v, w), with seed 1 and 6 restarts each. The worst normalised minimum was +6.6e-17, so it touches 0 without going negative. For the example (1/2)f1 + f53 + f12 + f7 + f18 the minimum was +5.9e-14 (30 restarts). Phi_{f18} alone reaches -0.0674.

### Why "presentation"

The step is correct, and the missing ingredient is a finite certificate: eight cubic vector-valued node functions. Once they are printed, a reader can check the claim by a polynomial computation. As written, the document only points to a companion program and to "already independently verified" coordinates that it never defines, so it is not self-contained. The repair is to state O(p,q), the node functions (or the features and T) and relations R1-R11. That is what the files listed below do.

## Files

- `replay2/blind_moc4_T44_coherent_chain_capacity_20261022__gapcheck.py`: exact replay, 21 checks.
- `replay2/blind_moc4_T44_coherent_chain_capacity_20261022__gapcheck_output.txt`: its output, ending in ALL EXACT CHECKS PASS.
